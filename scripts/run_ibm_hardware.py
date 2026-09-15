#!/usr/bin/env python3
"""Hardware validation: run the H2 VQE estimation protocol on real IBM hardware.

Compares, on the SAME circuit and observable:
  1. exact energy (classical diagonalization)
  2. our simulator prediction (asymmetric calibrated noise model, density matrix)
  3. the real device, raw sampled
  4. the real device + TREX-style readout mitigation (calibrated ON DEVICE)

Usage:
  export QISKIT_IBM_TOKEN=...      # free account: https://quantum.cloud.ibm.com
  python scripts/run_ibm_hardware.py              # real device
  python scripts/run_ibm_hardware.py --dry-run    # identical protocol on Aer

Free Open Plan runtime suffices: the H2 protocol is 8 calibration + 15
estimation circuits of a small ansatz.
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np
from qiskit import QuantumCircuit

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from qec_vqe.hamiltonian import build_hamiltonian                     # noqa: E402
from qec_vqe.ansatz import uccsd_ansatz                               # noqa: E402
from qec_vqe.expectation import (ExpectationEvaluator, build_groups,  # noqa: E402
                                 expectation_from_counts)
from qec_vqe.mitig.readout import ReadoutMitigator                    # noqa: E402
from qec_vqe.noise_models import build_noise_model                    # noqa: E402
from qec_vqe.benchmark import run_m1, TODAY_PROFILE                   # noqa: E402


def measurement_circuits(bound: "QuantumCircuit", groups):
    """Per-group measurement circuits: ansatz + basis rotation + measure."""
    circuits = []
    for g in groups:
        qc = bound.copy()
        qc.compose(g.rotation, inplace=True)
        qc.measure_all()
        circuits.append(qc)
    return circuits


def calibration_circuits(num_qubits: int):
    """Per-qubit |0>/|1> preparation circuits (2*nq circuits, in order)."""
    circuits = []
    for q in range(num_qubits):
        for prepared in (0, 1):
            qc = QuantumCircuit(num_qubits)
            if prepared:
                qc.x(q)
            qc.measure_all()
            circuits.append(qc)
    return circuits


def mitigator_from_device_counts(counts_list, num_qubits: int) -> ReadoutMitigator:
    """Tensored assignment matrix from per-qubit |0>/|1> device counts.

    counts_list order must match calibration_circuits(): for q in range(n),
    [prep0, prep1]. Kron order matches ReadoutMitigator._calibrate_tensored.
    """
    matrices = []
    for q in range(num_qubits):
        m = np.zeros((2, 2))
        for prepared in (0, 1):
            counts = counts_list[2 * q + prepared]
            total = sum(counts.values())
            for bitstring, cnt in counts.items():
                measured = int(bitstring[num_qubits - 1 - q])
                m[measured, prepared] += cnt / total
        matrices.append(m)
    cal = matrices[-1]
    for mm in reversed(matrices[:-1]):
        cal = np.kron(cal, mm)
    return ReadoutMitigator(cal_matrix=cal, num_qubits=num_qubits)


def energy_from_count_list(counts_list, groups, num_qubits: int,
                           mitigator: ReadoutMitigator | None = None,
                           offset: float = 0.0) -> float:
    energy = offset  # nuclear repulsion + frozen-core constant (not in Pauli terms)
    for counts, g in zip(counts_list, groups):
        if mitigator is not None:
            counts = mitigator.apply(counts, num_qubits)
        energy += expectation_from_counts(counts, g, num_qubits)
    return float(energy)


def run_protocol(count_fn, groups, num_qubits: int, cal_circuits, est_circuits,
                 shots: int, offset: float = 0.0):
    """count_fn(circuits, shots) -> list of counts dicts. Returns raw/TREX energies."""
    cal_counts = count_fn(cal_circuits, shots)
    rm = mitigator_from_device_counts(cal_counts, num_qubits)
    est_counts = count_fn(est_circuits, shots)
    raw = energy_from_count_list(est_counts, groups, num_qubits, offset=offset)
    trex = energy_from_count_list(est_counts, groups, num_qubits,
                                 mitigator=rm, offset=offset)
    print(f"readout calibration: joint P(00..0|00..0)={rm.cal_matrix[0, 0]:.3f}")
    return raw, trex


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--shots", type=int, default=4000)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--backend", type=str, default=None)
    args = ap.parse_args(argv)

    H = build_hamiltonian("H2")
    exact = H.exact_ground_energy()
    ansatz = uccsd_ansatz(H)
    theta = np.asarray(run_m1("H2")["params"])
    bound = ansatz.bind(theta)

    groups = build_groups(H.pauli_op)
    num_qubits = H.num_qubits
    cal_circuits = calibration_circuits(num_qubits)
    est_circuits = measurement_circuits(bound, groups)

    if args.dry_run:
        from qiskit_aer import AerSimulator
        nm = build_noise_model(TODAY_PROFILE)
        sim = AerSimulator(noise_model=nm, seed_simulator=11)

        def count_fn(circuits, shots):
            return [sim.run(qc, shots=shots).result().get_counts()
                    for qc in circuits]
    else:
        token = os.environ.get("QISKIT_IBM_TOKEN")
        if not token:
            print("Set QISK_IBM_TOKEN env var (free at "
                  "https://quantum.cloud.ibm.com) or use --dry-run.")
            return 1
        from qiskit_ibm_runtime import QiskitRuntimeService, SamplerV2 as Sampler
        service = QiskitRuntimeService(
            channel="ibm_quantum_platform", token=token)
        if args.backend:
            backend = service.backend(args.backend)
        else:
            cands = [b for b in service.backends()
                     if b.num_qubits >= 5 and b.status().operational]
            backend = min(cands, key=lambda b: b.status().pending_jobs)
        print(f"[hardware] backend: {backend.name} ({backend.num_qubits} qubits)")
        sampler = Sampler(mode=backend)

        def count_fn(circuits, shots):
            res = sampler.run(circuits, shots=shots).result()
            return [res[i].data.meas.get_counts() for i in range(len(circuits))]

    print(f"exact (diagonalization) : {exact:.6f} Ha")
    ev = ExpectationEvaluator(H, mode="exact_noisy",
                              noise_model=build_noise_model(TODAY_PROFILE))
    sim_pred = ev.energy(ansatz.circuit, theta)
    print(f"simulator prediction    : {sim_pred:.6f} Ha "
          f"({(sim_pred - exact) * 1e3:+.1f} mHa)")

    offset = float(H.nuclear_repulsion_energy + H.energy_offset)
    raw, trex = run_protocol(count_fn, groups, num_qubits, cal_circuits,
                             est_circuits, args.shots, offset)
    label = "Aer sampled" if args.dry_run else "device raw"
    print(f"{label:23s}: {raw:.6f} Ha ({(raw - exact) * 1e3:+.1f} mHa)")
    print(f"{'TREX-mitigated':23s}: {trex:.6f} Ha ({(trex - exact) * 1e3:+.1f} mHa)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
