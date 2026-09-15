"""Twirled Readout Error eXtinction (TREX) style readout mitigation.

Measurement twirling (random Pauli conjugation of the readout, implemented in
``ExpectationEvaluator`` via ``twirls``) converts coherent readout errors into
stochastic ones; the stochastic assignment error is then removed by inverting
a calibration matrix estimated from prepared computational basis states.

Two calibration strategies:

* ``full``      - prepare all 2^n basis states (exact, for small n)
* ``tensored``  - per-qubit 2x2 assignment matrices, Kronecker product
                  (scales linearly, standard choice)

Correction is a non-negative least-squares inversion
(``scipy.optimize.nnls``) so the corrected distribution is a valid
probability vector.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional

import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from scipy.optimize import nnls


@dataclass
class ReadoutMitigator:
    """Inverts a calibration matrix to correct measured bitstring counts."""

    cal_matrix: np.ndarray  # C[i, j] = P(measured i | prepared j)
    num_qubits: int

    @classmethod
    def calibrate(
        cls,
        noise_model,
        num_qubits: int,
        shots: int = 8192,
        full: bool = False,
        seed: Optional[int] = None,
    ) -> "ReadoutMitigator":
        """Estimate the assignment-error matrix on the noisy simulator."""
        sim = AerSimulator(method="automatic", noise_model=noise_model)
        if full and num_qubits <= 6:
            return cls._calibrate_full(sim, num_qubits, shots, seed)
        return cls._calibrate_tensored(sim, num_qubits, shots, seed)

    @staticmethod
    def _calibrate_tensored(sim, n: int, shots: int, seed) -> "ReadoutMitigator":
        matrices = []
        for q in range(n):
            m = np.zeros((2, 2))
            for prepared in (0, 1):
                qc = QuantumCircuit(n)
                if prepared:
                    qc.x(q)
                qc.measure_all()
                counts = sim.run(qc, shots=shots, seed_simulator=seed).result().get_counts(0)
                total = sum(counts.values())
                for bitstring, cnt in counts.items():
                    measured = int(bitstring[n - 1 - q])
                    m[measured, prepared] += cnt / total
            matrices.append(m)
        # index i = sum_q bit_q * 2^q, so qubit 0 is the LOWEST-order bit:
        # cal = M_{n-1} (x) ... (x) M_1 (x) M_0
        cal = matrices[-1]
        for m in reversed(matrices[:-1]):
            cal = np.kron(cal, m)
        return ReadoutMitigator(cal_matrix=cal, num_qubits=n)

    @staticmethod
    def _calibrate_full(sim, n: int, shots: int, seed) -> "ReadoutMitigator":
        cal = np.zeros((2 ** n, 2 ** n))
        for prepared in range(2 ** n):
            qc = QuantumCircuit(n)
            for q in range(n):
                if (prepared >> q) & 1:
                    qc.x(q)
            qc.measure_all()
            counts = sim.run(qc, shots=shots, seed_simulator=seed).result().get_counts(0)
            total = sum(counts.values())
            for bitstring, cnt in counts.items():
                measured = int(bitstring, 2)
                cal[measured, prepared] += cnt / total
        return ReadoutMitigator(cal_matrix=cal, num_qubits=n)

    def apply(self, counts: Dict[str, int], num_qubits: int) -> Dict[str, float]:
        """Correct measured counts -> corrected probability distribution."""
        if num_qubits != self.num_qubits:
            raise ValueError(
                f"mitigator calibrated for {self.num_qubits} qubits, got {num_qubits}"
            )
        total = float(sum(counts.values()))
        if total == 0:
            return counts
        meas = np.zeros(2 ** num_qubits)
        for bitstring, cnt in counts.items():
            meas[int(bitstring, 2)] = cnt / total
        corrected, _ = nnls(self.cal_matrix, meas)
        corrected = corrected / corrected.sum()
        # Aer count strings are the binary representation of the state index
        return {
            format(i, f"0{num_qubits}b"): float(p)
            for i, p in enumerate(corrected)
        }