"""Low-distance repetition codes tailored to biased noise.

For phase-flip-dominated (Z-biased) channels -- the common NISQ situation
(T2 < 2*T1) -- the correct code family is the **X-basis repetition code**
(stabilizers X_i X_{i+1}): it corrects Z errors. For bit-flip-dominated
channels the Z-basis code (stabilizers Z_i Z_{i+1}) is used instead. The
hardware-aware thesis implemented here is that under Z bias, the X-basis code
achieves a lower logical error rate than the Z-basis code at identical
physical error rates, and that the logical rate improves with code distance d.

Layout: data qubits on even indices (0, 2, ..., 2d-2), syndrome ancillas on
odd indices (1, 3, ..., 2d-3). Total physical qubits = 2d - 1.

Syndrome extraction (X-basis code, measuring X_i X_{i+1}):
    H(anc); CNOT(data_j -> anc); CNOT(data_{j+1} -> anc); H(anc); measure(anc)

Decoding is exact minimum-weight matching on the 1D syndrome chain: a
syndrome pattern s (length d-1) is consistent with the two error vectors e and
e XOR 111...1; the lower-weight one is chosen (exact below d/2 physical
errors).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence

import numpy as np
from qiskit import QuantumCircuit
from qiskit_aer import AerSimulator
from qiskit_aer.noise import NoiseModel, ReadoutError, pauli_error

from ..noise_models import biased_depolarizing, uniform_depolarizing_2q


@dataclass
class RepetitionCode:
    """A distance-d repetition code in the X or Z basis."""

    distance: int
    basis: str = "X"  # 'X': stabilizers X_i X_{i+1} (corrects Z errors)
    # 'Z': stabilizers Z_i Z_{i+1} (corrects X errors)

    def __post_init__(self):
        if self.distance < 3 or self.distance % 2 == 0:
            raise ValueError(f"distance must be an odd integer >= 3, got {self.distance}")
        if self.basis not in ("X", "Z"):
            raise ValueError(f"basis must be 'X' or 'Z', got {self.basis!r}")

    @property
    def num_data(self) -> int:
        return self.distance

    @property
    def num_ancilla(self) -> int:
        return self.distance - 1

    @property
    def num_qubits(self) -> int:
        return 2 * self.distance - 1

    def data_qubits(self) -> List[int]:
        return list(range(0, 2 * self.distance - 1, 2))

    def ancilla_qubits(self) -> List[int]:
        return list(range(1, 2 * self.distance - 2, 2))

    # ------------------------------------------------------------------ #
    def encode(self, logical_value: int = 0) -> QuantumCircuit:
        """Prepare |0_L> (or |1_L>): all-zeros / all-ones in the code basis."""
        qc = QuantumCircuit(self.num_qubits)
        data = self.data_qubits()
        if self.basis == "Z":
            if logical_value:
                for q in data:
                    qc.x(q)
        else:  # X-basis: |+...+> / |-...->
            for q in data:
                qc.h(q)
            if logical_value:
                for q in data:
                    qc.x(q)
        return qc

    def syndrome_extraction(self, reset_ancillas: bool = True) -> QuantumCircuit:
        """One round of syndrome extraction on fresh (or reset) ancillas."""
        qc = QuantumCircuit(self.num_qubits, self.num_ancilla)
        data = self.data_qubits()
        anc = self.ancilla_qubits()
        for j, a in enumerate(anc):
            if self.basis == "X":
                qc.h(a)
                qc.cx(data[j], a)
                qc.cx(data[j + 1], a)
                qc.h(a)
            else:
                qc.cx(data[j], a)
                qc.cx(data[j + 1], a)
            qc.measure(a, j)
            if reset_ancillas:
                qc.reset(a)
        return qc

    def logical_measurement(self) -> QuantumCircuit:
        """Measure the logical value: data qubits in the code basis."""
        qc = QuantumCircuit(self.num_qubits, self.distance)
        data = self.data_qubits()
        if self.basis == "X":
            for j, q in enumerate(data):
                qc.h(q)
                qc.measure(q, j)
        else:
            for j, q in enumerate(data):
                qc.measure(q, j)
        return qc

    # ------------------------------------------------------------------ #
    def decode(self, syndrome: Sequence[int], data_bits: Sequence[int]) -> int:
        """Decode a syndrome + data measurement into the corrected logical bit.

        ``syndrome`` is the accumulated parity pattern (length d-1);
        ``data_bits`` are the raw logical-basis measurements of the data
        qubits. Returns the corrected logical value (0 or 1).
        """
        s = np.asarray(syndrome, dtype=int)
        if s.shape[0] != self.num_ancilla:
            raise ValueError(f"expected {self.num_ancilla} syndrome bits")
        bits = np.asarray(data_bits, dtype=int)
        if bits.shape[0] != self.num_data:
            raise ValueError(f"expected {self.num_data} data bits")

        # min-weight error e with e_j XOR e_{j+1} = s_j, e_0 = 0
        e = np.zeros(self.num_data, dtype=int)
        for j in range(self.num_ancilla):
            e[j + 1] = e[j] ^ s[j]
        alt = 1 - e  # e XOR 111...1 (logical-X equivalent)
        correction = e if e.sum() <= alt.sum() else alt

        corrected = bits ^ correction
        # majority vote over the corrected data
        return int(1 if corrected.sum() >= (self.num_data + 1) // 2 else 0)

    # ------------------------------------------------------------------ #
    def gates_per_cycle(self) -> Dict[str, int]:
        """Physical gate counts of one syndrome-extraction cycle (overhead)."""
        counts = {"h": 0, "cx": 0, "measure": 0, "reset": 0}
        for op in self.syndrome_extraction().data:
            name = op.operation.name
            if name in counts:
                counts[name] += 1
        return counts

    def physical_gates_per_cycle(self) -> int:
        return sum(self.gates_per_cycle().values())


# --------------------------------------------------------------------------- #
def _memory_experiment_circuit(
    code: RepetitionCode,
    rounds: int,
    ancilla_error_rate: float = 0.0,
    readout_error: float = 0.0,
    logical_value: int = 0,
) -> QuantumCircuit:
    """Full memory-experiment circuit: encode, noisy rounds, logical measure.

    Classical layout (single register):
        clbits 0..d-1            : final data measurement (bit j = data qubit 2j)
        clbits d-1 .. d-1+d-1    : syndrome of round 0, then rounds 1..
    """
    n_data, n_anc = code.num_data, code.num_ancilla
    qc = QuantumCircuit(code.num_qubits, n_data + rounds * n_anc)
    qc.compose(code.encode(logical_value), inplace=True)
    data = code.data_qubits()
    anc = code.ancilla_qubits()
    for r in range(rounds):
        # idle on every qubit: data carry the logical state; ancillas decay too
        for q in range(code.num_qubits):
            qc.id(q)
        for j, a in enumerate(anc):
            if code.basis == "X":
                qc.h(a)
                qc.cx(data[j], a)
                qc.cx(data[j + 1], a)
                qc.h(a)
            else:
                qc.cx(data[j], a)
                qc.cx(data[j + 1], a)
            qc.measure(a, n_data + r * n_anc + j)
            qc.reset(a)
    for j, q in enumerate(data):
        if code.basis == "X":
            qc.h(q)
        qc.measure(q, j)
    return qc


def _build_noise_model(
    code: "RepetitionCode",
    p_x: float,
    p_y: float,
    p_z: float,
    ancilla_error_rate: float,
    readout_error: float,
) -> NoiseModel:
    """Circuit-level noise for a memory experiment.

    Idle (T1/T2-style) noise is applied only to the *data* qubits, which are
    the ones holding the logical state across rounds; freshly-reset ancillas
    pick up errors through the CX gates of the syndrome extraction instead.
    This matches the physical picture (ancillas are reprepared each round) and
    avoids an artificial syndrome-readout floor.
    """
    model = NoiseModel(basis_gates=["id", "cx", "h", "x", "measure", "reset"])
    idle_err = pauli_error(
        [("X", p_x), ("Y", p_y), ("Z", p_z), ("I", 1 - p_x - p_y - p_z)]
    )
    if p_x + p_y + p_z > 0:
        for q in code.data_qubits():
            model.add_quantum_error(idle_err, ["id"], [q])
    if ancilla_error_rate > 0:
        model.add_all_qubit_quantum_error(
            uniform_depolarizing_2q(ancilla_error_rate), ["cx"]
        )
    if readout_error > 0:
        model.add_all_qubit_readout_error(
            ReadoutError(
                [[1 - readout_error, readout_error], [readout_error, 1 - readout_error]]
            )
        )
    return model


def memory_experiment(
    code: RepetitionCode,
    p_x: float = 0.0,
    p_y: float = 0.0,
    p_z: float = 0.05,
    rounds: int = 1,
    shots: int = 2000,
    ancilla_error_rate: float = 0.0,
    readout_error: float = 0.0,
    logical_value: int = 0,
    seed: Optional[int] = None,
) -> Dict:
    """Estimate the logical error rate of a repetition code memory.

    Returns a dict with ``logical_error_rate``, ``logical_shots`` and the raw
    per-shot decoded outcomes. The physical per-round error budget on data
    qubits is ``p_x + p_y + p_z``.
    """
    noise = _build_noise_model(code, p_x, p_y, p_z, ancilla_error_rate, readout_error)
    sim = AerSimulator(method="automatic", noise_model=noise)

    circuit = _memory_experiment_circuit(
        code, rounds, ancilla_error_rate, readout_error, logical_value
    )
    result = sim.run(circuit, shots=shots, seed_simulator=seed).result()
    counts = result.get_counts(0)

    n_data = code.num_data
    n_clbits = n_data + rounds * code.num_ancilla
    errors = 0
    total = 0
    for bitstring, cnt in counts.items():
        # bitstring is little-endian: char position p = clbit (n_clbits - 1 - p)
        data_bits = [int(bitstring[n_clbits - 1 - j]) for j in range(n_data)]
        acc = [0] * code.num_ancilla
        for r in range(rounds):
            for j in range(code.num_ancilla):
                cl = n_data + r * code.num_ancilla + j
                acc[j] ^= int(bitstring[n_clbits - 1 - cl])
        decoded = code.decode(acc, data_bits)
        if decoded != logical_value:
            errors += cnt
        total += cnt

    return {
        "logical_error_rate": errors / total if total else 0.0,
        "logical_shots": total,
        "physical_error_per_round": p_x + p_y + p_z,
        "p_z_over_p_x": p_z / p_x if p_x > 0 else float("inf"),
    }