"""Expectation-value evaluation for VQE circuits.

Three evaluation modes:

* ``exact``        - noiseless statevector simulation (exact, no shot noise)
* ``exact_noisy``  - density-matrix simulation with a noise model
                     (exact noisy expectation, deterministic, no sampling)
* ``sampled``      - shot-based sampling with measurement twirling and an
                     optional readout-error filter (TREX)

Terms are grouped by measurement basis (see ``hamiltonian.pauli_groups``) so a
single circuit run estimates all terms sharing the same axis pattern.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable, Dict, List, Optional, Sequence, Union

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.quantum_info import SparsePauliOp
from qiskit_aer import AerSimulator

from .hamiltonian import MolecularHamiltonian, pauli_groups

#: callable(counts, num_qubits) -> corrected counts (e.g. TREX filter)
ReadoutFilter = Callable[[Dict[str, int], int], Dict[str, int]]


@dataclass
class Term:
    """A single Pauli term: P = coeff * X^x Z^z (label phases folded into coeff)."""

    coeff: complex
    label: str
    support: List[int]  # qubit indices with non-identity Pauli


@dataclass
class Group:
    """Terms sharing one measurement basis (same per-qubit axis pattern)."""

    axis_pattern: str
    terms: List[Term]
    rotation: QuantumCircuit  # appended before measurement (H / Sdg-H)


def _axis(char: str) -> str:
    return char if char in "XYZ" else "I"


def _qubit_from_label_pos(num_qubits: int, pos: int) -> int:
    """qiskit Pauli labels are little-endian: the rightmost char is qubit 0."""
    return num_qubits - 1 - pos


def build_groups(pauli_op: SparsePauliOp) -> List[Group]:
    """Group Pauli terms by measurement basis and precompute basis rotations."""
    num_qubits = pauli_op.num_qubits
    by_axis: Dict[str, List[Term]] = {}
    for label, coeff in pauli_op.to_list():
        # axis per qubit index (little-endian label -> qubit index)
        axis_chars = ["I"] * num_qubits
        support = []
        for pos, c in enumerate(label):
            if c != "I":
                q = _qubit_from_label_pos(num_qubits, pos)
                axis_chars[q] = c if c in "XYZ" else "I"
                support.append(q)
        axis = "".join(axis_chars)
        by_axis.setdefault(axis, []).append(
            Term(coeff=coeff, label=label, support=support)
        )

    groups = []
    for axis, terms in by_axis.items():
        rot = QuantumCircuit(num_qubits)
        for q, ch in enumerate(axis):
            if ch == "X":
                rot.h(q)
            elif ch == "Y":
                rot.sdg(q)
                rot.h(q)
        groups.append(Group(axis_pattern=axis, terms=terms, rotation=rot))
    return groups


def _z_signs(num_qubits: int, support: Sequence[int]) -> np.ndarray:
    """Array s[b] = product over q in support of (-1)^{bit_q(b)} for basis state b."""
    n = 2 ** num_qubits
    signs = np.ones(n, dtype=float)
    for q in support:
        bit = (np.arange(n) >> q) & 1
        signs *= 1 - 2 * bit
    return signs


def expectation_from_density_matrix(rho: np.ndarray, group: Group) -> float:
    """Energy contribution of a group from the rotated density matrix.

    After the basis rotation, each term is ``coeff * Z^support`` so its
    expectation is a signed sum over the diagonal of ``rho``.
    """
    diag = np.real(np.diag(rho))
    num_qubits = int(round(np.log2(rho.shape[0])))
    total = 0.0
    for term in group.terms:
        if not term.support:
            total += float(np.real(term.coeff))
            continue
        signs = _z_signs(num_qubits, term.support)
        expect_z = float(np.dot(diag, signs))
        total += float(np.real(term.coeff * expect_z))
    return total


def expectation_from_statevector(psi: np.ndarray, term: Term) -> float:
    """⟨ψ|P|ψ⟩ for one term via direct Pauli action on the statevector."""
    if not term.support:
        return float(np.real(term.coeff))
    num_qubits = int(round(np.log2(psi.shape[0])))
    xbits = []
    zbits = []
    n_y = 0
    for i, c in enumerate(term.label):
        q = _qubit_from_label_pos(num_qubits, i)
        if c == "X":
            xbits.append(q)
        elif c == "Z":
            zbits.append(q)
        elif c == "Y":  # Y = i * X * Z
            n_y += 1
            xbits.append(q)
            zbits.append(q)

    phases = np.ones(2 ** num_qubits, dtype=complex)
    for q in zbits:
        bit = (np.arange(2 ** num_qubits) >> q) & 1
        phases *= 1 - 2 * bit
    psi_z = psi * phases

    if xbits:
        idx = np.arange(2 ** num_qubits)
        for q in xbits:
            idx ^= 1 << q
        psi_x = psi_z[idx]
    else:
        psi_x = psi_z
    overlap = np.vdot(psi, psi_x)
    if n_y:
        overlap *= (1j) ** n_y
    return float(np.real(term.coeff * overlap))


def expectation_from_counts(
    counts: Dict[str, int], group: Group, num_qubits: int, total: Optional[int] = None
) -> float:
    """Energy contribution of a group from measured counts."""
    if total is None:
        total = sum(counts.values())
    if total == 0:
        return 0.0
    total = float(total)
    energy = 0.0
    for term in group.terms:
        if not term.support:
            energy += float(np.real(term.coeff))
            continue
        # sign of term on basis state b: (-1)^{popcount(b & support_mask)}
        mask = 0
        for q in term.support:
            mask |= 1 << q
        expect = 0.0
        for bitstring, cnt in counts.items():
            # Aer count strings are the binary representation of the state
            # index (leftmost char = highest qubit), so int(bs, 2) is the index
            b = int(bitstring, 2)
            parity = bin(b & mask).count("1")
            expect += cnt * (1.0 if parity % 2 == 0 else -1.0)
        expect /= total
        energy += float(np.real(term.coeff * expect))
    return energy


def _param_key(p) -> tuple:
    """Sort Parameter objects by (name, numeric index): theta[2] < theta[10]."""
    name = p.name
    if "[" in name and name.endswith("]"):
        head, idx = name.rsplit("[", 1)
        return (head, int(idx[:-1]))
    return (name, 0)


def _flip_bitstring(bitstring: str, qubits: Sequence[int]) -> str:
    # char position p (from the left) corresponds to qubit n-1-p
    bits = list(bitstring)
    n = len(bitstring)
    for q in qubits:
        p = n - 1 - q
        bits[p] = "1" if bits[p] == "0" else "0"
    return "".join(bits)


class ExpectationEvaluator:
    """Computes expectation values of a Hamiltonian for parameterized circuits."""

    def __init__(
        self,
        hamiltonian: MolecularHamiltonian,
        mode: str = "exact",
        noise_model=None,
        shots: Optional[int] = None,
        readout_filter: Optional[ReadoutFilter] = None,
        twirls: int = 1,
        seed: Optional[int] = None,
        method: str = "automatic",
    ):
        if mode not in ("exact", "exact_noisy", "sampled"):
            raise ValueError(f"Unknown mode {mode!r}")
        if mode == "exact" and noise_model is not None:
            raise ValueError("exact mode cannot be used with a noise model")
        if mode == "exact" and method == "density_matrix":
            raise ValueError("exact mode requires statevector simulation")
        if mode != "sampled" and shots is not None:
            raise ValueError("shots only apply to sampled mode")

        self.hamiltonian = hamiltonian
        self.mode = mode
        self.noise_model = noise_model
        self.shots = shots
        self.readout_filter = readout_filter
        self.twirls = twirls
        self.seed = seed
        self._rng = np.random.default_rng(seed)

        self.groups = build_groups(hamiltonian.pauli_op)
        self.num_qubits = hamiltonian.num_qubits
        self._offset = hamiltonian.nuclear_repulsion_energy + hamiltonian.energy_offset
        # circuits are transpiled to the noise model's basis gates so that
        # every emitted gate carries its calibrated error rate
        self._basis = None
        if noise_model is not None:
            self._basis = list(noise_model.basis_gates)

        sim_method = method
        if mode == "exact":
            sim_method = "statevector"
        elif mode == "exact_noisy":
            sim_method = "density_matrix"
        self._sim = AerSimulator(method=sim_method, noise_model=noise_model)
        self._qasm = AerSimulator(method="automatic", noise_model=noise_model)

    # ------------------------------------------------------------------ #
    def energy(self, circuit: QuantumCircuit, params=None) -> float:
        """Total energy (incl. NRE) of ``circuit`` with parameters bound."""
        bound = self._bind(circuit, params)
        if self.mode == "exact":
            return self._energy_exact(bound)
        if self.mode == "exact_noisy":
            return self._energy_exact_noisy(bound)
        return self._energy_sampled(bound)

    # ------------------------------------------------------------------ #
    def _bind(self, circuit: QuantumCircuit, params):
        if params is None:
            return circuit
        if isinstance(params, dict):
            return circuit.assign_parameters(params)
        pv = sorted(circuit.parameters, key=_param_key)
        if len(pv) != len(params):
            raise ValueError(
                f"expected {len(pv)} parameters, got {len(params)}"
            )
        return circuit.assign_parameters(dict(zip(pv, params)))

    def _energy_exact(self, circuit: QuantumCircuit) -> float:
        qc = circuit.copy()
        qc.save_statevector(label="psi")
        res = self._sim.run(qc, shots=1).result()
        psi = res.data(0)["psi"]
        if hasattr(psi, "data"):
            psi = psi.data
        energy = self._offset
        for group in self.groups:
            for term in group.terms:
                energy += expectation_from_statevector(np.asarray(psi), term)
        return float(energy)

    def _prepare(self, circuit: QuantumCircuit, suffix=None) -> QuantumCircuit:
        """Copy + compose rotation/measurement suffix + transpile to noise basis."""
        qc = circuit.copy()
        if suffix is not None:
            qc.compose(suffix, inplace=True)
        if self._basis is not None:
            qc = transpile(qc, basis_gates=self._basis, optimization_level=1)
        return qc

    def _energy_exact_noisy(self, circuit: QuantumCircuit) -> float:
        """Density-matrix simulation: one run, direct Pauli traces.

        A single noisy density-matrix simulation of the whole circuit yields
        ``rho``; each Pauli term is then evaluated analytically as
        ``Tr(P rho)`` using the sparsity pattern of P (a Pauli string has one
        nonzero entry per row/column along a single off-diagonal line). This
        avoids one simulation *per measurement group*, which is essential for
        larger molecules (LiH has ~276 terms vs H2's 6 groups).
        """
        qc = self._prepare(circuit, None)
        qc.save_density_matrix(label="rho")
        res = self._sim.run(qc, shots=1).result()
        rho = res.data(0)["rho"]
        if hasattr(rho, "data"):
            rho = np.asarray(rho.data)
        else:
            rho = np.asarray(rho)
        n = self.num_qubits
        idx = np.arange(2 ** n)
        energy = 0.0
        for group in self.groups:
            for term in group.terms:
                if not term.support:
                    energy += float(np.real(term.coeff))
                    continue
                xmask = 0
                zmask = 0
                n_y = 0
                for pos, c in enumerate(term.label):
                    q = _qubit_from_label_pos(n, pos)
                    if c == "X":
                        xmask |= 1 << q
                    elif c == "Z":
                        zmask |= 1 << q
                    elif c == "Y":  # Y = i*X*Z
                        xmask |= 1 << q
                        zmask |= 1 << q
                        n_y += 1
                # P|b> = i^{n_y} (-1)^{b.z} |b (+) x>  (Y = i X Z per qubit), so
                #   Tr(P rho) = i^{n_y} sum_b (-1)^{b.z} rho[b, b (+) x]
                flipped = idx ^ xmask
                sign = np.ones(2 ** n)
                z = zmask
                while z:
                    q = (z & -z).bit_length() - 1
                    bit = (idx >> q) & 1
                    sign *= 1 - 2 * bit
                    z &= z - 1
                trace = complex(np.sum(sign * rho[idx, flipped]))
                phase = (1j) ** n_y
                energy += float(np.real(term.coeff * phase * trace))
        return float(energy + self._offset)

    def _energy_sampled(self, circuit: QuantumCircuit) -> float:
        energy = self._offset
        shots = self.shots or 1024
        for group in self.groups:
            group_energy = 0.0
            for _ in range(self.twirls):
                qc = self._prepare(circuit, group.rotation)
                twirled = []
                if self.twirls > 1:
                    for q in range(self.num_qubits):
                        if self._rng.random() < 0.5:
                            qc.x(q)
                            twirled.append(q)
                qc.measure_all()
                counts = self._qasm.run(qc, shots=shots, seed_simulator=self.seed).result().get_counts(0)
                if twirled:
                    counts = {
                        _flip_bitstring(bs, twirled): c for bs, c in counts.items()
                    }
                if self.readout_filter is not None:
                    counts = self.readout_filter(counts, self.num_qubits)
                group_energy += expectation_from_counts(
                    counts, group, self.num_qubits
                )
            energy += group_energy / self.twirls
        return float(energy)