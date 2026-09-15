"""Molecular Hamiltonian construction for VQE targets.

Builds qubit Hamiltonians for small molecules (H2, LiH) in the STO-3G basis
using qiskit-nature's PySCF driver. The default fermion-to-qubit mapping is
Jordan-Wigner (qubit j = spin-orbital j), which keeps qiskit-nature's
Hartree-Fock / UCCSD circuit library consistent with our Hamiltonians with no
convention ambiguity. Optional parity mapping gives a 2-qubit reduction but is
not used by default because qiskit 2.x dropped the Z2-tapering helpers needed
to align ansatz circuits with parity-tapered Hamiltonians.

* H2  (bond 0.735 A, STO-3G)  ->  4 qubits,  15 Pauli terms (JW)
* LiH (bond 1.494 A, STO-3G)  -> 12 qubits, 631 Pauli terms (JW)
                                   (10 qubits / ~275 terms with Li 1s frozen)

All energies reported by this package include the nuclear repulsion energy,
i.e. they are total molecular energies in Hartree.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional

import numpy as np
from qiskit.quantum_info import SparsePauliOp

from qiskit_nature.second_q.drivers import PySCFDriver
from qiskit_nature.second_q.mappers import JordanWignerMapper, ParityMapper
from qiskit_nature.second_q.transformers import FreezeCoreTransformer

#: Reference ground-state energies (STO-3G, total energy incl. NRE, Hartree)
REFERENCE_ENERGIES: Dict[str, float] = {
    "H2": -1.137306,  # equilibrium bond 0.735 A
    "LiH": -7.882117,  # equilibrium bond 1.494 A
}

DEFAULT_GEOMETRIES: Dict[str, str] = {
    "H2": "H 0 0 0; H 0 0 0.735",
    "LiH": "Li 0 0 0; H 0 0 1.494",
}


@dataclass
class MolecularHamiltonian:
    """A qubit Hamiltonian plus the metadata needed to interpret it."""

    system: str
    pauli_op: SparsePauliOp
    num_particles: tuple
    num_spatial_orbitals: int
    nuclear_repulsion_energy: float
    #: constant energy offset not represented by Pauli terms (e.g. frozen-core
    #: energy kept by qiskit-nature's FreezeCoreTransformer in ``constants``)
    energy_offset: float = 0.0
    #: informational description of the fermion-to-qubit mapping used
    mapper_info: str = ""

    @property
    def num_qubits(self) -> int:
        return self.pauli_op.num_qubits

    @property
    def num_terms(self) -> int:
        return self.pauli_op.size

    def exact_eigenvalues(self) -> np.ndarray:
        """Eigenvalues of the *qubit* Hamiltonian (Ha, electronic part only)."""
        mat = self.pauli_op.to_matrix(sparse=True).toarray()
        return np.linalg.eigvalsh(mat)

    def exact_ground_energy(self, include_nre: bool = True) -> float:
        """Lowest eigenvalue of the qubit Hamiltonian, plus NRE + offsets."""
        e = float(self.exact_eigenvalues()[0])
        if include_nre:
            e += self.nuclear_repulsion_energy + self.energy_offset
        return e

    def reference_energy(self) -> Optional[float]:
        """Published STO-3G reference (total energy, Ha) if known."""
        return REFERENCE_ENERGIES.get(self.system)

    def sorted_terms(self) -> SparsePauliOp:
        """Pauli terms sorted by coefficient magnitude (descending)."""
        idx = np.argsort(-np.abs(self.pauli_op.coeffs))
        return self.pauli_op[idx]


def build_hamiltonian(
    system: str = "H2",
    geometry: Optional[str] = None,
    basis: str = "sto3g",
    freeze_core: Optional[bool] = None,
    mapping: str = "jw",
) -> MolecularHamiltonian:
    """Build the qubit Hamiltonian for a molecular VQE target.

    Parameters
    ----------
    system:
        "H2" or "LiH".
    geometry:
        PySCF-style atom string; defaults to the equilibrium geometry.
    basis:
        Gaussian basis set (default STO-3G).
    freeze_core:
        Freeze core orbitals (Li 1s for LiH). Default: True for LiH,
        False otherwise (H has no frozen core in STO-3G).
    mapping:
        "jw" (Jordan-Wigner, default) or "parity".
    """
    system = {"h2": "H2", "lih": "LiH"}.get(system.strip().lower())
    if system is None:
        raise ValueError(f"Unknown system {system!r}; expected H2 or LiH")
    geometry = geometry or DEFAULT_GEOMETRIES[system]
    if freeze_core is None:
        freeze_core = system == "LiH"

    driver = PySCFDriver(atom=geometry, basis=basis, charge=0, spin=0)
    problem = driver.run()

    reduced = problem
    energy_offset = 0.0
    if freeze_core:
        reduced = FreezeCoreTransformer(freeze_core=True).transform(reduced)
        energy_offset = float(
            reduced.hamiltonian.constants.get("FreezeCoreTransformer", 0.0)
        )

    if mapping == "jw":
        mapper = JordanWignerMapper()
    elif mapping == "parity":
        mapper = ParityMapper(num_particles=reduced.num_particles)
    else:
        raise ValueError(f"unknown mapping {mapping!r}; expected 'jw' or 'parity'")
    pauli_op = mapper.map(reduced.second_q_ops()[0]).simplify()

    return MolecularHamiltonian(
        system=system,
        pauli_op=pauli_op,
        num_particles=tuple(reduced.num_particles),
        num_spatial_orbitals=int(reduced.num_spatial_orbitals),
        nuclear_repulsion_energy=float(problem.nuclear_repulsion_energy),
        energy_offset=energy_offset,
        mapper_info=f"{mapper.__class__.__name__}, freeze_core={freeze_core}",
    )


def pauli_groups(pauli_op: SparsePauliOp) -> List[SparsePauliOp]:
    """Split a SparsePauliOp into groups sharing a measurement basis.

    Each returned operator has the same per-qubit Pauli *axis* pattern
    (I/X/Y/Z), so its terms can be estimated from a single measurement basis
    after the appropriate single-qubit rotations (see ``expectation.py``).
    """
    labels = pauli_op.paulis.to_labels()
    groups: Dict[str, List[int]] = {}
    for i, label in enumerate(labels):
        axis = "".join(ch if ch in "XYZ" else "I" for ch in label)
        groups.setdefault(axis, []).append(i)

    result = []
    for _, idxs in groups.items():
        result.append(pauli_op[idxs])
    return result