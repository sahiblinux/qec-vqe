"""M1 tests: Hamiltonian construction and exact-diagonalization validation."""

import numpy as np
import pytest

from qec_vqe import (
    REFERENCE_ENERGIES,
    MolecularHamiltonian,
    build_hamiltonian,
    pauli_groups,
)


def test_reference_energies_known():
    assert set(REFERENCE_ENERGIES) == {"H2", "LiH"}
    assert abs(REFERENCE_ENERGIES["H2"] - (-1.137306)) < 1e-4
    assert abs(REFERENCE_ENERGIES["LiH"] - (-7.882117)) < 1e-3


def test_h2_structure():
    H = build_hamiltonian("H2")
    assert isinstance(H, MolecularHamiltonian)
    assert H.num_qubits == 4  # Jordan-Wigner, 2 spatial orbitals
    assert H.system == "H2"


def test_h2_matches_reference():
    H = build_hamiltonian("H2")
    # cross-validated against exact classical diagonalization
    assert H.exact_ground_energy() == pytest.approx(H.exact_eigenvalues()[0]
        + H.nuclear_repulsion_energy + H.energy_offset, abs=1e-10)
    assert abs(H.exact_ground_energy() - REFERENCE_ENERGIES["H2"]) < 1e-4
    # reference ground energy reported in the module docstring
    assert abs(H.exact_ground_energy() - (-1.137306)) < 1e-4


def test_lih_structure_and_freeze_core():
    H = build_hamiltonian("LiH")
    # 12 qubits unfrozen? default freeze_core=True -> 10 active qubits (JW)
    # LiH STO-3G: 6 spatial orbitals -> 12 spin orbitals; freeze Li 1s -> 10
    assert H.num_qubits == 10
    assert abs(H.exact_ground_energy() - REFERENCE_ENERGIES["LiH"]) < 1e-3


def test_energy_offset_carries_frozen_core():
    """Frozen-core LiH energy must equal reference (offset constant included)."""
    H = build_hamiltonian("LiH")
    assert abs(H.exact_ground_energy() - REFERENCE_ENERGIES["LiH"]) < 1e-3
    assert H.energy_offset < 0  # frozen-core constant is a negative offset


def test_parity_mapping_reduces_qubits():
    H = build_hamiltonian("H2", mapping="parity")
    assert H.num_qubits == 2  # parity reduces H2 from 4 -> 2 qubits
    assert abs(H.exact_ground_energy() - REFERENCE_ENERGIES["H2"]) < 1e-4


def test_geometry_sensitivity():
    """Stretched H2 must sit above the equilibrium energy (dissociation)."""
    H = build_hamiltonian("H2", geometry="H 0 0 0; H 0 0 2.5")
    H_eq = build_hamiltonian("H2")
    assert H.exact_ground_energy() > H_eq.exact_ground_energy()


def test_pauli_groups_partition_terms():
    H = build_hamiltonian("H2")
    groups = pauli_groups(H.pauli_op)
    n_terms = sum(g.size for g in groups)
    assert n_terms == H.num_terms
    # every group must share one per-qubit axis pattern
    for g in groups:
        labels = g.paulis.to_labels()
        axes = {
            "".join(ch if ch in "XYZ" else "I" for ch in label)
            for label in labels
        }
        assert len(axes) == 1
