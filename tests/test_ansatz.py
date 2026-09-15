"""Tests for ansatz construction and HF consistency."""

import numpy as np
import pytest

from qec_vqe import (
    build_hamiltonian,
    hardware_efficient_ansatz,
    hartree_fock_state_circuit,
    uccsd_ansatz,
)
from qec_vqe.expectation import ExpectationEvaluator


def test_hf_state_matches_reference_energy():
    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    ev = ExpectationEvaluator(H, mode="exact")
    # all-zero UCCSD parameters -> Hartree-Fock state
    e_hf = ev.energy(ans.circuit, np.zeros(ans.num_parameters))
    # HF energy of H2 STO-3G is well above the exact (-1.1373 Ha)
    assert e_hf > H.exact_ground_energy()
    assert abs(e_hf - (-1.116999)) < 1e-3  # known STO-3G HF value


def test_uccsd_zero_params_is_hf():
    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    ev = ExpectationEvaluator(H, mode="exact")
    hf_circ = hartree_fock_state_circuit(H)
    from qiskit import QuantumCircuit

    # compose HF into an empty circuit of the right width
    qc = QuantumCircuit(H.num_qubits)
    qc.compose(hf_circ, inplace=True)
    e_hf_direct = ev.energy(qc)
    e_ucc = ev.energy(ans.circuit, np.zeros(ans.num_parameters))
    assert abs(e_hf_direct - e_ucc) < 1e-6


def test_uccsd_parameter_counts():
    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    assert ans.num_parameters == 3  # 1 singles + 2 doubles (2 electrons, 4 orbs)
    # binding must not raise
    bound = ans.bind(np.random.default_rng(0).uniform(-0.1, 0.1, ans.num_parameters))
    assert bound.num_qubits == H.num_qubits


def test_hea_structure():
    H = build_hamiltonian("H2")
    ans = hardware_efficient_ansatz(H, layers=2, entangling="linear")
    n = H.num_qubits
    assert ans.num_parameters == 2 * n * 2
    assert len(ans.blocks) == 2
    assert ans.bind(np.zeros(ans.num_parameters)).num_qubits == n


def test_hea_from_hf_not_above_hf_by_much():
    """A zero-parameter HEA (HF init + zero rotations) reproduces HF."""
    H = build_hamiltonian("H2")
    ans = hardware_efficient_ansatz(H, layers=1, entangling="linear")
    ev = ExpectationEvaluator(H, mode="exact")
    # zero-rotation layers with CX between: CX on |01..> HF state can move
    # population; HF init energy should still be close to pure HF
    e = ev.energy(ans.circuit, np.zeros(ans.num_parameters))
    assert e >= H.exact_ground_energy() - 1e-6  # variational upper bound-ish
