"""Tests for the expectation-value engine (exact / exact_noisy / sampled)."""

import numpy as np
import pytest

from qec_vqe import build_hamiltonian, uccsd_ansatz
from qec_vqe.expectation import ExpectationEvaluator
from qec_vqe.noise_models import build_noise_model, make_profile


def test_exact_mode_basic():
    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    ev = ExpectationEvaluator(H, mode="exact")
    e = ev.energy(ans.circuit, np.zeros(ans.num_parameters))
    # HF energy -1.117 Ha; exact diag -1.1373 Ha
    assert abs(e - (-1.116999)) < 1e-4


def test_exact_noisy_agrees_with_exact_without_noise():
    """With a zero-error profile, density-matrix mode reproduces statevector."""
    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    theta = np.zeros(ans.num_parameters)
    ev_exact = ExpectationEvaluator(H, mode="exact")
    e_exact = ev_exact.energy(ans.circuit, theta)

    prof = make_profile("zero", single_qubit_error_rate=0.0, two_qubit_error_rate=0.0,
                        bias_eta=1.0, readout_error=0.0)
    nm = build_noise_model(prof)
    ev_noisy = ExpectationEvaluator(H, mode="exact_noisy", noise_model=nm)
    e_noisy = ev_noisy.energy(ans.circuit, theta)
    assert abs(e_noisy - e_exact) < 1e-6


def test_noise_degrades_energy():
    """Under a noisy profile the noisy energy must move off the exact value."""
    H = build_hamiltonian("H2")
    exact = H.exact_ground_energy()
    ans = uccsd_ansatz(H)
    prof = make_profile("n", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    nm = build_noise_model(prof)
    ev = ExpectationEvaluator(H, mode="exact_noisy", noise_model=nm)
    e = ev.energy(ans.circuit, np.zeros(ans.num_parameters))
    assert abs(e - exact) > 1e-3  # measurably degraded


def test_sampled_matches_exact_within_shot_noise():
    """Shot-based sampling without noise converges to the exact energy."""
    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    theta = np.zeros(ans.num_parameters)
    ev_exact = ExpectationEvaluator(H, mode="exact")
    e_exact = ev_exact.energy(ans.circuit, theta)
    ev_sampled = ExpectationEvaluator(H, mode="sampled", shots=60000, seed=42)
    e_sampled = ev_sampled.energy(ans.circuit, theta)
    # several mHa of shot noise on a 4-qubit molecule at 60k shots/group
    assert abs(e_sampled - e_exact) < 0.01


def test_single_sim_density_matches_per_group_formula():
    """Direct Pauli-trace path == rotated-diagonal reference on a tiny system."""
    from qec_vqe.expectation import expectation_from_density_matrix

    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    theta = np.zeros(ans.num_parameters)
    prof = make_profile("n", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    nm = build_noise_model(prof)
    ev = ExpectationEvaluator(H, mode="exact_noisy", noise_model=nm)
    e_new = ev.energy(ans.circuit, theta)

    # reference: per-group rotation on the same density matrix
    bound = ev._bind(ans.circuit, theta)
    total = H.nuclear_repulsion_energy + H.energy_offset
    for group in ev.groups:
        qc = ev._prepare(bound, group.rotation)
        qc.save_density_matrix(label="rho")
        rho = ev._sim.run(qc, shots=1).result().data(0)["rho"]
        total += expectation_from_density_matrix(np.asarray(rho.data), group)
    assert abs(e_new - total) < 1e-6
