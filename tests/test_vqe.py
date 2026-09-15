"""M1/M3 tests: VQE convergence to chemical accuracy under each protocol."""

import numpy as np
import pytest

from qec_vqe import build_hamiltonian, uccsd_ansatz
from qec_vqe.expectation import ExpectationEvaluator
from qec_vqe.noise_models import build_noise_model, make_profile
from qec_vqe.vqe import SPSAOptions, minimize_layers, run_vqe, scipy_minimize_vqe, spsa_minimize

CHEM_ACC_MHA = 1.6  # chemical accuracy in millihartree


def test_spsa_perturbation_statistics():
    opts = SPSAOptions(seed=3)
    draws = np.array([opts.perturbation(50) for _ in range(200)])
    assert set(np.unique(draws)) <= {-1.0, 1.0}
    assert abs(draws.mean()) < 0.3  # roughly balanced signs


def test_spsa_decreases_on_convex_quadratic():
    def obj(x):
        return float(np.sum((x - np.array([1.0, -2.0, 3.0])) ** 2))

    res = spsa_minimize(obj, 3, SPSAOptions(maxiter=150, a=0.8, c=0.3, seed=1))
    assert res.energy < 2.0  # started near 14, should end well below


def test_h2_noiseless_reaches_chemical_accuracy():
    H = build_hamiltonian("H2")
    exact = H.exact_ground_energy()
    ans = uccsd_ansatz(H)
    ev = ExpectationEvaluator(H, mode="exact")
    res = scipy_minimize_vqe(
        ans, lambda th: ev.energy(ans.circuit, th), method="COBYLA", maxiter=400
    )
    assert (res.energy - exact) * 1e3 < CHEM_ACC_MHA
    assert res.energy <= H.exact_ground_energy() + 1e-3  # near-optimal


def test_noiseless_h2_spsa_converges():
    H = build_hamiltonian("H2")
    exact = H.exact_ground_energy()
    ans = uccsd_ansatz(H)
    ev = ExpectationEvaluator(H, mode="exact")
    res = run_vqe(ans, lambda th: ev.energy(ans.circuit, th), maxiter=150, seed=7)
    assert (res.energy - exact) * 1e3 < CHEM_ACC_MHA
    assert len(res.history) > 10


def test_noisy_vqe_degrades_below_chemical_accuracy():
    """Raw noisy optimization must fail chemical accuracy (motivates M3)."""
    H = build_hamiltonian("H2")
    exact = H.exact_ground_energy()
    ans = uccsd_ansatz(H)
    prof = make_profile("n", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    nm = build_noise_model(prof)
    ev = ExpectationEvaluator(H, mode="exact_noisy", noise_model=nm)
    res = run_vqe(ans, lambda th: ev.energy(ans.circuit, th), maxiter=40, seed=7)
    assert (res.energy - exact) * 1e3 > CHEM_ACC_MHA * 3


def test_minimize_layers_runs():
    H = build_hamiltonian("H2")
    ans = uccsd_ansatz(H)
    ev = ExpectationEvaluator(H, mode="exact")

    def obj(th):
        return ev.energy(ans.circuit, th)

    res = minimize_layers(obj, [list(range(ans.num_parameters))],
                          SPSAOptions(maxiter=30, seed=1), rounds=1)
    assert res.energy < -1.0
