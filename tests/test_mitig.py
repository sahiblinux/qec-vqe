"""M3 tests: ZNE extrapolation and TREX readout mitigation."""

import numpy as np
import pytest

from qec_vqe import build_hamiltonian, uccsd_ansatz
from qec_vqe.expectation import ExpectationEvaluator
from qec_vqe.mitig.readout import ReadoutMitigator
from qec_vqe.mitig.zne import (
    fold_circuit,
    run_zne,
    scale_profile,
    zne_extrapolate,
)
from qec_vqe.noise_models import build_noise_model, make_profile


def test_fold_circuit_preserves_qubits_and_clbits():
    from qiskit import QuantumCircuit

    qc = QuantumCircuit(2, 2)
    qc.h(0)
    qc.cx(0, 1)
    qc.measure(0, 0)
    qc.measure(1, 1)
    f3 = fold_circuit(qc, 3)
    assert f3.num_qubits == 2
    assert f3.num_clbits == 2
    assert f3.depth() >= qc.depth()
    # folded circuit still executable
    from qiskit_aer import AerSimulator

    out = AerSimulator().run(f3, shots=10).result().get_counts(0)
    assert sum(out.values()) == 10


def test_zne_extrapolation_linear_exact():
    # E(lambda) = 5 + 2 lambda -> zero-noise limit = 5
    assert zne_extrapolate([1.0, 3.0], [7.0, 11.0], "linear") == pytest.approx(5.0)


def test_zne_extrapolation_recovers_intercept():
    x = np.array([1.0, 2.0, 3.0])
    y = 1.7 + 0.3 * x - 0.02 * x ** 2
    assert zne_extrapolate(x, y, "richardson") == pytest.approx(1.7, abs=1e-8)
    # exponential on pure linear-ish data also lands near intercept
    assert abs(zne_extrapolate(x, y, "exponential") - 1.7) < 0.05


def test_scale_profile_multiplies_rates():
    prof = make_profile("p", single_qubit_error_rate=1e-4, two_qubit_error_rate=1e-3,
                        readout_error=0.01)
    s3 = scale_profile(prof, 3.0)
    assert s3["single_qubit_error_rate"] == pytest.approx(3e-4)
    assert s3["two_qubit_error_rate"] == pytest.approx(3e-3)
    assert s3["readout_error"] == pytest.approx(0.03)


def test_zne_reduces_noise_bias_on_h2():
    """ZNE (noise-rate scaling) must pull the noisy energy toward exact."""
    H = build_hamiltonian("H2")
    exact = H.exact_ground_energy()
    ans = uccsd_ansatz(H)
    ev = ExpectationEvaluator(H, mode="exact")
    # converge noiselessly (fast on H2)
    from qec_vqe.vqe import scipy_minimize_vqe

    res = scipy_minimize_vqe(ans, lambda th: ev.energy(ans.circuit, th),
                             method="COBYLA", maxiter=400)
    theta = res.params
    prof = make_profile("n", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    z = run_zne(H, ans.circuit, theta, profile=prof,
                scales=(1, 2, 3, 4, 5), method="richardson", mode="exact_noisy")
    raw_err = abs(z["energies"][0] - exact) * 1e3
    zne_err = abs(z["extrapolated"] - exact) * 1e3
    assert zne_err < raw_err  # mitigation reduces the bias
    # and recovers within a small fraction of the raw bias
    assert zne_err < max(10.0, 0.1 * raw_err)


def test_readout_mitigator_tensored_calibration():
    prof = make_profile("n", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    nm = build_noise_model(prof)
    rm = ReadoutMitigator.calibrate(nm, 4, shots=8192, seed=7)
    diag = np.diag(rm.cal_matrix)
    # per-qubit assignment fidelity ~0.98, so diagonal ~ 0.92
    assert np.all(diag > 0.85)
    # columns must be valid probability vectors
    assert np.allclose(rm.cal_matrix.sum(axis=0), 1.0, atol=1e-6)


def test_readout_mitigator_full_calibration():
    prof = make_profile("n", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    nm = build_noise_model(prof)
    rm = ReadoutMitigator.calibrate(nm, 3, shots=8192, full=True, seed=7)
    assert rm.cal_matrix.shape == (8, 8)
    assert np.all(np.diag(rm.cal_matrix) > 0.8)


def test_trex_corrects_ground_state_counts():
    """On prepared |0000>, TREX must push P(0000) toward 1."""
    from qiskit import QuantumCircuit
    from qiskit_aer import AerSimulator

    prof = make_profile("n", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    nm = build_noise_model(prof)
    rm = ReadoutMitigator.calibrate(nm, 4, shots=8192, seed=7)
    sim = AerSimulator(method="automatic", noise_model=nm)
    qc = QuantumCircuit(4)
    qc.measure_all()
    counts = sim.run(qc, shots=20000, seed_simulator=7).result().get_counts(0)
    raw_p = counts.get("0000", 0) / sum(counts.values())
    corrected = rm.apply(counts, 4)
    corr_p = corrected.get("0000", 0.0)
    assert raw_p < 0.99  # readout error present
    assert corr_p > raw_p  # mitigation recovers
    assert corr_p > 0.99
