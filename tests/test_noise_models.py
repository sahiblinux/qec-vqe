"""M2 tests: asymmetric noise profiles and calibration helpers."""

import numpy as np
import pytest

from qec_vqe.noise_models import (
    DEFAULT_BASIS_GATES,
    biased_depolarizing,
    build_noise_model,
    default_profiles,
    load_profile,
    make_profile,
    profile_from_t1t2,
    t1t2_derived_rates,
    uniform_depolarizing_2q,
)


def _pauldi_probs(err):
    """Map Pauli name -> probability from an Aer pauli_error."""
    d = err.to_dict()
    out = {}
    for inst, prob in zip(d["instructions"], d["probabilities"]):
        if inst[0]["name"] == "pauli":
            pauli = inst[0]["params"][0]
        else:
            pauli = "".join(g["name"].upper() for g in inst)
        out[pauli] = float(prob)
    return out


def test_biased_depolarizing_normalization():
    err = biased_depolarizing(0.02, eta=5.0)
    probs = _pauldi_probs(err)
    assert abs(sum(probs.values()) - 1.0) < 1e-9
    # eta = p_z / p_x ; p_y == p_x
    px, py, pz = probs["X"], probs["Y"], probs["Z"]
    assert abs(py - px) < 1e-12
    assert abs(pz / px - 5.0) < 1e-9
    assert abs(px + py + pz - 0.02) < 1e-9


def test_uniform_depolarizing_2q():
    err = uniform_depolarizing_2q(0.01)
    probs = _pauldi_probs(err)
    assert abs(probs["II"] - 0.99) < 1e-9
    non_id = {k: v for k, v in probs.items() if k != "II"}
    assert len(non_id) == 15
    for v in non_id.values():
        assert abs(v - 0.01 / 15) < 1e-12
    assert abs(sum(probs.values()) - 1.0) < 1e-9


def test_t1t2_derived_bias():
    rates = t1t2_derived_rates(t1_us=300.0, t2_us=120.0, gate_time_ns=50.0)
    # T2 < 2*T1 dephasing => phase-flip (Z) dominates bit-flip (X)
    assert rates["p_z"] > rates["p_x"]
    assert abs(rates["p_x"] - rates["p_y"]) < 1e-12
    assert rates["p_z"] >= 0


def test_profile_from_t1t2():
    prof = profile_from_t1t2(300.0, 120.0)
    assert prof["bias_eta"] > 1.0
    assert 0 < prof["single_qubit_error_rate"] < 1
    assert prof["two_qubit_error_rate"] > 0


def test_make_and_build_noise_model():
    prof = make_profile("t", single_qubit_error_rate=3e-4, two_qubit_error_rate=8e-3,
                        bias_eta=5.0, readout_error=0.02)
    nm = build_noise_model(prof)
    assert set(nm.basis_gates) <= set(DEFAULT_BASIS_GATES + ["cx"])
    assert not nm.is_ideal()
    assert "cx" in nm.noise_instructions  # 2q errors attached
    assert "measure" in nm.noise_instructions  # readout errors attached
    assert "id" in nm.noise_instructions


def test_shipped_profiles_load():
    profs = default_profiles()
    assert "ibm_brisbane_style" in profs
    assert "ibm_brisbane_biased" in profs
    p = load_profile("ibm_brisbane_style")
    assert p["bias_eta"] >= 1.0
    # load by path also works
    from qec_vqe.noise_models import PROFILES_DIR
    p2 = load_profile(str(PROFILES_DIR / "ibm_brisbane_biased.json"))
    assert p2["name"] == "ibm_brisbane_biased"
