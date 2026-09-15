"""Asymmetric Pauli noise models for NISQ simulation.

Builds Qiskit Aer ``NoiseModel`` objects from calibration-style JSON profiles.
The key hardware-aware feature is *noise bias*: real devices (IBM Eagle/Heron,
Google Sycamore) exhibit phase-flip error rates several times larger than
bit-flip rates (T2 < 2*T1), i.e. p_z >> p_x. We model this with asymmetric
Pauli channels ``(p_x, p_y, p_z)`` with bias ``eta = p_z / p_x``.

A profile looks like::

    {
      "name": "ibm_brisbane_style",
      "basis_gates": ["id", "rz", "x", "sx", "cx"],
      "single_qubit_error_rate": 3e-4,   # total per-gate error probability
      "two_qubit_error_rate": 8e-3,
      "bias_eta": 5.0,                    # p_z / p_x on single-qubit gates
      "readout_error": 0.02,
      "readout_asymmetry": 1.2            # P(0|1) / P(1|0)
    }

Circuits must be transpiled to the profile's basis gates for every gate to
carry its calibrated error rate (the evaluator does this automatically).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from qiskit_aer.noise import NoiseModel, ReadoutError, pauli_error

PROFILES_DIR = Path(__file__).resolve().parent.parent / "profiles"

#: gate names assumed by the canonical transpile basis of this package
DEFAULT_BASIS_GATES = ["id", "rz", "x", "sx", "cx"]


def biased_depolarizing(p_total: float, eta: float) -> "pauli_error":
    """Asymmetric single-qubit Pauli channel with bias ``eta = p_z / p_x``.

    With p_x = p_y = p/(2+eta) and p_z = eta*p/(2+eta), eta=1 reduces to
    uniform depolarizing (each Pauli with p/3). The returned object is a
    Qiskit ``pauli_error`` (a mixture of X/Y/Z errors).
    """
    p_total = float(p_total)
    if not 0.0 <= p_total <= 1.0:
        raise ValueError(f"error rate must be in [0, 1], got {p_total}")
    eta = float(eta)
    if eta <= 0:
        raise ValueError(f"bias eta must be positive, got {eta}")
    px = p_total / (2.0 + eta)
    pz = eta * px
    return pauli_error(
        [
            ("X", px),
            ("Y", px),
            ("Z", pz),
            ("I", max(0.0, 1.0 - p_total)),
        ]
    )


def uniform_depolarizing_2q(p_total: float) -> "pauli_error":
    """Two-qubit depolarizing channel (15 non-identity Paulis, equal weights)."""
    p = float(p_total)
    errs = [("II", 1.0 - p)]
    for i in range(1, 16):
        pauli = "".join(
            ["I", "X", "Y", "Z"][(i >> shift) & 3] for shift in range(0, 4, 2)
        )
        errs.append((pauli, p / 15.0))
    return pauli_error(errs)


def t1t2_derived_rates(t1_us: float, t2_us: float, gate_time_ns: float) -> Dict[str, float]:
    """Single-qubit Pauli rates from T1/T2 coherence times.

    Standard mapping (relaxation + dephasing approximated as Pauli channel):

        p_x = p_y = (1 - exp(-t/T1)) / 4
        p_z = (1 - exp(-t/T2)) / 2 - (1 - exp(-t/T1)) / 4

    Returns ``{"p_x": ..., "p_y": ..., "p_z": ...}``.
    """
    t = gate_time_ns * 1e-3  # ns -> us
    p_therm = (1.0 - np.exp(-t / t1_us)) / 4.0
    p_deph = (1.0 - np.exp(-t / t2_us)) / 2.0
    p_z = p_deph - p_therm
    return {"p_x": float(p_therm), "p_y": float(p_therm), "p_z": float(max(0.0, p_z))}


def profile_from_t1t2(
    t1_us: float,
    t2_us: float,
    gate_time_ns: float = 50.0,
    two_qubit_error_rate: float = 8e-3,
    readout_error: float = 0.02,
) -> Dict:
    """Synthesize a calibration profile from coherence times (T1, T2)."""
    rates = t1t2_derived_rates(t1_us, t2_us, gate_time_ns)
    p_total = rates["p_x"] + rates["p_y"] + rates["p_z"]
    eta = rates["p_z"] / rates["p_x"] if rates["p_x"] > 0 else 1.0
    return {
        "name": f"t1t2_derived_T1={t1_us}us_T2={t2_us}us",
        "basis_gates": DEFAULT_BASIS_GATES,
        "single_qubit_error_rate": float(p_total),
        "two_qubit_error_rate": two_qubit_error_rate,
        "bias_eta": float(eta),
        "readout_error": readout_error,
        "readout_asymmetry": 1.2,
        "t1_us": t1_us,
        "t2_us": t2_us,
        "gate_time_ns": gate_time_ns,
    }


def build_noise_model(profile: Dict) -> NoiseModel:
    """Construct an Aer NoiseModel from a calibration profile dict.

    Errors are keyed on the profile's basis gates; transpile circuits to
    ``profile["basis_gates"]`` so every gate carries its error (the
    :class:`qec_vqe.expectation.ExpectationEvaluator` does this).
    """
    basis = list(profile.get("basis_gates", DEFAULT_BASIS_GATES))
    p1 = float(profile["single_qubit_error_rate"])
    p2 = float(profile["two_qubit_error_rate"])
    eta = float(profile.get("bias_eta", 1.0))
    r_err = float(profile.get("readout_error", 0.0))
    r_asym = float(profile.get("readout_asymmetry", 1.0))

    model = NoiseModel(basis_gates=basis)

    # single-qubit gates: asymmetric (biased) Pauli channel
    one_q_gates = [g for g in basis if g != "cx"]
    for gate in one_q_gates:
        model.add_all_qubit_quantum_error(biased_depolarizing(p1, eta), [gate])

    # two-qubit gates: uniform depolarizing
    if "cx" in basis:
        model.add_all_qubit_quantum_error(uniform_depolarizing_2q(p2), ["cx"])

    # readout: asymmetric assignment error per qubit
    if r_err > 0:
        p_01 = r_err
        p_10 = r_err * r_asym
        readout = ReadoutError([[1.0 - p_01, p_01], [p_10, 1.0 - p_10]])
        model.add_all_qubit_readout_error(readout)

    return model


def default_profiles() -> Dict[str, Dict]:
    """Load the shipped calibration profiles (see ``profiles/``)."""
    profiles = {}
    for path in sorted(PROFILES_DIR.glob("*.json")):
        with open(path) as fh:
            profiles[path.stem] = json.load(fh)
    return profiles


def load_profile(name_or_path: str) -> Dict:
    """Load a profile by shipped name (e.g. "ibm_brisbane_style") or file path."""
    p = Path(name_or_path)
    if p.exists():
        with open(p) as fh:
            return json.load(fh)
    path = PROFILES_DIR / f"{name_or_path}.json"
    if path.exists():
        with open(path) as fh:
            return json.load(fh)
    raise FileNotFoundError(
        f"no profile named {name_or_path!r}; shipped profiles: "
        f"{', '.join(sorted(p.stem for p in PROFILES_DIR.glob('*.json')))}"
    )


def make_profile(
    name: str = "custom",
    single_qubit_error_rate: float = 3e-4,
    two_qubit_error_rate: float = 8e-3,
    bias_eta: float = 5.0,
    readout_error: float = 0.02,
    readout_asymmetry: float = 1.2,
    basis_gates: Optional[list] = None,
) -> Dict:
    """Helper to build a profile dict from explicit parameters."""
    return {
        "name": name,
        "basis_gates": basis_gates or DEFAULT_BASIS_GATES,
        "single_qubit_error_rate": single_qubit_error_rate,
        "two_qubit_error_rate": two_qubit_error_rate,
        "bias_eta": bias_eta,
        "readout_error": readout_error,
        "readout_asymmetry": readout_asymmetry,
    }