"""Zero-Noise Extrapolation (ZNE).

Noise scaling is performed by *unitary folding*: the full circuit unitary U is
replaced by ``U (U^dag U)^n``, multiplying the effective noise by roughly the
odd scale factor ``lambda = 2n + 1``. Noiseless expectation values are
invariant under folding, so any observed change is pure noise signature that
can be extrapolated away at ``lambda -> 0``.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Sequence

import numpy as np
from qiskit import QuantumCircuit
from scipy.optimize import curve_fit

SCALE_FACTORS = (1, 3, 5)


def fold_circuit(circuit: QuantumCircuit, scale_factor: float) -> QuantumCircuit:
    """Return a unitary-folded copy of ``circuit`` (odd scale factor).

    Only the unitary part is folded; measurements attached to the original
    circuit are re-appended (by classical-bit index) so the folded circuit can
    be measured directly.
    """
    if scale_factor < 1:
        raise ValueError(f"scale factor must be >= 1, got {scale_factor}")
    n = int(round((scale_factor - 1) / 2))
    if 2 * n + 1 != scale_factor:
        raise ValueError(f"scale factor must be odd (2n+1), got {scale_factor}")
    if n == 0:
        return circuit.copy()

    base = QuantumCircuit(circuit.num_qubits)
    meas = []  # (qubit index, clbit index) pairs, in order
    for op in circuit.data:
        if op.operation.name == "measure":
            q = circuit.find_bit(op.qubits[0]).index
            c = circuit.find_bit(op.clbits[0]).index
            meas.append((q, c))
        elif op.operation.name in ("reset", "barrier"):
            continue
        else:
            base.append(op.operation, op.qubits)
    inverse = base.inverse()

    folded = QuantumCircuit(circuit.num_qubits, circuit.num_clbits)
    folded.compose(base, inplace=True)
    for _ in range(n):
        folded.compose(inverse, inplace=True)
        folded.compose(base, inplace=True)

    # re-attach measurements by (qubit, clbit) index
    for q, c in meas:
        folded.measure(folded.qubits[q], folded.clbits[c])
    return folded


def zne_extrapolate(
    scales: Sequence[float], energies: Sequence[float], method: str = "exponential"
) -> float:
    """Extrapolate energies measured at ``scales`` to the zero-noise limit.

    Methods
    -------
    linear:
        first-order polynomial through the two lowest scales.
    richardson:
        polynomial of degree ``len(scales)-1`` evaluated at scale 0.
    exponential:
        fit ``E(lambda) = a + b * exp(c * lambda)`` (physical for decaying
        noise signatures), falling back to Richardson if the fit fails.
    """
    scales = np.asarray(scales, dtype=float)
    energies = np.asarray(energies, dtype=float)
    order = np.argsort(scales)
    scales, energies = scales[order], energies[order]

    if method == "linear":
        poly = np.polyfit(scales, energies, 1)
        return float(np.polyval(poly, 0.0))
    if method == "richardson":
        degree = max(1, len(scales) - 1)
        poly = np.polyfit(scales, energies, degree)
        return float(np.polyval(poly, 0.0))
    if method == "exponential":
        try:
            def model(x, a, b, c):
                return a + b * np.exp(c * x)

            p0 = [float(energies[-1]), float(energies[0] - energies[-1]), -0.5]
            popt, _ = curve_fit(model, scales, energies, p0=p0, maxfev=20000)
            return float(model(0.0, *popt))
        except Exception:
            return zne_extrapolate(scales, energies, "richardson")
    raise ValueError(f"unknown extrapolation method {method!r}")


def scale_profile(profile: Dict, scale: float) -> Dict:
    """Return a copy of ``profile`` with all error rates multiplied by ``scale``."""
    scaled = dict(profile)
    scaled["single_qubit_error_rate"] = float(profile["single_qubit_error_rate"]) * scale
    scaled["two_qubit_error_rate"] = float(profile["two_qubit_error_rate"]) * scale
    if "readout_error" in profile:
        scaled["readout_error"] = float(profile["readout_error"]) * scale
    return scaled


def run_zne(
    hamiltonian,
    circuit: QuantumCircuit,
    params=None,
    profile: Dict = None,
    scales: Sequence[float] = SCALE_FACTORS,
    method: str = "exponential",
    mode: str = "exact_noisy",
    shots: Optional[int] = None,
    seed: Optional[int] = None,
    readout_filter=None,
) -> Dict:
    """Zero-noise extrapolation by error-rate scaling (noise amplification).

    The noise model is rebuilt at each scale factor with all rates multiplied
    by ``lambda`` (single-qubit, two-qubit and readout). This is the canonical
    ZNE variant that avoids unitary folding altogether: folding is also
    available (:func:`fold_circuit`) but is fragile against transpiler gate
    cancellation. ``hamiltonian`` is a :class:`MolecularHamiltonian`;
    returns ``{"scales", "energies", "extrapolated", "method"}``.
    """
    from ..expectation import ExpectationEvaluator
    from ..noise_models import build_noise_model

    if profile is None:
        raise ValueError("run_zne requires a noise profile dict")
    energies = []
    for scale in scales:
        nm = build_noise_model(scale_profile(profile, scale))
        evaluator = ExpectationEvaluator(
            hamiltonian, mode=mode, noise_model=nm, shots=shots, seed=seed,
            readout_filter=readout_filter,
        )
        energies.append(evaluator.energy(circuit, params))
    extrapolated = zne_extrapolate(scales, energies, method)
    return {
        "scales": list(scales),
        "energies": energies,
        "extrapolated": extrapolated,
        "method": method,
    }