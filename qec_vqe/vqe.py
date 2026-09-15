"""Variational Quantum Eigensolver loop.

Uses SPSA (simultaneous perturbation stochastic approximation), which is
robust to the noisy cost landscapes typical of NISQ devices. Supports
layer-by-layer initialization to counter vanishing gradients: parameters are
optimized in blocks (one ansatz layer at a time) over several passes.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, List, Optional, Sequence

import numpy as np

from .ansatz import ParameterizedAnsatz


@dataclass
class SPSAOptions:
    """Gain-sequence hyperparameters: a_k = a/(k+1+A)^alpha, c_k = c/(k+1)^gamma."""

    maxiter: int = 100
    a: float = 0.5
    c: float = 0.1
    A: int = 0
    alpha: float = 0.602
    gamma: float = 0.101
    seed: Optional[int] = None

    def __post_init__(self):
        self._rng = np.random.default_rng(self.seed)

    def gain_a(self, k: int) -> float:
        return self.a / (k + 1 + self.A) ** self.alpha

    def gain_c(self, k: int) -> float:
        return self.c / (k + 1) ** self.gamma

    def perturbation(self, n_params: int) -> np.ndarray:
        """Bernoulli +/-1 perturbation vector (per iteration)."""
        return self._rng.choice([-1.0, 1.0], size=n_params)


@dataclass
class VQEResult:
    params: np.ndarray
    energy: float
    history: List[float] = field(default_factory=list)
    n_evaluations: int = 0
    converged: bool = False


def spsa_minimize(
    objective: Callable[[np.ndarray], float],
    n_params: int,
    options: Optional[SPSAOptions] = None,
    initial_params: Optional[np.ndarray] = None,
    callback: Optional[Callable[[int, np.ndarray, float], None]] = None,
) -> VQEResult:
    """Minimize ``objective`` with SPSA, returning params + energy history."""
    opts = options or SPSAOptions()
    theta = (
        np.zeros(n_params)
        if initial_params is None
        else np.asarray(initial_params, dtype=float).copy()
    )
    history: List[float] = []
    n_eval = 0
    best_theta, best_energy = theta.copy(), objective(theta)
    n_eval += 1
    history.append(best_energy)

    for k in range(opts.maxiter):
        delta = opts.perturbation(n_params)
        e_plus = objective(theta + opts.gain_c(k) * delta)
        e_minus = objective(theta - opts.gain_c(k) * delta)
        n_eval += 2
        grad = (e_plus - e_minus) / (2 * opts.gain_c(k) * delta)
        theta = theta - opts.gain_a(k) * grad
        energy = objective(theta)
        n_eval += 1
        history.append(energy)
        if energy < best_energy:
            best_theta, best_energy = theta.copy(), energy
        if callback is not None:
            callback(k + 1, theta, energy)

    return VQEResult(
        params=best_theta,
        energy=best_energy,
        history=history,
        n_evaluations=n_eval,
        converged=bool(n_eval > 3),
    )


def minimize_layers(
    objective: Callable[[np.ndarray], float],
    blocks: Sequence[Sequence[int]],
    options: Optional[SPSAOptions] = None,
    rounds: int = 2,
    initial_params: Optional[np.ndarray] = None,
    callback: Optional[Callable[[int, np.ndarray, float], None]] = None,
) -> VQEResult:
    """Layer-wise SPSA: optimize each parameter block in turn over several rounds.

    Holding all other parameters fixed, each block is optimized with a short
    SPSA run, then the next block is tackled. This prevents early layers from
    dominating the landscape (barren-plateau mitigation).
    """
    opts = options or SPSAOptions()
    n_params = sum(len(b) for b in blocks)
    theta = (
        np.zeros(n_params)
        if initial_params is None
        else np.asarray(initial_params, dtype=float).copy()
    )
    history: List[float] = []
    n_eval = 0
    energy = objective(theta)
    n_eval += 1
    history.append(energy)

    def blocked_objective(sub: np.ndarray, active: Sequence[int]) -> float:
        nonlocal n_eval
        full = theta.copy()
        for j, idx in enumerate(active):
            full[idx] = sub[j]
        e = objective(full)
        n_eval += 1
        return e

    for r in range(rounds):
        for block in blocks:
            sub_opts = SPSAOptions(
                maxiter=max(20, opts.maxiter // max(1, len(blocks))),
                a=opts.a,
                c=opts.c,
                alpha=opts.alpha,
                gamma=opts.gamma,
                seed=opts.seed,
            )
            result = spsa_minimize(
                lambda sub, b=block: blocked_objective(sub, b),
                len(block),
                sub_opts,
                initial_params=theta[list(block)],
            )
            for j, idx in enumerate(block):
                theta[idx] = result.params[j]
            energy = objective(theta)
            n_eval += 1
            history.append(energy)
            if callback is not None:
                callback(r * len(blocks) + len(history), theta, energy)

    return VQEResult(
        params=theta,
        energy=energy,
        history=history,
        n_evaluations=n_eval,
        converged=True,
    )


def scipy_minimize_vqe(
    ansatz: ParameterizedAnsatz,
    objective: Callable[[np.ndarray], float],
    method: str = "COBYLA",
    maxiter: int = 800,
    tol: float = 1e-6,
    initial_params: Optional[np.ndarray] = None,
    seed: Optional[int] = None,
    callback: Optional[Callable[[int, np.ndarray, float], None]] = None,
) -> VQEResult:
    """Noiseless-baseline VQE using a SciPy optimizer (exact gradients-free).

    COBYLA/BFGS on an exact (statevector) cost function converge reliably to
    chemical accuracy for the small molecules in this package; use SPSA
    (:func:`run_vqe`) when the cost function itself is noisy.
    """
    from scipy.optimize import minimize

    n_params = ansatz.num_parameters
    if initial_params is None:
        rng = np.random.default_rng(seed if seed is not None else 0)
        initial_params = rng.uniform(-0.3, 0.3, n_params)

    history: List[float] = []
    history.append(float(objective(np.asarray(initial_params))))

    def wrapped(theta):
        e = objective(np.asarray(theta))
        history.append(float(e))
        if callback is not None:
            callback(len(history), np.asarray(theta), e)
        return e

    res = minimize(
        wrapped,
        np.asarray(initial_params),
        method=method,
        options={"maxiter": maxiter, "tol": tol},
    )
    return VQEResult(
        params=np.asarray(res.x),
        energy=float(res.fun),
        history=history,
        n_evaluations=len(history),
        converged=bool(res.success or abs(res.fun - history[-1]) < 1e-8),
    )


def run_vqe(
    ansatz: ParameterizedAnsatz,
    objective: Callable[[np.ndarray], float],
    maxiter: int = 100,
    layer_wise: bool = False,
    rounds: int = 2,
    seed: Optional[int] = None,
    initial_params: Optional[np.ndarray] = None,
    callback: Optional[Callable[[int, np.ndarray, float], None]] = None,
    spsa: Optional[SPSAOptions] = None,
) -> VQEResult:
    """Run the VQE optimization loop with the given objective function.

    ``objective`` must map a parameter vector to an energy (Hartree).
    """
    opts = spsa or SPSAOptions(maxiter=maxiter, seed=seed)
    if layer_wise and len(ansatz.blocks) > 1:
        return minimize_layers(
            objective,
            [idx for _, idx in ansatz.blocks],
            options=opts,
            rounds=rounds,
            initial_params=initial_params,
            callback=callback,
        )
    return spsa_minimize(
        objective, ansatz.num_parameters, opts, initial_params, callback
    )