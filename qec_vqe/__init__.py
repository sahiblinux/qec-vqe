"""Hardware-aware QEC and noise mitigation for NISQ-era VQE."""

__version__ = "0.1.0"

from .hamiltonian import (
    MolecularHamiltonian,
    REFERENCE_ENERGIES,
    build_hamiltonian,
    pauli_groups,
)
from .ansatz import (
    ParameterizedAnsatz,
    hardware_efficient_ansatz,
    hartree_fock_state_circuit,
    uccsd_ansatz,
)
from .expectation import ExpectationEvaluator
from .vqe import SPSAOptions, VQEResult, minimize_layers, run_vqe, scipy_minimize_vqe, spsa_minimize
from .noise_models import (
    biased_depolarizing,
    build_noise_model,
    default_profiles,
    load_profile,
    make_profile,
    profile_from_t1t2,
    t1t2_derived_rates,
)

__all__ = [
    "MolecularHamiltonian",
    "REFERENCE_ENERGIES",
    "build_hamiltonian",
    "pauli_groups",
    "ParameterizedAnsatz",
    "hardware_efficient_ansatz",
    "hartree_fock_state_circuit",
    "uccsd_ansatz",
    "ExpectationEvaluator",
    "SPSAOptions",
    "VQEResult",
    "minimize_layers",
    "run_vqe",
    "scipy_minimize_vqe",
    "spsa_minimize",
    "biased_depolarizing",
    "build_noise_model",
    "default_profiles",
    "load_profile",
    "make_profile",
    "profile_from_t1t2",
    "t1t2_derived_rates",
]