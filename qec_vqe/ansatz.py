"""Parameterized ansatz circuits for VQE.

Two ansatz families:

* ``hardware_efficient_ansatz`` - layered Ry/Rz + nearest-neighbor CX, with
  an optional Hartree-Fock initial state. Parameter blocks per layer enable
  layer-by-layer optimization (risk mitigation for vanishing gradients).
* ``uccsd_ansatz`` - chemistry-grade UCCSD from qiskit-nature, starting from
  the Hartree-Fock state (all-zero parameters = HF, physics-informed init).

Both return a :class:`ParameterizedAnsatz` whose ``blocks`` attribute lists
``(name, parameter_indices)`` for layer-wise optimization.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Tuple

import numpy as np
from qiskit import QuantumCircuit, transpile
from qiskit.circuit import ParameterVector

from .hamiltonian import MolecularHamiltonian

#: canonical gate set emitted by ansatz builders (Aer can execute all of these)
_CANONICAL_BASIS = ["rx", "rz", "cx", "x", "h", "s", "sdg", "id"]

try:  # qiskit-nature circuit library (optional but expected)
    from qiskit_nature.second_q.circuit.library import HartreeFock, UCCSD
    from qiskit_nature.second_q.mappers import JordanWignerMapper

    _HAS_NATURE = True
except Exception:  # pragma: no cover - qiskit-nature is a hard dependency
    _HAS_NATURE = False


@dataclass
class ParameterizedAnsatz:
    """A parameterized circuit with named parameter blocks for layer-wise runs."""

    circuit: QuantumCircuit
    parameters: ParameterVector
    blocks: List[Tuple[str, List[int]]] = field(default_factory=list)

    @property
    def num_parameters(self) -> int:
        return len(self.parameters)

    def bind(self, theta: np.ndarray) -> QuantumCircuit:
        return self.circuit.assign_parameters(
            {self.parameters[i]: theta[i] for i in range(len(theta))}
        )

    def block_indices(self, name: str) -> List[int]:
        for n, idx in self.blocks:
            if n == name:
                return idx
        raise KeyError(f"no parameter block named {name!r}")


def hartree_fock_state_circuit(hamiltonian: MolecularHamiltonian) -> QuantumCircuit:
    """Prepare the Hartree-Fock state in the Jordan-Wigner qubit space.

    With Jordan-Wigner mapping (the default used by ``build_hamiltonian``)
    qubit j = spin-orbital j and the HF determinant is a computational basis
    state, so this circuit is exactly consistent with our Hamiltonians.
    """
    if not _HAS_NATURE:
        raise RuntimeError("qiskit-nature is required for Hartree-Fock preparation")
    mapper = JordanWignerMapper()
    hf = HartreeFock(
        num_spatial_orbitals=hamiltonian.num_spatial_orbitals,
        num_particles=tuple(hamiltonian.num_particles),
        qubit_mapper=mapper,
    )
    return hf


def hardware_efficient_ansatz(
    hamiltonian: MolecularHamiltonian,
    layers: int = 2,
    entangling: str = "linear",
    hf_init: bool = True,
    param_prefix: str = "theta",
) -> ParameterizedAnsatz:
    """Build a layered Ry/Rz + CX hardware-efficient ansatz.

    Parameters
    ----------
    layers:
        Number of entangling layers.
    entangling:
        "linear" (nearest-neighbor CX) or "full" (all-to-all).
    hf_init:
        Prepend the Hartree-Fock state (physics-informed initialization).
    """
    n = hamiltonian.num_qubits
    pv = ParameterVector(param_prefix, 2 * n * layers)
    qc = QuantumCircuit(n)
    if hf_init:
        qc.compose(hartree_fock_state_circuit(hamiltonian), inplace=True)

    blocks = []
    for layer in range(layers):
        idx = 2 * n * layer
        for q in range(n):
            qc.ry(pv[idx + 2 * q], q)
            qc.rz(pv[idx + 2 * q + 1], q)
        if entangling == "linear":
            for q in range(n - 1):
                qc.cx(q, q + 1)
        elif entangling == "full":
            for q in range(n):
                for q2 in range(q + 1, n):
                    qc.cx(q, q2)
        else:
            raise ValueError(f"unknown entangling pattern {entangling!r}")
        blocks.append((f"layer{layer}", list(range(idx, idx + 2 * n))))

    return ParameterizedAnsatz(circuit=qc, parameters=pv, blocks=blocks)


def uccsd_ansatz(
    hamiltonian: MolecularHamiltonian, reps: int = 1, include_imaginary: bool = False
) -> ParameterizedAnsatz:
    """Build a UCCSD ansatz (qiskit-nature) on the Jordan-Wigner space.

    With all-zero parameters the circuit reduces to the Hartree-Fock state.
    """
    if not _HAS_NATURE:
        raise RuntimeError("qiskit-nature is required for UCCSD")
    mapper = JordanWignerMapper()
    initial = HartreeFock(
        num_spatial_orbitals=hamiltonian.num_spatial_orbitals,
        num_particles=tuple(hamiltonian.num_particles),
        qubit_mapper=mapper,
    )
    uccsd = UCCSD(
        num_spatial_orbitals=hamiltonian.num_spatial_orbitals,
        num_particles=tuple(hamiltonian.num_particles),
        qubit_mapper=mapper,
        reps=reps,
        initial_state=initial,
        include_imaginary=include_imaginary,
    )
    # qiskit-nature emits EvolvedOps/PauliEvolution custom instructions that
    # Aer cannot execute directly; transpile to the canonical gate set used
    # throughout this package (rx/rz/cx/x/h/s/sdg).
    pv = ParameterVector("theta", uccsd.num_parameters)
    circuit = uccsd.assign_parameters(
        {uccsd.parameters[i]: pv[i] for i in range(uccsd.num_parameters)}
    )
    circuit = transpile(
        circuit,
        basis_gates=_CANONICAL_BASIS,
        optimization_level=0,
    )
    return ParameterizedAnsatz(
        circuit=circuit,
        parameters=pv,
        blocks=[("uccsd", list(range(uccsd.num_parameters)))],
    )