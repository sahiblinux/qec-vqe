---
title: 'qec-vqe: A Tested Python Framework for Controlled Bias–Variance Comparison of Quantum Error-Mitigation Stacks in NISQ VQE'
tags:
  - Python
  - quantum computing
  - quantum error mitigation
  - zero-noise extrapolation
  - readout mitigation
  - quantum error correction
  - variational quantum eigensolver
authors:
  - name: Sahibjot Singh
    affiliation: 1
affiliations:
  - name: Independent researcher
    index: 1
date: 15 September 2026
bibliography: paper.bib
---

# Summary

`qec-vqe` is a Python framework for studying how quantum error-mitigation (EM) and
error-correction (QEC) techniques behave on device-calibrated noise models, built
around a controlled comparison of mitigation stacks on the variational quantum
eigensolver (VQE) [@peruzzo2014variational]. The package provides:

- Molecular Hamiltonian construction (H₂, LiH) via `qiskit-nature` with verified
  frozen-core offsets, and UCCSD or hardware-efficient ansätze with exact,
  density-matrix (noisy), and shot-sampled evaluation modes;
- Configurable asymmetric depolarizing noise models derived from calibration-style
  profiles (1q/2q gate error, readout error, T₁/T₂ decay);
- Zero-noise extrapolation (ZNE) [@temme2017error] by both gate folding and
  noise-rate scaling with Richardson and least-squares estimators, and measurement
  readout mitigation ("TREX"-style calibration inversion);
- A repetition-code QEC layer with syndrome-extraction circuits and software
  decoding, reproducing textbook distance scaling (3p², 10p³, 35p⁴);
- Benchmark drivers that reproduce every figure and table of the accompanying
  research paper from a single command, with multi-seed statistics.

The framework is designed for controlled experiments: the same circuit, shot
budget, and noise instance can be evaluated through different estimator stacks
(raw, ZNE, readout-mitigated, combined), enabling bias–variance analyses that are
difficult to run ad hoc. It is intended as a research and teaching instrument for
the NISQ-era EM/QEC community, and every numerical claim in the documentation is
regenerable from the repository.

# Statement of need

Published comparisons of EM techniques are often method-specific: a paper
introducing an extrapolation scheme evaluates it on its own terms, making
controlled cross-technique comparison — equal circuit, equal shot budget, equal
noise instance — surprisingly rare. General-purpose frameworks (e.g. Mitiq
[@mitiq2022], Qiskit Runtime error mitigation) prioritize breadth of technique
support over experimental control and pedagogical transparency.

`qec-vqe` fills this niche with a deliberately small, fully tested surface
(49 unit tests) in which:

1. **Every stage is inspectable.** Grouping, expectation evaluation, transpilation,
   folding, calibration, and decoding are separate, unit-tested functions —
   including the endianness and parameter-ordering conventions that are the most
   common source of silent errors in hand-rolled VQE stacks (documented in the
   repository's method notes).
2. **Estimators are composable under identical conditions.** The benchmark suite
   evaluates raw, ZNE (several estimator orders), and readout-mitigated arms on
   the same circuit and shot budget with seeded statistics, exposing the
   bias–variance tradeoff of ZNE: e.g. on H₂ under device-grade noise, ZNE reduces
   deterministic error from 303.6 to 0.094 mHa while amplifying estimator variance
   ~25× at a 20k-shot budget, whereas readout mitigation removes readout bias at
   essentially zero variance cost.
3. **Scaling boundaries are reproducible.** Recovery-vs-depth studies locate where
   EM stops working (H₂ 56→168 CX, LiH 1616 CX) and where QEC must take over, and
   the QEC layer reproduces repetition-code distance scaling with zero free
   parameters.

The package has been used to produce a complete research manuscript (included in
the repository) whose figures and tables regenerate from
`python scripts/run_milestones.py`.

# Key features and example usage

```python
from qec_vqe import build_hamiltonian, uccsd_ansatz, SPSA

h = build_hamiltonian("H2")            # verified against exact diagonalization
ansatz = uccsd_ansatz(h)                # transpiled, simulation-ready
vqe = SPSA(h, ansatz)                   # stochastic optimizer
result = vqe.minimize(shots=None)       # noiseless; noisy modes available
```

Benchmarks (mitigation ladder, ZNE bias–variance frontier, recovery-vs-depth,
QEC distance/rounds scaling):

```bash
python scripts/run_milestones.py        # regenerates out/*.png + out/results.json
```

# Related work

EM surveys and technique papers [@cai2023quantum; @temme2017error;
@gururani2022-experimental] establish the individual methods; Mitiq [@mitiq2022]
provides broad EM tooling. `qec-vqe` differs in focus: a compact, fully tested
experimental harness for controlled cross-estimator comparison with end-to-end
reproducibility, plus an integrated QEC layer for mapping where EM ends and QEC
must begin.

# AI usage disclosure

The author used an AI coding assistant (Codebuff, an AI agent for software
development) during the preparation of this submission. Scope of assistance:
scaffolding and refactoring of the `qec_vqe` package, test generation,
implementation of benchmark drivers and figure generation, drafting and
copy-editing of documentation, CI configuration, and assistance with manuscript
preparation. All AI-assisted code and text were reviewed, edited, and validated
by the author, who takes full responsibility for the content of this submission
and for the accuracy of all reported results.

# References
