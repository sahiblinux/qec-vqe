# qec-vqe — Hardware-Aware Quantum Error Correction & Noise Mitigation for NISQ VQE

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![CI](https://github.com/sahiblinux/qec-vqe/actions/workflows/ci.yml/badge.svg)](https://github.com/sahiblinux/qec-vqe/actions/workflows/ci.yml)

A complete, tested simulation pipeline for studying how quantum error mitigation
and error correction behave on device-calibrated noise — built around a
controlled **bias–variance comparison** of mitigation stacks.

If you use this software, please cite it (see [Citation](#citation)).

**Headline results** (all regenerable with one command):

| Finding | Number |
|---|---|
| Noiseless VQE (H₂, LiH) vs exact diagonalization | **0.000 mHa** both |
| H₂ raw energy error under device-grade noise | 303.6 mHa (190× chemical accuracy) |
| ZNE-corrected (deterministic) | **0.094 mHa** — 3200× reduction |
| ZNE statistical cost at 20k shots | variance amplified 25× (σ 5.6 → 140 mHa) |
| TREX readout mitigation | −42 mHa bias, zero variance penalty |
| ZNE recovery vs depth (56 → 168 CX) | 99.97% → 95.5%; exhausted at 1616 CX |
| Repetition-code logical error vs theory | 3p², 10p³, 35p⁴ (no free parameters) |
| QEC gate overhead | exactly linear: 6d − 6 per cycle |

## Installation

Requires Python 3.9+.

```bash
git clone https://github.com/sahiblinux/qec-vqe
cd qec-vqe
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
```

## Quickstart

```bash
.venv/bin/pytest tests/ -q                 # 49 unit tests, < 1 min, CPU only
.venv/bin/python - <<'EOF'
from qec_vqe import build_hamiltonian, uccsd_ansatz, SPSA
h = build_hamiltonian("H2")
vqe = SPSA(h, uccsd_ansatz(h))
print(vqe.minimize(shots=None))            # noiseless ground-state search
EOF
.venv/bin/python scripts/run_milestones.py # all benchmarks + figures -> out/
```

Results, raw data: `out/results.json`. Figures: `out/fig_*.png`.

## Layout

```
qec_vqe/
  hamiltonian.py      # H2 / LiH STO-3G -> qubit Hamiltonians (JW), validated
  ansatz.py           # UCCSD + hardware-efficient ansatz
  expectation.py      # exact / density-matrix / sampled evaluators, Pauli grouping
  noise_models.py     # asymmetric Pauli channels (p_z != p_x) + readout error
  vqe.py              # SPSA + SciPy optimizers
  mitig/zne.py        # ZNE: noise-rate scaling (production) + unitary folding
  mitig/readout.py    # calibrated assignment-matrix inversion (TREX-style)
  qec/repetition_code.py  # distance-d codes, syndrome extraction, decoding
  benchmark.py        # every study in the paper; incremental caching
profiles/             # calibration-style noise profiles (JSON)
scripts/
  run_milestones.py   # end-to-end driver (--m1 / --mitigation / --qec)
  run_ibm_hardware.py # hardware validation (real IBM device or --dry-run)
  make_paper_html.py  # paper/paper.md -> print-ready HTML
paper/                # JOSS paper (paper.md) + full research paper (research_paper.md)
```

## Reproducing the paper

```bash
.venv/bin/python scripts/run_milestones.py     # regenerates every number + figure
.venv/bin/python scripts/make_paper_html.py    # build paper/paper.html
# open paper/paper.html in a browser -> Print -> Save as PDF
```

## Running on real IBM hardware

```bash
export QISKIT_IBM_TOKEN=...   # free account: https://quantum.cloud.ibm.com
.venv/bin/python scripts/run_ibm_hardware.py            # least-busy backend
.venv/bin/python scripts/run_ibm_hardware.py --dry-run  # identical protocol on Aer
```

The protocol (8 readout-calibration + 15 estimation circuits, ~4k shots each)
fits comfortably in the free Open Plan quota and closes the loop: device
results vs the simulator prediction vs exact diagonalization, with TREX
calibrated **on the device itself**.

## Method notes (interview defense)

- **Noise-rate scaling ZNE** was chosen over unitary folding after discovering
  the transpiler cancels folded gate pairs (`mitig/zne.py` keeps both).
- **Count bitstrings** are parsed as state indices (`int(bs, 2)`), matching
  Aer's big-endian printing; Pauli labels are little-endian — both conventions
  are unit-tested.
- **Energy offsets** (nuclear repulsion, LiH frozen-core constant) live outside
  the Pauli operator and are added explicitly (`MolecularHamiltonian.energy_offset`).
- All statistics use fixed seeds; the head-to-head holds circuit, noise and
  shot budget constant across estimator arms.

## Citation

If you use `qec-vqe` in your work, please cite:

> Singh, S. (2026). *qec-vqe: A Tested Python Framework for Controlled
> Bias–Variance Comparison of Quantum Error-Mitigation Stacks in NISQ VQE*.
> Journal of Open Source Software (in review). DOI: *(assigned on acceptance —
> also minted via Zenodo for the repository archive; see CITATION.cff)*.

The accompanying research paper (`paper/research_paper.md`) and all raw
benchmark data (`out/results.json`) are archived with the repository.
