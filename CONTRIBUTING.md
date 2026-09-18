# Contributing to qec-vqe

Thank you for considering contributing! This document explains how to set up,
test, and submit changes.

## Setup

```bash
git clone https://github.com/sahiblinux/qec-vqe
cd qec-vqe
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Running the tests

```bash
python -m pytest tests/ -q
```

All tests must pass before any PR is merged. Tests run on CPU only and complete
in under a minute; no IBM Quantum account is needed for the test suite.

## Reproducing the paper's figures

```bash
python scripts/run_milestones.py
```

writes `out/results.json` and `out/fig_*.png`. Benchmark results are cached
incrementally; delete `out/results.json` for a fully fresh run (this takes
longer — the LiH runs are the expensive part).

## Ground rules

- **Correctness conventions are load-bearing.** Qiskit uses little-endian qubit
  ordering; several silent bugs in this project's history came from endianness,
  parameter sorting, and energy-offset handling. If you touch `expectation.py`,
  `readout.py`, or `hamiltonian.py`, run the full suite and state in your PR what
  convention you preserved.
- **Every number in the paper must regenerate.** If you change a module that
  feeds `benchmark.py`, re-run the affected stage and confirm `out/results.json`
  is unchanged (or explain why it should change).
- Keep the dependency set minimal; new dependencies need a justification in the PR.

## Reporting bugs

Open a GitHub issue with: the command you ran, the full traceback, your
`qiskit`/`qiskit-aer`/`qiskit-nature` versions (`pip show qiskit qiskit-aer qiskit-nature`),
and whether the failure reproduces with a fresh `out/` directory.
