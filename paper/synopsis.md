# Project Synopsis — IRIS National Fair 2026–27

**Title:** When Does Zero-Noise Extrapolation Help? A Controlled Bias–Variance Analysis of Error Mitigation for Variational Quantum Chemistry

**Author:** Sahibjotsingh (Class 12, independent research)
**Category:** Physics & Astronomy / Mathematics & Computational Sciences

## The Problem

Quantum computers could simulate molecules that classical computers cannot, but today's processors are too noisy: every gate fails with ~0.1–1% probability, and readout misreads outcomes ~2% of the time. The Variational Quantum Eigensolver (VQE) estimates molecular ground-state energies despite this noise, but under realistic device noise the hydrogen molecule's energy error inflates to ~300 millihartree — about 190× worse than "chemical accuracy" (1.6 mHa), the threshold chemists require. Software fixes exist — zero-noise extrapolation (ZNE), readout mitigation (TREX) — but research papers report each fix separately, so a practitioner with a limited measurement budget cannot know which fix to buy.

## What's New

This project is a **controlled experimental comparison**: same circuit, same calibrated noise model, same shot budget, many random seeds — measuring the **bias and variance** of four mitigation stacks head-to-head. The key results are quantitative exchange rates, not qualitative claims:

1. **ZNE is a bias–variance tradeoff.** A 7-point Richardson extrapolant cuts the deterministic energy bias 3200-fold (303.6 → 0.094 mHa), but amplifies per-run statistical noise 25× (σ: 5.6 → 140 mHa). Lower-order extrapolants (2-point, 3-point) sit at measured points along a bias–variance frontier that matches theoretical amplification factors (1.6, 5.2, 5214). TREX removes readout bias with zero variance penalty.
2. **ZNE has a finite operating window in circuit depth.** Recovery is 99.97% at 56 entangling gates, 95.5% at 168, and exhausted entirely at 1616 gates (LiH) at near-term error rates — recovering only at roadmap-grade (100× improved) error rates. This locates the boundary where error *correction* must take over.
3. **Bias-tailored error correction works at linear cost.** A repetition code matched to the dominant (phase-flip) error axis suppresses logical error rates exactly as theory predicts (3p², 10p³, 35p⁴ for distances 3–7) with gate overhead of exactly 6d − 6 per cycle — the sub-quadratic scaling required for practical correction.
4. **In-loop mitigation rescues optimization.** An optimizer running on the raw noisy objective stalls at 20.4 mHa; the same optimizer with a ZNE-corrected cost function converges to 0.09 mHa.

## Method

Built a complete open-source simulation pipeline (Python, Qiskit): molecular Hamiltonians (H₂, LiH) validated against exact diagonalization to 0.000 mHa; hardware-calibrated asymmetric noise models (phase-flip bias η = 5, per gate-dependent error rates); three evaluation modes (exact, density-matrix, shot-sampled at 20,000 shots); implementations of ZNE (noise-rate scaling and unitary folding) and calibrated readout inversion; distance-3/5/7 repetition codes with full syndrome extraction and minimum-weight decoding; and a multi-seed statistics protocol. 49 unit tests pass; every figure regenerates from one command.

## Impact

The measured bias–variance frontier is a decision tool: it tells experimentalists how to split a fixed measurement budget between noise amplification, calibration, and sampling — the same budgeting problem faced by recent industrial experiments demonstrating quantum utility. The depth-window result gives device engineers a concrete requirement (≈100× gate improvement for LiH-scale chemistry under mitigation alone). All code, data and a hardware-validation pipeline for IBM Quantum processors are open-sourced for reproduction.

## Future Work

Run the head-to-head on real IBM hardware via the included validation pipeline; add measurement-basis twirling (full TREX) and probabilistic error cancellation; extend to surface-code patches and tensor-network simulation beyond 20 qubits.
