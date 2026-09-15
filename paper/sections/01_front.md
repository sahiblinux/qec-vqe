# When Does Zero-Noise Extrapolation Help? A Controlled Bias–Variance Analysis of Error Mitigation for Variational Quantum Chemistry

**Author:** Sahibjot Singh
**Category:** Computational Biology and Bioinformatics / Physics and Astronomy (quantum computing)
**Platform:** Simulation study with hardware-calibrated noise models; validation pipeline for IBM Quantum hardware included
**Code:** `qec_vqe` package, 49 passing unit tests; all figures reproducible via `scripts/run_milestones.py`

---

## Abstract

Variational Quantum Eigensolvers (VQE) promise molecular ground-state energies on today's noisy intermediate-scale quantum (NISQ) processors, but device noise biases the estimated energy severely — under a device-calibrated noise model, the hydrogen molecule's energy error grows from 0.000 mHa (noiseless) to 303.6 mHa (noisy), roughly 190× the "chemical accuracy" threshold of 1.6 mHa that computational chemistry requires. Error-mitigation techniques repair such biases in software, but published results usually report each technique in isolation, making it impossible to know *which* mitigation to buy with a *fixed* shot budget. This work answers that question with a controlled comparison. Under identical circuits, noise profiles and shot budgets, I measure the **bias and variance** of four estimator stacks — raw sampling, twirled-readout extinction (TREX-style calibrated readout inversion), zero-noise extrapolation (ZNE), and ZNE+TREX — across independent random seeds. Three findings emerge. (1) ZNE is a bias–variance tradeoff, not a free lunch: a 7-point Richardson extrapolant reduces the deterministic bias 3200-fold (303.6 → 0.094 mHa) but amplifies statistical noise by a measured factor of up to 25× per run (σ: 5.6 → 140 mHa), while TREX removes readout bias with no variance penalty. (2) ZNE has a finite operating window in circuit depth: error recovery falls 99.97% → 99.1% → 95.5% as CX count grows 56 → 112 → 168, and is exhausted entirely for a 1616-CX lithium-hydride circuit at near-term error rates — recovering only at roadmap-grade (100× improved) error rates. (3) A noise-bias-tailored repetition code corrects the dominant error axis with logical error rates matching theory (3p², 10p³, 35p⁴ for distances 3, 5, 7) at strictly linear gate overhead (gates = 6d − 6, R² ≈ 1), satisfying the sub-quadratic overhead requirement for practical error correction. The results reframe mitigation selection as budget allocation on an empirically measured bias–variance frontier and provide open-source tooling to reproduce every number.

---

## Research Questions

**RQ1.** At a fixed shot budget on a fixed circuit and noise profile, which error-mitigation stack — raw, TREX-only, ZNE-only, or ZNE+TREX — achieves the best bias–variance tradeoff for VQE energy estimation?

**RQ2.** How does ZNE's effectiveness degrade as circuit depth (entangling-gate count) increases, and where exactly does it stop working?

**RQ3.** Does a repetition code tailored to a biased noise channel (phase-flip rate ≫ bit-flip rate) suppress the logical error rate as theory predicts, at gate overhead that scales sub-quadratically with code distance?

These questions target a gap in the literature: mitigation methods are standardly validated one at a time against a raw baseline [1,2,3], so a practitioner with 20,000 shots per observable has no principled answer to "which mitigation should I run?"
