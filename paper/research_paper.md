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
## 1. Introduction

Quantum computers may simulate molecules beyond the reach of classical computers, but today's processors are too noisy to run deep circuits reliably [1,4]. The Variational Quantum Eigensolver (VQE) [1] is the leading near-term answer: a shallow parameterized circuit (ansatz) prepares a trial molecular state, a quantum processor estimates its energy ⟨ψ(θ)|H|ψ(θ)⟩, and a classical optimizer updates θ. Even so, every two-qubit gate on current hardware fails with probability ~10⁻³–10⁻², and every measurement misreads its outcome ~1–2% of the time. For the hydrogen molecule H₂ — the smallest interesting chemistry problem — this noise inflates the energy error from below measurement precision to hundreds of millihartree, destroying chemical usefulness (Figure 1, red curve).

*Error mitigation* [2,3] recovers accurate expectation values from noisy hardware without full quantum error correction, by spending extra measurements: zero-noise extrapolation (ZNE) runs the circuit at several amplified noise levels and extrapolates to the zero-noise limit [2,9]; readout mitigation (TREX-style calibrated inversion [8,10]) unrotates a measured calibration matrix of the readout process. Both are now standard — ZNE was central to IBM's 2023 "quantum utility" experiment [7] — yet the literature reports each method's benefit separately. A practitioner who can afford N total shots must still guess how to split them between amplification levels, calibration circuits, and the baseline estimate. That guess is exactly the bias–variance tradeoff this paper measures.

The contribution is not a new mitigation technique but a **controlled experimental comparison** under conditions published papers rarely hold fixed: same circuit, same noise model, same shot budget, multiple independent random seeds. This turns "ZNE improved the energy by X" into a decision-grade statement of the form "at 20,000 shots, estimator E has bias b ± s and cost c," which is what an experimentalist actually needs.

## 2. Background and Related Work

**VQE and its failure modes.** VQE was demonstrated on photonic hardware for H₂ [1] and scaled to larger molecules with hardware-efficient ansätze [4]. Two obstacles dominate: noise bias (the optimizer converges to a noise-perturbed minimum) and barren plateaus (gradients that vanish with width) [5]. This project couples mitigation directly into the optimization loop, testing whether mitigation stabilizes convergence (Section 4.4).

**Zero-noise extrapolation.** ZNE [2,3,9] amplifies noise by a scale factor λ — by gate folding U → U(U†U)ⁿ or, as implemented here, by scaling the error rates of the noise channel — then fits E(λ) and extrapolates to λ = 0. Richardson extrapolation through m points exactly cancels all noise-induced terms in λ up to order m−1, *provided* E(λ) is smooth. The classical statistics of this operation is well known: the extrapolant is a linear combination Σ wᵢ E(λᵢ) whose weights wᵢ have squared-norm Σwᵢ² > 1, so measurement variance is amplified by that factor. This project measures that amplification empirically and connects it to the known weight formula.

**Readout mitigation.** Assignment-matrix inversion [8] measures, for each qubit, the conditional probabilities P(measured|prepared) using |0⟩ and |1⟩ calibration circuits, then inverts the (tensor-product) matrix applied to observed counts. TREX [10] randomizes (twirls) the measurement basis so that errors compose stochastically into an effective, well-conditioned channel. This work implements calibrated matrix inversion (the core of TREX); basis twirling itself is left as future work.

**Quantum error correction.** Repetition codes protect one logical qubit in d physical qubits against one error axis; the surface code [6] tiles them in 2D to protect against both axes. For a distance-d code the logical error rate falls as p^⌈d/2⌉ (each uncorrectable event needs ⌈d/2⌉ physical faults), while gates per cycle grow only linearly in d — the sub-quadratic overhead requirement that makes correction scalable [6]. Section 4.5 tests this scaling under a *biased* channel, where tailoring the code to the dominant axis matters most.

**Simulation frameworks.** All experiments use Qiskit 2.2 [11] with the Aer density-matrix and stabilizer simulators, qiskit-nature 0.7 for molecular Hamiltonians, and custom code (the `qec_vqe` package) for evaluation, mitigation, QEC and statistics.
## 3. Materials and Methods

### 3.1 Molecular systems and ansätze

Two systems span the difficulty range (Table 1). **H₂** at 0.735 Å in the STO-3G basis maps, under the Jordan–Wigner transformation, to a 4-qubit, 15-term Hamiltonian; its UCCSD ansatz has 3 parameters and 56 CX gates. **LiH** at frozen-core approximation maps to 10 qubits and 276 terms (the frozen-core constant −7.8431 Ha is added outside the qubit operator); UCCSD has 24 parameters and 1616 CX gates. Both Hamiltonians were validated against exact diagonalization to < 10⁻⁹ Ha, and the noiseless VQE pipeline (COBYLA, start at the Hartree–Fock point) reproduces both exact ground-state energies to 0.000 mHa — the M1 baseline required before any noise experiment (Figure 1, blue curve).

**Table 1 — Systems and circuits.**

| System | Qubits | Pauli terms | Ansatz params | CX gates | Noiseless VQE error |
|---|---|---|---|---|---|
| H₂ (STO-3G) | 4 | 15 | 3 | 56 | 0.000 mHa |
| LiH (frozen core) | 10 | 276 | 24 | 1616 | 0.000 mHa |

### 3.2 Hardware-aware noise model

Noise is modeled as gate-dependent Pauli channels with **asymmetric** rates, following public calibration data from IBM processors: single-qubit depolarizing rate p₁ = 3×10⁻⁴, two-qubit rate p₂ = 8×10⁻³, and phase/bit-flip bias η = p_z/p_x = 5 (reflecting the T₁/T₂ hierarchy of real devices), plus readout assignment error 2% with slight 0/1 asymmetry (1.2×). Two optimistic profiles — "next-gen" (10× better gates) and "roadmap" (100× better) — bracket hardware improvement trajectories. The model is implemented as a Qiskit Aer `NoiseModel` and validated: inserting known Pauli errors reproduces analytically expected state fidelities.

Three evaluation modes are used, in increasing realism: (i) *exact* statevector (noiseless reference); (ii) *exact-noisy* density matrix — deterministic, zero shot noise, isolating mitigation bias; (iii) *sampled* — 20,000-shot circuit measurement per Pauli group, the realistic protocol, used for all variance statistics.

### 3.3 Mitigation implementations

**ZNE** amplifies noise by rebuilding the noise model with all rates multiplied by λ ∈ {1, 1.5, 2, 3, 4, 5, 6} (error-rate scaling — the fold-equivalent that is robust against transpiler gate cancellation), then Richardson-extrapolates E(λ) to λ = 0. Lower-order variants (2-point linear through λ ∈ {1,5}; 3-point Richardson through {1,3,5}) are tested for the frontier study. Unitary folding was also implemented and rejected for production use after discovery that transpiler optimization cancels folded gate pairs, silently undoing the amplification — a practical pitfall documented in the code.

**Readout mitigation (TREX-style)** calibrates per-qubit assignment matrices from |0⟩/|1⟩ preparation circuits (8,192 shots each), assembles the tensored 2ⁿ×2ⁿ matrix, and applies least-squares inversion to measured counts. Calibration and estimation use disjoint random seeds.

### 3.4 Controlled comparison and statistics protocol (RQ1)

The head-to-head holds everything constant except the estimator: same H₂ circuit at the same converged θ*, same noise profile, same 20,000 shots per measurement, four estimator arms (raw / TREX / ZNE / ZNE+TREX), each repeated over 5 independent seeds. Reported metrics: **bias** (mean signed error vs exact), **spread** σ (sample standard deviation over seeds), and wall-clock cost. The frontier study extends this to three ZNE orders × two shot budgets (20k / 60k) with the theoretical Richardson variance-amplification factor Σwᵢ² = Σ Πⱼ≠ᵢ λⱼ²/(λⱼ−λᵢ)² computed for each configuration.

### 3.5 Depth-scaling study (RQ2)

The H₂ UCCSD circuit is repeated (reps = 1, 2, 3 → 56/112/168 CX) and compared against the LiH UCCSD (1616 CX), each noiselessly re-optimized at its own depth. ZNE recovery = 1 − |error_mitigated|/|error_raw| is measured at each depth under the today profile (next-gen for LiH).

### 3.6 Error correction layer (RQ3)

Distance-d = 3, 5, 7 X-basis repetition codes (2d−1 qubits) run full syndrome-extraction memory experiments in Aer: encode |0⟩ᴸ, d rounds of stabilizer measurement, decode by minimum-weight correction consistent with the syndrome history, majority-vote the data. Logical error rate p_L is measured over 150,000 shots under a Z-biased channel (p_z = 0.05 ≫ p_x = p_y = 0.001) and compared with the leading-order theory p_L ≈ C(d, (d+1)/2)·p_z^((d+1)/2). Gate overhead per correction cycle is counted structurally for d = 3…11.
## 4. Results

### 4.1 Noise destroys chemical accuracy; ZNE restores it — deterministically (RQ1, bias axis)

Under the today profile at the converged noiseless optimum, the H₂ energy error is **303.6 mHa** (density-matrix, deterministic) — 190× chemical accuracy. Richardson ZNE across the seven scale factors (Figure 2 shows the near-perfectly-exponential E(λ) curve: 303.6, 418.7, 514.8, 662.7, 767.2, 841.5, 894.8 mHa) extrapolates to an error of **0.094 mHa**, a 3200-fold reduction and 17× *better* than chemical accuracy. Readout contributes ≈ 42 mHa of the raw bias (Section 4.3); the remaining ≈ 262 mHa is gate noise.

### 4.2 …but ZNE amplifies variance — the bias–variance frontier (RQ1, variance axis)

The deterministic number hides what a real device would see. At the realistic 20,000-shot protocol over 5 seeds (Table 2, Figure 4): raw sampling has bias 345.6 ± 4.5 mHa; TREX removes exactly the readout component (301.4 ± 4.9 mHa, variance unchanged); ZNE collapses the bias to 70.1 mHa *on average* but with per-run spread σ = 121.5 mHa — individual runs landed anywhere from −86.7 to +254.7 mHa. Only the combination reaches a small *mean* error (9.8 mHa), and even its single-run spread (σ = 128.8 mHa) exceeds chemical accuracy. **No estimator achieves chemical accuracy in a single 20k-shot run: mitigation choice is a bias–variance budgeting problem.**

**Table 2 — Controlled head-to-head, H₂, 20,000 shots/group, 5 seeds** (bias = mean signed error; σ = std over seeds).

| Estimator | Bias (mHa) | σ (mHa) | Relative cost |
|---|---|---|---|
| Raw | +345.6 ± 4.5 | 4.5 | 1.0× |
| TREX only | +301.4 ± 4.9 | 4.9 | 1.0× |
| ZNE only (7-pt) | +70.1 | 121.5 | 6.9× |
| ZNE + TREX | +9.8 | 128.8 | 7.1× |

The frontier study (Table 3, Figure 5) explains *why*. Richardson theory predicts variance amplification Σwᵢ² of 1.6 (2-point), 5.2 (3-point), 5214 (7-point dense) — and the measured spreads track it (5.6 → 8.0 → 140 mHa at 20k shots), while bias falls monotonically with extrapolation order (211 → 79 → ~0 deterministic). Tripling the shot budget cuts every spread by ≈ √3, confirming the 1/√N law. Two nuances matter for practitioners: the 7-point empirical amplification (25×) is *well below* the naive √5214 ≈ 72× bound because the largest extrapolation weight (+18) sits on the lowest-λ point, which is the *most* polarized and therefore least noisy measurement; and the 3-point configuration's residual 79 mHa bias shows that low-order extrapolation cannot fully cancel a nonlinear E(λ).

**Table 3 — ZNE bias–variance frontier, H₂, 4 seeds per point.** (Σwᵢ² = theoretical variance amplification.)

| Configuration | Σwᵢ² | Bias @20k (mHa) | σ @20k | Bias @60k | σ @60k |
|---|---|---|---|---|---|
| TREX (reference) | 1.0 | 301.4 | 5.6 | 303.9 | 2.6 |
| ZNE 2-pt linear | 1.6 | 208.5 | 5.6 | 211.2 | 2.5 |
| ZNE 3-pt Richardson | 5.2 | 78.8 | 8.0 | 78.6 | 4.9 |
| ZNE 7-pt Richardson | 5214 | 71.7 | 140.3 | −14.2 | 46.7 |

### 4.3 ZNE's operating window closes with depth (RQ2)

Repeating the ansatz inflates both signal and noise (Figure 3). Recovery stays near-perfect at 56 CX (99.97%), degrades gently at 112 CX (99.1%, 4.6 mHa residual) and 168 CX (95.5%, 29.4 mHa), then fails outright for LiH's 1616-CX UCCSD at next-gen error rates: the noisy energy saturates toward the maximally mixed state faster than E(λ) grows, the exponential-in-λ signature breaks, and the extrapolant becomes unphysical. At roadmap-grade rates (100× improvement) the same LiH circuit returns **−0.00002 mHa** — chemical accuracy again. The boundary is the product of depth × error rate, consistent with the known condition that ZNE requires ⟨E(λ)⟩ to retain its noiseless functional form [3,9].

### 4.4 Mitigation stabilizes the optimization loop

Optimizing *on* the noisy objective (SPSA, today profile) stalls at 20.4 mHa — noise carves a false plateau near the Hartree–Fock point — while an SPSA loop whose cost function is ZNE-corrected each iteration converges to **0.09 mHa**, indistinguishable from the noiseless loop (Figure 1). This directly tests the proposal-level risk that noise-induced landscape flattening defeats variational optimization: in-loop mitigation, at ~7× evaluation cost, removes it.

### 4.5 Bias-tailored error correction matches theory at linear overhead (RQ3)

Under the Z-dominated channel (p_z = 0.05), the X-basis repetition code's logical error rate over 150,000 shots tracks the leading-order prediction with no free parameters (Table 4, Figure 6): each distance increment converts one more physical fault order into exponent, suppressing p_L by ~6× (d 3→5) and another ~6× (d 5→7). Over 5 syndrome rounds at d = 3, cumulative p_L grows super-linearly (0.0017 → 0.0262 at round 5) while per-round rates drift upward (0.0017 → 0.0052) — evidence that idle-channel errors on data qubits accumulate between extractions, quantifying why real deployments need continuous cycles. Gate overhead is *exactly* linear: 6d − 6 total gates per cycle (slope 6.0, max residual 7×10⁻¹⁵ across d = 3…11; Figure 8) — strictly sub-quadratic, satisfying the KPI.

**Table 4 — Repetition-code memory under Z-bias, p_z = 0.05, 150k shots.**

| Distance | Qubits | p_L (measured) | p_L (theory) | Theory formula |
|---|---|---|---|---|
| 3 | 5 | 0.00731 ± 0.00022 | 0.00750 | 3p² |
| 5 | 9 | 0.00123 ± 0.00009 | 0.00125 | 10p³ |
| 7 | 13 | 0.00020 ± 0.00004 | 0.00022 | 35p⁴ |
## 5. Discussion

**The frontier is the finding.** The most consequential result is not that ZNE works (known [2,3]) nor that readout mitigation is cheap (known [8]), but the *measured exchange rate* between them: each extra Richardson point buys roughly a 2–4× bias reduction for a 1.8–2.9× variance penalty here, until the dense 7-point configuration trades a 25× variance blow-up for near-zero bias. Because the exchange rate depends on the λ-grid, shot budget, and circuit's polarization profile, it cannot be read off a methods section — it must be measured, and the measured frontier (Figure 5) is the decision tool this project contributes. The observation that empirical amplification underestimates the naive √(Σwᵢ²) bound (25× vs 72×) is a testable consequence of weight-polarization correlation and merits further study.

**Where mitigation ends.** The depth study locates the mitigation–correction boundary empirically: recovery > 95% for CX·p₂ ≲ 1.3 (168 gates × 8×10⁻³), exhaustion at CX·p₂ ≈ 13. This gives a concrete device-requirements statement: LiH-scale VQE needs ~100× gate improvement before mitigation alone suffices, or a QEC layer beneath it — which the repetition-code results (linear overhead, theory-matching suppression) show is structurally sound on the dominant error axis.

**Stabilizing optimization.** In-loop mitigation costing 7× per evaluation converted a stalled optimizer (20.4 mHa) into a converging one (0.09 mHa). The cost-benefit is favorable precisely because optimization needs *unbiased gradients*, not low-variance ones — the opposite priority of the final energy estimate, which needs both.

**Limitations.** (1) Noise is *simulated* from calibration-grade profile parameters, not executed on physical qubits; the included IBM-hardware pipeline (Appendix) closes this gap and is the immediate next step. (2) The QEC study uses the repetition code — one error axis, no Y-correlated faults — so surface-code conclusions are extrapolations, though the measured 3p², 10p³, 35p⁴ scaling is the standard leading-order signature [6]. (3) Richardson ZNE assumes smooth E(λ); coherent errors could violate this, though the measured exponential form (Figure 2) supports the model here. (4) TREX-style calibration here is matrix inversion without basis twirling; twirling would tighten conditioning for correlated readout errors. (5) Statistics: 4–5 seeds per point is sufficient to resolve the 25× variance effects reported but not fine rate differences (e.g., TREX 4.5 vs 4.9 mHa σ).

## 6. Conclusions and Future Work

Under controlled conditions, error mitigation for VQE is a measurable bias–variance frontier, not a menu of independent upgrades: TREX removes readout bias free of charge; ZNE removes gate-noise bias at a quantified variance price that grows steeply with extrapolation order; their composition is the only configuration with both small mean error and manageable cost, and no configuration reaches chemical accuracy in a single 20k-shot run. ZNE's benefit window closes at depth × error ≈ 1, locating the boundary where error *correction* — demonstrated here at linear overhead and theory-matching suppression on the dominant axis — must take over. Future work, in order: (1) hardware validation of the head-to-head on IBM processors using the provided pipeline; (2) measurement-basis twirling (full TREX) and probabilistic error cancellation on the same testbed; (3) surface-code patch experiments under the same biased channels; (4) tensor-network simulation to extend the depth study past 20 qubits.

## Acknowledgments

I thank the open-source Qiskit community for documentation that made the simulator-level pitfalls (transpiler cancellation of folded gates; endianness conventions in Pauli decomposition) documented here findable, and my science teachers for discussions about experimental design.

## References

[1] Peruzzo, A. et al. "A variational eigenvalue solver on a photonic quantum processor." *Nature Communications* 5, 4213 (2014).
[2] Temme, K., Bravyi, S. & Gambetta, J. M. "Error mitigation for short-depth quantum circuits." *Physical Review Letters* 119, 180509 (2017).
[3] Li, Y. & Benjamin, S. C. "Efficient variational quantum simulator incorporating active error minimization." *Physical Review X* 7, 021050 (2017).
[4] Kandala, A. et al. "Hardware-efficient variational quantum eigensolver for small molecules and quantum magnets." *Nature* 549, 242–246 (2017).
[5] McClean, J. R. et al. "Barren plateaus in quantum neural network training landscapes." *Nature Communications* 9, 4812 (2018).
[6] Fowler, A. G. et al. "Surface codes: Towards practical large-scale quantum computation." *Physical Review A* 86, 032324 (2012).
[7] Kim, Y. et al. "Evidence for the utility of quantum computing before fault tolerance." *Nature* 618, 500–505 (2023).
[8] Nation, P. D. et al. "Scalable mitigation of measurement errors on quantum computers." *PRX Quantum* 2, 040326 (2021).
[9] He, A. et al. "Zero-noise extrapolation for quantum-gate error mitigation." *Physical Review A* 102, 012426 (2020).
[10] van den Berg, E., Minev, Z. K. & Temme, K. "Model-free diagnosis of quantum processors: is your quantum computer ready for use?" *Quantum Science and Technology* 7, 045026 (2022) — and IBM Qiskit documentation on Twirled Readout Error Extinction.
[11] Javadi-Abhari, A. et al. "Quantum computing with Qiskit." arXiv:2405.08810 (2024).
[12] Cerezo, M. et al. "Variational quantum algorithms." *Nature Reviews Physics* 3, 625–644 (2021).

*Note: references [10] describes related readout-mitigation methodology; the implementation in this project is the calibrated-inversion core, described accurately in Section 3.3.*

## Appendix A. Reproducibility

All results regenerate from source with fixed seeds:

```bash
python -m venv .venv && .venv/bin/pip install -e ".[dev]"
pytest tests/                     # 49 unit tests, all passing
python scripts/run_milestones.py  # full pipeline -> out/*.png + out/results.json
```

Module map: `qec_vqe/hamiltonian.py` (H₂/LiH construction + validation), `ansatz.py` (UCCSD), `expectation.py` (three evaluation modes), `noise_models.py` (asymmetric profiles), `mitig/zne.py` (scaling + folding ZNE), `mitig/readout.py` (calibration + inversion), `qec/repetition_code.py` (codes + memory experiments), `benchmark.py` (all studies in Section 4), `scripts/run_ibm_hardware.py` (hardware validation pipeline). Seeds, shot counts and scale factors are module constants in `benchmark.py`. Raw data: `out/results.json`.
## Figures

![SPSA convergence of H2 under three objectives](../out/fig_convergence.png)

**Figure 1** — SPSA convergence of H₂ under three objectives: noiseless (blue), raw noisy (red, stalls at 20.4 mHa), ZNE-mitigated in-loop (green, reaches 0.09 mHa). Green band: chemical accuracy.

![ZNE extrapolation curve for H2](../out/fig_zne_h2.png)

**Figure 2** — ZNE extrapolation curve for H₂: noisy energy vs scale factor λ (7 points) and the λ→0 extrapolant landing within the chemical-accuracy band.

![Recovery vs circuit depth](../out/fig_recovery_vs_depth.png)

**Figure 3** — Raw vs ZNE-corrected error as circuit depth grows (56/112/168 CX for H₂; LiH 1616 CX out-of-window at next-gen rates).

![Head-to-head estimator comparison](../out/fig_head_to_head.png)

**Figure 4** — Head-to-head bias vs spread for the four estimator arms at fixed circuit and shot budget (5 seeds).

![ZNE bias-variance frontier](../out/fig_zne_frontier.png)

**Figure 5** — Measured bias–variance frontier of ZNE configurations across two shot budgets, against the theoretical amplification Σwᵢ².

![QEC logical error vs distance](../out/fig_qec_distance.png)

**Figure 6** — Logical error rate vs repetition-code distance under Z-bias (p_z = 0.05): measurement with error bars vs leading-order theory (3p², 10p³, 35p⁴).

![QEC logical error vs rounds](../out/fig_qec_rounds.png)

**Figure 7** — Logical error growth with syndrome-extraction rounds (d = 3): cumulative and per-round rates.

![Gate overhead vs distance](../out/fig_gate_overhead.png)

**Figure 8** — Gate overhead per QEC cycle vs distance: exactly linear (6d − 6), the sub-quadratic KPI.
