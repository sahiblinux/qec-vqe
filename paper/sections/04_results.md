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
