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
