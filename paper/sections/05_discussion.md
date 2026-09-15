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
