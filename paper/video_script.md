# 90-Second Video Script — IRIS Submission

**Total: ~230 words ≈ 90 seconds at natural pace. On-screen text suggestions in [brackets].**

---

**[0:00–0:12] HOOK** *(close-up of molecule animation)*

"Chemists dream of simulating molecules exactly — new drugs, better batteries, catalysts that pull carbon from air. Quantum computers could do it. But there's a catch: today's quantum processors are *noisy*."

**[0:12–0:30] PROBLEM** *(graph: clean blue curve vs corrupted red curve)*

"When my simulation of the hydrogen molecule runs on a calibrated model of real IBM hardware noise, the energy error explodes from zero to over 300 millihartree — 190 times worse than what chemistry actually needs. Quantum error *mitigation* promises to fix this in software — but which technique, and at what price?"

**[0:30–0:58] WHAT I DID** *(code scrolling, then frontier plot appears)*

"So I built a full open-source quantum simulation pipeline in Python and ran a controlled experiment: same circuit, same noise, same shot budget — comparing zero-noise extrapolation, readout mitigation, and their combination over many random seeds. The result is a measured *bias–variance frontier*: extrapolation kills 99.97% of the bias, but amplifies random noise up to 25 times — so at realistic shot budgets, no single method achieves chemical accuracy alone. I also mapped exactly where extrapolation *breaks* as circuits get deeper — and built error-correcting codes that suppress logical errors exactly as theory predicts, at strictly linear hardware cost."

**[0:58–1:15] RESULTS** *(bar chart: 345 → 9.8 mHa)*

"Combining techniques cut the final error from 345 to under 10 millihartree — and a corrected optimizer that previously stalled now converges to the exact molecular energy."

**[1:15–1:30] WHY IT MATTERS** *(IBM chip photo, then title card)*

"This gives experimentalists a practical rulebook for spending their measurement budget — the same decision facing industrial teams running quantum utility experiments today. Everything is open-source and reproducible in one command. I'm Sahibjot Singh — thank you for watching."

---

**Production notes:** keep every number on screen while spoken; show 2–3 s of the actual code and the reproducibility command (`python scripts/run_milestones.py`) for credibility; end card with repo URL and figure montage (Figures 1–8).
