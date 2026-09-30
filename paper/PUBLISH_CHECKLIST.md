# Publishing checklist — sahiblinux/qec-vqe

**Done for you (2026-09-30):**
- Code pushed → https://github.com/sahiblinux/qec-vqe
- Release tagged → https://github.com/sahiblinux/qec-vqe/releases/tag/v0.1.0 (this is what Zenodo mints the software DOI from)
- Note: `.github/workflows/ci.yml` was removed from the push because the token saved on this Mac lacks the `workflow` scope (backup: `.freebuff/ci.yml.bak`). The CI badge in the README will show "no status" until the file is restored. Fix later in 2 minutes: open https://github.com/sahiblinux/qec-vqe/new/main?filename=.github/workflows/ci.yml → paste the contents of `.freebuff/ci.yml.bak` → commit (web edits don't need the workflow scope). Then run `git pull` locally.

---

## Step 1 — Zenodo: mint the research-paper DOI ✅ DONE 2026-09-30

**DOI: https://doi.org/10.5281/zenodo.23061923**

1. Open https://zenodo.org → **Sign up** (or **Log in with GitHub** — easiest).
2. Top-right **Upload** → **New upload**.
3. Paste the abstract (copy from Step 1b below).
4. **Files → Upload**: drag in
   - `paper/research_paper.pdf` (the 15-page paper)
   - `out/results.json` (raw data backing every number)
5. **License**: `Creative Commons Attribution 4.0 International` (cc-by-4.0).
6. **Resource type**: `Text / Preprint`.
7. Scroll to **Related identifiers**: *Identifier* `https://github.com/sahiblinux/qec-vqe` · *Relationship* `is supplemented by`.
8. **Publish**. Copy the DOI it gives you (looks like `10.5281/zenodo.NNNNNNN`).

### 1b — Abstract to paste (verbatim)

Variational Quantum Eigensolvers (VQE) promise molecular ground-state energies on today's noisy intermediate-scale quantum (NISQ) processors, but device noise biases the estimated energy severely — under a device-calibrated noise model, the hydrogen molecule's energy error grows from 0.000 mHa (noiseless) to 303.6 mHa (noisy), roughly 190× the "chemical accuracy" threshold of 1.6 mHa that computational chemistry requires. Error-mitigation techniques repair such biases in software, but published results usually report each technique in isolation, making it impossible to know which mitigation to buy with a fixed shot budget. This work answers that question with a controlled comparison. Under identical circuits, noise profiles and shot budgets, I measure the bias and variance of four estimator stacks — raw sampling, twirled-readout extinction (TREX-style calibrated readout inversion), zero-noise extrapolation (ZNE), and ZNE+TREX — across independent random seeds. Three findings emerge. (1) ZNE is a bias–variance tradeoff, not a free lunch: a 7-point Richardson extrapolant reduces the deterministic bias 3200-fold (303.6 → 0.094 mHa) but amplifies statistical noise by a measured factor of up to 25× per run, while TREX removes readout bias with no variance penalty. (2) ZNE has a finite operating window in circuit depth: error recovery falls 99.97% → 99.1% → 95.5% as CX count grows 56 → 112 → 168, and is exhausted entirely for a 1616-CX lithium-hydride circuit at near-term error rates. (3) A noise-bias-tailored repetition code corrects the dominant error axis with logical error rates matching theory (3p², 10p³, 35p⁴ for distances 3, 5, 7) at strictly linear gate overhead (gates = 6d − 6). The results reframe mitigation selection as budget allocation on an empirically measured bias–variance frontier and provide open-source tooling to reproduce every number.

### 1c — Also: get the GitHub↔Zenodo auto-link (optional, 2 min)

1. Zenodo top-right menu → **GitHub** → **Connect** / Log in with GitHub → grant access.
2. On the GitHub side accept, then at https://github.com/sahiblinux/qec-vqe/releases edit v0.1.0 nothing more is needed — the existing release syncs and Zenodo mints a **software version DOI** automatically.

## Step 2 — JOSS submission (~5 min)

1. Open https://joss.theoj.org → **Submit / Login via GitHub** (your sahiblinux account).
2. Start a new submission; paste the manuscript (copy from `paper/paper.md`, everything including the `---` YAML header).
3. Attach `paper/paper.bib` when the form asks for the bibliography.
4. Paper URL / repo field: `https://github.com/sahiblinux/qec-vqe`.
5. Archive DOI field: paste the Zenodo **software** DOI from Step 1c (or the one from Step 1 if you didn't do 1c).
6. Submit → an editor will assign an eic and reviewers open issues on the repo. When someone comments you'll get an email; just come back here and I'll draft every reply.

## Step 3 — Wire DOI into repo ✅ DONE 2026-09-30

- `CITATION.cff` now carries DOI `10.5281/zenodo.23061923`
- `README.md` has a Zenodo DOI badge + full research-paper citation block

## Login status on this machine (checked)

- **GitHub**: already logged in via saved token — push and release both succeeded without you touching anything. Repo: https://github.com/sahiblinux/qec-vqe
- **Zenodo**: no login found → do Step 1 in your browser.
- **JOSS**: no account yet → use "Log in via GitHub" in Step 2, no new password needed.
