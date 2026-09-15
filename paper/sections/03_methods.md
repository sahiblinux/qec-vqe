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
