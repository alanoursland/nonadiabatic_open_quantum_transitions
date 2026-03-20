# WP7 Progress Log

Reference plan: `dev/plan_wp7_260318.md`

---

## Task 1: Single-Spin Analytical Benchmark

**Status: Done**

**Files created:**
- `src/nmr_analytical.py` — analytical NMR relaxation rate functions
- `src/test_nmr_analytical.py` — 39 tests, all passing

**What it does:**
Computes ¹⁵N backbone amide relaxation rates from the Solomon-Bloembergen equations using a spectral density function J(ω). Serves as the "standard Redfield" benchmark for comparison with the numerical simulator and the nonadiabatic CGME method.

**Physical constants defined:**
- µ₀/(4π) = 1×10⁻⁷ T·m/A
- γ(¹⁵N) = −2.7116×10⁷ rad/s/T
- γ(¹H) = 2.6752×10⁸ rad/s/T
- r(NH) = 1.02 Å (default)
- Δσ(¹⁵N) = −170 ppm (default)

**Standalone functions:**
- `dipolar_coupling_constant(gamma_I, gamma_S, r)` → d = (µ₀/4π)|γI||γS|ℏ/r³ ≈ 72 krad/s for ¹⁵N-¹H
- `csa_coupling_constant(omega_0, delta_sigma)` → c = ω₀|Δσ|/√3 ≈ 31 krad/s at 11.7 T
- `compute_R1(d, c, J_fn, omega_N, omega_H)` → (d²/4)[J(ωH−ωN) + 3J(ωN) + 6J(ωH+ωN)] + c²J(ωN)
- `compute_R2(d, c, J_fn, omega_N, omega_H)` → (d²/8)[4J(0) + J(ωH−ωN) + 3J(ωN) + 6J(ωH) + 6J(ωH+ωN)] + (c²/6)[4J(0) + 3J(ωN)]
- `compute_NOE(d, J_fn, omega_N, omega_H, R1, gamma_N, gamma_H)` → 1 + (γH/γN)·σNH/R₁
- `compute_R1rho(R1, R2, omega1, delta_omega, Rex)` → R₁cos²θ + (R₂+Rex)sin²θ

**Convenience class `NMRRelaxation`:**
- Constructor takes B₀, J_fn, and optional physical parameters (γN, γH, rNH, Δσ)
- Computes ωN, ωH, d, c at construction
- Methods: `R1()`, `R2()`, `NOE()`, `R1rho(omega1_hz, offset_hz)`, `R1rho_dispersion(omega1_array_hz, offset_hz)`
- User-facing inputs in Hz; internal conversion to rad/s

**Design decisions:**
- Pure Python (`math` module only) — no torch/numpy dependency, matching `spectral_densities.py` style
- Reuses `LipariSzabo` from `spectral_densities.py` (J includes 2/5 orientational averaging factor)
- Imports `HBAR`, `KB` from `bath.py` (no constant redefinition)
- All angular frequencies in rad/s; NMRRelaxation accepts Hz at the user interface
- d always positive (absolute values of γ); signed γ used only in NOE formula

**Convention note:** The dipolar coupling constant d ≈ 72 krad/s differs from the plan's estimate of ~24 krad/s. The plan's value appears to have a numerical error. The correct value (verified against published R₁ ranges) is d = (µ₀/4π)|γN||γH|ℏ/rNH³ ≈ 72,100 rad/s. With this d, the Solomon equations produce R₁ ≈ 1.5–2.5 s⁻¹ for ubiquitin-like parameters at 500 MHz, matching published data.

**Tests cover (39 tests across 9 classes):**

*Coupling constants (5 tests):*
- d magnitude ~72 krad/s, positivity, 1/r³ scaling
- c magnitude ~31 krad/s, proportional to B₀

*Spectral density sanity (2 tests):*
- J(0) > J(ωN) > J(ωH) in motional narrowing regime
- J(ω) ≥ 0 for all ω

*Relaxation rate ranges (14 tests, including 9 parametrized):*
- R₁ ∈ [0.5, 5] s⁻¹, R₂ ∈ [3, 20] s⁻¹ for ubiquitin-like params
- R₂ > R₁ always
- Both positive for 9 combinations of τc and S²

*R₁ρ limiting cases (5 tests):*
- On resonance: R₁ρ = R₂ exactly
- Large offset: R₁ρ → R₁ (< 1% error at Δω/ω₁ = 10⁴)
- Monotonic decrease with offset
- Bounded between R₁ and R₂
- On-resonance dispersion flat for all ω₁

*NOE (2 tests):*
- Range check (−5 < NOE < 1.5 for ¹⁵N)
- NOE < 1 for fast tumbling (τc = 1 ns)

*Dipolar vs CSA (2 tests):*
- Dipolar dominates R₁ at low field (7 T)
- c scales linearly with B₀

*Published values (3 tests):*
- R₁ ∈ [1.0, 3.5] for ubiquitin at 500 MHz
- R₂ ∈ [4.0, 18.0]
- R₂/R₁ ratio ∈ [2, 10]

*NMRRelaxation class (4 tests):*
- Construction, Larmor frequencies, dispersion length
- Rigid tumbler (S²=1, τe=0) matches SimpleLorentzian exactly

*Phase 0 parameters (3 tests):*
- R₁, R₂ in expected ranges for τc=5ns, S²=0.85, B₀=11.7T
- On-resonance R₁ρ dispersion flat

**Full test suite: 183/183 passing (7.0s).**

---

## Task 2: Nonadiabatic Correction Estimate (Phase 0 Decision Gate)

**Status: Done — VERDICT: NEGATIVE**

**Files created:**
- `src/nonadiabatic_estimate.py` — analytical correction estimates and Phase 0 report
- `src/test_nonadiabatic_estimate.py` — 29 tests, all passing

**What it does:**
Estimates the magnitude of the nonadiabatic correction to NMR R₁ρ through two independent mechanisms, and compares the secular Lindblad and nonadiabatic CGME methods numerically on a toy 2-level system.

**Two correction mechanisms analyzed:**

1. **Perturbative adiabatic mixing:** |a₁|² = (ω₁/(2ωN))²
   - For NMR: ω₁ ~ 10² Hz → ~600 rad/s, ωN ~ 50 MHz → 3.2×10⁸ rad/s
   - |a₁|² ~ 10⁻¹³ for all spin-lock powers 25–1000 Hz
   - Maximum relative correction to R₁ρ: **6.3×10⁻¹¹**

2. **CGME vs secular sinc damping:** sinc((ωab − ωcd)·δτ/2)
   - For all NMR Bohr frequency pairs with δτ = 1 µs: |sinc| < 6.3×10⁻³
   - Non-secular Redfield elements are effectively zero → CGME agrees with secular Lindblad

**Phase 0 report output (B₀=11.7 T, τc=5 ns, S²=0.85, τe=50 ps):**
- R₁ = 2.449 s⁻¹ (T₁ = 0.408 s)
- R₂ = 6.789 s⁻¹ (T₂ = 0.147 s)
- NOE = 0.790
- Max relative correction: 6.27×10⁻¹¹ (threshold: 10⁻²)
- **Verdict: NEGATIVE** — NMR R₁ρ cannot distinguish the nonadiabatic framework

**Functions implemented:**
- `adiabatic_mixing_fraction(omega1, omega0)` → (ω₁/(2ω₀))²
- `cgme_sinc_retention(bohr_freq_diff, delta_tau)` → sinc(x) = sin(x)/x
- `nonadiabatic_correction_table(nmr, omega1_values_hz)` → sweep of corrections
- `cgme_sinc_table(nmr, delta_tau_s)` → sinc factors for all NMR Bohr frequency pairs
- `phase0_report(B0, tau_c, S2, tau_e, delta_tau_s)` → structured report dict
- `print_phase0_report(**kwargs)` → formatted output

**Tests cover (29 tests across 8 classes):**

*Adiabatic mixing fraction (6 tests):*
- ω₁² scaling, 1/ω₀² scaling, exact value
- NMR magnitude < 10⁻¹⁰
- Rejects zero/negative ω₀

*CGME sinc retention (6 tests):*
- Zero frequency diff → 1, small diff → ~1, large diff → ~0
- Exact zero at sinc(π)
- Bounded by 1, NMR ωN gives |sinc| < 0.01

*Correction table (4 tests):*
- Correct length, monotonic increase with ω₁
- All corrections < 10⁻⁸, correct keys

*Sinc table (4 tests):*
- Returns 6 frequency pairs, all |sinc| < 0.01
- Includes ωN and ωH labels

*Phase 0 report (6 tests):*
- Verdict is NEGATIVE, max correction < 10⁻⁸
- Sinc factors < 0.01, rates in physical ranges
- Report structure complete, threshold = 0.01

*Numerical simulator agreement (3 tests):*
- Secular Lindblad vs CGME populations agree (atol=10⁻⁸) with zero drive
- Coherences agree approximately (populations exact, off-diagonals within 10⁻²)
  - Genuine difference due to non-secular sinc terms at toy parameters where sinc ≈ O(1)
- Both methods decay from excited state and agree on equilibrium (atol=10⁻⁶)

**Phase 0 conclusion:** The nonadiabatic correction to NMR R₁ρ is ~10⁻¹¹, which is 9 orders of magnitude below any measurable threshold (~10⁻³). The correction is so small because the spin-lock frequency ω₁ ~ 10² Hz is negligible compared to the Larmor frequency ωN ~ 10⁸ rad/s. Recommend: pivot to HITRAN/FTMW experimental candidates where the perturbation/splitting ratio can be much larger.

**Full test suite: 212/212 passing (7.9s).**

---

## Task 3: Rotating-Frame Simulator Adapter — Not started

## Task 4: Full NMR R₁ρ Sweep — Not started

## Task 5: Validation and Comparison — Not started

## Task 6: Comparison Against Published Data — Not started
