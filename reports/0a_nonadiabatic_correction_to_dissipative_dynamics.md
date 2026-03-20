# Phase 0 Feasibility Report: Experimental Detection of the Nonadiabatic Correction to Dissipative Dynamics

## Summary

We computed order-of-magnitude estimates for two candidate experiments — NMR rotating-frame relaxation (R₁ρ) and HITRAN pressure-broadened line mixing — to determine whether the nonadiabatic correction to the Redfield relaxation tensor is large enough to be detected with existing published data. The nonadiabatic correction arises from the Landau-Lifshitz decomposition of Dirac's time-dependent perturbation theory coefficients, which separates reversible adiabatic polarization from genuine nonadiabatic transitions. When a bath is present during driving, the standard (Dirac) and nonadiabatic frameworks couple the bath to different components of the density matrix, producing different relaxation rates.

The correction scales as (V/ΔE)², where V is the off-diagonal perturbation strength and ΔE is the relevant energy gap.

**NMR R₁ρ (Phase 0): NEGATIVE.** For a ¹⁵N backbone amide under a spin-lock field at 11.7 T, V/ΔE ~ 10⁻⁷, giving corrections of order 10⁻¹³. This is thirteen orders of magnitude below experimental precision.

**HITRAN line mixing (Phase 0A): POSITIVE for CO₂ at ≥5 atm.** For CO₂ (B = 0.39 cm⁻¹) with N₂ buffer gas at 5 atm, V/ΔE ~ 0.45 at J = 0, giving corrections of order 20%. This exceeds the 1% detection threshold by more than an order of magnitude.

The recommendation is to abandon NMR R₁ρ as a test system and proceed to a full comparison of the nonadiabatic coarse-grained master equation against published HITRAN CO₂ line-mixing data at elevated pressures.


---

## 1. Background

### 1.1 The Problem

The standard treatment of driven, dissipative quantum systems uses Dirac's transition coefficients c_k(t) to define state populations. However, c_k(t) conflates two physically distinct processes: the reversible polarization of the quantum state in response to the driving field (adiabatic), and genuine irreversible transitions between energy levels (nonadiabatic). When a thermal bath is coupled to the system during driving, this conflation causes the bath to relax the adiabatic polarization as if it were a real excitation, producing incorrect relaxation rates, fictitious heat flow, and the wrong thermal steady state.

The Landau-Lifshitz integration-by-parts decomposition separates c_k(t) = a_k(t) + b_k(t), where a_k is the adiabatic (polarization) term and b_k is the nonadiabatic (transition) term. Mandal, Hunt, and collaborators have shown that |b_k|² — not |c_k|² — gives physically correct transition probabilities that satisfy gauge invariance, the energy-work theorem, and thermodynamic consistency.

The theoretical case is established. What is missing is a comparison against independent experimental data using independently known parameters with no adjustable quantities.

### 1.2 Detection Criterion

The nonadiabatic correction modifies the Redfield relaxation tensor by changing which density matrix elements the bath acts on. The leading-order correction to any relaxation rate scales as:

|a_k|² = (V_k0 / ω_k0)²

where V_k0 = ⟨k|V|0⟩ is the perturbation matrix element and ω_k0 = (E_k - E_0)/ℏ is the Bohr frequency. This is the fraction of the state amplitude that is adiabatic polarization rather than genuine transition.

For an experiment to distinguish the two frameworks, the correction must exceed the experimental uncertainty. We adopt a threshold of 1%: if (V/ΔE)² < 0.01, the experiment cannot resolve the correction. If (V/ΔE)² ≥ 0.01, it can in principle, and detailed numerical comparison is warranted.

### 1.3 The Three Master Equations

Our simulator implements three dissipator variants derived from the same microscopic Redfield tensor R_{ab,cd}:

- **Full Redfield**: Uses R_{ab,cd} directly. Captures all non-secular coupling between populations and coherences. May violate positivity.
- **Secular Lindblad (GKSL)**: Zeros all R_{ab,cd} where |ω_ab - ω_cd| > 0. Guarantees positivity but discards coherence transfer between non-degenerate transitions.
- **Coarse-grained master equation (CGME)**: Damps each R_{ab,cd} by sinc((ω_ab - ω_cd)Δτ/2). Interpolates between full Redfield (Δτ → 0) and secular Lindblad (Δτ → ∞). Preserves positivity for appropriate Δτ.

The nonadiabatic correction enters through *which density matrix* the dissipator acts on (σ_nad vs σ_Dirac), while the secular/non-secular distinction enters through *which elements of the Redfield tensor* are retained. Both matter; the Phase 0 estimates address both.


---

## 2. Phase 0: NMR R₁ρ

### 2.1 System

A ¹⁵N backbone amide in a protein (standard parameters):

| Parameter | Value |
|---|---|
| Static field B₀ | 11.7 T (500 MHz ¹H) |
| ¹⁵N Larmor frequency ω_N | 3.17 × 10⁸ rad/s (50.5 MHz) |
| ¹H Larmor frequency ω_H | 3.13 × 10⁹ rad/s (497.8 MHz) |
| Dipolar coupling d | ~24,350 rad/s |
| CSA coupling c | ~10,200 rad/s |
| Correlation time τ_c | 5.0 ns |
| Order parameter S² | 0.85 |
| Internal correlation time τ_e | 50 ps |
| Temperature | 298 K |

The perturbation is a transverse spin-lock field of strength ω₁/(2π) = 25–1000 Hz applied to the ¹⁵N spin in the rotating frame. The bath is molecular tumbling, characterized by the Lipari-Szabo spectral density.

### 2.2 Correction 1: Adiabatic Mixing

The perturbation matrix element is V₁₀ = ω₁/2 (the Ix matrix element between Iz eigenstates). The unperturbed splitting is the ¹⁵N Larmor frequency ω_N. The adiabatic mixing fraction is:

|a₁|² = (ω₁ / 2ω_N)²

| ω₁/(2π) (Hz) | ω₁ (rad/s) | \|a₁\|² | δR₁ρ / R₁ρ |
|---|---|---|---|
| 25 | 157 | 6.2 × 10⁻¹⁴ | ~10⁻¹³ |
| 50 | 314 | 2.5 × 10⁻¹³ | ~10⁻¹³ |
| 100 | 628 | 9.9 × 10⁻¹³ | ~10⁻¹² |
| 250 | 1,571 | 6.2 × 10⁻¹² | ~10⁻¹¹ |
| 500 | 3,142 | 2.5 × 10⁻¹¹ | ~10⁻¹¹ |
| 1000 | 6,283 | 9.9 × 10⁻¹¹ | ~10⁻¹⁰ |

The correction is bounded by |a₁|² × |R₂ - R₁|, giving δR₁ρ/R₁ρ ~ 10⁻¹³ to 10⁻¹⁰. Experimental precision for R₁ρ is ~1–3%.

### 2.3 Correction 2: CGME vs Secular

The CGME retains non-secular Redfield elements damped by sinc((ω_ab - ω_cd)Δτ/2). For the ¹⁵N-¹H system, the non-degenerate Bohr frequency pairs are separated by ω_N, ω_H, ω_H ± ω_N, all of order 10⁸–10⁹ rad/s. With a coarse-graining timescale Δτ ~ 1 µs:

| Frequency pair | \|Δω\| (rad/s) | \|sinc\| |
|---|---|---|
| ω_N | 3.17 × 10⁸ | ~10⁻³ |
| ω_H | 3.13 × 10⁹ | ~10⁻⁴ |
| ω_H - ω_N | 2.81 × 10⁹ | ~10⁻⁴ |
| ω_H + ω_N | 3.44 × 10⁹ | ~10⁻⁴ |

All sinc factors are ≪ 1. The CGME is indistinguishable from the secular Lindblad for this system.

### 2.4 Verdict

**NEGATIVE.** Both correction mechanisms produce effects far below any achievable measurement precision. NMR R₁ρ cannot distinguish the nonadiabatic framework from the standard treatment. The fundamental reason is that the spin-lock field strength (10²–10³ rad/s) is negligible compared to the Larmor frequency splitting (10⁸ rad/s), making the adiabatic polarization immeasurably small.


---

## 3. Phase 0A: HITRAN Line Mixing

### 3.1 Motivation

Pressure-broadened line mixing in gas-phase molecular spectroscopy is a realization of the same driven-dissipative physics: the radiation field drives transitions between rotational states, and buffer-gas collisions provide the thermal bath. The standard theoretical treatment (Anderson-Tsao-Curnutte formalism) uses the Redfield relaxation matrix to predict off-diagonal line-coupling coefficients. This is precisely the Redfield tensor that the nonadiabatic framework modifies.

The key difference from NMR: the perturbation strength V (collisional coupling, proportional to pressure × broadening coefficient) can be made comparable to the energy gap ΔE (rotational level spacing, set by the rotational constant B) by using molecules with small B at moderate pressures. This places the system in a regime where the adiabatic mixing fraction is non-negligible.

### 3.2 System

A linear rigid rotor with energy levels E_J = B·J(J+1) in cm⁻¹. Adjacent levels (ΔJ = 1) are separated by ΔE = 2B(J+1) cm⁻¹. The collisional perturbation strength is V = γ_L × P, where γ_L is the pressure-broadening coefficient (cm⁻¹/atm) and P is the buffer gas pressure.

| Molecule | B (cm⁻¹) | γ_L (cm⁻¹/atm) | Smallest ΔE (cm⁻¹) |
|---|---|---|---|
| CO | 1.9313 | 0.065 | 3.86 |
| CO₂ | 0.3902 | 0.070 | 0.78 |
| N₂O | 0.4190 | 0.073 | 0.84 |
| HCl | 10.440 | 0.050 | 20.9 |

### 3.3 Correction 1: Adiabatic Mixing

The perturbation ratio V/ΔE at J = 0 (the most favorable case, smallest gap):

| Molecule | P (atm) | V (cm⁻¹) | V/ΔE | \|a\|² | Verdict |
|---|---|---|---|---|---|
| CO | 1 | 0.065 | 0.017 | 2.8 × 10⁻⁴ | NEGATIVE |
| CO | 10 | 0.650 | 0.168 | 2.8 × 10⁻² | POSITIVE |
| CO₂ | 1 | 0.070 | 0.090 | 8.0 × 10⁻³ | NEGATIVE (borderline) |
| CO₂ | 5 | 0.350 | 0.449 | 2.0 × 10⁻¹ | **POSITIVE** |
| CO₂ | 10 | 0.700 | 0.897 | 8.0 × 10⁻¹ | POSITIVE |
| N₂O | 5 | 0.365 | 0.436 | 1.9 × 10⁻¹ | POSITIVE |
| HCl | 1 | 0.050 | 0.002 | 5.7 × 10⁻⁶ | NEGATIVE |
| HCl | 50 | 2.500 | 0.120 | 1.4 × 10⁻² | POSITIVE |

CO₂ at 5 atm is the strongest candidate: a 20% correction at J = 0, with the first four J levels (J = 0 through J = 3) all above the 1% threshold.

The correction decreases with increasing J because the energy gap 2B(J+1) grows while V is J-independent (in the rigid-rotor approximation). This produces a J-dependent signature: the nonadiabatic correction is largest for low-J lines and falls off as ~1/(J+1)². This J-dependence is a distinctive prediction that could be tested line by line in published HITRAN data.

### 3.4 Correction 2: CGME Sinc Factors

For gas-phase collisions, the bath correlation time is τ_c ~ 1 ps (kinetic collision duration at room temperature). The CGME sinc argument for adjacent rotational lines separated by 2B is:

x = (2B × CM_TO_RAD_S) × τ_c / 2

For CO₂: x = (0.78 cm⁻¹)(1.88 × 10¹¹ rad/s/cm⁻¹)(10⁻¹² s)/2 = 0.073, giving sinc(x) ≈ 0.999.

| Molecule | Pair | Δω (cm⁻¹) | \|sinc\| |
|---|---|---|---|
| CO | Adjacent (2B) | 3.86 | 0.978 |
| CO | Next-nearest (4B) | 7.73 | 0.914 |
| CO₂ | Adjacent (2B) | 0.78 | 0.999 |
| CO₂ | Next-nearest (4B) | 1.56 | 0.996 |

All sinc factors are close to 1. This means the CGME retains essentially all non-secular Redfield elements for this system. The coarse-grained dissipator is nearly identical to the full Redfield dissipator, and both differ substantially from the secular Lindblad (which discards all non-secular coupling). This is the opposite of the NMR result, where all three methods collapsed to the same answer.

For HITRAN line mixing, the three-way comparison (Redfield / CGME / secular Lindblad) is genuine: secular Lindblad corresponds to ignoring line mixing entirely, full Redfield and CGME include it, and the nonadiabatic correction modifies the magnitude of the off-diagonal relaxation matrix elements by ~20% at low J.

### 3.5 Perturbative Validity

The estimate |a|² = (V/ΔE)² uses first-order perturbation theory. This is valid when |a|² ≪ 1. For CO₂ at 5 atm, |a|² ≈ 0.20 at J = 0 — large, but still below unity. At 10 atm (|a|² ≈ 0.80) the perturbative estimate is being stretched. At 50 atm (|a|² ≈ 20) it has broken down completely and the number is unphysical.

The useful operating regime for a perturbative comparison is CO₂ at 1–10 atm: large enough to detect (above the 1% threshold for P ≥ 5 atm), small enough that first-order perturbation theory remains valid (|a|² < 1). Higher pressures require non-perturbative treatment.

### 3.6 Verdict

**POSITIVE.** CO₂ at 5 atm gives a 20% nonadiabatic correction to the relaxation matrix at J = 0. N₂O at 5 atm gives a comparable 19% correction. CO requires ≥10 atm for a detectable signal. HCl requires impractically high pressures due to its large rotational constant.


---

## 4. Why NMR Fails and HITRAN Succeeds

The contrast between the two systems is entirely explained by the perturbation ratio V/ΔE:

| System | V | ΔE | V/ΔE | \|a\|² |
|---|---|---|---|---|
| NMR ¹⁵N spin-lock | ω₁ ~ 150 rad/s | ω_N ~ 3 × 10⁸ rad/s | 5 × 10⁻⁷ | 2.5 × 10⁻¹³ |
| CO₂ + N₂ at 5 atm | 0.35 cm⁻¹ | 0.78 cm⁻¹ | 0.45 | 0.20 |

The ratio differs by six orders of magnitude. NMR has an intrinsically unfavorable geometry for this test: the spin-lock field is deliberately weak (to probe slow dynamics), while the Larmor frequency splitting set by a multi-Tesla magnet is enormous. The perturbation-to-splitting ratio is built into the experiment to be negligible.

Molecular rotational spectroscopy inverts this: the collisional perturbation (proportional to gas pressure, an easily tunable knob) can be made comparable to the rotational level spacing (set by the molecular geometry, which for heavy triatomic molecules like CO₂ is intrinsically small). The experimenter controls V/ΔE directly.

A secondary factor is the bath timescale. In NMR, the bath correlation time (molecular tumbling, τ_c ~ ns) is much longer than the inverse Larmor frequency (~ns), making the secular approximation excellent. In gas-phase spectroscopy, the collision duration (~1 ps) is much shorter than the inverse rotational spacing (~10 ps for CO₂), so non-secular Redfield elements survive and the CGME/secular/Redfield distinction is physically meaningful.


---

## 5. Integration Test: Redfield Tensor with Rotational System

To verify that the existing codebase (WP1–WP6) can handle the HITRAN system class, we constructed a SpectralDensityBath using:

- A 6-level rigid rotor (CO, J_max = 5) with energy levels from the RotationalSystem class.
- A CollisionalLorentzian spectral density (γ₀ = 0.065 × CM_TO_RAD_S ≈ 1.2 × 10¹⁰ rad/s, τ_c = 1 ps).
- A nearest-neighbor (ΔJ = ±1) dimensionless coupling operator.

The Redfield tensor was built via bath.precompute() using the lab-frame eigenvalues and identity eigenvectors (the rigid rotor Hamiltonian is already diagonal in the |J⟩ basis).

All tests passed:

| Test | Result |
|---|---|
| Redfield dissipator trace preservation (Tr[D[ρ]] = 0) | PASS |
| Secular Lindblad dissipator trace preservation | PASS |
| CGME dissipator trace preservation (Δτ = 1 ps) | PASS |
| Redfield dissipator Hermiticity (D[ρ] = D[ρ]†) | PASS |
| Population decay from J = 5 (excited state loses population) | PASS |
| Population gain at J = 4 (adjacent state gains population) | PASS |

The existing bath machinery works without modification for the rotational system. No new physics code is needed for the dissipator — only the system class (RotationalSystem) and spectral density model (CollisionalLorentzian) are new.


---

## 6. Recommendations

### 6.1 Immediate Next Step

Proceed to a full numerical comparison of standard Redfield vs. nonadiabatic CGME predictions for CO₂ R-branch line-mixing parameters against published HITRAN data at 5–10 atm. The target observable is the off-diagonal relaxation matrix element W_{J,J'} that governs line coupling between adjacent rotational transitions. HITRAN tabulates first-order line-mixing coefficients (the Y-parameters) derived from these matrix elements. The nonadiabatic framework predicts that W_{J,J'} is reduced by a J-dependent factor of order (V/2B(J+1))² relative to the standard Redfield value.

### 6.2 Specific Data Target

CO₂ in N₂ buffer gas at 296 K, 1–10 atm, R-branch transitions from J = 0 to J = 20. Published data sources:

- HITRAN2020 database (Gordon et al., JQSRT 2022): line-mixing parameters for CO₂ ν₃ band.
- Hartmann et al., *Collisional Effects on Molecular Spectra* (Elsevier, 2021): comprehensive tables of relaxation matrix elements.
- Gamache et al.: computed broadening and line-mixing parameters from ab initio potentials.

Known discrepancies between standard Redfield predictions and measured HITRAN line shapes at elevated pressures have been documented, particularly for CO₂ Q-branches. The nonadiabatic correction provides a new candidate explanation for these discrepancies.

### 6.3 Theoretical Work Required

The connection between the nonadiabatic framework (formulated in the time domain for pulsed excitation) and the line-shape formalism (frequency domain, steady-state absorption) must be established explicitly. They are related by Fourier transform, but the nonadiabatic decomposition has not previously been developed in the frequency domain. This is a theoretical prerequisite for the HITRAN comparison, not merely an implementation task.

### 6.4 What Has Been Ruled Out

NMR R₁ρ is definitively excluded as a test system for the nonadiabatic framework. The correction is thirteen orders of magnitude too small. No change in spin-lock power, magnetic field strength, or molecular system can overcome the fundamental (ω₁/ω_N)² scaling. The NMR codebase (WP1–WP6) remains useful as a validated simulator for driven-dissipative dynamics, but the NMR experiment plan should not be pursued further.

HCl is excluded at practical pressures (≤10 atm) due to its large rotational constant. CO requires elevated pressures (≥10 atm) and is a secondary target. N₂O is comparable to CO₂ and is a useful cross-check.


---

## 7. Codebase Status

| Component | Status | Tests |
|---|---|---|
| SpinSystem (WP1) | Complete | 9 passing |
| SpinLockDrive (WP2) | Complete | 13 passing |
| SpectralDensityBath (WP3) | Complete | ~61 passing |
| NonadiabaticDecomposition (WP4) | Complete | 11 passing |
| Simulator (WP5) | Complete | 25 passing |
| R₁ρ Extraction (WP6) | Complete | 13 passing |
| NMR analytical rates | Complete | — |
| NMR Phase 0 estimate | Complete (NEGATIVE) | — |
| RotationalSystem | Complete | 17 passing |
| CollisionalLorentzian | Complete | 6 passing |
| HITRAN Phase 0A estimate | Complete (POSITIVE) | 27 passing |
| **Total** | | **~182 tests passing** |

All code is in PyTorch with float64/complex128 precision. The simulator uses RK4 with automatic sub-stepping. The Redfield tensor is built element-wise (O(dim⁴)) and stored as a superoperator for fast dissipator evaluation.


---

## References

- Mandal, A. & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
- Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields: Dephasing and population relaxation due to contact with a bath. *J. Chem. Phys.*, 158, 164107.
- Schaller, G. & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Phys. Rev. A*, 78, 022106.
- Hartmann, J.-M., Boulet, C., & Robert, D. (2021). *Collisional Effects on Molecular Spectra* (2nd ed.). Elsevier.
- Gordon, I. E. et al. (2022). The HITRAN2020 molecular spectroscopic database. *JQSRT*, 277, 107949.
- Lipari, G. & Szabo, A. (1982). Model-free approach to the interpretation of nuclear magnetic resonance relaxation in macromolecules. *JACS*, 104, 4546–4559.
- Korzhnev, D. M., Orekhov, V. Y., & Kay, L. E. (2005). Off-resonance R₁ρ NMR studies of exchange dynamics in proteins with low spin-lock fields. *JACS*, 127, 713–721.
- Massi, F., Johnson, E., Wang, C., Rance, M., & Palmer, A. G. (2004). NMR R₁ρ rotating-frame relaxation with weak radio frequency fields. *JACS*, 126, 2247–2256.