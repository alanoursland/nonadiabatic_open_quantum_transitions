# WP7 Implementation Plan: NMR R₁ρ Main Script

## Goal

Replace the transmon demo `main.py` with a script that sets up a ¹⁵N-¹H spin pair with physical NMR parameters, runs all three master equation methods (Redfield, secular Lindblad, nonadiabatic CGME) across a sweep of spin-lock powers, extracts R₁ρ dispersion curves, and produces a quantitative comparison. This is the payoff step: everything in WP1–WP6 exists to answer one question — do the standard and nonadiabatic frameworks predict measurably different R₁ρ?

The script implements Phase 0 (feasibility estimate) and Phase 3 (numerical comparison) of the experiment plan, using the codebase built in WP1–WP6.

---

## Problem 1: Unit Consistency

This is the single hardest part of WP7, and it must be solved before anything else runs.

The codebase currently has three different frequency scales that need to be reconciled:

- **SpinSystem** computes Larmor frequencies in **rad/s** (γ·B₀ ~ 10⁸ rad/s for ¹⁵N at 11.7 T).
- **SpinLockDrive** works in **rad/µs** internally (converts from Hz at construction).
- **SpectralDensityBath** evaluates J(ω) at the system's Bohr frequencies and computes γ(ω) using ℏ and k_B in SI. The Redfield tensor R has units of inverse time in whatever time unit the Bohr frequencies use.

The simulator propagates `drho/dt = -i[H,ρ] + D[ρ]`. For this equation to be dimensionally consistent, H and D must use the same time unit. Right now:

- H₀ from `SpinSystem.get_H0()` is in rad/s.
- V(t) from `SpinLockDrive.__call__(t)` is in rad/µs (because t is in µs).
- The Redfield tensor from `SpectralDensityBath` is in rad/s (because it's built from Bohr frequencies in rad/s and γ(ω) in s⁻¹).
- The RK4 integrator steps in whatever units `times` uses.

**These are incompatible.** The Hamiltonian has two terms in different time units, and the dissipator is in yet another.

### Resolution: Work in µs throughout the propagation

The choice of µs as the propagation time unit is forced by the `SpinLockDrive` (its envelope and derivative are defined in µs). The conversions needed:

1. **H₀**: Convert from rad/s to rad/µs by multiplying by 10⁻⁶. This means the Bohr frequencies passed to `bath.precompute()` must also be in rad/µs. Alternatively, convert the eigenvalues after `SpinSystem` computes them.

2. **Spectral density and bath rates**: J(ω) expects ω in rad/s (the NMR spectral density models are parameterized in seconds). But the Redfield tensor must be in µs⁻¹ for consistency with the propagator. Two options:
   - (a) Build the Redfield tensor at the physical (rad/s) Bohr frequencies, then scale the entire tensor by 10⁻⁶ to convert from s⁻¹ to µs⁻¹.
   - (b) Pass Bohr frequencies in rad/s to `precompute()` (so J(ω) is called at the right frequencies), but divide the resulting R-tensor by 10⁶ before use.

   Option (a) is cleaner. The bath's `_gamma(ω)` returns rates in s⁻¹ when J is in seconds. The full R-tensor has units of s⁻¹. Multiplying by 10⁻⁶ converts to µs⁻¹.

3. **Drive V(t)**: Already in rad/µs. No conversion needed.

4. **Output times**: In µs. The R₁ρ fit returns rates in µs⁻¹, which must be converted to s⁻¹ (multiply by 10⁶) for comparison with published data.

### Implementation

Create a helper function that:
- Takes the SpinSystem's H₀ eigenvalues (rad/s) and eigenvectors.
- Calls `bath.precompute(energies_rad_s, eigvecs)` with the original rad/s eigenvalues (so J(ω) is evaluated at correct physical frequencies).
- Stores a unit conversion factor `R_SCALE = 1e-6` to apply to the Redfield tensor.

**This requires a small modification to the simulator or a wrapper.** The cleanest approach: after `bath.precompute()`, scale `bath.R_full`, `bath.R_secular`, `bath.R_full_super`, and `bath.R_secular_super` by 1e-6 in-place. The `_bohr_diff` matrix (used for the sinc in CGME) should NOT be scaled — it's a frequency difference, not a rate. And the H₀ diagonal used in the eigenbasis propagation must be in rad/µs.

Concretely, the main script will:
```
# After bath.precompute(energies_rad_s, eigvecs):
bath.R_full *= 1e-6
bath.R_secular *= 1e-6
bath.R_full_super *= 1e-6
bath.R_secular_super *= 1e-6
# And override bath.energies with the µs-scaled values:
bath.energies = energies_rad_s * 1e-6
```

This is admittedly a hack — the right long-term fix is a unit system in the bath class. But for WP7's purposes, scaling the precomputed tensors is correct, transparent, and testable.

**Verification**: Run with zero drive. The secular Lindblad dissipator should produce population decay with T₁ values consistent with the known NMR R₁ for the ¹⁵N-¹H system (~1–2 s⁻¹ at 11.7 T with τ_c = 5 ns). If the rates come out as 10⁶ × too large or too small, the scaling is wrong.

---

## Problem 2: Coupling Operators

The `SpectralDensityBath` takes a list of coupling operators S_α that represent how the bath couples to the system. For an amide ¹⁵N-¹H pair, there are two dominant relaxation mechanisms:

### 2a: ¹⁵N-¹H Dipolar Coupling

The dipolar coupling between the ¹⁵N and ¹H spins, modulated by molecular tumbling, is the primary source of ¹⁵N relaxation. In the high-field secular limit (valid for heteronuclear pairs where ω_N ≠ ω_H), the dipolar relaxation is driven by terms of the form:

```
S_dipolar = d_NH * Iz_N
```

where `d_NH` is the dipolar coupling constant:
```
d_NH = -(µ₀/4π) * γ_N * γ_H * ℏ / r_NH³
```

with r_NH = 1.02 Å (the N-H bond length). This gives d_NH ≈ -24,350 rad/s (the sign is from γ_N < 0).

However, the full dipolar relaxation involves spectral density sampled at multiple frequencies: J(0), J(ω_N), J(ω_H), J(ω_H - ω_N), J(ω_H + ω_N). In the Redfield tensor formalism as implemented, this is automatically captured if the coupling operator has the correct matrix elements between eigenstates. The issue is that the dipolar Hamiltonian has both secular and non-secular parts, and the bath's Redfield tensor construction already handles this by evaluating γ(ω) at all Bohr frequencies.

For the 4×4 two-spin system, the dipolar coupling operator in the product basis is:

```
H_DD ∝ 2*Iz_N⊗Iz_H - Ix_N⊗Ix_H - Iy_N⊗Iy_H
     = 2*Iz_N⊗Iz_H - 0.5*(I+_N⊗I-_H + I-_N⊗I+_H)
```

The full operator (with the prefactor) is:
```
S_DD = d_NH * (2*Iz_N*Iz_H - 0.5*(Ip_N*Im_H + Im_N*Ip_H))
```

where d_NH includes the geometric and distance factors. The (2/5) factor in the Lipari-Szabo spectral density already accounts for the orientational averaging of the rank-2 spherical harmonics.

### 2b: ¹⁵N Chemical Shift Anisotropy (CSA)

The CSA mechanism arises from the anisotropic chemical shielding tensor of ¹⁵N. In the secular approximation, the CSA coupling operator is:

```
S_CSA = δ_CSA * ω_N * Iz_N
```

where δ_CSA is the CSA anisotropy (typically -170 ppm for backbone amide ¹⁵N) and ω_N is the ¹⁵N Larmor frequency. The factor ω_N·δ_CSA gives the CSA interaction strength in rad/s.

More precisely:
```
Δσ = σ_∥ - σ_⊥ ≈ -170 ppm
S_CSA = (2/3)^(1/2) * γ_N * B₀ * Δσ * 1e-6 * Iz_N
```

The (2/5) in the spectral density again handles the orientational averaging.

### 2c: Cross-Correlated Relaxation

The dipolar and CSA mechanisms are correlated (the same molecular tumbling modulates both), producing cross-correlated relaxation. For the Phase 0 feasibility check, this can be neglected — it affects the absolute R₁ρ values but not the relative difference between Dirac and nonadiabatic predictions, which is what we're testing. Add a note to revisit if the Phase 0 difference is borderline.

### Implementation

Build the coupling operators from SpinSystem's operator accessors:

```python
Ix_N = system.get_spin_operator(0, 'x')  # ¹⁵N is spin 0
Iy_N = system.get_spin_operator(0, 'y')
Iz_N = system.get_spin_operator(0, 'z')
Iz_H = system.get_spin_operator(1, 'z')
Ix_H = system.get_spin_operator(1, 'x')
Iy_H = system.get_spin_operator(1, 'y')

# Dipolar coupling constant
MU0_4PI = 1e-7  # T·m/A
HBAR = 1.054571817e-34
gamma_N = -27.116e6  # rad/s/T
gamma_H = 267.522e6
r_NH = 1.02e-10  # meters
d_NH = -MU0_4PI * gamma_N * gamma_H * HBAR / r_NH**3

S_DD = d_NH * (2*Iz_N @ Iz_H - 0.5*(Ix_N @ Ix_H + Iy_N @ Iy_H)
               - 0.5*(-1j*Iy_N @ (1j*Iy_H)))  
# Simplify: the I+I- terms. Actually:
# Ix⊗Ix + Iy⊗Iy = 0.5*(I+⊗I- + I-⊗I+)
# So: S_DD = d_NH * (2*IzN*IzH - IxN*IxH - IyN*IyH)

S_DD = d_NH * (2 * (Iz_N @ Iz_H) - (Ix_N @ Ix_H) - (Iy_N @ Iy_H))

# CSA coupling
delta_sigma = -170e-6  # convert ppm to dimensionless
omega_N = abs(gamma_N * B0)  # Larmor frequency magnitude in rad/s
S_CSA = (2/3)**0.5 * omega_N * delta_sigma * Iz_N

coupling_operators = [S_DD, S_CSA]
```

**Note on Iz operators**: In the product basis, `Iz_N @ Iz_H` means Iz_N ⊗ Iz_H. But `SpinSystem.get_spin_operator(0, 'z')` already returns the full 4×4 operator (Iz_N ⊗ I_H), so the product `Iz_N @ Iz_H` in the code is actually (Iz_N ⊗ I_H)(I_N ⊗ Iz_H) = Iz_N ⊗ Iz_H. This is correct via the matrix multiplication of the Kronecker product representations. Verify by checking matrix elements.

---

## Problem 3: Rotating Frame vs Lab Frame

The SpinSystem constructs H₀ in the lab frame (Zeeman + J-coupling). But the spin-lock experiment operates in the rotating frame at (or near) the ¹⁵N Larmor frequency. In the rotating frame:

- The ¹⁵N Zeeman term disappears (or is replaced by the offset Δω·Iz_N).
- The ¹H Zeeman term is unchanged (the frame rotates at the ¹⁵N frequency, not the ¹H frequency).
- The J-coupling is unchanged (it's a scalar interaction, frame-independent).

For the on-resonance case (Δω = 0), the rotating-frame Hamiltonian is:
```
H₀_rot = ω_H * Iz_H + 2π * J_NH * (I_N · I_H)
```

where ω_H is the ¹H Larmor frequency. The ¹⁵N Zeeman term is gone.

**This is important.** The Bohr frequencies of H₀_rot determine the Redfield tensor — they control at which frequencies the spectral density is sampled. In the lab frame, the Bohr frequencies include ω_N ~ 10⁸ rad/s. In the rotating frame, the relevant frequencies are ω_H (~10⁹ rad/s), J_NH (~600 rad/s), and small differences. The relaxation physics depends on both (J is evaluated at ω_N, ω_H, etc.), but the propagation frame affects which terms are secular and which are not.

### Resolution

The standard NMR relaxation theory (Solomon equations, Redfield for dipolar/CSA) derives relaxation rates from the spectral density evaluated at the **lab-frame** Larmor frequencies: J(0), J(ω_N), J(ω_H), J(ω_H ± ω_N). These are the frequencies at which the dipolar interaction fluctuates due to molecular tumbling. The choice of rotating frame affects the coherent evolution but not the bath coupling.

The cleanest approach:

1. **Construct H₀ in the rotating frame** (remove the ¹⁵N Zeeman term, keep ¹H Zeeman and J-coupling). This is the Hamiltonian whose eigenstates define the propagation basis.

2. **Evaluate the spectral density at the lab-frame Bohr frequencies** for constructing the Redfield tensor. This means the bath's `precompute()` should receive eigenvalues that reflect the lab-frame energy splittings, not the rotating-frame ones.

This creates a tension: the Redfield tensor must be built at lab-frame frequencies, but the propagation (and the eigenbasis) uses rotating-frame eigenvalues. The resolution is that the Redfield tensor elements R_{abcd} involve γ(ω_ca) where ω_ca = E_c - E_a. If we build R using **lab-frame** eigenvalues, the rates are correct. The secular approximation (which terms to keep) should use the **rotating-frame** Bohr frequencies (because those determine which elements of ρ oscillate slowly enough to survive coarse-graining).

**Practical implementation**: This is actually how the standard NMR R₁/R₂ calculation works. The approach for WP7:

- Do NOT use SpinSystem.get_H0() directly for precompute. Instead, build a modified Hamiltonian or pass custom eigenvalues.
- Better: since the 4×4 system is small, compute the relaxation rates semi-analytically using the standard Solomon-Bloembergen-Morgan expressions as a cross-check, and use the full numerical Redfield machinery for the actual comparison.

**Simplification for Phase 0**: For the feasibility estimate, use a **single-spin** model first. A single ¹⁵N spin in the rotating frame has H₀ = Δω·Iz (2×2). The Redfield tensor can be constructed from J(ω_N), J(0), etc., using the known NMR relaxation expressions. The coupling operators reduce to Iz (for both dipolar and CSA, in the secular high-field limit). This sidesteps the rotating-frame subtlety entirely and gives a direct estimate of whether the nonadiabatic correction exceeds 1%.

If the single-spin feasibility is positive, then implement the full two-spin rotating-frame version.

---

## Problem 4: Time Grid and Integration Stability

At physical NMR parameters, the frequency scales in the problem are:

- ¹H Larmor: ω_H ≈ 3.13 × 10⁹ rad/s = 3.13 × 10³ rad/µs
- ¹⁵N Larmor (in rotating frame): Δω ≈ 0 (on resonance) or up to ~10⁴ rad/s = 0.01 rad/µs
- J-coupling: 2π × 92 Hz ≈ 578 rad/s = 5.78 × 10⁻⁴ rad/µs
- Spin-lock power: ω₁ = 2π × 25–1000 Hz ≈ 157–6283 rad/s = 1.57×10⁻⁴ – 6.28×10⁻³ rad/µs
- Relaxation rates: R₁, R₂ ~ 1–10 s⁻¹ = 10⁻⁶–10⁻⁵ µs⁻¹

If we work in the ¹⁵N rotating frame and the system is on-resonance, the fastest frequency in the Hamiltonian is ω_H ≈ 3000 rad/µs. The RK4 stability limit is dt < 1.5/3000 ≈ 0.0005 µs. A 50 µs plateau requires 100,000 substeps. This is feasible but slow.

**The ¹H frequency is the problem.** In the ¹⁵N rotating frame, the ¹H Zeeman term still oscillates at ~500 MHz. But the proton spin is not being driven or observed — it's just there as the other half of the dipolar pair. Its fast precession does not affect the ¹⁵N relaxation physics (the relevant spectral density frequencies are fixed).

### Resolution: Truncate to the ¹⁵N Subspace

For a single-spin Phase 0 calculation: work with a 2×2 system (just ¹⁵N). The ¹H is part of the bath. The coupling to ¹H enters through the spectral density (dipolar relaxation rates), not through the system Hamiltonian. This is the standard NMR relaxation approach: the "system" is the observed spin, and the "bath" includes both the lattice and the heteronuclear partner.

For the two-spin version: implement the propagation in the **doubly rotating frame** (rotating at ω_N for ¹⁵N and ω_H for ¹H). In this frame, both Zeeman terms vanish (replaced by offsets), and the fastest frequency is J_NH ≈ 600 rad/s = 6×10⁻⁴ rad/µs. The RK4 stability limit becomes dt < 1.5/(6×10⁻⁴) ≈ 2500 µs — essentially unlimited. The propagation would then be trivially fast.

The doubly rotating frame requires transforming the coupling operators accordingly. The J-coupling Hamiltonian is unchanged (scalar), but the dipolar coupling operator picks up oscillating factors at ω_H - ω_N and ω_H + ω_N — these become the non-secular dipolar terms that contribute to relaxation but not to coherent evolution in the doubly rotating frame. Since we handle relaxation through the Redfield tensor (evaluated at lab-frame frequencies), the coherent Hamiltonian in the doubly rotating frame is just:
```
H₀_drf = Δω_N * Iz_N + Δω_H * Iz_H + 2π * J_NH * Iz_N * Iz_H
```
(keeping only the secular part of J for heteronuclear pairs, which is the Iz·Iz term).

---

## Staged Implementation Plan

### Stage 1: Single-Spin Feasibility Check (Phase 0)

**Goal**: Estimate |ΔR₁ρ|/R₁ρ for a single ¹⁵N spin. If < 1%, stop and pivot.

**System**: 2×2, H₀ = Δω·Iz in the ¹⁵N rotating frame.

**Bath**: Lipari-Szabo spectral density with τ_c = 5 ns, S² = 0.85, τ_e = 50 ps. Coupling operator: `S = d_eff * Iz` where d_eff encodes the effective dipolar + CSA coupling strength. The combined R₁ and R₂ rates for ¹⁵N are:
```
R₁ = (d²/4)[J(ωH-ωN) + 3J(ωN) + 6J(ωH+ωN)] + c²J(ωN)
R₂ = (d²/8)[4J(0) + J(ωH-ωN) + 3J(ωN) + 6J(ωH) + 6J(ωH+ωN)] 
    + (c²/6)[4J(0) + 3J(ωN)]
```
where d = dipolar coupling constant, c = CSA coupling constant.

However, for the Redfield tensor approach with a single coupling operator, we can't reproduce the full multi-frequency structure of NMR relaxation (which involves J evaluated at 5 different frequencies). The single-operator Redfield tensor evaluates γ at the system's Bohr frequencies only.

**Key insight**: For a 2×2 system with eigenvalues {0, Δω}, the Bohr frequencies are {0, ±Δω}. If Δω is the spin-lock offset (small, possibly zero for on-resonance), then the Redfield tensor samples J only at J(0) and J(Δω) — missing J(ω_N), J(ω_H), etc. This is wrong for NMR relaxation.

**This means the simple Redfield tensor construction doesn't work for NMR.** The NMR relaxation rates involve spectral densities at lab-frame frequencies (hundreds of MHz), but the rotating-frame Bohr frequencies are near-zero (Hz to kHz). The mismatch is fundamental to the rotating-frame treatment.

### Resolution: Pre-compute NMR rates, inject as effective Redfield tensor

The standard NMR approach is to compute R₁, R₂ (and cross-correlated rates) analytically from the Solomon equations and inject them directly. For the Redfield tensor of the 2×2 system, this means:

```
R_{00,00} = -R₁ * p_eq,1    (ground state gains population from excited)
R_{11,11} = -R₁ * p_eq,0    (excited state loses population)
R_{01,01} = -R₂             (coherence decays)
R_{10,10} = -R₂
etc.
```

But this defeats the purpose — it makes all three methods use the same rates, which is exactly what WP3 was designed to avoid.

**The real resolution**: The Redfield tensor must be built in the **lab frame** (where the Bohr frequencies are ω_N, ω_H, etc.) and then **transformed** to the rotating frame. The lab-frame Redfield tensor has the correct frequency dependence. The transformation to the rotating frame changes which elements are secular (oscillate slowly) but preserves the rates.

For the 4×4 two-spin system in the lab frame, the eigenvalues are dominated by the Zeeman interaction: {½(ω_N+ω_H), ½(ω_N-ω_H), ½(-ω_N+ω_H), ½(-ω_N-ω_H)} plus small J-coupling corrections. The Bohr frequencies include ω_N, ω_H, ω_H±ω_N, and 0 — exactly the frequencies where J should be sampled.

So the procedure is:

1. Build the **lab-frame** Hamiltonian H₀_lab from SpinSystem. Diagonalize it.
2. Call `bath.precompute(energies_lab, eigvecs_lab)` — this builds R_{abcd} with γ evaluated at lab-frame Bohr frequencies. Correct.
3. Transform to the **rotating frame** for propagation. The rotating-frame transformation is a unitary U_rot(t) = exp(-iω_N·Iz_N·t - iω_H·Iz_H·t). In the interaction picture, the density matrix transforms as ρ_rot = U†·ρ·U. The Redfield tensor in the rotating frame picks up oscillating phase factors exp(i·Δω·t) on the non-secular elements, which are the same factors that the coarse-graining (sinc damping) or secular approximation zeros out.

**Key realization**: the Redfield tensor `R_{ab,cd}` as built by `bath.precompute()` is already the correct lab-frame tensor. When the simulator propagates in the eigenbasis of H₀, it's effectively in the lab frame's energy eigenbasis. The secular/non-secular distinction is automatic: secular elements have ω_ab = ω_cd, non-secular have ω_ab ≠ ω_cd.

The issue is not the Redfield tensor — it's the **propagation frame**. The simulator currently propagates in the eigenbasis of whatever H₀ it's given. If we give it the lab-frame H₀ (with huge Zeeman splittings), the propagation has fast oscillations (10⁸ rad/s) requiring tiny time steps. If we give it the rotating-frame H₀, the Bohr frequencies are wrong for the Redfield tensor.

### The Correct Architecture

Use the lab-frame eigenbasis for the **Redfield tensor** but propagate the density matrix in the **interaction picture** where the fast Zeeman evolution has been removed.

Concretely:

1. Diagonalize H₀_lab → eigenvalues E_a, eigenvectors U.
2. Build Redfield tensor R_{ab,cd} using these eigenvalues (bath.precompute). This is correct.
3. In the interaction picture, ρ̃_{ab}(t) = ρ_{ab}(t) · exp(iω_{ab}t). The master equation becomes:
   ```
   dρ̃_{ab}/dt = sum_{cd} R_{ab,cd} · exp(i(ω_{ab}-ω_{cd})t) · ρ̃_{cd} 
                + drive terms in interaction picture
   ```
4. The secular approximation: keep only terms where ω_{ab} = ω_{cd}, so the exponentials are 1. This is time-independent → easy to propagate.
5. The CGME: the sinc factor damps terms by sinc((ω_{ab}-ω_{cd})·Δτ/2), which for large Δτ kills non-secular terms. This is also time-independent after coarse-graining.
6. The full Redfield: keeps all terms, but now the non-secular ones oscillate at ω_{ab}-ω_{cd}. For heteronuclear pairs, these differences are ~ω_H ≈ 10⁹ rad/s, making full Redfield impractical to propagate (and physically, these rapidly oscillating terms average to zero anyway — which is exactly what the secular approximation does).

**Practical upshot**: For the NMR R₁ρ problem, the secular Lindblad and CGME are the physically relevant methods. The full Redfield method will either require impossibly small time steps or will need the non-secular terms dropped manually (making it equivalent to secular). The three-way comparison (Redfield / secular / CGME) reduces to a two-way comparison (secular Lindblad / CGME) for this system, with the full Redfield serving only as a cross-check in a simpler test system.

### Revised Implementation Architecture

Given the above analysis, the main script should be structured as:

**Step A: Compute lab-frame Redfield tensor.**
```python
system = SpinSystem(...)  # lab-frame H₀
E_lab, U_lab = torch.linalg.eigh(system.get_H0().real.double())
bath = SpectralDensityBath(J_fn, T, coupling_ops)
bath.precompute(E_lab, U_lab.to(torch.complex128))
```

**Step B: Transform to rotating frame for propagation.**

Build the rotating-frame Hamiltonian:
```python
H0_rot = system.get_H0() - omega_N_carrier * Iz_N - omega_H_carrier * Iz_H
```
where the carrier frequencies are the on-resonance Larmor frequencies. This H₀_rot has eigenvalues of order J_NH (~100 Hz), not ω_N or ω_H.

**Step C: Modify the simulator to use the lab-frame Redfield tensor with the rotating-frame Hamiltonian.**

This requires decoupling the eigenbasis used for the Redfield tensor (lab-frame) from the propagation Hamiltonian (rotating-frame). The current simulator assumes they're the same (`self._U = self.bath.U`). This needs to change.

The dissipator D[ρ] = Σ R_{ab,cd} ρ_{cd} is defined in the lab-frame eigenbasis. If ρ is in the rotating frame, it must be transformed to the lab-frame eigenbasis, the dissipator applied, and the result transformed back:
```python
rho_eig = U_lab† @ U_rot @ rho_rot @ U_rot† @ U_lab
D_eig = R @ rho_eig  # in lab eigenbasis
D_rot = U_rot† @ U_lab @ D_eig @ U_lab† @ U_rot
```

Alternatively (and more efficiently): precompute R in the rotating-frame basis by a basis change of the 4-index tensor. This is a one-time cost.

**This is a nontrivial modification to the simulator.** It's the price of getting the physics right for NMR.

---

## Revised Task Breakdown

### Task 1: Single-Spin Analytical Benchmark (1–2 days)

Before touching the propagator, compute R₁ρ analytically for a single ¹⁵N spin using the standard NMR expressions. This establishes the target values that the numerical code must reproduce.

**Deliverables:**
- `nmr_analytical.py`: Functions for R₁, R₂, R₁ρ(ω₁, Δω) from Lipari-Szabo parameters.
- Compute R₁ρ for the Phase 0 parameter set (τ_c=5ns, S²=0.85, B₀=11.7T, ω₁=25–1000 Hz).
- Print results. These are the "standard Redfield" predictions.
- `test_nmr_analytical.py`: Verify R₁ and R₂ against published values for ubiquitin backbone ¹⁵N at 500 MHz.

This task has no dependency on the simulator at all. It uses only the spectral density models and the known NMR relaxation formulae.

### Task 2: Estimate Nonadiabatic Correction Magnitude (1 day)

Using the analytical single-spin model, estimate the nonadiabatic correction.

For a single spin under a constant spin-lock, the nonadiabatic correction enters through the difference between:
- **Dirac**: The bath relaxes the full density matrix including the adiabatic polarization. R₁ρ = R₁cos²θ + R₂sin²θ.
- **Nonadiabatic**: The bath relaxes only the nonadiabatic density matrix σ_nad. During the plateau, σ_nad is constant (dV/dt = 0), so the relaxation acts only on the deviation from the adiabatic state.

The magnitude of the correction is proportional to (ω₁/ω_N)², which for ω₁ = 25 Hz and ω_N = 50 MHz is ~(25/5×10⁷)² ≈ 2.5×10⁻¹³. This is absurdly small.

**This is the Phase 0 decision point.** If the correction is this small, NMR R₁ρ cannot distinguish the frameworks, and the answer from Phase 0 is: proceed to a different experimental candidate.

**However**, the correction structure may be different from the naive (ω₁/ω₀)² estimate. The nonadiabatic framework changes *which bath matrix elements act on which density matrix elements*, not just the overall scale. The correction could be larger if the Redfield tensor has specific structure that amplifies the difference. This needs to be checked numerically.

**Deliverable**: A clear numerical estimate of ΔR₁ρ/R₁ρ for the Phase 0 parameters. If < 0.01, document the negative result and stop. If ≥ 0.01, proceed to Task 3.

### Task 3: Rotating-Frame Simulator Adapter (2–3 days)

Implement the rotating-frame / lab-frame decoupling described above.

**Option A (minimal changes):** Create a `RotatingFrameAdapter` class that:
- Accepts the lab-frame bath (with precomputed Redfield tensor) and the rotating-frame Hamiltonian.
- Provides `get_redfield_dissipator(rho_rot)` etc. that internally transform ρ to the lab eigenbasis, apply the Redfield tensor, and transform back.
- Passes through to the existing Simulator unchanged.

**Option B (Redfield tensor basis change):** Transform R_{ab,cd} from the lab eigenbasis to the rotating-frame eigenbasis as a one-time precomputation. This is more efficient (no per-step basis change) but requires implementing the 4-index tensor transformation.

Option A is simpler and sufficient for a 4×4 system. Go with Option A unless profiling shows the per-step transformation is a bottleneck.

**Deliverable**: `rotating_frame.py` with the adapter, plus tests verifying that the secular dissipator in the rotating frame produces the correct R₁ and R₂.

### Task 4: Full NMR R₁ρ Sweep (1–2 days)

Wire everything together:

```python
# 1. Lab-frame system and bath
system = SpinSystem(...)
E_lab, U_lab = torch.linalg.eigh(system.H0.real.double())
J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
bath = SpectralDensityBath(J_fn, 298.0, [S_DD, S_CSA])
bath.precompute(E_lab, U_lab.to(torch.complex128))

# 2. Rotating-frame Hamiltonian
H0_rot = build_rotating_frame_H0(system, carrier_N, carrier_H)
adapter = RotatingFrameAdapter(bath, U_lab, H0_rot)

# 3. Sweep spin-lock powers
dispersion = R1rhoDispersion()
for omega1_hz in [25, 50, 100, 250, 500, 1000]:
    Ix_N = system.get_spin_operator(0, 'x')
    Iz_N = system.get_spin_operator(0, 'z')
    drive = SpinLockDrive(omega1_hz, Ix_N, Iz_N, offset_hz=0.0,
                          ramp_time_us=1.0, plateau_time_us=50.0)
    
    sim = Simulator(system_rot, drive, adapter)
    
    for method in ['secular_lindblad', 'nonadiabatic_cgme']:
        # Prepare initial state: Ix eigenstate (spin locked along x)
        # In rotating frame, this is the +x eigenstate of ¹⁵N
        rho_0 = build_spin_locked_initial_state(system, spin_index=0)
        
        times = np.linspace(0, 52.0, 5000)  # µs, covers ramp + plateau
        result = sim.run(times, method=method, initial_state=rho_0,
                         delta_tau=delta_tau_cgme)
        
        # Extract R₁ρ from plateau region
        observable = Ix_N  # or the appropriate tilted-frame operator
        extractor = R1rhoExtractor(result, observable)
        fit = extractor.fit(t_start=1.0, t_end=51.0)  # plateau only
        
        dispersion.add_point(method, omega1_hz, fit['R1rho'] * 1e6)  # µs⁻¹ → s⁻¹

# 4. Plot
fig, ax = plt.subplots()
dispersion.plot(ax)
# Add analytical prediction
ax.plot(omega1_array, R1rho_analytical, 'k--', label='Analytical (standard)')
```

**Deliverable**: `main.py` that runs the full sweep and produces a dispersion plot.

### Task 5: Validation and Comparison (1 day)

1. **Cross-check**: Secular Lindblad R₁ρ must match the analytical R₁ρ from Task 1 to within numerical tolerance. If it doesn't, there's a bug in the rotating-frame adapter or the coupling operators.

2. **Nonadiabatic correction**: Compare CGME R₁ρ against secular Lindblad. Quantify ΔR₁ρ/R₁ρ at each spin-lock power.

3. **Convergence checks**: Verify that results don't change with halved time step, doubled plateau time, or doubled number of output points.

4. **δτ sensitivity**: CGME results should be insensitive to the choice of coarse-graining timescale δτ within the window τ_B ≪ δτ ≪ T₁. Sweep δτ from 10 ns to 100 µs and verify R₁ρ is stable.

**Deliverable**: Validation report with convergence plots.

### Task 6: Comparison Against Published Data (1–2 days)

If Phase 0 is positive (correction > 1%):
- Hardcode Korzhnev-Kay (JACS 2005) Fyn SH3 R₁ρ dispersion data.
- Run predictions with published parameters for 3–5 non-exchanging residues.
- Compute χ² for each method.
- Produce the comparison figure.

If Phase 0 is negative (correction < 1%):
- Document the result.
- Compute the precision that would be needed to detect the correction.
- Assess whether HITRAN or FTMW candidates are more promising.

**Deliverable**: Final comparison figure and summary statistics.

---

## Files Created or Modified

| File | Action | Purpose |
|---|---|---|
| `nmr_analytical.py` | New | Analytical R₁, R₂, R₁ρ from NMR formulae |
| `test_nmr_analytical.py` | New | Verify against published values |
| `rotating_frame.py` | New | Lab↔rotating frame adapter for bath |
| `test_rotating_frame.py` | New | Verify dissipator produces correct R₁, R₂ |
| `main.py` | Replace | Full NMR R₁ρ sweep |
| `simulator.py` | Minor edits | Accept adapter in place of bare bath (may be unnecessary if adapter mimics the bath interface) |

---

## Estimated Effort

| Task | Days | Dependency |
|---|---|---|
| Task 1: Analytical benchmark | 1–2 | None |
| Task 2: Nonadiabatic correction estimate | 1 | Task 1 |
| **Decision gate** | — | Task 2 result |
| Task 3: Rotating-frame adapter | 2–3 | Task 2 positive |
| Task 4: Full R₁ρ sweep | 1–2 | Task 3 |
| Task 5: Validation | 1 | Task 4 |
| Task 6: Data comparison | 1–2 | Task 5 |
| **Total** | **7–11 days** | |

The critical path is Task 1 → Task 2 → decision. If the nonadiabatic correction is too small, Tasks 3–6 are not needed, and the effort is 2–3 days to document a negative Phase 0 result.

---

## Risk Register

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Nonadiabatic correction < 1% for NMR R₁ρ | High | Project pivot | The (ω₁/ω₀)² scaling suggests the correction is ~10⁻¹³. Task 2 will confirm. Have HITRAN/FTMW fallback identified. |
| Unit inconsistency producing wrong rates | Medium | Silent errors | Task 1's analytical benchmark provides a known-good reference. Cross-check secular Lindblad against analytical R₁ρ before trusting any numerical results. |
| Rotating-frame adapter introduces bugs | Medium | Wrong physics | Test against single-spin analytical R₁, R₂. These have known values (~1–2 s⁻¹ for R₁, ~5–15 s⁻¹ for R₂ at 500 MHz). |
| RK4 stability with ¹H frequency in lab frame | High if using lab frame | Unusable | Doubly rotating frame eliminates the problem. Fall back to single-spin (2×2) if the two-spin rotating-frame implementation is intractable. |
| Published data not precise enough to distinguish | Medium | Inconclusive result | This is Outcome B from the experiment plan — a valid scientific result. Document the predicted correction magnitude and the precision needed. |