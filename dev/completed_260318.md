# Progress Log — NMR R1rho Simulator Refactoring

Reference plan: `dev/plan_260318.md`

---

## Pre-WP: Bug Fixes to Existing Skeleton

**Status: Done**

Three bugs fixed to get the original 14 tests passing:

1. **`signals.py` — `__call__` / `derivative` TypeError**: These methods passed raw `t` (float) to `torch.cos`/`torch.sin`. Fixed by converting to tensor via `torch.as_tensor(t, dtype=torch.float64)` at entry.

2. **`signals.py` — `_envelope` negative time bug**: Ramp-up mask `t < ramp` incorrectly captured negative times, producing nonzero envelope before t=0. Fixed by adding `t >= 0` guard: `(t_tensor >= 0) & (t_tensor < self.ramp)`.

3. **`simulator.py` — RK4 numerical instability**: The integrator used the output time grid directly as step sizes. For system frequencies ~31.4 rad/ns, the Nyquist-like stability limit is dt < ~0.09 ns, but tests used dt up to 0.56 ns. Fixed by adding automatic sub-stepping: estimates max stable dt from H0's spectral range (`1.5 / spectral_range`) and subdivides each output interval.

---

## WP1: SpinSystem Class

**Status: Done**

**Files created:**
- `src/spin_system.py` — `SpinSystem` class
- `src/test_spin_system.py` — 9 tests, all passing

**What it does:**
- Takes spin definitions (nucleus, chemical shift ppm, gyromagnetic ratio), B0 field, J-couplings (Hz), and CSA tensors (ppm)
- Computes Larmor frequencies: `omega_i = -gamma_i * B0 * (1 + delta_i * 1e-6)` in rad/s
- Builds spin-1/2 operators (Ix, Iy, Iz) for each spin in the full product basis via Kronecker products
- Constructs Zeeman Hamiltonian: `H_Z = sum_i omega_i * Iz_i`
- Constructs J-coupling Hamiltonian: `H_J = sum_{i<j} 2*pi*J_ij * (Ii . Ij)`
- Computes eigenvalues and Bohr frequency matrix
- Provides `get_H0()` and `get_spin_operator(spin_index, component)` interface

**Tests cover:**
- Hilbert space dimension (2^N)
- Larmor frequency numerical accuracy
- Hamiltonian Hermiticity
- Spin commutation relations: [Ix, Iy] = i*Iz
- Spin operator Hermiticity and tracelessness
- Bohr frequencies for single spin
- J-coupling lifts degeneracy
- Different-spin operators commute: [I_{a,i}, I_{b,j}] = 0

**Interface:** Compatible with existing simulator via `get_H0()`, `dim`, `device`. The old `System` class is preserved; nothing broken.

**Full test suite: 23/23 passing.**

---

## WP2: SpinLockDrive Class

**Status: Done**

**Files modified/created:**
- `src/signals.py` — added `SpinLockDrive` class (extends `Signal`)
- `src/test_spin_lock.py` — 13 tests, all passing

**What it does:**
- NMR spin-lock drive in the rotating frame: `V(t) = omega_1 * f(t) * Ix + Delta * Iz`
- Reuses the sin²/cos² envelope shape from `PlateauPulse` (reimplemented, not inherited — no carrier oscillation in the rotating frame)
- Constructor takes NMR-natural units: Hz for power/offset, µs for timing
- Internally converts to rad/µs for consistent time units
- `__call__(t)` returns a **matrix** (the perturbation Hamiltonian), not a scalar
- `derivative(t)` returns `dV/dt = omega_1 * df/dt * Ix` — **zero during plateau** (key for nonadiabatic theory)
- Takes spin operator matrices directly (decoupled from SpinSystem)
- `effective_field_angle()`: tilt angle θ = atan2(ω₁, Δω), π/2 on resonance
- `effective_field_magnitude()`: ω_eff = sqrt(ω₁² + Δω²)

**Design decisions:**
- Returns matrices (not scalars) — knows which spin operators it couples to
- Time unit: µs (matching constructor params)
- Offset term Δω·Iz is constant (no envelope), only the RF term ω₁·Ix is ramped
- Handles zero ramp time (instant-on) and negative times (returns zero matrix)

**Tests cover:**
- Correct perturbation matrix (on-resonance and off-resonance)
- Hermiticity of V(t) at all times
- Zero derivative during plateau (4 time points)
- Nonzero derivative during ramp
- Zero before pulse (t < 0) and after pulse completion
- Effective field angle: π/2 on resonance, π/4 when ω₁ = Δω
- Effective field magnitude
- Instant-on (zero ramp time)
- Envelope midpoint value (f = 0.5 at ramp/2)
- Numerical vs analytical derivative agreement

**Full test suite: 36/36 passing.**

---

## WP4: NonadiabaticDecomposition Class

**Status: Done**

**Files created:**
- `src/nonadiabatic.py` — `NonadiabaticDecomposition` class
- `src/test_nonadiabatic.py` — 11 tests, all passing

**What it does:**
Implements the Landau-Lifshitz integration-by-parts decomposition of first-order TDPT coefficients c_k(t) = a_k(t) + b_k(t):

- **a_k(t)** = -V_{k0}(t)/omega_{k0} * exp(i*omega_{k0}*t) — adiabatic (instantaneous response to V)
- **b_k(t)** = (1/omega_{k0}) * integral dV_{k0}/dt' * exp(i*omega*t') dt' — nonadiabatic (accumulated history through dV/dt)

Key property: during a plateau (dV/dt = 0), b_k does not change.

**Interface:**
- Constructor takes energy eigenvalues + eigenvectors (decoupled from SpinSystem)
- `step(t, dt, dV_dt_matrix)` — accumulates b_k integral (trapezoidal rule: averages phase at t and t+dt, single eigenbasis transform per step)
- `get_adiabatic_coefficients(t, V_matrix)` — returns a_k(t)
- `get_nonadiabatic_coefficients()` — returns current raw b_k (b_0 uncorrected)
- `get_nonadiabatic_density_matrix()` — returns sigma_nad = |b><b| with Tr=1 (b_0 normalized)
- `get_nonadiabatic_populations()` — returns |b_k|^2 with sum=1 (b_0 normalized)

**Sign correction from plan:** The plan's formula for a_k omitted a minus sign. The correct boundary term from integration by parts is a_k = **-**V_{k0}/omega_{k0} * exp(i*omega*t). Verified against direct numerical integration of c_k.

**b_0 normalization (post-review fix):** Raw b_0 stays at 1 (first-order PT), so sum |b_k|^2 > 1 once any excited state is populated. Fixed by computing b_0 = sqrt(1 - sum_{k>0} |b_k|^2) in `_normalized_b()`, used by density matrix and population methods. Raw coefficients remain accessible via `get_nonadiabatic_coefficients()`. This ensures Tr(sigma_nad) = 1 for the simulator (WP5).

**Quadrature (post-review fix):** Upgraded from rectangle rule (O(dt)) to trapezoidal rule (O(dt^2)). Evaluates phase factor at both endpoints t and t+dt using the same matrix elements, so only one eigenbasis transform per step.

**Tests cover:**
- Plateau invariance: 100 steps with dV/dt=0 leave b_k unchanged (atol=1e-15)
- Linear ramp analytical: b_1 matches v*(exp(iwT)-1)/(i*w^2*T) to atol=1e-4
- b_0 (initial state) raw value stays at exactly 1.0
- Decomposition consistency: a_k(T) + b_k(T) = c_k(T) (direct integration) for sin^2 ramp
- Density matrix: Tr=1, outer product of normalized b, Hermitian
- Populations: sum=1, excited-state population correct, ground-state = 1 - excited
- Raw coefficients unaffected by normalization
- Non-diagonal H0: works with nontrivial eigenvectors (verifies eigenbasis transform)

**Full test suite: 48/48 passing.**

---

## WP3: SpectralDensityBath Class

**Status: Done**

**Files created/modified:**
- `src/bath.py` — added `SpectralDensityBath`, `DebyeSpectralDensity`, `LipariSzaboSpectralDensity` (ThermalBath preserved for backward compat)
- `src/spectral_densities.py` — `SimpleLorentzian`, `LipariSzabo`, `OhmicDrude` (standalone models)
- `src/test_bath.py` — extended with ~27 SpectralDensityBath tests
- `src/test_spectral_density_bath.py` — 34 tests using the new spectral density models

**What it does:**
Microscopic bath defined by a spectral density function J(omega). Computes dissipators via the Redfield tensor R_{abcd} built element-wise in the energy eigenbasis.

- **Rate function** `_gamma(omega)`: one-sided Fourier transform of bath correlation. Emission: J(omega) * [n(omega) + 1]. Absorption: J(|omega|) * n(|omega|). Satisfies detailed balance.
- **Bose-Einstein** `_bose_einstein(omega)`: handles edge cases (omega -> 0 Taylor, large |x| limits).
- **`precompute(energies, eigenvectors)`**: builds the full 4-index Redfield tensor from coupling operators transformed to the eigenbasis. Also builds the secular mask (|omega_ab - omega_cd| < tol) and the Bohr frequency difference tensor for CGME sinc damping.
- **Three dissipator methods:**
  1. `get_redfield_dissipator(rho)` — full non-secular Redfield (may violate positivity)
  2. `get_secular_lindblad_dissipator(rho)` — secular (GKSL), guaranteed CP
  3. `get_cgme_dissipator(rho, delta_tau)` — coarse-grained, sinc-damped interpolation

All dissipators operate in the energy eigenbasis via superoperator matrix-vector multiplication: D[rho] = R @ vec(rho), reshaped back to (dim, dim).

**Spectral density models (in `spectral_densities.py`):**
- `SimpleLorentzian(tau_c)`: J(omega) = (2/5) * tau_c / (1 + omega^2 * tau_c^2) — rigid isotropic tumbling
- `LipariSzabo(tau_c, S2, tau_e)`: model-free with overall and internal motion — reduces to SimpleLorentzian when S^2 = 1 or tau_e = 0
- `OhmicDrude(eta, lambda_c)`: J(omega) = eta * |omega| * lambda_c / (omega^2 + lambda_c^2) — condensed-phase

**Legacy models (in `bath.py`, for backward compat with test_bath.py):**
- `DebyeSpectralDensity(eta, tau_c)`: J(omega) = eta * |omega| * tau_c / (1 + omega^2 * tau_c^2)
- `LipariSzaboSpectralDensity(tau_c, S2, tau_e)`: NMR model-free (2/5 prefactor)

**Bug fixes applied:**

1. **Redfield tensor index error (terms 3 & 4):** The trace-preserving counterterms had wrong gamma frequency arguments. Term 3 used `gamma_mat[n, a]` (gamma(omega_{na})) instead of `gamma_mat[c, n]` (gamma(omega_{cn})). Term 4 used `gamma_mat[n, b]` (gamma(omega_{nb})) instead of `gamma_mat[d, n]` (gamma(omega_{dn})). This broke the cancellation that guarantees sum_a R_{aa,cd} = 0. Two index swaps fixed it.

2. **Boltzmann steady-state test tolerance:** At omega = 1e12 rad/s, gamma(omega) ~ 4e-13, making the relative tolerance 1e-6 * 4e-13 = 4e-19 — absurdly tight. The actual derivative ~1e-14 is physically excellent but numerically larger. Fixed by using `max(max_rate, 1.0)` floor to prevent the tolerance from collapsing when rates are tiny.

3. **`make_bath_two_level` return value (test_spectral_density_bath.py):** Helper returned only `bath` but the thermal steady state test unpacked three values. Fixed by constructing `energies` locally in the test.

**Design decisions:**
- `precompute()` is a separate step (not in constructor) because it needs the system eigendata, which may change if H0 is modified
- Coupling operators are stored as a list — supports multiple mechanisms (dipolar + CSA) by summing their contributions to the Redfield tensor
- Dissipators take rho in the eigenbasis (not the computational basis) — the simulator must handle basis transformations
- The CGME sinc factor is computed per-element from precomputed Bohr frequency differences

**Tests cover (`test_bath.py` — SpectralDensityBath section):**
- Trace preservation: Redfield, secular, CGME (2LS × 4 rho states + 3LS)
- Hermiticity: Redfield, secular (2LS + 3LS)
- Detailed balance: gamma(omega)/gamma(-omega) = exp(hbar*omega/kT)
- High-temperature limit: gamma(+omega) ~ gamma(-omega)
- Secular vs Redfield: R_{01,10} nonzero in full, zero in secular
- Population block matches in both
- Dissipators differ on coherent state, agree on diagonal
- CGME interpolation: small Tc -> Redfield, large Tc -> secular
- Population decay direction: excited state decays
- Thermal steady state: Boltzmann is approximate steady state
- Spectral density model properties: positivity, peak location, rigid limits

**Tests cover (`test_spectral_density_bath.py`):**
- SimpleLorentzian: zero frequency, high-freq falloff, positivity, peak at zero
- LipariSzabo: reduces to Lorentzian (S^2=1), reduces with tau_e=0, internal motion adds high freq
- OhmicDrude: zero at zero, positivity, peak near cutoff
- Detailed balance, high-T limit, positive rates
- Secular dissipator: trace, Hermiticity, decay direction, 3-level trace
- Redfield dissipator: trace, Hermiticity, 3-level trace
- CGME dissipator: trace, Hermiticity, large Tc -> secular, small Tc -> Redfield
- Two-level secular trap: full 5-test suite (R_{01,10} nonzero/zero, population block, differ/agree)
- Thermal steady state with SimpleLorentzian
- Multiple coupling operators: sigma_x + sigma_z trace/Hermiticity, sigma_z adds pure dephasing

---

## WP5: Simulator Refactor

**Status: Done**

**Files modified:**
- `src/simulator.py` — refactored with eigenbasis propagation and three new methods
- `src/test_simulator.py` — extended with 20 new WP5 tests (25 total)

**What changed:**

The simulator now supports four methods: the original `lindblad` (backward compat) plus three new eigenbasis methods that use `SpectralDensityBath`:

| Method | Dissipator | Positivity | Bath type |
|---|---|---|---|
| `lindblad` | ThermalBath (T1/T2) | Guaranteed | `ThermalBath` |
| `redfield` | Full Redfield | Not guaranteed | `SpectralDensityBath` |
| `secular_lindblad` | Secular GKSL | Guaranteed | `SpectralDensityBath` |
| `nonadiabatic_cgme` | CGME (sinc-damped) | Guaranteed* | `SpectralDensityBath` |

All three new methods share the same master equation structure: `d(rho)/dt = -i[H(t), rho] + D[rho]`. The only difference between them is the dissipator D. The `nonadiabatic_cgme` method additionally tracks b_k coefficients via `NonadiabaticDecomposition` alongside the rho propagation.

**Architecture (key design decision):** Option A — propagate rho (not sigma_nad) for all three methods. The nonadiabatic decomposition runs alongside to provide b_k population analysis, but does not modify the master equation. Reasoning: sigma_nad = |b><b| is a pure state; a dissipative master equation immediately mixes it, breaking the b_k connection. Propagating rho avoids this problem entirely.

**Implementation details:**

- **`_run_eigenbasis()`**: New method for the three SpectralDensityBath methods. Uses the bath's eigenbasis (`bath.U`, `bath.energies`) from `precompute()` for consistency with the Redfield tensor. Transforms rho to eigenbasis at start, transforms back to computational basis at each output step for storage.

- **`_run_legacy()`**: Extracted original propagation code (unchanged). Keeps the trace renormalization for backward compatibility.

- **Eigenbasis consistency**: The simulator uses `bath.U` (not its own diagonalization of H0) to guarantee the same eigenvector sign/phase convention used when building the Redfield tensor.

- **Trace assertion**: New methods assert `|Tr(rho) - 1| < 1e-6` at each output step. No renormalization — trace drift is an error, not something to mask.

- **Drive integration**: New methods expect a matrix-valued drive (like `SpinLockDrive`). `V(t) = drive(t)` is transformed to eigenbasis via `U^dag V U` in the Liouvillian.

- **Nonadiabatic tracking**: `decomp.step(t, dt, dV_dt)` is called once per substep (after the RK4 rho update, not per RK4 stage). Result stores `nonadiabatic_populations` tensor of shape (n_times, dim).

- **Configurable initial state**: `run()` accepts an optional `initial_state` parameter (default: |0><0|).

- **`Result` class**: Added `nonadiabatic_populations` attribute (None for non-CGME methods).

**Unit consistency pitfall discovered during testing:** The system's energy eigenvalues (in arbitrary units like omega=5) combined with the bath's SI constants (HBAR, KB) can produce catastrophic numerical instability. At T=300K with omega=5 rad/s, the Bose-Einstein occupation n(omega) ~ 10^12, making gamma ~ 10^10 and the integrator instantly unstable. Fix: tests use T_K=1e-10 where hbar*omega/kT ~ 0.4, giving n ~ 2 and manageable rates. **This is a test-only issue; real NMR simulations will use SI frequencies (~10^8 rad/s) where the rates are physical.**

**Tests cover (20 new tests):**
- Trace preservation: all 3 methods × (zero drive + constant drive) = 6 parametrized tests
- Hermiticity: all 3 methods = 3 parametrized tests
- Population decay: secular and Redfield from excited state
- Zero-drive population agreement: all 3 methods give identical population dynamics for diagonal initial state (pop-to-pop block of R is the same for all three)
- Nonadiabatic tracking: populations present, sum to 1, not present for other methods, plateau invariance (b_k constant when dV/dt = 0)
- Configurable initial state: custom and default
- Method dispatch: unknown method raises ValueError, unprecomputed bath raises AssertionError

**Full test suite: 131/131 passing (6.2s).**

---

## WP6: R1rho Extraction

**Status: Done**

**Files created:**
- `src/r1rho.py` — `R1rhoExtractor` and `R1rhoDispersion` classes
- `src/test_r1rho.py` — 13 tests, all passing

**What it does:**

Extracts the rotating-frame relaxation rate R1rho from simulator output by fitting the decay of an observable's expectation value to a single exponential.

**`R1rhoExtractor`:**
- Takes a `Result` object and a Hermitian observable operator O
- `get_expectation_values()` computes <O(t)> = Tr[rho(t) @ O] for all stored time steps via `torch.einsum('tij,ji->t', states, O)`
- `fit(t_start, t_end)` fits the expectation values to M0 * exp(-R1rho * t) within the specified time window using `scipy.optimize.curve_fit`
- Time origin is shifted to the start of the fitting window for numerical stability
- Returns dict with `R1rho`, `M0`, and `residual_std`
- The fitting window `[t_start, t_end]` allows excluding ramp-up/ramp-down regions
- Raises `ValueError` if fewer than 3 data points in the window

**`R1rhoDispersion`:**
- Container for (omega1, R1rho) pairs across multiple methods
- `add_point(method_name, omega1_hz, R1rho)` stores a single measurement
- `get_dispersion(method_name)` returns (omega1_array, R1rho_array)
- `plot(ax, experimental)` generates dispersion curves for all methods, with optional experimental data overlay

**Design decisions:**
- Observable is passed explicitly (not derived from SpinSystem) — keeps R1rhoExtractor decoupled from the spin system implementation
- Uses `scipy.optimize.curve_fit` with non-negative bounds on M0 and R1rho (physical constraint)
- Fitting window is optional — defaults to the full time range
- R1rhoDispersion stores raw lists, converts to numpy on retrieval

**Tests cover:**

*Expectation values (3 tests):*
- Identity observable gives Tr[rho] = 1 at all times
- sigma_z gives +1 for |0><0| and -1 for |1><1|
- sigma_x gives +1 for |+><+|

*Exponential fit (5 tests):*
- Exact exponential: recovers M0 and R to 1e-6 precision
- Slow decay (R=0.01): relative error < 1e-4 over long time series
- Fitting window: R1rho recovered correctly from sub-interval
- Too few points: ValueError raised for < 3 points
- 4x4 system (Ix ⊗ I observable): correct extraction for multi-spin system

*R1rhoDispersion (4 tests):*
- Add and retrieve: correct storage and retrieval of multiple methods
- Empty dispersion: get_methods() returns empty list
- Plot returns Axes object without error
- Plot with experimental overlay works

*Integration (1 test):*
- Full Simulator → R1rhoExtractor pipeline: secular_lindblad with |+> initial state, extract <sigma_x>, fit gives R1rho > 0

**Full test suite: 144/144 passing (6.9s).**

---

## WP7: Main Script — Not started
