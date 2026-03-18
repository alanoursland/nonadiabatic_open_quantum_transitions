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

## WP3: SpectralDensityBath — Not started
## WP4: NonadiabaticDecomposition — Not started
## WP5: Simulator Refactor — Not started
## WP6: R1rho Extraction — Not started
## WP7: Main Script — Not started
