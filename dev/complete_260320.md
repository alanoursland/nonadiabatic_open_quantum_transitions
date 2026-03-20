# Phase 0A: HITRAN Line Mixing Feasibility Estimate — Progress

## Status: COMPLETE (53 new tests, 265/265 total passing)

## Files Created/Modified

| File | Action | Tests |
|------|--------|-------|
| `src/rotational_system.py` | New | 20 |
| `src/test_rotational_system.py` | New | — |
| `src/spectral_densities.py` | Appended `CollisionalLorentzian` | 6 |
| `src/test_collisional_spectral_density.py` | New | — |
| `src/hitran_estimate.py` | New | 20 + 7 integration |
| `src/test_hitran_estimate.py` | New | — |

## Results

### Correction 1: Perturbative mixing |a|² = (V/ΔE)²

| System | Pressure | V/ΔE (J=0) | |a|² | Verdict |
|--------|----------|-------------|------|---------|
| CO 1 atm | 1.0 | 0.017 | 2.8e-4 | NEGATIVE |
| CO 10 atm | 10.0 | 0.168 | 0.028 | **POSITIVE** |
| CO₂ 1 atm | 1.0 | 0.090 | 8.0e-3 | NEGATIVE (borderline) |
| CO₂ 5 atm | 5.0 | 0.449 | 0.201 | **POSITIVE** |
| N₂O 5 atm | 5.0 | 0.436 | 0.190 | **POSITIVE** |
| HCl 50 atm | 50.0 | 0.120 | 0.014 | **POSITIVE** |

### Correction 2: CGME sinc factors (τ_c = 1 ps)

Unlike NMR (where sinc factors were ~10⁻³), HITRAN sinc factors are near 1.
Adjacent R-branch pairs: sinc ≈ 0.98 (CO) to 0.999 (CO₂).
This means the CGME retains nearly all non-secular Redfield elements —
the secular approximation is NOT valid for these timescales.

### Key Finding

The nonadiabatic correction to HITRAN line mixing is **POSITIVE for CO₂ and
N₂O at moderate pressures (≥5 atm)**, and for CO at higher pressures (≥10 atm).
CO₂ at 1 atm is borderline (|a|² = 0.8%, just below the 1% threshold).

**Best candidate: CO₂ Q-branch at 5-10 atm in N₂ buffer gas.**
- V/ΔE ~ 0.45-0.90 for the lowest rotational transitions
- Expected correction to relaxation matrix: 20-80%
- Well-characterized HITRAN data available
- Standard Redfield theory known to have discrepancies for CO₂ Q-branches

### Important Physics Note

At V/ΔE > 0.5, first-order perturbation theory breaks down and |a|² > 1 is
unphysical. The values for CO₂ at 10+ atm and CO at 50 atm are in this
regime. The full nonadiabatic framework (not just the perturbative estimate)
is needed there. But the key point is: the correction is clearly large enough
to be measurable for CO₂ at 5 atm, where V/ΔE ≈ 0.45 and perturbation
theory is still marginally valid.
