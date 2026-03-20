"""
Phase 0A feasibility estimate: nonadiabatic correction to HITRAN line mixing.

Estimates the magnitude of the nonadiabatic correction to pressure-broadened
rotational line mixing for linear molecules (CO, CO2, N2O, HCl). Mirrors
the structure of nonadiabatic_estimate.py (the NMR Phase 0 that returned
NEGATIVE).

Two independent correction mechanisms are analyzed:

1. Perturbative mixing: |a|^2 = (V / DeltaE)^2
   The fraction of the relaxation matrix that is adiabatically modified by
   collisional coupling. V is the off-diagonal coupling strength (~pressure
   broadening coefficient * pressure) and DeltaE is the energy gap between
   adjacent rotational levels.

2. CGME vs secular: sinc((omega_ab - omega_cd) * delta_tau / 2)
   The coarse-grained dissipator retains non-secular terms damped by sinc
   factors. For widely spaced rotational lines with fast collisions
   (tau_c ~ 1 ps), the sinc factors are small. But for closely spaced
   lines or slow collisions, they can be O(1).

The correction scales as (V/DeltaE)^2. For CO at 1 atm, V/DeltaE ~ 0.017
(correction ~ 3e-4, NEGATIVE). For CO2 at 5 atm, V/DeltaE ~ 0.45
(correction ~ 0.20, strongly POSITIVE).
"""
import math
from rotational_system import RotationalSystem, CM_TO_RAD_S
from nonadiabatic_estimate import cgme_sinc_retention


def perturbation_ratio_table(B_cm, J_max, gamma_pressure_cm, pressure_atm=1.0):
    """Compute V/DeltaE and mixing fraction for all adjacent rotational levels.

    For each pair (J, J+1), computes:
    - DeltaE = 2 * B * (J+1) in cm^-1
    - V = gamma_pressure * pressure in cm^-1
    - ratio = V / DeltaE
    - mixing = ratio^2 (the nonadiabatic correction scale)

    Args:
        B_cm: rotational constant in cm^-1
        J_max: maximum J quantum number
        gamma_pressure_cm: pressure-broadening coefficient in cm^-1/atm
        pressure_atm: buffer gas pressure in atm

    Returns:
        list of dicts with keys: J, delta_E_cm, V_cm, ratio, mixing_fraction
    """
    V_cm = gamma_pressure_cm * pressure_atm
    results = []
    for J in range(J_max):
        delta_E = 2.0 * B_cm * (J + 1)
        ratio = V_cm / delta_E
        results.append({
            'J': J,
            'delta_E_cm': delta_E,
            'V_cm': V_cm,
            'ratio': ratio,
            'mixing_fraction': ratio ** 2,
        })
    return results


def cgme_sinc_table_hitran(B_cm, J_max, delta_tau_s=1e-12):
    """Compute CGME sinc factors for adjacent Bohr frequency pairs.

    For a rigid rotor, the non-secular Redfield elements couple different
    pairs of transitions. The CGME damps these by sinc((omega_ab - omega_cd)
    * delta_tau / 2). We compute the sinc factor for the most relevant pairs:
    adjacent R-branch transitions separated by 2*B.

    Args:
        B_cm: rotational constant in cm^-1
        J_max: maximum J quantum number
        delta_tau_s: coarse-graining timescale in seconds (collision duration)

    Returns:
        list of dicts with keys: label, freq_diff_rad_s, freq_diff_cm, sinc_factor
    """
    system = RotationalSystem(B_cm=B_cm, J_max=J_max)
    results = []

    # Adjacent level pairs: (J, J+1) vs (J+1, J+2) differ by 2*B cm^-1
    spacing_cm = 2.0 * B_cm
    spacing_rad = spacing_cm * CM_TO_RAD_S
    sinc_val = cgme_sinc_retention(spacing_rad, delta_tau_s)
    results.append({
        'label': 'Adjacent R-branch (2B)',
        'freq_diff_cm': spacing_cm,
        'freq_diff_rad_s': spacing_rad,
        'sinc_factor': sinc_val,
        'sinc_abs': abs(sinc_val),
    })

    # Next-nearest: differ by 4*B
    spacing2_cm = 4.0 * B_cm
    spacing2_rad = spacing2_cm * CM_TO_RAD_S
    sinc_val2 = cgme_sinc_retention(spacing2_rad, delta_tau_s)
    results.append({
        'label': 'Next-nearest R-branch (4B)',
        'freq_diff_cm': spacing2_cm,
        'freq_diff_rad_s': spacing2_rad,
        'sinc_factor': sinc_val2,
        'sinc_abs': abs(sinc_val2),
    })

    # Smallest Bohr frequency: J=0 -> J=1, DeltaE = 2B
    smallest_cm = 2.0 * B_cm
    smallest_rad = smallest_cm * CM_TO_RAD_S
    sinc_smallest = cgme_sinc_retention(smallest_rad, delta_tau_s)
    results.append({
        'label': 'Smallest Bohr freq (2B)',
        'freq_diff_cm': smallest_cm,
        'freq_diff_rad_s': smallest_rad,
        'sinc_factor': sinc_smallest,
        'sinc_abs': abs(sinc_smallest),
    })

    # Largest Bohr frequency: J=0 -> J=J_max
    largest_cm = B_cm * J_max * (J_max + 1)
    largest_rad = largest_cm * CM_TO_RAD_S
    sinc_largest = cgme_sinc_retention(largest_rad, delta_tau_s)
    results.append({
        'label': f'Largest Bohr freq (J=0->J={J_max})',
        'freq_diff_cm': largest_cm,
        'freq_diff_rad_s': largest_rad,
        'sinc_factor': sinc_largest,
        'sinc_abs': abs(sinc_largest),
    })

    return results


def phase0a_report(B_cm=1.9313, J_max=10, gamma_pressure_cm=0.065,
                   pressure_atm=1.0, T_K=296.0, tau_collision_s=1e-12):
    """Generate the Phase 0A feasibility report for HITRAN line mixing.

    Computes the nonadiabatic correction magnitude for collisional line
    mixing of a linear rigid rotor molecule and returns a structured report.

    Args:
        B_cm: rotational constant in cm^-1 (default: CO)
        J_max: maximum J quantum number
        gamma_pressure_cm: pressure-broadening coefficient in cm^-1/atm
        pressure_atm: buffer gas pressure in atm
        T_K: temperature in Kelvin
        tau_collision_s: collision correlation time in seconds

    Returns:
        dict with full analysis including verdict
    """
    system = RotationalSystem(B_cm=B_cm, J_max=J_max)

    # Correction 1: perturbation ratios
    ratio_table = perturbation_ratio_table(B_cm, J_max, gamma_pressure_cm,
                                           pressure_atm)
    max_ratio = max(r['ratio'] for r in ratio_table)
    max_mixing = max(r['mixing_fraction'] for r in ratio_table)

    # Correction 2: CGME sinc factors
    sinc_table = cgme_sinc_table_hitran(B_cm, J_max, tau_collision_s)
    max_sinc = max(s['sinc_abs'] for s in sinc_table)

    return {
        'parameters': {
            'B_cm': B_cm,
            'J_max': J_max,
            'gamma_pressure_cm': gamma_pressure_cm,
            'pressure_atm': pressure_atm,
            'T_K': T_K,
            'tau_collision_s': tau_collision_s,
            'V_cm': gamma_pressure_cm * pressure_atm,
            'V_rad_s': gamma_pressure_cm * pressure_atm * CM_TO_RAD_S,
            'smallest_gap_cm': 2.0 * B_cm,
            'smallest_gap_rad_s': 2.0 * B_cm * CM_TO_RAD_S,
        },
        'perturbation_table': ratio_table,
        'cgme_sinc_table': sinc_table,
        'max_perturbation_ratio': max_ratio,
        'max_mixing_fraction': max_mixing,
        'max_sinc_retention': max_sinc,
        'verdict': 'POSITIVE' if max_mixing >= 0.01 else 'NEGATIVE',
        'threshold': 0.01,
    }


def print_phase0a_report(**kwargs):
    """Print a formatted Phase 0A feasibility report to stdout."""
    report = phase0a_report(**kwargs)
    p = report['parameters']

    print("=" * 72)
    print("PHASE 0A FEASIBILITY REPORT: Nonadiabatic Correction to HITRAN Line Mixing")
    print("=" * 72)
    print()
    print("System parameters:")
    print(f"  B = {p['B_cm']:.4f} cm^-1")
    print(f"  J_max = {p['J_max']}")
    print(f"  Pressure broadening: gamma = {p['gamma_pressure_cm']:.4f} cm^-1/atm")
    print(f"  Pressure: P = {p['pressure_atm']:.1f} atm")
    print(f"  Temperature: T = {p['T_K']:.0f} K")
    print(f"  Collision time: tau_c = {p['tau_collision_s']*1e12:.1f} ps")
    print(f"  V = gamma * P = {p['V_cm']:.4f} cm^-1  ({p['V_rad_s']:.3e} rad/s)")
    print(f"  Smallest gap: 2B = {p['smallest_gap_cm']:.4f} cm^-1  ({p['smallest_gap_rad_s']:.3e} rad/s)")
    print()

    print("-" * 72)
    print("Correction 1: Perturbative mixing  |a|^2 = (V / DeltaE)^2")
    print("-" * 72)
    print(f"  {'J':>3s}  {'DeltaE (cm^-1)':>14s}  {'V (cm^-1)':>10s}  "
          f"{'V/DeltaE':>10s}  {'|a|^2':>12s}")
    for r in report['perturbation_table']:
        print(f"  {r['J']:3d}  {r['delta_E_cm']:14.4f}  {r['V_cm']:10.4f}  "
              f"{r['ratio']:10.4e}  {r['mixing_fraction']:12.4e}")
    print()

    print("-" * 72)
    print(f"Correction 2: CGME sinc factors  (tau_c = {p['tau_collision_s']*1e12:.0f} ps)")
    print("-" * 72)
    print(f"  {'Pair':>32s}  {'Dw (cm^-1)':>12s}  {'Dw (rad/s)':>12s}  {'|sinc|':>10s}")
    for s in report['cgme_sinc_table']:
        print(f"  {s['label']:>32s}  {s['freq_diff_cm']:12.4f}  "
              f"{s['freq_diff_rad_s']:12.3e}  {s['sinc_abs']:10.2e}")
    print()

    print("=" * 72)
    print(f"VERDICT: {report['verdict']}")
    print(f"  Max perturbation ratio V/DeltaE: {report['max_perturbation_ratio']:.4e}")
    print(f"  Max mixing fraction |a|^2:       {report['max_mixing_fraction']:.4e}")
    print(f"  Max sinc retention:              {report['max_sinc_retention']:.4e}")
    print(f"  Threshold for proceeding:        {report['threshold']:.0e}")
    if report['verdict'] == 'POSITIVE':
        print()
        print("  The nonadiabatic correction exceeds 1% for this system.")
        print("  Proceed to full comparison against HITRAN data.")
    elif report['verdict'] == 'NEGATIVE':
        print()
        print("  The nonadiabatic correction is below 1% for this system.")
        print("  Consider higher pressure or a molecule with smaller B.")
    print("=" * 72)


def print_molecule_pressure_sweep():
    """Print Phase 0A results for multiple molecules and pressures."""
    molecules = [
        ('CO',   1.9313, 0.065),
        ('CO2',  0.3902, 0.070),
        ('N2O',  0.4190, 0.073),
        ('HCl', 10.4400, 0.050),
    ]
    pressures = [0.1, 1.0, 5.0, 10.0, 50.0]

    print()
    print("=" * 80)
    print("MOLECULE-PRESSURE SWEEP: Nonadiabatic Correction Landscape")
    print("=" * 80)
    print()
    print(f"  {'Molecule':>8s}  {'B (cm^-1)':>10s}  {'P (atm)':>8s}  "
          f"{'V/DeltaE':>10s}  {'|a|^2':>10s}  {'Verdict':>8s}")
    print("-" * 80)

    for name, B, gamma in molecules:
        for P in pressures:
            report = phase0a_report(B_cm=B, J_max=15,
                                    gamma_pressure_cm=gamma,
                                    pressure_atm=P)
            print(f"  {name:>8s}  {B:10.4f}  {P:8.1f}  "
                  f"{report['max_perturbation_ratio']:10.4e}  "
                  f"{report['max_mixing_fraction']:10.4e}  "
                  f"{report['verdict']:>8s}")
        print()

    print("=" * 80)


if __name__ == '__main__':
    # Default: CO at 1 atm
    print_phase0a_report()
    print()

    # CO2 at 5 atm (expected POSITIVE)
    print_phase0a_report(B_cm=0.3902, gamma_pressure_cm=0.070, pressure_atm=5.0)
    print()

    # Full sweep
    print_molecule_pressure_sweep()
