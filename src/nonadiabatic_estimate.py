"""
Nonadiabatic correction estimates for NMR R1rho.

Phase 0 feasibility check: estimates the magnitude of the nonadiabatic
correction to R1rho for a 15N backbone amide under a spin-lock field.

Two independent correction mechanisms are analyzed:

1. Perturbative mixing: |a_1|^2 = (omega1 / (2*omega_N))^2
   The fraction of the density matrix that is adiabatically polarized
   by the spin-lock field. For NMR: omega1 ~ 10^2 rad/s,
   omega_N ~ 10^8 rad/s, giving |a_1|^2 ~ 10^-13.

2. CGME vs secular: sinc((omega_ab - omega_cd) * delta_tau / 2)
   The coarse-grained dissipator retains slightly more non-secular terms
   than the secular approximation. For NMR with delta_tau ~ 1 us,
   the sinc factors for non-degenerate Bohr frequency pairs are ~ 10^-3,
   and these non-secular terms contribute negligibly to the total rate.

Both corrections are far below any measurable threshold.
Phase 0 verdict: NMR R1rho CANNOT distinguish the nonadiabatic framework.
"""
import math
from nmr_analytical import NMRRelaxation, GAMMA_15N, GAMMA_1H
from spectral_densities import LipariSzabo


def adiabatic_mixing_fraction(omega1, omega0):
    """Leading-order adiabatic population in the excited state.

    For a 2-level system with unperturbed splitting omega0 and a transverse
    perturbation of strength omega1 (e.g. spin-lock), the first-order TDPT
    adiabatic coefficient gives:

        |a_1|^2 = (V_10 / omega0)^2 = (omega1 / (2 * omega0))^2

    where V_10 = omega1/2 is the matrix element of omega1*Ix between the
    two spin-1/2 eigenstates of Iz.

    This fraction represents the part of the density matrix that responds
    instantaneously to the perturbation (adiabatic polarization). In the
    nonadiabatic framework, this part is excluded from bath-induced relaxation.

    Args:
        omega1: perturbation strength in rad/s
        omega0: unperturbed Bohr frequency in rad/s (must be > 0)

    Returns:
        |a_1|^2 (dimensionless)
    """
    if omega0 <= 0:
        raise ValueError(f"omega0 must be positive, got {omega0}")
    return (omega1 / (2.0 * omega0)) ** 2


def cgme_sinc_retention(bohr_freq_diff, delta_tau):
    """CGME sinc damping factor for a non-secular Redfield tensor element.

    The CGME dissipator multiplies each Redfield element R_{ab,cd} by:
        sinc((omega_ab - omega_cd) * delta_tau / 2)

    where sinc(x) = sin(x)/x (unnormalized sinc).

    For secular elements (omega_ab = omega_cd): factor = 1.
    For non-secular elements with large frequency difference: factor -> 0.

    Args:
        bohr_freq_diff: |omega_ab - omega_cd| in rad/s
        delta_tau: coarse-graining timescale in seconds

    Returns:
        sinc factor (dimensionless, between -1 and 1)
    """
    x = bohr_freq_diff * delta_tau / 2.0
    if abs(x) < 1e-15:
        return 1.0
    return math.sin(x) / x


def nonadiabatic_correction_table(nmr, omega1_values_hz):
    """Compute nonadiabatic correction estimates for a sweep of spin-lock powers.

    For each omega1, computes:
    - The adiabatic mixing fraction |a_1|^2 = (omega1/(2*omega_N))^2
    - The resulting correction to R1rho: delta ~ |a_1|^2 * |R2 - R1|
    - The relative correction delta / R1rho

    Args:
        nmr: NMRRelaxation instance with precomputed R1, R2, omega_N
        omega1_values_hz: list of spin-lock powers in Hz

    Returns:
        list of dicts with keys: omega1_hz, omega1_rad_s, mixing_fraction,
        delta_R1rho, relative_correction
    """
    R1 = nmr.R1()
    R2 = nmr.R2()
    results = []
    for w1_hz in omega1_values_hz:
        omega1 = 2.0 * math.pi * w1_hz
        mixing = adiabatic_mixing_fraction(omega1, nmr.omega_N)
        # On resonance, R1rho = R2. The correction is bounded by
        # |a_1|^2 * |R2 - R1| (the adiabatic part sees R1 instead of R2).
        delta = mixing * abs(R2 - R1)
        relative = delta / R2
        results.append({
            'omega1_hz': w1_hz,
            'omega1_rad_s': omega1,
            'mixing_fraction': mixing,
            'delta_R1rho': delta,
            'relative_correction': relative,
        })
    return results


def cgme_sinc_table(nmr, delta_tau_s=1e-6):
    """Compute CGME sinc retention factors for all NMR Bohr frequency pairs.

    For a 15N-1H system, the lab-frame Bohr frequencies include:
    omega_N, omega_H, omega_H +/- omega_N, and 0.

    The CGME retains non-secular elements by sinc((w_ab - w_cd) * delta_tau / 2).
    For typical NMR frequencies and delta_tau ~ 1 us, all non-degenerate
    pairs have sinc factors << 1.

    Args:
        nmr: NMRRelaxation instance
        delta_tau_s: coarse-graining timescale in seconds (default: 1 us)

    Returns:
        list of dicts with keys: label, freq_diff_rad_s, sinc_factor
    """
    wN = nmr.omega_N
    wH = nmr.omega_H
    pairs = [
        ('omega_N', wN),
        ('omega_H', wH),
        ('omega_H - omega_N', wH - wN),
        ('omega_H + omega_N', wH + wN),
        ('2*omega_N', 2 * wN),
        ('2*omega_H', 2 * wH),
    ]
    results = []
    for label, freq_diff in pairs:
        sinc = cgme_sinc_retention(freq_diff, delta_tau_s)
        results.append({
            'label': label,
            'freq_diff_rad_s': freq_diff,
            'sinc_factor': sinc,
            'sinc_abs': abs(sinc),
        })
    return results


def phase0_report(B0=11.7, tau_c=5e-9, S2=0.85, tau_e=50e-12, delta_tau_s=1e-6):
    """Generate the Phase 0 feasibility report.

    Computes the nonadiabatic correction magnitude for the NMR R1rho
    experiment and returns a structured report.

    Args:
        B0: static field in Tesla
        tau_c: overall correlation time in seconds
        S2: order parameter
        tau_e: internal motion correlation time in seconds
        delta_tau_s: CGME coarse-graining timescale in seconds

    Returns:
        dict with full analysis including verdict
    """
    J_fn = LipariSzabo(tau_c=tau_c, S2=S2, tau_e=tau_e)
    nmr = NMRRelaxation(B0=B0, J_fn=J_fn)

    omega1_values = [25, 50, 100, 250, 500, 1000]
    corrections = nonadiabatic_correction_table(nmr, omega1_values)
    sinc_table = cgme_sinc_table(nmr, delta_tau_s)

    max_relative = max(r['relative_correction'] for r in corrections)
    max_sinc = max(r['sinc_abs'] for r in sinc_table)

    return {
        'parameters': {
            'B0_T': B0,
            'tau_c_s': tau_c,
            'S2': S2,
            'tau_e_s': tau_e,
            'delta_tau_s': delta_tau_s,
            'omega_N_rad_s': nmr.omega_N,
            'omega_H_rad_s': nmr.omega_H,
            'd_rad_s': nmr.d,
            'c_rad_s': nmr.c,
        },
        'rates': {
            'R1_per_s': nmr.R1(),
            'R2_per_s': nmr.R2(),
            'NOE': nmr.NOE(),
        },
        'perturbative_corrections': corrections,
        'cgme_sinc_factors': sinc_table,
        'max_relative_correction': max_relative,
        'max_sinc_retention': max_sinc,
        'verdict': 'NEGATIVE' if max_relative < 0.01 else 'POSITIVE',
        'threshold': 0.01,
    }


def print_phase0_report(**kwargs):
    """Print a formatted Phase 0 feasibility report to stdout."""
    report = phase0_report(**kwargs)
    p = report['parameters']
    r = report['rates']

    print("=" * 70)
    print("PHASE 0 FEASIBILITY REPORT: Nonadiabatic Correction to NMR R1rho")
    print("=" * 70)
    print()
    print("System parameters:")
    print(f"  B0 = {p['B0_T']:.2f} T")
    print(f"  tau_c = {p['tau_c_s']*1e9:.1f} ns")
    print(f"  S^2 = {p['S2']:.2f}")
    print(f"  tau_e = {p['tau_e_s']*1e12:.0f} ps")
    print(f"  omega_N = {p['omega_N_rad_s']:.3e} rad/s  ({p['omega_N_rad_s']/2/math.pi/1e6:.2f} MHz)")
    print(f"  omega_H = {p['omega_H_rad_s']:.3e} rad/s  ({p['omega_H_rad_s']/2/math.pi/1e6:.2f} MHz)")
    print(f"  d = {p['d_rad_s']:.0f} rad/s")
    print(f"  c = {p['c_rad_s']:.0f} rad/s")
    print()
    print("Relaxation rates (analytical, standard Redfield):")
    print(f"  R1 = {r['R1_per_s']:.4f} s^-1   (T1 = {1/r['R1_per_s']:.3f} s)")
    print(f"  R2 = {r['R2_per_s']:.4f} s^-1   (T2 = {1/r['R2_per_s']:.3f} s)")
    print(f"  NOE = {r['NOE']:.4f}")
    print()

    print("-" * 70)
    print("Correction 1: Perturbative adiabatic mixing  |a_1|^2 = (w1/(2*wN))^2")
    print("-" * 70)
    print(f"  {'omega1 (Hz)':>12s}  {'|a_1|^2':>12s}  {'dR1rho (s^-1)':>14s}  {'dR1rho/R1rho':>14s}")
    for c in report['perturbative_corrections']:
        print(f"  {c['omega1_hz']:12.0f}  {c['mixing_fraction']:12.2e}  "
              f"{c['delta_R1rho']:14.2e}  {c['relative_correction']:14.2e}")
    print()

    print("-" * 70)
    print(f"Correction 2: CGME vs secular sinc factors  (delta_tau = {p['delta_tau_s']*1e6:.0f} us)")
    print("-" * 70)
    print(f"  {'Frequency pair':>22s}  {'|Dw| (rad/s)':>14s}  {'|sinc|':>10s}")
    for s in report['cgme_sinc_factors']:
        print(f"  {s['label']:>22s}  {s['freq_diff_rad_s']:14.3e}  {s['sinc_abs']:10.2e}")
    print()
    print("  Non-secular Redfield elements are suppressed by these sinc factors.")
    print("  Since all |sinc| << 1, CGME agrees with secular Lindblad to high precision.")
    print()

    print("=" * 70)
    print(f"VERDICT: {report['verdict']}")
    print(f"  Max relative correction: {report['max_relative_correction']:.2e}")
    print(f"  Threshold for proceeding: {report['threshold']:.0e}")
    if report['verdict'] == 'NEGATIVE':
        print()
        print("  The nonadiabatic correction to NMR R1rho is ~10^-13, far below")
        print("  any achievable measurement precision (~10^-3 at best).")
        print("  NMR R1rho CANNOT distinguish the nonadiabatic framework.")
        print("  Recommend: pivot to HITRAN/FTMW experimental candidates.")
    print("=" * 70)


if __name__ == '__main__':
    print_phase0_report()
