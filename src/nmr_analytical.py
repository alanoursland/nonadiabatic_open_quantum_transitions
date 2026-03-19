"""
Analytical NMR relaxation rates for 15N backbone amides.

Implements the Solomon-Bloembergen equations for R1, R2, NOE, and R1rho
from spectral density functions. These serve as:
  - Benchmark targets for validating the numerical Redfield simulator
  - "Standard Redfield" predictions for comparison with CGME
  - Phase 0 feasibility check for the NMR R1rho experiment

All angular frequencies are in rad/s. Rates (R1, R2, R1rho) are in s^-1.
The spectral density J(omega) must include the (2/5) orientational averaging
factor, as in the LipariSzabo class from spectral_densities.py.

Reference: Cavanagh et al., "Protein NMR Spectroscopy", 2nd ed., Ch. 5;
           Farrow et al. (1994) Biochemistry 33, 5984.
"""
import math
from bath import HBAR, KB

# --- Physical constants for 15N-1H NMR ---
MU0_OVER_4PI = 1e-7          # mu_0 / (4*pi) in T*m/A
GAMMA_15N = -2.7116e7         # 15N gyromagnetic ratio, rad/s/T
GAMMA_1H = 2.6752e8           # 1H gyromagnetic ratio, rad/s/T
R_NH_DEFAULT = 1.02e-10       # N-H bond length in meters
DELTA_SIGMA_DEFAULT = -170e-6 # 15N CSA anisotropy (dimensionless, ppm -> fraction)


def dipolar_coupling_constant(gamma_I, gamma_S, r):
    """Dipolar coupling constant d for a heteronuclear spin pair.

    d = (mu_0 / 4*pi) * |gamma_I| * |gamma_S| * hbar / r^3

    Args:
        gamma_I: gyromagnetic ratio of spin I in rad/s/T
        gamma_S: gyromagnetic ratio of spin S in rad/s/T
        r: internuclear distance in meters

    Returns:
        d: dipolar coupling constant in rad/s (always positive)
    """
    return MU0_OVER_4PI * abs(gamma_I) * abs(gamma_S) * HBAR / r**3


def csa_coupling_constant(omega_0, delta_sigma):
    """CSA coupling constant c for a spin with chemical shift anisotropy.

    c = omega_0 * |delta_sigma| / sqrt(3)

    Args:
        omega_0: Larmor frequency magnitude in rad/s (positive)
        delta_sigma: CSA anisotropy (dimensionless, e.g. -170e-6 for 15N)

    Returns:
        c: CSA coupling constant in rad/s (always positive)
    """
    return omega_0 * abs(delta_sigma) / math.sqrt(3)


def compute_R1(d, c, J_fn, omega_N, omega_H):
    """Longitudinal relaxation rate R1 for a 15N-1H pair.

    R1 = (d^2/4) * [J(omega_H - omega_N) + 3*J(omega_N) + 6*J(omega_H + omega_N)]
         + c^2 * J(omega_N)

    All frequencies in rad/s. J_fn must include the (2/5) factor.

    Args:
        d: dipolar coupling constant (rad/s)
        c: CSA coupling constant (rad/s)
        J_fn: callable spectral density J(omega) -> float
        omega_N: 15N Larmor frequency magnitude (rad/s, positive)
        omega_H: 1H Larmor frequency magnitude (rad/s, positive)

    Returns:
        R1 in s^-1
    """
    d2_over_4 = d**2 / 4.0
    R1_dipolar = d2_over_4 * (
        J_fn(omega_H - omega_N) + 3.0 * J_fn(omega_N) + 6.0 * J_fn(omega_H + omega_N)
    )
    R1_csa = c**2 * J_fn(omega_N)
    return R1_dipolar + R1_csa


def compute_R2(d, c, J_fn, omega_N, omega_H):
    """Transverse relaxation rate R2 for a 15N-1H pair.

    R2 = (d^2/8) * [4*J(0) + J(omega_H - omega_N) + 3*J(omega_N)
                     + 6*J(omega_H) + 6*J(omega_H + omega_N)]
         + (c^2/6) * [4*J(0) + 3*J(omega_N)]

    Args:
        d: dipolar coupling constant (rad/s)
        c: CSA coupling constant (rad/s)
        J_fn: callable spectral density J(omega) -> float
        omega_N: 15N Larmor frequency magnitude (rad/s, positive)
        omega_H: 1H Larmor frequency magnitude (rad/s, positive)

    Returns:
        R2 in s^-1
    """
    d2_over_8 = d**2 / 8.0
    R2_dipolar = d2_over_8 * (
        4.0 * J_fn(0.0)
        + J_fn(omega_H - omega_N)
        + 3.0 * J_fn(omega_N)
        + 6.0 * J_fn(omega_H)
        + 6.0 * J_fn(omega_H + omega_N)
    )
    c2_over_6 = c**2 / 6.0
    R2_csa = c2_over_6 * (4.0 * J_fn(0.0) + 3.0 * J_fn(omega_N))
    return R2_dipolar + R2_csa


def compute_NOE(d, J_fn, omega_N, omega_H, R1, gamma_N=GAMMA_15N, gamma_H=GAMMA_1H):
    """Heteronuclear NOE enhancement for a 15N-1H pair.

    sigma_NH = (d^2/4) * [6*J(omega_H + omega_N) - J(omega_H - omega_N)]
    NOE = 1 + (gamma_H / gamma_N) * sigma_NH / R1

    Args:
        d: dipolar coupling constant (rad/s)
        J_fn: callable spectral density J(omega) -> float
        omega_N: 15N Larmor frequency magnitude (rad/s, positive)
        omega_H: 1H Larmor frequency magnitude (rad/s, positive)
        R1: longitudinal relaxation rate (s^-1)
        gamma_N: 15N gyromagnetic ratio (rad/s/T, signed)
        gamma_H: 1H gyromagnetic ratio (rad/s/T, signed)

    Returns:
        NOE enhancement (dimensionless)
    """
    d2_over_4 = d**2 / 4.0
    sigma_NH = d2_over_4 * (
        6.0 * J_fn(omega_H + omega_N) - J_fn(omega_H - omega_N)
    )
    return 1.0 + (gamma_H / gamma_N) * sigma_NH / R1


def compute_R1rho(R1, R2, omega1, delta_omega, Rex=0.0):
    """R1rho relaxation rate in the rotating frame.

    theta = arctan2(omega1, delta_omega)
    R1rho = R1 * cos^2(theta) + (R2 + Rex) * sin^2(theta)

    On resonance (delta_omega = 0): theta = pi/2, R1rho = R2 + Rex.
    Far off resonance (|delta_omega| >> omega1): theta -> 0, R1rho -> R1.

    Args:
        R1: longitudinal relaxation rate (s^-1)
        R2: transverse relaxation rate (s^-1)
        omega1: spin-lock field strength in rad/s
        delta_omega: resonance offset in rad/s
        Rex: exchange contribution to R2 (s^-1, default 0)

    Returns:
        R1rho in s^-1
    """
    theta = math.atan2(omega1, delta_omega)
    cos2 = math.cos(theta) ** 2
    sin2 = math.sin(theta) ** 2
    return R1 * cos2 + (R2 + Rex) * sin2


class NMRRelaxation:
    """Convenience wrapper for computing all 15N NMR relaxation rates.

    Encapsulates physical parameters (field, nuclei, spectral density) and
    computes coupling constants, Larmor frequencies, and all standard
    relaxation rates (R1, R2, NOE, R1rho).

    Usage:
        from spectral_densities import LipariSzabo
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        nmr = NMRRelaxation(B0=11.74, J_fn=J_fn)
        print(nmr.R1(), nmr.R2(), nmr.NOE())
        print(nmr.R1rho(omega1_hz=500.0))
    """

    def __init__(self, B0, J_fn,
                 gamma_N=GAMMA_15N, gamma_H=GAMMA_1H,
                 r_NH=R_NH_DEFAULT, delta_sigma=DELTA_SIGMA_DEFAULT):
        """
        Args:
            B0: static magnetic field in Tesla
            J_fn: callable spectral density J(omega) -> float, with (2/5) factor
            gamma_N: 15N gyromagnetic ratio (rad/s/T, signed)
            gamma_H: 1H gyromagnetic ratio (rad/s/T, signed)
            r_NH: N-H bond length in meters
            delta_sigma: 15N CSA anisotropy (dimensionless, e.g. -170e-6)
        """
        self.B0 = B0
        self.J_fn = J_fn
        self.gamma_N = gamma_N
        self.gamma_H = gamma_H

        # Larmor frequencies (positive magnitudes)
        self.omega_N = abs(gamma_N) * B0
        self.omega_H = abs(gamma_H) * B0

        # Coupling constants
        self.d = dipolar_coupling_constant(gamma_N, gamma_H, r_NH)
        self.c = csa_coupling_constant(self.omega_N, delta_sigma)

    def R1(self):
        """Longitudinal relaxation rate in s^-1."""
        return compute_R1(self.d, self.c, self.J_fn, self.omega_N, self.omega_H)

    def R2(self):
        """Transverse relaxation rate in s^-1."""
        return compute_R2(self.d, self.c, self.J_fn, self.omega_N, self.omega_H)

    def NOE(self):
        """Heteronuclear NOE enhancement (dimensionless)."""
        return compute_NOE(
            self.d, self.J_fn, self.omega_N, self.omega_H,
            self.R1(), self.gamma_N, self.gamma_H,
        )

    def R1rho(self, omega1_hz, offset_hz=0.0, Rex=0.0):
        """R1rho at a given spin-lock power and offset.

        Args:
            omega1_hz: spin-lock field strength in Hz
            offset_hz: resonance offset in Hz (default: on-resonance)
            Rex: exchange contribution in s^-1

        Returns:
            R1rho in s^-1
        """
        omega1 = 2.0 * math.pi * omega1_hz
        delta_omega = 2.0 * math.pi * offset_hz
        return compute_R1rho(self.R1(), self.R2(), omega1, delta_omega, Rex)

    def R1rho_dispersion(self, omega1_array_hz, offset_hz=0.0):
        """Compute R1rho dispersion curve over an array of spin-lock powers.

        Args:
            omega1_array_hz: iterable of spin-lock powers in Hz
            offset_hz: resonance offset in Hz

        Returns:
            list of R1rho values in s^-1
        """
        return [self.R1rho(w1, offset_hz) for w1 in omega1_array_hz]
