"""
Spectral density models J(omega) for NMR and condensed-phase quantum dynamics.

A spectral density J(omega) characterizes the frequency-dependent coupling
strength between a quantum system and its thermal environment. It enters the
Redfield tensor through the rate function gamma(omega), which is computed
from J and the Bose-Einstein distribution in SpectralDensityBath.

All spectral densities are callable: J(omega) -> float.
Convention: omega in rad/s, J(omega) in seconds (or rad^{-1} depending on
the prefactor convention — the (2/5) in the NMR forms absorbs the geometric
factors from the dipolar/CSA interaction).

Physical requirement: J(omega) >= 0 for all omega >= 0.
"""


class SimpleLorentzian:
    """
    Single Lorentzian spectral density (rigid isotropic tumbling):

        J(omega) = (2/5) * tau_c / (1 + omega^2 * tau_c^2)

    This is the rigid-body limit (S^2 = 1) of the Lipari-Szabo model.
    It corresponds to exponentially decaying bath correlation functions
    with a single correlation time tau_c.

    Args:
        tau_c: rotational correlation time in seconds
    """

    def __init__(self, tau_c):
        self.tau_c = tau_c

    def __call__(self, omega):
        return (2.0 / 5.0) * self.tau_c / (1.0 + omega ** 2 * self.tau_c ** 2)


class LipariSzabo:
    """
    Lipari-Szabo model-free spectral density:

        J(omega) = (2/5) * [ S^2 * tau_c / (1 + omega^2 * tau_c^2)
                           + (1 - S^2) * tau' / (1 + omega^2 * tau'^2) ]

    where 1/tau' = 1/tau_c + 1/tau_e is the effective internal correlation time.

    In the rigid limit (S^2 = 1 or tau_e = 0), the second term vanishes
    and this reduces to SimpleLorentzian.

    Reference: Lipari & Szabo, JACS 104, 4546 (1982).

    Args:
        tau_c: overall rotational correlation time in seconds
        S2: generalized order parameter (0 <= S^2 <= 1)
        tau_e: internal motion correlation time in seconds (0 = rigid limit)
    """

    def __init__(self, tau_c, S2, tau_e=0.0):
        self.tau_c = tau_c
        self.S2 = S2
        self.tau_e = tau_e

        # Effective internal correlation time
        if tau_e > 0 and tau_c > 0:
            self.tau_eff = 1.0 / (1.0 / tau_c + 1.0 / tau_e)
        else:
            self.tau_eff = 0.0

    def __call__(self, omega):
        w2 = omega ** 2
        term1 = self.S2 * self.tau_c / (1.0 + w2 * self.tau_c ** 2)
        term2 = 0.0
        if self.tau_eff > 0:
            term2 = (1.0 - self.S2) * self.tau_eff / (1.0 + w2 * self.tau_eff ** 2)
        return (2.0 / 5.0) * (term1 + term2)


class OhmicDrude:
    """
    Ohmic spectral density with Drude (Lorentzian) cutoff:

        J(omega) = eta * |omega| * lambda_c / (omega^2 + lambda_c^2)

    Common in condensed-phase quantum dynamics (spin-boson model,
    excitation energy transfer). The |omega| ensures J >= 0 for all omega.

    Note: some references define this without the absolute value, using
    the convention that J is only evaluated for omega > 0. We use |omega|
    so the function is safe to call at any frequency.

    Args:
        eta: dimensionless coupling strength
        lambda_c: cutoff frequency in rad/s
    """

    def __init__(self, eta, lambda_c):
        self.eta = eta
        self.lambda_c = lambda_c

    def __call__(self, omega):
        return self.eta * abs(omega) * self.lambda_c / (omega ** 2 + self.lambda_c ** 2)