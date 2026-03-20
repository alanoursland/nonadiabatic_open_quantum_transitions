"""Tests for CollisionalLorentzian spectral density."""
import math
import pytest

from spectral_densities import CollisionalLorentzian


@pytest.fixture
def co_n2_collision():
    """CO-N2 collisional spectral density at 296 K.

    gamma_L ~ 0.065 cm^-1/atm -> gamma_0 ~ 1.22e10 rad/s
    tau_c ~ 1e-12 s (picosecond collision duration)
    """
    CM_TO_RAD_S = 2.0 * math.pi * 2.99792458e10
    gamma_0 = 0.065 * CM_TO_RAD_S  # ~1.22e10 rad/s
    return CollisionalLorentzian(gamma_0=gamma_0, tau_c=1e-12)


class TestCollisionalLorentzian:

    def test_positive(self, co_n2_collision):
        """J(omega) >= 0 for all omega."""
        J = co_n2_collision
        for omega in [0, 1e8, 1e10, 1e11, 1e12, 1e13]:
            assert J(omega) >= 0.0

    def test_peak_at_zero(self, co_n2_collision):
        """J(0) > J(omega) for any omega > 0."""
        J = co_n2_collision
        J0 = J(0.0)
        for omega in [1e8, 1e10, 1e12]:
            assert J0 > J(omega)

    def test_zero_frequency_value(self, co_n2_collision):
        """J(0) = gamma_0 * tau_c."""
        J = co_n2_collision
        expected = J.gamma_0 * J.tau_c
        assert J(0.0) == pytest.approx(expected, rel=1e-10)

    def test_high_frequency_falloff(self):
        """J(omega) ~ gamma_0 / (omega^2 * tau_c) for omega * tau_c >> 1."""
        gamma_0 = 1e10
        tau_c = 1e-12
        J = CollisionalLorentzian(gamma_0=gamma_0, tau_c=tau_c)
        omega = 1e14  # omega * tau_c = 100 >> 1
        expected = gamma_0 / (omega ** 2 * tau_c)
        assert J(omega) == pytest.approx(expected, rel=1e-2)

    def test_no_nmr_prefactor(self):
        """CollisionalLorentzian has no (2/5) factor (not NMR)."""
        gamma_0 = 1.0
        tau_c = 1.0
        J = CollisionalLorentzian(gamma_0=gamma_0, tau_c=tau_c)
        # J(0) = gamma_0 * tau_c = 1.0, not (2/5) * 1.0
        assert J(0.0) == pytest.approx(1.0, rel=1e-10)
        assert J(0.0) != pytest.approx(2.0 / 5.0, rel=1e-2)

    def test_callable_interface(self, co_n2_collision):
        """Can be called with a float and returns a float."""
        result = co_n2_collision(1e10)
        assert isinstance(result, float)
