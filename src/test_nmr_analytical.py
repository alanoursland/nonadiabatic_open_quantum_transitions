"""
Tests for nmr_analytical.py: analytical NMR relaxation rates for 15N.

Validates coupling constants, relaxation rates, and physical consistency
against known properties and published values for backbone amide 15N.

Reference values: ubiquitin-like parameters at 500 MHz (11.74 T),
tau_c = 4.1 ns, S^2 = 0.85, tau_e = 50 ps.
Expected: R1 ~ 1.5-2.5 s^-1, R2 ~ 5-12 s^-1.
"""
import pytest
import math
from spectral_densities import LipariSzabo, SimpleLorentzian
from nmr_analytical import (
    MU0_OVER_4PI, GAMMA_15N, GAMMA_1H, R_NH_DEFAULT, DELTA_SIGMA_DEFAULT,
    dipolar_coupling_constant, csa_coupling_constant,
    compute_R1, compute_R2, compute_NOE, compute_R1rho,
    NMRRelaxation,
)
from bath import HBAR


# --- Fixtures ---

@pytest.fixture
def ubiquitin_params():
    """Ubiquitin-like backbone 15N at 500 MHz / 11.74 T."""
    J_fn = LipariSzabo(tau_c=4.1e-9, S2=0.85, tau_e=50e-12)
    return NMRRelaxation(B0=11.74, J_fn=J_fn)


@pytest.fixture
def phase0_params():
    """Phase 0 feasibility parameters: tau_c=5ns, S^2=0.85, 11.7 T."""
    J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
    return NMRRelaxation(B0=11.7, J_fn=J_fn)


# --- Test classes ---

class TestCouplingConstants:
    """Verify dipolar and CSA coupling constant magnitudes."""

    def test_dipolar_coupling_magnitude(self):
        """d for 15N-1H at 1.02 A should be ~72,000 rad/s."""
        d = dipolar_coupling_constant(GAMMA_15N, GAMMA_1H, R_NH_DEFAULT)
        assert d > 0, "d must be positive (absolute values used)"
        assert 60000 < d < 90000, f"d = {d:.1f} rad/s, expected ~72,000"

    def test_dipolar_coupling_positive(self):
        """d is always positive regardless of gamma signs."""
        d = dipolar_coupling_constant(-2.7116e7, 2.6752e8, 1.02e-10)
        assert d > 0

    def test_dipolar_coupling_distance_dependence(self):
        """d scales as 1/r^3: doubling r should reduce d by factor 8."""
        d1 = dipolar_coupling_constant(GAMMA_15N, GAMMA_1H, 1.0e-10)
        d2 = dipolar_coupling_constant(GAMMA_15N, GAMMA_1H, 2.0e-10)
        assert abs(d1 / d2 - 8.0) < 0.01

    def test_csa_coupling_magnitude(self):
        """c for 15N at 11.74 T with delta_sigma = -170 ppm should be ~31,000 rad/s."""
        omega_N = abs(GAMMA_15N) * 11.74
        c = csa_coupling_constant(omega_N, DELTA_SIGMA_DEFAULT)
        assert c > 0, "c must be positive"
        assert 25000 < c < 40000, f"c = {c:.1f} rad/s, expected ~31,000"

    def test_csa_coupling_field_dependence(self):
        """c is proportional to B0 (through omega_0)."""
        c1 = csa_coupling_constant(abs(GAMMA_15N) * 11.74, -170e-6)
        c2 = csa_coupling_constant(abs(GAMMA_15N) * 23.48, -170e-6)
        assert abs(c2 / c1 - 2.0) < 1e-10


class TestSpectralDensitySanity:
    """Verify J(omega) ordering for motional narrowing regime."""

    def test_monotone_decrease(self):
        """J(0) > J(omega_N) > J(omega_H) for tau_c ~ ns."""
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        omega_N = abs(GAMMA_15N) * 11.7
        omega_H = abs(GAMMA_1H) * 11.7
        assert J_fn(0.0) > J_fn(omega_N) > J_fn(omega_H) > 0

    def test_positive_all_frequencies(self):
        """J(omega) >= 0 for all omega >= 0."""
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        for omega in [0.0, 1e6, 1e8, 1e9, 1e10, 1e12]:
            assert J_fn(omega) >= 0


class TestRelaxationRateRanges:
    """R1 and R2 must be in published ranges for backbone 15N."""

    def test_R1_range(self, ubiquitin_params):
        """R1 should be 1-3 s^-1 at 500 MHz for tau_c ~ 4 ns."""
        R1 = ubiquitin_params.R1()
        assert 0.5 < R1 < 5.0, f"R1 = {R1:.3f} s^-1, expected 1-3"

    def test_R2_range(self, ubiquitin_params):
        """R2 should be 5-15 s^-1 at 500 MHz for tau_c ~ 4 ns."""
        R2 = ubiquitin_params.R2()
        assert 3.0 < R2 < 20.0, f"R2 = {R2:.3f} s^-1, expected 5-15"

    def test_R2_greater_than_R1(self, ubiquitin_params):
        """R2 > R1 always in the slow-tumbling regime (omega*tau_c > 1)."""
        assert ubiquitin_params.R2() > ubiquitin_params.R1()

    def test_R1_R2_positive(self, ubiquitin_params):
        """Both rates must be positive."""
        assert ubiquitin_params.R1() > 0
        assert ubiquitin_params.R2() > 0

    @pytest.mark.parametrize("tau_c,S2", [
        (1e-9, 0.5), (1e-9, 0.85), (1e-9, 1.0),
        (5e-9, 0.5), (5e-9, 0.85), (5e-9, 1.0),
        (20e-9, 0.5), (20e-9, 0.85), (20e-9, 1.0),
    ])
    def test_R1_R2_positive_various_conditions(self, tau_c, S2):
        """R1, R2 > 0 for various physically valid parameters."""
        J_fn = LipariSzabo(tau_c=tau_c, S2=S2, tau_e=50e-12)
        nmr = NMRRelaxation(B0=11.74, J_fn=J_fn)
        assert nmr.R1() > 0, f"R1 < 0 for tau_c={tau_c}, S2={S2}"
        assert nmr.R2() > 0, f"R2 < 0 for tau_c={tau_c}, S2={S2}"


class TestR1rhoLimits:
    """R1rho must interpolate between R1 and R2 with correct limits."""

    def test_on_resonance_equals_R2(self, ubiquitin_params):
        """On resonance (delta_omega = 0): R1rho = R2 (no Rex)."""
        R2 = ubiquitin_params.R2()
        R1rho = ubiquitin_params.R1rho(omega1_hz=500.0, offset_hz=0.0)
        assert abs(R1rho - R2) < 1e-10, f"R1rho = {R1rho}, R2 = {R2}"

    def test_large_offset_approaches_R1(self, ubiquitin_params):
        """Large offset (|delta_omega| >> omega1): R1rho -> R1."""
        R1 = ubiquitin_params.R1()
        R1rho = ubiquitin_params.R1rho(omega1_hz=100.0, offset_hz=1e6)
        assert abs(R1rho - R1) / R1 < 0.01, \
            f"R1rho = {R1rho}, R1 = {R1}, should converge"

    def test_monotonic_offset_dependence(self, ubiquitin_params):
        """R1rho decreases monotonically as offset increases (R2 > R1)."""
        R1rho_values = [
            ubiquitin_params.R1rho(omega1_hz=500.0, offset_hz=offset)
            for offset in [0, 100, 500, 1000, 5000, 50000]
        ]
        for i in range(len(R1rho_values) - 1):
            assert R1rho_values[i] >= R1rho_values[i + 1] - 1e-15, \
                f"R1rho not monotonically decreasing with offset"

    def test_R1rho_between_R1_and_R2(self, ubiquitin_params):
        """R1rho must be between R1 and R2 for any omega1, offset (no Rex)."""
        R1 = ubiquitin_params.R1()
        R2 = ubiquitin_params.R2()
        for omega1_hz in [25, 100, 500, 1000]:
            for offset_hz in [0, 100, 1000, 10000]:
                R1rho = ubiquitin_params.R1rho(omega1_hz, offset_hz)
                assert R1 - 1e-10 <= R1rho <= R2 + 1e-10, \
                    f"R1rho={R1rho:.4f} outside [{R1:.4f}, {R2:.4f}]"

    def test_R1rho_dispersion_on_resonance_flat(self, ubiquitin_params):
        """On resonance, R1rho = R2 independent of omega1 (no Rex)."""
        R2 = ubiquitin_params.R2()
        omega1_array = [25, 50, 100, 250, 500, 1000]
        R1rho_array = ubiquitin_params.R1rho_dispersion(omega1_array, offset_hz=0.0)
        for R1rho in R1rho_array:
            assert abs(R1rho - R2) < 1e-10


class TestNOE:
    """NOE enhancement must be in physical range for 15N."""

    def test_NOE_range(self, ubiquitin_params):
        """15N NOE should be between -5 and +1 (gamma_N < 0)."""
        noe = ubiquitin_params.NOE()
        assert -5.0 < noe < 1.5, f"NOE = {noe}"

    def test_NOE_less_than_one_for_fast_tumbling(self):
        """For fast tumbling (tau_c ~ 1 ns), 15N NOE < 1."""
        J_fn = LipariSzabo(tau_c=1e-9, S2=0.85, tau_e=50e-12)
        nmr = NMRRelaxation(B0=11.74, J_fn=J_fn)
        noe = nmr.NOE()
        assert noe < 1.0, f"NOE = {noe}, expected < 1 for 15N fast tumbling"


class TestDipolarVsCSA:
    """Verify relative contributions of dipolar and CSA to relaxation."""

    def test_dipolar_dominates_R1_at_low_field(self):
        """At low field (~7 T), dipolar contribution to R1 should dominate."""
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        nmr = NMRRelaxation(B0=7.05, J_fn=J_fn)
        d, c = nmr.d, nmr.c
        omega_N, omega_H = nmr.omega_N, nmr.omega_H
        R1_dip = (d**2 / 4) * (
            J_fn(omega_H - omega_N) + 3 * J_fn(omega_N) + 6 * J_fn(omega_H + omega_N)
        )
        R1_csa = c**2 * J_fn(omega_N)
        assert R1_dip > R1_csa, \
            f"Dipolar ({R1_dip:.4f}) should dominate R1 over CSA ({R1_csa:.4f})"

    def test_csa_scales_with_field(self):
        """CSA coupling c is proportional to B0."""
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        nmr_low = NMRRelaxation(B0=7.05, J_fn=J_fn)
        nmr_high = NMRRelaxation(B0=18.8, J_fn=J_fn)
        assert abs(nmr_high.c / nmr_low.c - 18.8 / 7.05) < 0.01


class TestPublishedValues:
    """Cross-check against known values for backbone 15N relaxation.

    Reference: typical backbone amide 15N at 500 MHz (11.74 T) for a
    small globular protein (MW ~ 8.5 kDa, tau_c ~ 4 ns, S^2 ~ 0.85).
    """

    def test_ubiquitin_R1(self, ubiquitin_params):
        """R1 for ubiquitin-like parameters at 500 MHz."""
        R1 = ubiquitin_params.R1()
        assert 1.0 < R1 < 3.5, f"R1 = {R1:.3f} s^-1"

    def test_ubiquitin_R2(self, ubiquitin_params):
        """R2 for ubiquitin-like parameters at 500 MHz."""
        R2 = ubiquitin_params.R2()
        assert 4.0 < R2 < 18.0, f"R2 = {R2:.3f} s^-1"

    def test_R2_over_R1_ratio(self, ubiquitin_params):
        """R2/R1 ratio ~ 3-8 for small proteins at 500 MHz."""
        ratio = ubiquitin_params.R2() / ubiquitin_params.R1()
        assert 2.0 < ratio < 10.0, f"R2/R1 = {ratio:.2f}"


class TestNMRRelaxationClass:
    """Test the convenience class interface."""

    def test_construction(self):
        """Class should construct without error."""
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        nmr = NMRRelaxation(B0=11.74, J_fn=J_fn)
        assert nmr.omega_N > 0
        assert nmr.omega_H > 0
        assert nmr.d > 0
        assert nmr.c > 0

    def test_larmor_frequencies(self):
        """Larmor frequencies should match |gamma| * B0."""
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        nmr = NMRRelaxation(B0=11.74, J_fn=J_fn)
        assert abs(nmr.omega_N - abs(GAMMA_15N) * 11.74) < 1.0
        assert abs(nmr.omega_H - abs(GAMMA_1H) * 11.74) < 1.0

    def test_R1rho_dispersion_length(self):
        """R1rho_dispersion returns correct number of values."""
        J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        nmr = NMRRelaxation(B0=11.74, J_fn=J_fn)
        omega1_array = [25, 50, 100, 250, 500, 1000]
        result = nmr.R1rho_dispersion(omega1_array)
        assert len(result) == 6

    def test_rigid_tumbler_matches_simple_lorentzian(self):
        """With S^2=1 and tau_e=0, results should match SimpleLorentzian."""
        tau_c = 5e-9
        J_ls = LipariSzabo(tau_c, S2=1.0, tau_e=0.0)
        J_simple = SimpleLorentzian(tau_c)
        nmr_ls = NMRRelaxation(B0=11.74, J_fn=J_ls)
        nmr_simple = NMRRelaxation(B0=11.74, J_fn=J_simple)
        assert abs(nmr_ls.R1() - nmr_simple.R1()) < 1e-10
        assert abs(nmr_ls.R2() - nmr_simple.R2()) < 1e-10


class TestPhase0Parameters:
    """Compute relaxation rates for Phase 0 feasibility parameters.

    tau_c = 5 ns, S^2 = 0.85, tau_e = 50 ps, B0 = 11.7 T.
    """

    def test_R1_phase0(self, phase0_params):
        R1 = phase0_params.R1()
        assert 0.5 < R1 < 4.0, f"Phase 0 R1 = {R1:.3f}"

    def test_R2_phase0(self, phase0_params):
        R2 = phase0_params.R2()
        assert 3.0 < R2 < 20.0, f"Phase 0 R2 = {R2:.3f}"

    def test_R1rho_on_resonance_dispersion_flat(self, phase0_params):
        """On-resonance R1rho dispersion should be flat (= R2) for all omega1."""
        R2 = phase0_params.R2()
        for omega1_hz in [25, 50, 100, 250, 500, 1000]:
            R1rho = phase0_params.R1rho(omega1_hz, offset_hz=0.0)
            assert abs(R1rho - R2) < 1e-10, \
                f"omega1={omega1_hz} Hz: R1rho={R1rho:.6f}, R2={R2:.6f}"
