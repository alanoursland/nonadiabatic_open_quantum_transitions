"""
Tests for SpectralDensityBath using the spectral density models from spectral_densities.py.

These tests verify:
    - Spectral density model properties (positivity, limits, reductions)
    - Rate function detailed balance and positivity
    - Dissipator trace preservation, Hermiticity, and decay direction
    - CGME interpolation between Redfield and secular limits
    - The two-level secular trap: Redfield != secular for a 2LS with sigma_x coupling
    - Multiple coupling operators
    - Thermal steady state

IMPORTANT: For a two-level system with off-diagonal coupling (sigma_x),
full Redfield and secular Lindblad give DIFFERENT results because R_{01,10}
is nonzero in the full Redfield tensor but zeroed by the secular approximation.
This is correct physics, not a bug. See Gemini's "Two-Level Secular Trap" analysis.
"""
import pytest
import torch
import numpy as np
from bath import SpectralDensityBath, HBAR, KB
from spectral_densities import SimpleLorentzian, LipariSzabo, OhmicDrude


# ============================================================
# Helpers
# ============================================================

def make_two_level_system(omega=5.0):
    """Two-level system with gap omega, already in eigenbasis."""
    energies = torch.tensor([0.0, omega], dtype=torch.float64)
    eigvecs = torch.eye(2, dtype=torch.complex128)
    return energies, eigvecs


def sigma_x():
    return torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)


def sigma_z():
    return torch.tensor([[1, 0], [0, -1]], dtype=torch.complex128)


def make_bath_two_level(omega=5.0, tau_c=1e-9, T=300.0):
    """SpectralDensityBath with Lorentzian J(w), sigma_x coupling, two-level system."""
    J = SimpleLorentzian(tau_c)
    bath = SpectralDensityBath(J, T, [sigma_x()])
    energies, eigvecs = make_two_level_system(omega)
    bath.precompute(energies, eigvecs)
    return bath


def make_three_level_bath(omega1=5.0, omega2=12.0, tau_c=1e-9, T=300.0):
    """Three-level system with nearest-neighbor coupling."""
    energies = torch.tensor([0.0, omega1, omega2], dtype=torch.float64)
    eigvecs = torch.eye(3, dtype=torch.complex128)
    S = torch.tensor(
        [[0.0, 1.0, 0.0],
         [1.0, 0.0, 1.0],
         [0.0, 1.0, 0.0]], dtype=torch.complex128)
    J = SimpleLorentzian(tau_c)
    bath = SpectralDensityBath(J, T, [S])
    bath.precompute(energies, eigvecs)
    return bath


# ============================================================
# Spectral density model tests
# ============================================================

class TestSimpleLorentzian:

    def test_zero_frequency(self):
        J = SimpleLorentzian(tau_c=5e-9)
        expected = (2.0 / 5.0) * 5e-9
        assert abs(J(0.0) - expected) < 1e-25

    def test_high_frequency_falloff(self):
        """J(omega) ~ (2/5) / (omega^2 * tau_c) for omega*tau >> 1."""
        J = SimpleLorentzian(tau_c=1e-9)
        omega = 1e12
        expected = (2.0 / 5.0) * 1e-9 / (1.0 + (1e12 * 1e-9) ** 2)
        assert abs(J(omega) - expected) / expected < 1e-10

    def test_positive(self):
        J = SimpleLorentzian(tau_c=1e-9)
        for omega in [0.0, 1e6, 1e9, 1e12]:
            assert J(omega) >= 0

    def test_peak_at_zero(self):
        """Lorentzian J(omega) is maximum at omega=0."""
        J = SimpleLorentzian(tau_c=1e-9)
        assert J(0.0) > J(1e9)
        assert J(0.0) > J(1e12)


class TestLipariSzabo:

    def test_reduces_to_lorentzian_s2_one(self):
        """With S^2 = 1, Lipari-Szabo = SimpleLorentzian regardless of tau_e."""
        tau_c = 5e-9
        J_simple = SimpleLorentzian(tau_c)
        J_ls = LipariSzabo(tau_c, S2=1.0, tau_e=50e-12)
        for omega in [0.0, 1e8, 1e9, 1e10]:
            assert abs(J_simple(omega) - J_ls(omega)) < 1e-20

    def test_reduces_to_lorentzian_tau_e_zero(self):
        """With tau_e = 0 (rigid limit), second term vanishes."""
        tau_c = 5e-9
        J_simple = SimpleLorentzian(tau_c)
        J_ls = LipariSzabo(tau_c, S2=0.85, tau_e=0.0)
        for omega in [0.0, 1e8, 1e9, 1e10]:
            # Only the S^2 * tau_c term survives
            expected = 0.85 * J_simple(omega)
            assert abs(J_ls(omega) - expected) < 1e-20

    def test_internal_motion_adds_high_freq(self):
        """S^2 < 1 with finite tau_e adds spectral density at high frequencies."""
        tau_c = 5e-9
        J_rigid = SimpleLorentzian(tau_c)
        J_flex = LipariSzabo(tau_c, S2=0.85, tau_e=50e-12)
        omega_high = 1e11
        assert J_flex(omega_high) > J_rigid(omega_high)

    def test_positive(self):
        J = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        for omega in [0.0, 1e6, 1e9, 1e12]:
            assert J(omega) >= 0


class TestOhmicDrude:

    def test_zero_at_zero(self):
        J = OhmicDrude(eta=0.01, lambda_c=1e12)
        assert J(0.0) == 0.0

    def test_positive_for_positive_omega(self):
        J = OhmicDrude(eta=0.01, lambda_c=1e12)
        for omega in [1e6, 1e9, 1e12, 1e15]:
            assert J(omega) > 0

    def test_positive_for_negative_omega(self):
        """Uses |omega|, so J >= 0 for all omega."""
        J = OhmicDrude(eta=0.01, lambda_c=1e12)
        for omega in [-1e6, -1e9, -1e12]:
            assert J(omega) >= 0

    def test_peak_near_cutoff(self):
        """J(omega) peaks at omega = lambda_c."""
        lc = 1e12
        J = OhmicDrude(eta=1.0, lambda_c=lc)
        J_at_peak = J(lc)
        assert J_at_peak > J(0.1 * lc)
        assert J_at_peak > J(10.0 * lc)


# ============================================================
# Rate function tests
# ============================================================

class TestDetailedBalance:

    def test_ratio(self):
        """gamma(omega) / gamma(-omega) = exp(hbar*omega/kT)."""
        omega = 1e10
        T = 300.0
        J = SimpleLorentzian(tau_c=1e-9)
        bath = SpectralDensityBath(J, T, [sigma_x()])

        g_plus = bath._gamma(omega)
        g_minus = bath._gamma(-omega)
        expected = np.exp(HBAR * omega / (KB * T))
        actual = g_plus / g_minus

        assert abs(actual - expected) / expected < 1e-10

    def test_high_temperature_limit(self):
        """At high T, gamma(+omega) ~ gamma(-omega)."""
        omega = 1e6
        T = 300.0
        J = SimpleLorentzian(tau_c=1e-9)
        bath = SpectralDensityBath(J, T, [sigma_x()])

        g_plus = bath._gamma(omega)
        g_minus = bath._gamma(-omega)
        assert abs(g_plus - g_minus) / g_plus < 0.01

    def test_positive_rates(self):
        """gamma(omega) >= 0 for all omega."""
        J = SimpleLorentzian(tau_c=1e-9)
        bath = SpectralDensityBath(J, 300.0, [sigma_x()])
        for omega in [-1e12, -1e9, -1e6, 0.0, 1e6, 1e9, 1e12]:
            assert bath._gamma(omega) >= 0


# ============================================================
# Dissipator property tests — Secular Lindblad
# ============================================================

class TestSecularDissipator:

    def test_trace_preservation(self):
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho)
        assert torch.abs(torch.trace(D)) < 1e-12

    def test_hermiticity(self):
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho)
        assert torch.allclose(D, D.mH, atol=1e-12)

    def test_excited_state_decays(self):
        bath = make_bath_two_level()
        rho = torch.tensor([[0.0, 0.0], [0.0, 1.0]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho)
        assert D[1, 1].real < 0
        assert D[0, 0].real > 0

    def test_three_level_trace(self):
        bath = make_three_level_bath()
        rho = torch.tensor(
            [[0.5, 0.1, 0.05],
             [0.1, 0.3, 0.1],
             [0.05, 0.1, 0.2]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho)
        assert torch.abs(torch.trace(D)) < 1e-12


# ============================================================
# Dissipator property tests — Full Redfield
# ============================================================

class TestRedfieldDissipator:

    def test_trace_preservation(self):
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D = bath.get_redfield_dissipator(rho)
        assert torch.abs(torch.trace(D)) < 1e-12

    def test_hermiticity(self):
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D = bath.get_redfield_dissipator(rho)
        assert torch.allclose(D, D.mH, atol=1e-12)

    def test_three_level_trace(self):
        bath = make_three_level_bath()
        rho = torch.tensor(
            [[0.5, 0.1, 0.05],
             [0.1, 0.3, 0.1],
             [0.05, 0.1, 0.2]], dtype=torch.complex128)
        D = bath.get_redfield_dissipator(rho)
        assert torch.abs(torch.trace(D)) < 1e-12


# ============================================================
# Dissipator property tests — CGME
# ============================================================

class TestCGMEDissipator:

    def test_trace_preservation(self):
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D = bath.get_cgme_dissipator(rho, delta_tau=1.0)
        assert torch.abs(torch.trace(D)) < 1e-12

    def test_hermiticity(self):
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D = bath.get_cgme_dissipator(rho, delta_tau=1.0)
        assert torch.allclose(D, D.mH, atol=1e-12)

    def test_large_tau_recovers_secular(self):
        """CGME with very large Tc -> secular Lindblad."""
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D_sec = bath.get_secular_lindblad_dissipator(rho)
        D_cgme = bath.get_cgme_dissipator(rho, delta_tau=1e10)
        assert torch.allclose(D_sec, D_cgme, atol=1e-8)

    def test_small_tau_recovers_redfield(self):
        """CGME with very small Tc -> full Redfield."""
        bath = make_bath_two_level()
        rho = torch.tensor([[0.3, 0.2 + 0.1j], [0.2 - 0.1j, 0.7]], dtype=torch.complex128)
        D_red = bath.get_redfield_dissipator(rho)
        D_cgme = bath.get_cgme_dissipator(rho, delta_tau=1e-15)
        assert torch.allclose(D_red, D_cgme, atol=1e-8)


# ============================================================
# The Two-Level Secular Trap: Redfield != Secular for a 2LS
# ============================================================

class TestTwoLevelSecularTrap:
    """
    For a two-level system with off-diagonal coupling (sigma_x),
    the full Redfield tensor has R_{01,10} != 0, which couples
    rho_{01} (at freq -omega) to rho_{10} (at freq +omega).
    The secular approximation zeros this term because the frequencies
    don't match. This means full Redfield and secular Lindblad
    give DIFFERENT results for any rho with nonzero coherences.

    This is correct physics, NOT a bug. See:
    - Gemini's "Two-Level Secular Trap" analysis
    - Breuer & Petruccione, §3.3.1
    """

    def test_R_01_10_nonzero_in_full_redfield(self):
        bath = make_bath_two_level(omega=10.0)
        R0110 = bath.R_full[0, 1, 1, 0]
        assert torch.abs(R0110) > 1e-15, \
            f"R_{{01,10}} should be nonzero in full Redfield, got {R0110}"

    def test_R_01_10_zero_in_secular(self):
        bath = make_bath_two_level(omega=10.0)
        R0110_sec = bath.R_secular[0, 1, 1, 0]
        assert torch.abs(R0110_sec) < 1e-15, \
            f"R_{{01,10}} should be zero in secular, got {R0110_sec}"

    def test_population_block_matches(self):
        """The population-to-population block (all at freq 0) is the same
        in both Redfield and secular."""
        bath = make_bath_two_level(omega=10.0)
        for a in [0, 1]:
            for c in [0, 1]:
                R_full = bath.R_full[a, a, c, c]
                R_sec = bath.R_secular[a, a, c, c]
                assert torch.allclose(R_full, R_sec, atol=1e-10), \
                    f"R_{{{a}{a},{c}{c}}}: full={R_full}, secular={R_sec}"

    def test_dissipators_differ_for_coherent_rho(self):
        """With nonzero coherences, Redfield and secular give different D[rho]."""
        bath = make_bath_two_level(omega=10.0)
        rho = torch.tensor([[0.5, 0.3], [0.3, 0.5]], dtype=torch.complex128)
        D_red = bath.get_redfield_dissipator(rho)
        D_sec = bath.get_secular_lindblad_dissipator(rho)
        diff = torch.max(torch.abs(D_red - D_sec)).item()
        assert diff > 1e-10, \
            f"Redfield and secular MUST differ for coherent rho, but max diff = {diff}"

    def test_dissipators_agree_for_diagonal_rho(self):
        """For diagonal rho (no coherences), both agree on the population dynamics."""
        bath = make_bath_two_level(omega=10.0)
        rho = torch.tensor([[0.3, 0.0], [0.0, 0.7]], dtype=torch.complex128)
        D_red = bath.get_redfield_dissipator(rho)
        D_sec = bath.get_secular_lindblad_dissipator(rho)
        # Diagonal elements should match (populations only couple to populations)
        assert torch.allclose(D_red[0, 0], D_sec[0, 0], atol=1e-10)
        assert torch.allclose(D_red[1, 1], D_sec[1, 1], atol=1e-10)


# ============================================================
# Thermal steady state
# ============================================================

class TestThermalSteadyState:

    def test_boltzmann_is_approximate_steady_state(self):
        """At the Boltzmann distribution, the secular dissipator should give
        approximately zero population derivatives."""
        omega = 1e12  # 1 THz — makes hbar*omega/kT nontrivial at 300 K
        T = 300.0
        bath = make_bath_two_level(omega=omega, T=T)
        energies = torch.tensor([0.0, omega], dtype=torch.float64)

        beta = 1.0 / (KB * T)
        Z = np.exp(-HBAR * energies[0].item() * beta) + np.exp(-HBAR * energies[1].item() * beta)
        p0 = np.exp(-HBAR * energies[0].item() * beta) / Z
        p1 = np.exp(-HBAR * energies[1].item() * beta) / Z

        rho_boltz = torch.tensor([[p0, 0.0], [0.0, p1]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho_boltz)

        max_rate = max(abs(bath._gamma(omega)), abs(bath._gamma(-omega)))
        assert torch.abs(D[0, 0].real) < 1e-6 * max(max_rate, 1.0)
        assert torch.abs(D[1, 1].real) < 1e-6 * max(max_rate, 1.0)


# ============================================================
# Multiple coupling operators
# ============================================================

class TestMultipleCouplingOperators:

    def test_trace_and_hermiticity_with_two_ops(self):
        """Bath with sigma_x + sigma_z coupling: basic properties hold."""
        J = SimpleLorentzian(tau_c=1e-9)
        bath = SpectralDensityBath(J, 300.0, [sigma_x(), sigma_z()])
        energies, eigvecs = make_two_level_system(5.0)
        bath.precompute(energies, eigvecs)

        rho = torch.tensor([[0.5, 0.3j], [-0.3j, 0.5]], dtype=torch.complex128)

        for name, method in [('secular', bath.get_secular_lindblad_dissipator),
                             ('redfield', bath.get_redfield_dissipator)]:
            D = method(rho)
            assert torch.abs(torch.trace(D)) < 1e-12, f"{name}: trace != 0"
            assert torch.allclose(D, D.mH, atol=1e-12), f"{name}: not Hermitian"

    def test_sigma_z_adds_pure_dephasing(self):
        """sigma_z coupling contributes pure dephasing: faster coherence decay
        but no population transfer (for diagonal sigma_z in the eigenbasis)."""
        J = SimpleLorentzian(tau_c=1e-9)

        bath_x = SpectralDensityBath(J, 300.0, [sigma_x()])
        bath_xz = SpectralDensityBath(J, 300.0, [sigma_x(), sigma_z()])

        energies, eigvecs = make_two_level_system(5.0)
        bath_x.precompute(energies, eigvecs)
        bath_xz.precompute(energies, eigvecs)

        # sigma_z coupling in the eigenbasis is diagonal, so it only adds
        # to the coherence decay rate, not the population transfer rates.
        # R_{00,11} should be the same (population transfer from sigma_x only)
        R_pop_x = bath_x.R_secular[0, 0, 1, 1]
        R_pop_xz = bath_xz.R_secular[0, 0, 1, 1]
        assert torch.allclose(R_pop_x, R_pop_xz, atol=1e-15), \
            "sigma_z coupling should not affect population transfer rates"

        # R_{01,01} should be larger in magnitude with sigma_z (more dephasing)
        R_coh_x = bath_x.R_secular[0, 1, 0, 1]
        R_coh_xz = bath_xz.R_secular[0, 1, 0, 1]
        assert R_coh_xz.real < R_coh_x.real, \
            "sigma_z coupling should increase coherence decay rate (more negative R_{01,01})"