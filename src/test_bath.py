"""
Tests for bath.py: ThermalBath (backward compat) + SpectralDensityBath (microscopic).

The SpectralDensityBath tests verify the Redfield tensor R_{abcd} and its
secular and coarse-grained variants against known properties:

    1. Trace preservation:   Tr[D[rho]] = 0 for any rho
    2. Hermiticity:          D[rho] is Hermitian when rho is Hermitian
    3. Detailed balance:     gamma(omega)/gamma(-omega) = exp(hbar*omega/kT)
    4. Secular vs Redfield:  they differ for a 2LS (R_{01,10} != 0 in full Redfield)
    5. Population equations: secular R decouples populations from coherences
    6. CGME interpolation:   CGME -> Redfield as Tc -> 0, CGME -> secular as Tc -> inf
    7. Thermal steady state: the secular dissipator drives toward Boltzmann populations
"""
import pytest
import torch
import numpy as np
from system import System
from bath import (
    ThermalBath, SpectralDensityBath,
    DebyeSpectralDensity, LipariSzaboSpectralDensity,
    HBAR, KB,
)


# ============================================================
# Backward-compatibility tests for ThermalBath (original tests)
# ============================================================

@pytest.fixture
def qubit_system():
    return System([0.0, 5.0], [[0.0, 1.0], [1.0, 0.0]])


@pytest.fixture
def thermal_bath():
    return ThermalBath(
        temperature=0.015,
        t1_times={(1, 0): 50.0},
        t2_times={(1, 0): 30.0},
    )


@pytest.fixture
def density_matrix():
    return torch.tensor([[0.5, 0.5], [0.5, 0.5]], dtype=torch.complex128)


def test_dissipator_trace_preservation(qubit_system, thermal_bath, density_matrix):
    D_rho = thermal_bath.get_lindblad_dissipator(qubit_system, density_matrix)
    assert torch.abs(torch.trace(D_rho)) < 1e-10


def test_dissipator_hermiticity(qubit_system, thermal_bath, density_matrix):
    D_rho = thermal_bath.get_lindblad_dissipator(qubit_system, density_matrix)
    assert torch.allclose(D_rho, D_rho.mH, atol=1e-10)


def test_population_decay(qubit_system, thermal_bath):
    rho = torch.tensor([[0.0, 0.0], [0.0, 1.0]], dtype=torch.complex128)
    D_rho = thermal_bath.get_lindblad_dissipator(qubit_system, rho)
    assert torch.real(D_rho[1, 1]) < 0
    assert torch.real(D_rho[0, 0]) > 0


# ============================================================
# SpectralDensityBath tests
# ============================================================

def _make_two_level_bath(omega=10.0, eta=0.01, tau_c=1.0, T_K=300.0):
    """Helper: two-level system with Debye bath coupled via sigma_x."""
    energies = torch.tensor([0.0, omega], dtype=torch.float64)
    eigvecs = torch.eye(2, dtype=torch.complex128)
    S = torch.tensor([[0.0, 1.0], [1.0, 0.0]], dtype=torch.complex128)

    J = DebyeSpectralDensity(eta=eta, tau_c=tau_c)
    bath = SpectralDensityBath(J, T_K, [S])
    bath.precompute(energies, eigvecs)
    return bath, energies, eigvecs


def _make_three_level_bath(omega1=5.0, omega2=12.0, eta=0.01, tau_c=1.0, T_K=300.0):
    """Helper: three-level system with Debye bath."""
    energies = torch.tensor([0.0, omega1, omega2], dtype=torch.float64)
    eigvecs = torch.eye(3, dtype=torch.complex128)
    S = torch.tensor(
        [[0.0, 1.0, 0.0],
         [1.0, 0.0, 1.0],
         [0.0, 1.0, 0.0]], dtype=torch.complex128,
    )
    J = DebyeSpectralDensity(eta=eta, tau_c=tau_c)
    bath = SpectralDensityBath(J, T_K, [S])
    bath.precompute(energies, eigvecs)
    return bath, energies, eigvecs


class TestRedfieldTracePreservation:
    """Tr[D[rho]] = 0 for all three dissipator types and arbitrary rho."""

    @pytest.mark.parametrize("rho_data", [
        [[1.0, 0.0], [0.0, 0.0]],
        [[0.0, 0.0], [0.0, 1.0]],
        [[0.5, 0.3 + 0.1j], [0.3 - 0.1j, 0.5]],
        [[0.7, 0.4j], [-0.4j, 0.3]],
    ])
    def test_redfield_trace(self, rho_data):
        bath, _, _ = _make_two_level_bath()
        rho = torch.tensor(rho_data, dtype=torch.complex128)
        D = bath.get_redfield_dissipator(rho)
        assert torch.abs(torch.trace(D)) < 1e-10, f"Tr[D] = {torch.trace(D)}"

    @pytest.mark.parametrize("rho_data", [
        [[1.0, 0.0], [0.0, 0.0]],
        [[0.5, 0.3 + 0.1j], [0.3 - 0.1j, 0.5]],
    ])
    def test_secular_trace(self, rho_data):
        bath, _, _ = _make_two_level_bath()
        rho = torch.tensor(rho_data, dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho)
        assert torch.abs(torch.trace(D)) < 1e-10

    @pytest.mark.parametrize("rho_data", [
        [[1.0, 0.0], [0.0, 0.0]],
        [[0.5, 0.3 + 0.1j], [0.3 - 0.1j, 0.5]],
    ])
    def test_cgme_trace(self, rho_data):
        bath, _, _ = _make_two_level_bath()
        rho = torch.tensor(rho_data, dtype=torch.complex128)
        D = bath.get_cgme_dissipator(rho, delta_tau=0.5)
        assert torch.abs(torch.trace(D)) < 1e-10

    def test_three_level_trace(self):
        bath, _, _ = _make_three_level_bath()
        rho = torch.tensor(
            [[0.5, 0.1, 0.05],
             [0.1, 0.3, 0.1],
             [0.05, 0.1, 0.2]], dtype=torch.complex128,
        )
        for method in ['redfield', 'secular', 'cgme']:
            if method == 'redfield':
                D = bath.get_redfield_dissipator(rho)
            elif method == 'secular':
                D = bath.get_secular_lindblad_dissipator(rho)
            else:
                D = bath.get_cgme_dissipator(rho, delta_tau=0.5)
            scale = torch.max(torch.abs(D)).item()
            assert torch.abs(torch.trace(D)) < 1e-10 * max(scale, 1.0), \
                f"{method}: Tr[D] = {torch.trace(D)}, scale = {scale}"

class TestRedfieldHermiticity:
    """D[rho] must be Hermitian when rho is Hermitian."""

    def test_redfield_hermitian(self):
        bath, _, _ = _make_two_level_bath()
        rho = torch.tensor([[0.6, 0.2 + 0.1j], [0.2 - 0.1j, 0.4]], dtype=torch.complex128)
        D = bath.get_redfield_dissipator(rho)
        assert torch.allclose(D, D.mH, atol=1e-10)

    def test_secular_hermitian(self):
        bath, _, _ = _make_two_level_bath()
        rho = torch.tensor([[0.6, 0.2 + 0.1j], [0.2 - 0.1j, 0.4]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho)
        assert torch.allclose(D, D.mH, atol=1e-10)

    def test_three_level_hermitian(self):
        bath, _, _ = _make_three_level_bath()
        rho = torch.tensor(
            [[0.5, 0.1 + 0.05j, 0.02],
             [0.1 - 0.05j, 0.3, 0.1j],
             [0.02, -0.1j, 0.2]], dtype=torch.complex128,
        )
        D = bath.get_redfield_dissipator(rho)
        assert torch.allclose(D, D.mH, atol=1e-10)


class TestDetailedBalance:
    """gamma(omega) / gamma(-omega) = exp(hbar * omega / kT)."""

    def test_detailed_balance_ratio(self):
        T_K = 300.0
        omega = 1e12  # 1 THz — a frequency where n(omega) is moderate
        J = DebyeSpectralDensity(eta=0.01, tau_c=1e-12)
        bath = SpectralDensityBath(J, T_K, [torch.eye(2, dtype=torch.complex128)])

        g_plus = bath._gamma(omega)
        g_minus = bath._gamma(-omega)
        expected_ratio = np.exp(HBAR * omega / (KB * T_K))

        assert g_minus > 0, "Absorption rate must be positive"
        assert g_plus > 0, "Emission rate must be positive"
        actual_ratio = g_plus / g_minus
        assert abs(actual_ratio - expected_ratio) / expected_ratio < 1e-6

    def test_high_temperature_limit(self):
        """At high T, gamma(omega) ~ gamma(-omega) ~ J(omega) * kT / (hbar * omega)."""
        T_K = 1e6  # very high temperature
        omega = 1.0
        J = DebyeSpectralDensity(eta=0.01, tau_c=1.0)
        bath = SpectralDensityBath(J, T_K, [torch.eye(2, dtype=torch.complex128)])

        g_plus = bath._gamma(omega)
        g_minus = bath._gamma(-omega)
        # At high T, n(omega) >> 1 and n+1 ~ n, so g+ ~ g-
        assert abs(g_plus - g_minus) / g_plus < 1e-4


class TestSecularVsRedfield:
    """Full Redfield and secular differ for a 2LS: R_{01,10} != 0 in full Redfield."""

    def test_R_01_10_nonzero_in_full(self):
        bath, _, _ = _make_two_level_bath(omega=10.0)
        # R_{01,10} couples rho_{01} (freq -omega) to rho_{10} (freq +omega)
        R0110 = bath.R_full[0, 1, 1, 0]
        assert torch.abs(R0110) > 1e-15, f"R_{{01,10}} should be nonzero, got {R0110}"

    def test_R_01_10_zero_in_secular(self):
        bath, _, _ = _make_two_level_bath(omega=10.0)
        R0110_sec = bath.R_secular[0, 1, 1, 0]
        assert torch.abs(R0110_sec) < 1e-15, f"Secular R_{{01,10}} should be zero, got {R0110_sec}"

    def test_population_block_matches(self):
        """For a 2LS with sigma_x coupling, the population-to-population block
        of the Redfield tensor should be the same in both full and secular,
        since populations all evolve at omega=0."""
        bath, _, _ = _make_two_level_bath(omega=10.0)
        # Population block: R_{aa,cc} for a,c in {0,1}
        for a in [0, 1]:
            for c in [0, 1]:
                R_full = bath.R_full[a, a, c, c]
                R_sec = bath.R_secular[a, a, c, c]
                assert torch.allclose(R_full, R_sec, atol=1e-10), \
                    f"R_{{{a}{a},{c}{c}}}: full={R_full}, secular={R_sec}"

    def test_dissipators_differ_on_coherent_state(self):
        """Redfield and secular give different D[rho] for rho with coherences."""
        bath, _, _ = _make_two_level_bath(omega=10.0)
        rho = torch.tensor([[0.5, 0.3], [0.3, 0.5]], dtype=torch.complex128)
        D_full = bath.get_redfield_dissipator(rho)
        D_sec = bath.get_secular_lindblad_dissipator(rho)
        # They should differ (non-secular terms contribute)
        diff = torch.max(torch.abs(D_full - D_sec)).item()
        assert diff > 1e-10, f"Redfield and secular should differ, max diff = {diff}"

    def test_dissipators_agree_on_diagonal(self):
        """For a diagonal rho (no coherences), full and secular should agree
        on the diagonal of D[rho], since coherence coupling terms don't contribute."""
        bath, _, _ = _make_two_level_bath(omega=10.0)
        rho = torch.tensor([[0.3, 0.0], [0.0, 0.7]], dtype=torch.complex128)
        D_full = bath.get_redfield_dissipator(rho)
        D_sec = bath.get_secular_lindblad_dissipator(rho)
        assert torch.allclose(D_full[0, 0], D_sec[0, 0], atol=1e-10)
        assert torch.allclose(D_full[1, 1], D_sec[1, 1], atol=1e-10)


class TestCGMEInterpolation:
    """CGME interpolates between full Redfield (Tc->0) and secular (Tc->inf)."""

    def test_small_tc_approaches_redfield(self):
        bath, _, _ = _make_two_level_bath(omega=10.0)
        rho = torch.tensor([[0.5, 0.3], [0.3, 0.5]], dtype=torch.complex128)
        D_full = bath.get_redfield_dissipator(rho)
        D_cgme = bath.get_cgme_dissipator(rho, delta_tau=1e-10)
        assert torch.allclose(D_full, D_cgme, atol=1e-8), \
            f"CGME(Tc->0) should match Redfield, max diff = {torch.max(torch.abs(D_full - D_cgme))}"

    def test_large_tc_approaches_secular(self):
        bath, _, _ = _make_two_level_bath(omega=10.0)
        rho = torch.tensor([[0.5, 0.3], [0.3, 0.5]], dtype=torch.complex128)
        D_sec = bath.get_secular_lindblad_dissipator(rho)
        D_cgme = bath.get_cgme_dissipator(rho, delta_tau=1e6)
        assert torch.allclose(D_sec, D_cgme, atol=1e-8), \
            f"CGME(Tc->inf) should match secular, max diff = {torch.max(torch.abs(D_sec - D_cgme))}"


class TestPopulationDecayDirection:
    """Excited-state population should decrease, ground state should increase."""

    def test_two_level_decay(self):
        bath, _, _ = _make_two_level_bath(omega=10.0, T_K=300.0)
        rho = torch.tensor([[0.0, 0.0], [0.0, 1.0]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho)
        assert torch.real(D[1, 1]) < 0, "Excited state should decay"
        assert torch.real(D[0, 0]) > 0, "Ground state should gain population"

    def test_three_level_decay(self):
        bath, _, _ = _make_three_level_bath()
        rho = torch.zeros(3, 3, dtype=torch.complex128)
        rho[2, 2] = 1.0  # all population in highest state
        D = bath.get_secular_lindblad_dissipator(rho)
        assert torch.real(D[2, 2]) < 0, "Highest state should decay"


class TestThermalSteadyState:
    """The secular dissipator should drive populations toward Boltzmann."""

    def test_boltzmann_is_approximate_steady_state(self):
        """At the Boltzmann distribution, population derivatives should be small."""
        omega = 1e12  # 1 THz
        T_K = 300.0
        bath, energies, _ = _make_two_level_bath(omega=omega, T_K=T_K)

        # Boltzmann populations
        beta = 1.0 / (KB * T_K)
        Z = np.exp(-HBAR * energies[0].item() * beta) + np.exp(-HBAR * energies[1].item() * beta)
        p0 = np.exp(-HBAR * energies[0].item() * beta) / Z
        p1 = np.exp(-HBAR * energies[1].item() * beta) / Z

        rho_boltz = torch.tensor([[p0, 0.0], [0.0, p1]], dtype=torch.complex128)
        D = bath.get_secular_lindblad_dissipator(rho_boltz)

        # Population derivatives should be approximately zero
        max_rate = max(abs(bath._gamma(omega)), abs(bath._gamma(-omega)))
        assert torch.abs(torch.real(D[0, 0])) < 1e-6 * max(max_rate, 1.0)
        assert torch.abs(torch.real(D[1, 1])) < 1e-6 * max(max_rate, 1.0)


class TestSpectralDensityModels:
    """Verify spectral density functions have correct properties."""

    def test_debye_positive(self):
        J = DebyeSpectralDensity(eta=0.01, tau_c=1.0)
        for omega in [0.01, 0.1, 1.0, 10.0, 100.0]:
            assert J(omega) >= 0

    def test_debye_zero_at_zero(self):
        J = DebyeSpectralDensity(eta=0.01, tau_c=1.0)
        assert J(0.0) == 0.0

    def test_debye_peak(self):
        """Debye J(omega) peaks near omega = 1/tau_c."""
        tau_c = 1.0
        J = DebyeSpectralDensity(eta=1.0, tau_c=tau_c)
        # J(omega) = eta * omega * tau_c / (1 + omega^2 * tau_c^2)
        # dJ/domega = 0 at omega = 1/tau_c
        J_at_peak = J(1.0 / tau_c)
        J_below = J(0.1 / tau_c)
        J_above = J(10.0 / tau_c)
        assert J_at_peak > J_below
        assert J_at_peak > J_above

    def test_lipari_szabo_rigid_limit(self):
        """With S^2 = 1.0, Lipari-Szabo reduces to simple Lorentzian."""
        tau_c = 5e-9
        J_ls = LipariSzaboSpectralDensity(tau_c=tau_c, S2=1.0, tau_e=0.0)
        omega = 1e9
        expected = 0.4 * tau_c / (1.0 + omega ** 2 * tau_c ** 2)
        assert abs(J_ls(omega) - expected) < 1e-20

    def test_lipari_szabo_positive(self):
        J = LipariSzaboSpectralDensity(tau_c=5e-9, S2=0.85, tau_e=50e-12)
        for omega in [1e6, 1e8, 1e10, 1e12]:
            assert J(omega) >= 0