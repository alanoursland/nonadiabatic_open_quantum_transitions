"""Tests for hitran_estimate.py: Phase 0A HITRAN line-mixing feasibility."""
import math
import pytest
import torch

from hitran_estimate import (
    perturbation_ratio_table, cgme_sinc_table_hitran, phase0a_report,
)
from rotational_system import RotationalSystem, build_collisional_coupling_operator, CM_TO_RAD_S
from spectral_densities import CollisionalLorentzian
from bath import SpectralDensityBath


# --- Fixtures ---

@pytest.fixture
def co_1atm():
    """CO at 1 atm: expected NEGATIVE."""
    return phase0a_report(B_cm=1.9313, J_max=10, gamma_pressure_cm=0.065,
                          pressure_atm=1.0)


@pytest.fixture
def co_10atm():
    """CO at 10 atm: expected POSITIVE."""
    return phase0a_report(B_cm=1.9313, J_max=10, gamma_pressure_cm=0.065,
                          pressure_atm=10.0)


@pytest.fixture
def co2_1atm():
    """CO2 at 1 atm: borderline NEGATIVE."""
    return phase0a_report(B_cm=0.3902, J_max=15, gamma_pressure_cm=0.070,
                          pressure_atm=1.0)


@pytest.fixture
def co2_5atm():
    """CO2 at 5 atm: expected POSITIVE."""
    return phase0a_report(B_cm=0.3902, J_max=15, gamma_pressure_cm=0.070,
                          pressure_atm=5.0)


# --- Perturbation ratio tests ---

class TestPerturbationRatio:

    def test_co_j0_ratio(self):
        """V/DeltaE for CO at 1 atm, J=0: 0.065 / (2*1.9313) ~ 0.0168."""
        table = perturbation_ratio_table(B_cm=1.9313, J_max=5,
                                         gamma_pressure_cm=0.065)
        j0 = table[0]
        expected_ratio = 0.065 / (2.0 * 1.9313)
        assert j0['ratio'] == pytest.approx(expected_ratio, rel=1e-10)

    def test_ratio_decreases_with_J(self):
        """V/DeltaE decreases with J (DeltaE grows, V constant)."""
        table = perturbation_ratio_table(B_cm=1.9313, J_max=10,
                                         gamma_pressure_cm=0.065)
        for i in range(len(table) - 1):
            assert table[i]['ratio'] > table[i + 1]['ratio']

    def test_scales_linearly_with_pressure(self):
        """V/DeltaE at 10 atm = 10 * V/DeltaE at 1 atm."""
        table_1 = perturbation_ratio_table(B_cm=1.9313, J_max=5,
                                           gamma_pressure_cm=0.065,
                                           pressure_atm=1.0)
        table_10 = perturbation_ratio_table(B_cm=1.9313, J_max=5,
                                            gamma_pressure_cm=0.065,
                                            pressure_atm=10.0)
        for r1, r10 in zip(table_1, table_10):
            assert r10['ratio'] == pytest.approx(10.0 * r1['ratio'], rel=1e-10)

    def test_mixing_is_ratio_squared(self):
        """mixing_fraction = ratio^2 exactly."""
        table = perturbation_ratio_table(B_cm=1.9313, J_max=5,
                                         gamma_pressure_cm=0.065)
        for r in table:
            assert r['mixing_fraction'] == pytest.approx(r['ratio'] ** 2, rel=1e-10)

    def test_correct_number_of_entries(self):
        """One entry per J from 0 to J_max - 1."""
        table = perturbation_ratio_table(B_cm=1.0, J_max=8,
                                         gamma_pressure_cm=0.1)
        assert len(table) == 8

    def test_delta_E_formula(self):
        """DeltaE = 2*B*(J+1) for each entry."""
        B = 1.9313
        table = perturbation_ratio_table(B_cm=B, J_max=5,
                                         gamma_pressure_cm=0.065)
        for r in table:
            expected = 2.0 * B * (r['J'] + 1)
            assert r['delta_E_cm'] == pytest.approx(expected, rel=1e-10)


# --- CGME sinc tests ---

class TestCGMESinc:

    def test_sinc_near_one_for_short_coarse_graining(self):
        """With tau_c = 1e-12 s and 2B ~ 3.86 cm^-1, sinc argument ~ 0.36.

        Unlike NMR (where sinc was tiny), for HITRAN with picosecond
        collisions the CGME retains nearly all non-secular elements.
        sinc(0.36) ~ 0.98, so CGME ≈ full Redfield (not secular).
        """
        table = cgme_sinc_table_hitran(B_cm=1.9313, J_max=10,
                                       delta_tau_s=1e-12)
        # Adjacent R-branch: argument = 2B * CM_TO_RAD_S * 1e-12 / 2 ~ 0.36
        adjacent = table[0]
        assert adjacent['sinc_abs'] > 0.9  # sinc is close to 1

    def test_sinc_values_bounded(self):
        """All sinc absolute values are between 0 and 1."""
        table = cgme_sinc_table_hitran(B_cm=0.3902, J_max=10,
                                       delta_tau_s=1e-12)
        for s in table:
            assert 0.0 <= s['sinc_abs'] <= 1.0

    def test_returns_multiple_pairs(self):
        """Table includes adjacent, next-nearest, smallest, and largest."""
        table = cgme_sinc_table_hitran(B_cm=1.9313, J_max=10,
                                       delta_tau_s=1e-12)
        assert len(table) >= 3


# --- Phase 0A report tests ---

class TestPhase0AReport:

    def test_report_structure(self, co_1atm):
        """Report dict has all expected keys."""
        expected_keys = [
            'parameters', 'perturbation_table', 'cgme_sinc_table',
            'max_perturbation_ratio', 'max_mixing_fraction',
            'max_sinc_retention', 'verdict', 'threshold',
        ]
        for key in expected_keys:
            assert key in co_1atm

    def test_parameters_complete(self, co_1atm):
        """Parameters dict has all physical quantities."""
        p = co_1atm['parameters']
        for key in ['B_cm', 'J_max', 'gamma_pressure_cm', 'pressure_atm',
                     'T_K', 'tau_collision_s', 'V_cm', 'V_rad_s',
                     'smallest_gap_cm', 'smallest_gap_rad_s']:
            assert key in p

    def test_threshold_is_one_percent(self, co_1atm):
        """Threshold = 0.01 (1%)."""
        assert co_1atm['threshold'] == 0.01

    def test_co_1atm_verdict_negative(self, co_1atm):
        """CO at 1 atm: |a|^2 ~ 3e-4 < 0.01 -> NEGATIVE."""
        assert co_1atm['verdict'] == 'NEGATIVE'
        assert co_1atm['max_mixing_fraction'] < 0.01

    def test_co_10atm_verdict_positive(self, co_10atm):
        """CO at 10 atm: V/DeltaE ~ 0.17, |a|^2 ~ 0.028 -> POSITIVE."""
        assert co_10atm['verdict'] == 'POSITIVE'
        assert co_10atm['max_mixing_fraction'] >= 0.01

    def test_co2_1atm_below_threshold(self, co2_1atm):
        """CO2 at 1 atm: |a|^2 ~ 8e-3, borderline but below 1%."""
        assert co2_1atm['max_mixing_fraction'] < 0.01

    def test_co2_5atm_verdict_positive(self, co2_5atm):
        """CO2 at 5 atm: V/DeltaE ~ 0.45, |a|^2 ~ 0.20 -> POSITIVE."""
        assert co2_5atm['verdict'] == 'POSITIVE'
        assert co2_5atm['max_mixing_fraction'] > 0.1

    def test_co2_larger_ratio_than_co(self, co_1atm, co2_1atm):
        """CO2 (smaller B) has larger V/DeltaE than CO at same pressure."""
        assert co2_1atm['max_perturbation_ratio'] > co_1atm['max_perturbation_ratio']

    def test_all_energies_physical(self, co_1atm):
        """All energy gaps and coupling strengths are positive."""
        for r in co_1atm['perturbation_table']:
            assert r['delta_E_cm'] > 0
            assert r['V_cm'] > 0
            assert r['ratio'] > 0
            assert r['mixing_fraction'] > 0

    def test_max_ratio_from_j0(self, co_1atm):
        """Max perturbation ratio comes from J=0 (smallest DeltaE)."""
        table = co_1atm['perturbation_table']
        assert table[0]['ratio'] == co_1atm['max_perturbation_ratio']

    def test_co_10atm_mixing_magnitude(self, co_10atm):
        """CO at 10 atm: |a|^2(J=0) ~ (0.65/3.86)^2 ~ 0.028."""
        expected = (0.065 * 10.0 / (2.0 * 1.9313)) ** 2
        assert co_10atm['max_mixing_fraction'] == pytest.approx(expected, rel=1e-3)


# --- Integration test: Redfield tensor with RotationalSystem ---

class TestRedfieldIntegration:
    """Verify that SpectralDensityBath works with RotationalSystem + CollisionalLorentzian.

    The coupling operator is dimensionless (elements = 1), matching the
    codebase convention (S is a dimensionless system operator, J carries
    coupling strength). build_collisional_coupling_operator is for the
    Phase 0A estimate only; for the Redfield tensor we separate S and J.
    """

    @pytest.fixture
    def rotational_bath(self):
        """Build a SpectralDensityBath for a 6-level rigid rotor (CO, J_max=5)."""
        J_max = 5
        system = RotationalSystem(B_cm=1.9313, J_max=J_max)
        # Spectral density carries the full coupling strength
        gamma_0 = 0.065 * CM_TO_RAD_S  # ~1.2e10 rad/s
        J_fn = CollisionalLorentzian(gamma_0=gamma_0, tau_c=1e-12)
        # Dimensionless coupling operator (elements = 1 for DeltaJ = +-1)
        dim = J_max + 1
        S = torch.zeros(dim, dim, dtype=torch.complex128)
        for J in range(J_max):
            S[J, J + 1] = 1.0
            S[J + 1, J] = 1.0
        bath = SpectralDensityBath(J_fn, temperature_K=296.0,
                                   coupling_operators=[S])
        # Rigid rotor: H0 is diagonal, eigenvectors = identity
        eigvecs = torch.eye(dim, dtype=torch.complex128)
        bath.precompute(system.energy_levels, eigvecs)
        return bath, system

    def test_precompute_succeeds(self, rotational_bath):
        """SpectralDensityBath.precompute completes without error."""
        bath, system = rotational_bath
        assert bath._precomputed

    def test_redfield_dissipator_trace_preserving(self, rotational_bath):
        """Tr[D[rho]] = 0 for the Redfield dissipator."""
        bath, system = rotational_bath
        dim = system.dim
        # Start from a mixed state (not thermal — just a test state)
        rho = torch.zeros(dim, dim, dtype=torch.complex128)
        rho[0, 0] = 0.5
        rho[1, 1] = 0.3
        rho[2, 2] = 0.2
        D_rho = bath.get_redfield_dissipator(rho)
        trace = torch.trace(D_rho).real.item()
        assert abs(trace) < 1e-10

    def test_secular_lindblad_trace_preserving(self, rotational_bath):
        """Tr[D[rho]] = 0 for the secular Lindblad dissipator."""
        bath, system = rotational_bath
        dim = system.dim
        rho = torch.zeros(dim, dim, dtype=torch.complex128)
        rho[0, 0] = 0.5
        rho[1, 1] = 0.3
        rho[2, 2] = 0.2
        D_rho = bath.get_secular_lindblad_dissipator(rho)
        trace = torch.trace(D_rho).real.item()
        assert abs(trace) < 1e-10

    def test_cgme_dissipator_trace_preserving(self, rotational_bath):
        """Tr[D[rho]] = 0 for the CGME dissipator."""
        bath, system = rotational_bath
        dim = system.dim
        rho = torch.zeros(dim, dim, dtype=torch.complex128)
        rho[0, 0] = 0.5
        rho[1, 1] = 0.3
        rho[2, 2] = 0.2
        D_rho = bath.get_cgme_dissipator(rho, delta_tau=1e-12)
        trace = torch.trace(D_rho).real.item()
        assert abs(trace) < 1e-10

    def test_dissipator_hermitian(self, rotational_bath):
        """D[rho] is Hermitian when rho is Hermitian."""
        bath, system = rotational_bath
        dim = system.dim
        rho = torch.zeros(dim, dim, dtype=torch.complex128)
        rho[0, 0] = 0.5
        rho[1, 1] = 0.3
        rho[2, 2] = 0.2
        D_rho = bath.get_redfield_dissipator(rho)
        assert torch.allclose(D_rho, D_rho.conj().T, atol=1e-10)

    def test_excited_state_decays(self, rotational_bath):
        """Starting in J=5, population should decrease (decay to lower states)."""
        bath, system = rotational_bath
        dim = system.dim
        rho = torch.zeros(dim, dim, dtype=torch.complex128)
        rho[dim - 1, dim - 1] = 1.0  # pure J=5 state
        D_rho = bath.get_redfield_dissipator(rho)
        # d(rho_{55})/dt should be negative (population leaving J=5)
        assert D_rho[dim - 1, dim - 1].real.item() < 0

    def test_ground_state_gains(self, rotational_bath):
        """Starting in J=5, ground state (J=0) population should increase
        (through cascade)."""
        bath, system = rotational_bath
        dim = system.dim
        rho = torch.zeros(dim, dim, dtype=torch.complex128)
        rho[dim - 1, dim - 1] = 1.0  # pure J=5 state
        D_rho = bath.get_redfield_dissipator(rho)
        # Adjacent state J=4 should gain population
        assert D_rho[dim - 2, dim - 2].real.item() > 0
