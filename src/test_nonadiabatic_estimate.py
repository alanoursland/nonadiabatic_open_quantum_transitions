"""
Tests for nonadiabatic correction estimates (Phase 0 feasibility check).

Verifies:
- Analytical correction formulas (mixing fraction, sinc damping)
- Phase 0 report outputs and verdict
- Numerical simulator comparison: secular_lindblad vs nonadiabatic_cgme
"""
import math
import pytest
import numpy as np
import torch

from nonadiabatic_estimate import (
    adiabatic_mixing_fraction,
    cgme_sinc_retention,
    nonadiabatic_correction_table,
    cgme_sinc_table,
    phase0_report,
)
from nmr_analytical import NMRRelaxation, GAMMA_15N, GAMMA_1H
from spectral_densities import LipariSzabo
from system import System
from bath import SpectralDensityBath, DebyeSpectralDensity
from simulator import Simulator


# ============================================================
# Fixtures
# ============================================================

@pytest.fixture
def phase0_nmr():
    """NMRRelaxation instance with Phase 0 parameters."""
    J_fn = LipariSzabo(tau_c=5e-9, S2=0.85, tau_e=50e-12)
    return NMRRelaxation(B0=11.7, J_fn=J_fn)


# ============================================================
# Test adiabatic_mixing_fraction
# ============================================================

class TestAdiabaticMixingFraction:

    def test_scales_as_omega1_squared(self):
        """Mixing fraction should scale as omega1^2."""
        omega0 = 1e8
        f1 = adiabatic_mixing_fraction(100.0, omega0)
        f2 = adiabatic_mixing_fraction(200.0, omega0)
        assert abs(f2 / f1 - 4.0) < 1e-10

    def test_scales_as_omega0_inverse_squared(self):
        """Mixing fraction should scale as 1/omega0^2."""
        omega1 = 1000.0
        f1 = adiabatic_mixing_fraction(omega1, 1e8)
        f2 = adiabatic_mixing_fraction(omega1, 2e8)
        assert abs(f1 / f2 - 4.0) < 1e-10

    def test_nmr_magnitude(self, phase0_nmr):
        """For NMR: omega1 ~ 100 Hz -> ~600 rad/s, omega_N ~ 50 MHz.
        |a_1|^2 ~ (600 / (2 * 3e8))^2 ~ 10^-12."""
        omega1 = 2.0 * math.pi * 100.0  # 100 Hz spin-lock
        mixing = adiabatic_mixing_fraction(omega1, phase0_nmr.omega_N)
        assert mixing < 1e-10, f"Expected < 1e-10, got {mixing}"
        assert mixing > 0

    def test_exact_value(self):
        """Check exact formula: (omega1 / (2 * omega0))^2."""
        omega1, omega0 = 300.0, 1e6
        expected = (omega1 / (2.0 * omega0)) ** 2
        assert abs(adiabatic_mixing_fraction(omega1, omega0) - expected) < 1e-20

    def test_rejects_zero_omega0(self):
        with pytest.raises(ValueError):
            adiabatic_mixing_fraction(100.0, 0.0)

    def test_rejects_negative_omega0(self):
        with pytest.raises(ValueError):
            adiabatic_mixing_fraction(100.0, -1e8)


# ============================================================
# Test cgme_sinc_retention
# ============================================================

class TestCGMESincRetention:

    def test_zero_frequency_diff(self):
        """Secular elements (zero frequency difference) -> sinc = 1."""
        assert abs(cgme_sinc_retention(0.0, 1e-6) - 1.0) < 1e-12

    def test_small_frequency_diff(self):
        """Small frequency difference -> sinc ~ 1."""
        sinc = cgme_sinc_retention(1.0, 1e-6)  # x = 0.5e-6, very small
        assert abs(sinc - 1.0) < 1e-6

    def test_large_frequency_diff(self):
        """Large frequency difference -> |sinc| << 1."""
        # omega_N ~ 3e8 rad/s, delta_tau = 1e-6 s -> x = 150 >> 1
        sinc = cgme_sinc_retention(3e8, 1e-6)
        assert abs(sinc) < 0.01

    def test_sinc_at_pi(self):
        """sinc(pi) = 0: x = pi -> bohr_freq_diff * delta_tau/2 = pi."""
        # x = bohr_freq_diff * delta_tau / 2 = pi
        # -> bohr_freq_diff = 2*pi / delta_tau
        delta_tau = 1e-6
        bohr_freq_diff = 2.0 * math.pi / delta_tau
        sinc = cgme_sinc_retention(bohr_freq_diff, delta_tau)
        assert abs(sinc) < 1e-10

    def test_bounded_by_one(self):
        """sinc(x) is always bounded by 1 in magnitude."""
        for freq in [1e3, 1e6, 1e8, 1e10]:
            sinc = cgme_sinc_retention(freq, 1e-6)
            assert abs(sinc) <= 1.0 + 1e-12

    def test_nmr_omega_n(self, phase0_nmr):
        """For omega_N ~ 50 MHz and delta_tau = 1 us, sinc should be very small."""
        sinc = cgme_sinc_retention(phase0_nmr.omega_N, 1e-6)
        assert abs(sinc) < 1e-2


# ============================================================
# Test nonadiabatic_correction_table
# ============================================================

class TestCorrectionTable:

    def test_returns_correct_length(self, phase0_nmr):
        omega1_values = [25, 50, 100, 500]
        table = nonadiabatic_correction_table(phase0_nmr, omega1_values)
        assert len(table) == 4

    def test_correction_increases_with_omega1(self, phase0_nmr):
        omega1_values = [25, 50, 100, 500, 1000]
        table = nonadiabatic_correction_table(phase0_nmr, omega1_values)
        for i in range(len(table) - 1):
            assert table[i + 1]['mixing_fraction'] > table[i]['mixing_fraction']

    def test_all_corrections_negligible(self, phase0_nmr):
        """All relative corrections should be << 1 for NMR parameters."""
        omega1_values = [25, 50, 100, 250, 500, 1000]
        table = nonadiabatic_correction_table(phase0_nmr, omega1_values)
        for row in table:
            assert row['relative_correction'] < 1e-8, \
                f"omega1={row['omega1_hz']} Hz: relative correction = {row['relative_correction']}"

    def test_keys_present(self, phase0_nmr):
        table = nonadiabatic_correction_table(phase0_nmr, [100])
        row = table[0]
        expected_keys = {'omega1_hz', 'omega1_rad_s', 'mixing_fraction',
                         'delta_R1rho', 'relative_correction'}
        assert set(row.keys()) == expected_keys


# ============================================================
# Test cgme_sinc_table
# ============================================================

class TestSincTable:

    def test_returns_six_pairs(self, phase0_nmr):
        table = cgme_sinc_table(phase0_nmr)
        assert len(table) == 6

    def test_all_sinc_small(self, phase0_nmr):
        """All sinc factors for NMR frequencies should be << 1."""
        table = cgme_sinc_table(phase0_nmr, delta_tau_s=1e-6)
        for row in table:
            assert row['sinc_abs'] < 0.01, \
                f"{row['label']}: |sinc| = {row['sinc_abs']}"

    def test_includes_omega_n(self, phase0_nmr):
        table = cgme_sinc_table(phase0_nmr)
        labels = [r['label'] for r in table]
        assert 'omega_N' in labels

    def test_includes_omega_h(self, phase0_nmr):
        table = cgme_sinc_table(phase0_nmr)
        labels = [r['label'] for r in table]
        assert 'omega_H' in labels


# ============================================================
# Test phase0_report
# ============================================================

class TestPhase0Report:

    def test_verdict_is_negative(self):
        """Phase 0 should be NEGATIVE for standard NMR parameters."""
        report = phase0_report()
        assert report['verdict'] == 'NEGATIVE'

    def test_max_relative_correction_tiny(self):
        report = phase0_report()
        assert report['max_relative_correction'] < 1e-8

    def test_max_sinc_small(self):
        report = phase0_report()
        assert report['max_sinc_retention'] < 0.01

    def test_rates_physical(self):
        report = phase0_report()
        R1 = report['rates']['R1_per_s']
        R2 = report['rates']['R2_per_s']
        assert 0.5 < R1 < 5.0, f"R1 = {R1} out of expected range"
        assert 3.0 < R2 < 20.0, f"R2 = {R2} out of expected range"
        assert R2 > R1

    def test_report_structure(self):
        report = phase0_report()
        assert 'parameters' in report
        assert 'rates' in report
        assert 'perturbative_corrections' in report
        assert 'cgme_sinc_factors' in report
        assert 'verdict' in report
        assert 'threshold' in report

    def test_threshold_is_one_percent(self):
        report = phase0_report()
        assert report['threshold'] == 0.01


# ============================================================
# Numerical simulator check: secular_lindblad vs nonadiabatic_cgme
# ============================================================

class ZeroDrive:
    """V(t) = 0 for all t."""
    def __init__(self, dim):
        self._zero = torch.zeros(dim, dim, dtype=torch.complex128)

    def __call__(self, t):
        return self._zero

    def derivative(self, t):
        return self._zero


class TestSimulatorAgreement:
    """Compare secular_lindblad and nonadiabatic_cgme decay rates on a
    toy 2-level system. With zero drive, both methods should produce
    identical dynamics (no adiabatic/nonadiabatic distinction when dV/dt = 0)."""

    def test_population_agreement_zero_drive(self):
        """With zero drive, secular and CGME give the same populations."""
        omega = 5.0
        sys = System([0.0, omega], [[0, 1], [1, 0]])
        S = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        J = DebyeSpectralDensity(eta=0.1, tau_c=1.0)
        bath = SpectralDensityBath(J, 1e-10, [S])

        energies = torch.tensor([0.0, omega], dtype=torch.float64)
        eigvecs = torch.eye(2, dtype=torch.complex128)
        bath.precompute(energies, eigvecs)

        drive = ZeroDrive(2)
        rho_0 = torch.tensor([[0.2, 0.0], [0.0, 0.8]], dtype=torch.complex128)
        times = np.linspace(0, 3.0, 30)

        sim_sec = Simulator(sys, drive, bath)
        res_sec = sim_sec.run(times, method='secular_lindblad', initial_state=rho_0)

        sim_cgme = Simulator(sys, drive, bath)
        res_cgme = sim_cgme.run(times, method='nonadiabatic_cgme',
                                initial_state=rho_0, delta_tau=0.5)

        for state_idx in [0, 1]:
            pop_sec = res_sec.get_populations(state_idx)
            pop_cgme = res_cgme.get_populations(state_idx)
            assert torch.allclose(pop_sec, pop_cgme, atol=1e-8), \
                f"State {state_idx}: secular and CGME populations differ\n" \
                f"  max diff = {(pop_sec - pop_cgme).abs().max().item():.2e}"

    def test_coherence_approximate_agreement_zero_drive(self):
        """With zero drive, off-diagonal elements agree approximately.

        The secular approximation drops all non-secular Redfield terms, while
        CGME retains them damped by sinc((omega_ab - omega_cd) * delta_tau/2).
        For toy parameters (omega=5, delta_tau=0.5), some sinc factors are O(1),
        so the methods differ slightly in the coherence sector. The populations
        still agree exactly (tested above). The coherence difference confirms
        the CGME retains more physics than the secular approximation.
        """
        omega = 5.0
        sys = System([0.0, omega], [[0, 1], [1, 0]])
        S = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        J = DebyeSpectralDensity(eta=0.1, tau_c=1.0)
        bath = SpectralDensityBath(J, 1e-10, [S])

        energies = torch.tensor([0.0, omega], dtype=torch.float64)
        eigvecs = torch.eye(2, dtype=torch.complex128)
        bath.precompute(energies, eigvecs)

        drive = ZeroDrive(2)
        rho_0 = torch.tensor([[0.5, 0.3 + 0.1j],
                               [0.3 - 0.1j, 0.5]], dtype=torch.complex128)
        times = np.linspace(0, 2.0, 20)

        sim_sec = Simulator(sys, drive, bath)
        res_sec = sim_sec.run(times, method='secular_lindblad', initial_state=rho_0)

        sim_cgme = Simulator(sys, drive, bath)
        res_cgme = sim_cgme.run(times, method='nonadiabatic_cgme',
                                initial_state=rho_0, delta_tau=0.5)

        for i in range(len(times)):
            rho_sec = res_sec.states[i]
            rho_cgme = res_cgme.states[i]
            # Populations agree exactly
            assert abs(rho_sec[0, 0].real.item() - rho_cgme[0, 0].real.item()) < 1e-8
            assert abs(rho_sec[1, 1].real.item() - rho_cgme[1, 1].real.item()) < 1e-8
            # Coherences agree approximately (non-secular CGME terms cause small diff)
            assert torch.allclose(rho_sec, rho_cgme, atol=1e-2), \
                f"States differ too much at t={times[i]:.2f}\n" \
                f"  max diff = {(rho_sec - rho_cgme).abs().max().item():.2e}"

    def test_both_decay_and_agree(self):
        """Both methods decay from pure excited state and agree on final populations."""
        omega = 5.0
        sys = System([0.0, omega], [[0, 1], [1, 0]])
        S = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        J = DebyeSpectralDensity(eta=0.5, tau_c=1.0)
        T_K = 1e-10
        bath = SpectralDensityBath(J, T_K, [S])

        energies = torch.tensor([0.0, omega], dtype=torch.float64)
        eigvecs = torch.eye(2, dtype=torch.complex128)
        bath.precompute(energies, eigvecs)

        drive = ZeroDrive(2)
        rho_0 = torch.tensor([[0.0, 0.0], [0.0, 1.0]], dtype=torch.complex128)
        times = np.linspace(0, 100.0, 100)

        final_pops = {}
        for method, kwargs in [('secular_lindblad', {}),
                                ('nonadiabatic_cgme', {'delta_tau': 0.5})]:
            sim = Simulator(sys, drive, bath)
            res = sim.run(times, method=method, initial_state=rho_0, **kwargs)

            pop_excited = res.get_populations(1)[-1].item()
            # Both methods should show decay from initial p1=1.0
            assert pop_excited < 0.95, \
                f"{method}: excited population barely decayed ({pop_excited:.3f})"
            final_pops[method] = pop_excited

        # Both methods should reach the same equilibrium
        assert abs(final_pops['secular_lindblad'] - final_pops['nonadiabatic_cgme']) < 1e-6, \
            f"Methods disagree on equilibrium: secular={final_pops['secular_lindblad']:.6f}, " \
            f"CGME={final_pops['nonadiabatic_cgme']:.6f}"
