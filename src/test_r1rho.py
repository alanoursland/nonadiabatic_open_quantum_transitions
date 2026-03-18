import pytest
import numpy as np
import torch
import matplotlib.pyplot as plt
from r1rho import R1rhoExtractor, R1rhoDispersion


# ============================================================
# Helpers
# ============================================================

class FakeResult:
    """Minimal Result-like object for unit testing R1rhoExtractor."""

    def __init__(self, times, states):
        self.times = times
        self.states = states
        self.nonadiabatic_populations = None


def _make_synthetic_decay(dim, observable, M0, R, times_np):
    """Build a FakeResult whose <O(t)> = M0 * exp(-R * t) exactly.

    Constructs rho(t) = (I/dim) + f(t) * O_normalized so that
    Tr[rho(t) @ O] = f(t) * Tr[O @ O] = M0 * exp(-R*t).
    """
    O = observable.to(torch.complex128)
    trOO = torch.trace(O @ O).real.item()
    assert trOO > 0, "Observable must have non-zero Frobenius norm"

    times = torch.tensor(times_np, dtype=torch.float64)
    identity = torch.eye(dim, dtype=torch.complex128) / dim

    states = []
    for t_val in times_np:
        f = M0 * np.exp(-R * t_val) / trOO
        rho = identity + f * O
        states.append(rho)

    return FakeResult(times, torch.stack(states))


# ============================================================
# R1rhoExtractor tests
# ============================================================

class TestExpectationValues:
    """Test that get_expectation_values computes Tr[rho @ O] correctly."""

    def test_identity_observable(self):
        """<I> = Tr[rho] = 1 at all times."""
        dim = 2
        O = torch.eye(dim, dtype=torch.complex128)
        times_np = np.linspace(0, 1, 10)
        # Any valid density matrices
        states = []
        for _ in times_np:
            rho = torch.tensor([[0.6, 0.1 + 0.2j],
                                [0.1 - 0.2j, 0.4]], dtype=torch.complex128)
            states.append(rho)
        result = FakeResult(torch.tensor(times_np, dtype=torch.float64),
                            torch.stack(states))
        ext = R1rhoExtractor(result, O)
        vals = ext.get_expectation_values()
        assert torch.allclose(vals, torch.ones(10, dtype=torch.float64), atol=1e-12)

    def test_sigma_z_expectation(self):
        """<sigma_z> for |0><0| is +1, for |1><1| is -1."""
        dim = 2
        sigma_z = torch.tensor([[1, 0], [0, -1]], dtype=torch.complex128)

        rho_up = torch.tensor([[1, 0], [0, 0]], dtype=torch.complex128)
        rho_down = torch.tensor([[0, 0], [0, 1]], dtype=torch.complex128)

        states = torch.stack([rho_up, rho_down])
        times = torch.tensor([0.0, 1.0], dtype=torch.float64)
        result = FakeResult(times, states)

        ext = R1rhoExtractor(result, sigma_z)
        vals = ext.get_expectation_values()
        assert abs(vals[0].item() - 1.0) < 1e-12
        assert abs(vals[1].item() - (-1.0)) < 1e-12

    def test_off_diagonal_observable(self):
        """<sigma_x> for a coherence state."""
        sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        # rho = |+><+| = [[0.5, 0.5], [0.5, 0.5]], so <sigma_x> = 1
        rho_plus = torch.tensor([[0.5, 0.5], [0.5, 0.5]], dtype=torch.complex128)
        states = torch.stack([rho_plus])
        times = torch.tensor([0.0], dtype=torch.float64)
        result = FakeResult(times, states)

        ext = R1rhoExtractor(result, sigma_x)
        vals = ext.get_expectation_values()
        assert abs(vals[0].item() - 1.0) < 1e-12


class TestExponentialFit:
    """Test that fit() recovers known exponential decay parameters."""

    def test_exact_exponential(self):
        """Perfect exponential data should be recovered exactly."""
        dim = 2
        M0_true = 0.5
        R_true = 2.0
        sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        times_np = np.linspace(0, 3.0, 50)

        result = _make_synthetic_decay(dim, sigma_x, M0_true, R_true, times_np)
        ext = R1rhoExtractor(result, sigma_x)
        fit = ext.fit()

        assert abs(fit['R1rho'] - R_true) < 1e-6, \
            f"R1rho = {fit['R1rho']}, expected {R_true}"
        assert abs(fit['M0'] - M0_true) < 1e-6, \
            f"M0 = {fit['M0']}, expected {M0_true}"
        assert fit['residual_std'] < 1e-10

    def test_slow_decay(self):
        """Very slow decay (R ~ 0.01) should still be recovered."""
        dim = 2
        M0_true = 0.3
        R_true = 0.01
        sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        times_np = np.linspace(0, 100.0, 200)

        result = _make_synthetic_decay(dim, sigma_x, M0_true, R_true, times_np)
        ext = R1rhoExtractor(result, sigma_x)
        fit = ext.fit()

        assert abs(fit['R1rho'] - R_true) / R_true < 1e-4
        assert abs(fit['M0'] - M0_true) / M0_true < 1e-4

    def test_fitting_window(self):
        """Fitting only over a sub-interval should still work."""
        dim = 2
        M0_true = 0.4
        R_true = 1.0
        sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        times_np = np.linspace(0, 5.0, 100)

        result = _make_synthetic_decay(dim, sigma_x, M0_true, R_true, times_np)
        ext = R1rhoExtractor(result, sigma_x)

        # Fit only over [1, 4] — the effective M0 at t=1 is M0*exp(-R*1)
        fit = ext.fit(t_start=1.0, t_end=4.0)

        # R1rho should still be recovered
        assert abs(fit['R1rho'] - R_true) < 1e-4
        # M0 is the amplitude at t_start, i.e. M0_true * exp(-R*1)
        # Tolerance is looser because the grid doesn't land exactly on t_start
        expected_M0 = M0_true * np.exp(-R_true * 1.0)
        assert abs(fit['M0'] - expected_M0) < 0.01

    def test_too_few_points_raises(self):
        """Fewer than 3 points in the window should raise ValueError."""
        dim = 2
        sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        times_np = np.array([0.0, 1.0])
        result = _make_synthetic_decay(dim, sigma_x, 0.5, 1.0, times_np)
        ext = R1rhoExtractor(result, sigma_x)
        with pytest.raises(ValueError, match="at least 3 points"):
            ext.fit()

    def test_larger_system(self):
        """4x4 system (2 spins) with Ix observable."""
        dim = 4
        # Ix for spin 0 in a 2-spin system: Ix ⊗ I
        Ix1 = torch.tensor([[0, 0.5], [0.5, 0]], dtype=torch.complex128)
        Id = torch.eye(2, dtype=torch.complex128)
        O = torch.kron(Ix1, Id)

        M0_true = 0.25
        R_true = 0.5
        times_np = np.linspace(0, 10.0, 80)

        result = _make_synthetic_decay(dim, O, M0_true, R_true, times_np)
        ext = R1rhoExtractor(result, O)
        fit = ext.fit()

        assert abs(fit['R1rho'] - R_true) < 1e-4
        assert abs(fit['M0'] - M0_true) < 1e-4


# ============================================================
# R1rhoDispersion tests
# ============================================================

class TestR1rhoDispersion:

    def test_add_and_retrieve(self):
        disp = R1rhoDispersion()
        disp.add_point('redfield', 100.0, 5.0)
        disp.add_point('redfield', 200.0, 3.0)
        disp.add_point('secular', 100.0, 4.8)

        assert set(disp.get_methods()) == {'redfield', 'secular'}

        omega1, r1rho = disp.get_dispersion('redfield')
        np.testing.assert_array_equal(omega1, [100.0, 200.0])
        np.testing.assert_array_equal(r1rho, [5.0, 3.0])

    def test_empty_dispersion(self):
        disp = R1rhoDispersion()
        assert disp.get_methods() == []

    def test_plot_returns_axes(self):
        """plot() should return a matplotlib Axes without error."""
        disp = R1rhoDispersion()
        disp.add_point('method_a', 50.0, 10.0)
        disp.add_point('method_a', 100.0, 8.0)
        ax = disp.plot()
        assert ax is not None
        plt.close('all')

    def test_plot_with_experimental(self):
        """plot() with experimental overlay should not raise."""
        disp = R1rhoDispersion()
        disp.add_point('theory', 100.0, 5.0)
        exp = {'omega1': np.array([100.0, 200.0]),
               'R1rho': np.array([5.1, 3.2])}
        ax = disp.plot(experimental=exp)
        assert ax is not None
        plt.close('all')


# ============================================================
# Integration: R1rhoExtractor with actual Simulator
# ============================================================

class TestR1rhoWithSimulator:
    """Smoke test using a real Simulator run + R1rhoExtractor."""

    def test_secular_lindblad_decay_gives_positive_rate(self):
        """Run secular_lindblad, extract <Ix(t)>, fit should give R1rho > 0."""
        from system import System
        from bath import SpectralDensityBath, DebyeSpectralDensity

        omega = 5.0
        sys = System([0.0, omega], [[0, 1], [1, 0]])
        S = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        J = DebyeSpectralDensity(eta=0.5, tau_c=1.0)
        bath = SpectralDensityBath(J, 1e-10, [S])

        energies = torch.tensor([0.0, omega], dtype=torch.float64)
        eigvecs = torch.eye(2, dtype=torch.complex128)
        bath.precompute(energies, eigvecs)

        # Zero drive
        zero = torch.zeros(2, 2, dtype=torch.complex128)

        class ZeroDrive:
            def __call__(self, t):
                return zero
            def derivative(self, t):
                return zero

        from simulator import Simulator
        sim = Simulator(sys, ZeroDrive(), bath)

        # Start in |+> state so <sigma_x> = 1 initially
        rho_0 = torch.tensor([[0.5, 0.5], [0.5, 0.5]], dtype=torch.complex128)
        times = np.linspace(0, 20.0, 100)
        result = sim.run(times, method='secular_lindblad', initial_state=rho_0)

        sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        ext = R1rhoExtractor(result, sigma_x, method_name='secular_lindblad')

        # Check expectation values are real and start near 1
        vals = ext.get_expectation_values()
        assert abs(vals[0].item() - 1.0) < 1e-6

        fit = ext.fit()
        assert fit['R1rho'] > 0, f"Expected positive R1rho, got {fit['R1rho']}"
        assert fit['M0'] > 0
