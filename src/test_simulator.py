import pytest
import numpy as np
import torch
from system import System
from signals import PlateauPulse, SpinLockDrive
from bath import ThermalBath, SpectralDensityBath, DebyeSpectralDensity
from simulator import Simulator

@pytest.fixture
def simple_setup():
    sys = System([0.0, 5.0], [[0.0, 1.0], [1.0, 0.0]])
    pulse = PlateauPulse(0.1, 5.0, 5.0, 20.0)
    bath = ThermalBath(0.015, {(1, 0): 50.0}, {(1, 0): 30.0})
    return Simulator(sys, pulse, bath)

def test_simulator_initialization(simple_setup):
    assert simple_setup.system.dim == 2

def test_commutator(simple_setup):
    A = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
    B = torch.tensor([[1, 0], [0, -1]], dtype=torch.complex128)
    # [X, Z] = XZ - ZX = -2iY
    comm = simple_setup.commutator(A, B)
    expected = torch.tensor([[0, -2], [2, 0]], dtype=torch.complex128)
    assert torch.allclose(comm, expected)

def test_run_trajectory_shape(simple_setup):
    times = np.linspace(0, 10, 50)
    res = simple_setup.run(times, method='lindblad')
    
    # Expect 50 time steps, 2x2 matrices
    assert res.states.shape == (50, 2, 2)

def test_run_trace_preservation(simple_setup):
    times = np.linspace(0, 15, 100)
    res = simple_setup.run(times, method='lindblad')
    
    # Check that trace is 1.0 at the final time step
    final_state = res.states[-1]
    assert torch.allclose(torch.trace(final_state), torch.tensor(1.0 + 0j, dtype=torch.complex128))

def test_get_populations(simple_setup):
    times = np.linspace(0, 5, 10)
    res = simple_setup.run(times, method='lindblad')
    
    pops_0 = res.get_populations(0)
    pops_1 = res.get_populations(1)
    
    # Populations should be 1D tensors of length len(times)
    assert pops_0.shape == (10,)
    # Total probability must equal 1 at all times
    total_prob = pops_0 + pops_1
    assert torch.allclose(total_prob, torch.ones_like(total_prob, dtype=torch.float64))


# ============================================================
# SpectralDensityBath eigenbasis method tests (WP5)
# ============================================================

class ZeroDrive:
    """V(t) = 0 for all t."""
    def __init__(self, dim):
        self._zero = torch.zeros(dim, dim, dtype=torch.complex128)

    def __call__(self, t):
        return self._zero

    def derivative(self, t):
        return self._zero


class ConstantMatrixDrive:
    """V(t) = strength * operator for all t. dV/dt = 0."""
    def __init__(self, strength, operator):
        self.strength = strength
        self.op = operator
        self._zero = torch.zeros_like(operator)

    def __call__(self, t):
        return self.strength * self.op

    def derivative(self, t):
        return self._zero


def _make_eigenbasis_simulator(omega=5.0, eta=0.1, tau_c=1.0, T_K=1e-10, drive=None):
    """2LS with diagonal H0, sigma_x coupling, Debye bath. Ready to run.

    Default T_K=1e-10 ensures hbar*omega/kT ~ 0.4 so that Bose-Einstein
    occupation is O(1) and bath rates are O(eta), not O(eta * kT/hbar*omega).
    At T=300K with omega=5 rad/s, n(omega) ~ 10^12, making rates enormous
    and the integrator unstable.
    """
    sys = System([0.0, omega], [[0, 1], [1, 0]])
    S = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
    J = DebyeSpectralDensity(eta=eta, tau_c=tau_c)
    bath = SpectralDensityBath(J, T_K, [S])

    energies = torch.tensor([0.0, omega], dtype=torch.float64)
    eigvecs = torch.eye(2, dtype=torch.complex128)
    bath.precompute(energies, eigvecs)

    if drive is None:
        drive = ZeroDrive(2)

    return Simulator(sys, drive, bath)


class TestEigenbasisTracePreservation:
    """Tr[rho(t)] = 1 at all times for all three new methods."""

    @pytest.mark.parametrize("method,kwargs", [
        ('redfield', {}),
        ('secular_lindblad', {}),
        ('nonadiabatic_cgme', {'delta_tau': 0.5}),
    ])
    def test_trace(self, method, kwargs):
        sim = _make_eigenbasis_simulator()
        rho_0 = torch.tensor([[0.3, 0.0], [0.0, 0.7]], dtype=torch.complex128)
        times = np.linspace(0, 2.0, 20)
        result = sim.run(times, method=method, initial_state=rho_0, **kwargs)
        for i in range(len(times)):
            trace = torch.trace(result.states[i]).real.item()
            assert abs(trace - 1.0) < 1e-8, \
                f"{method}: Tr[rho] = {trace} at t={times[i]}"

    @pytest.mark.parametrize("method,kwargs", [
        ('redfield', {}),
        ('secular_lindblad', {}),
        ('nonadiabatic_cgme', {'delta_tau': 0.5}),
    ])
    def test_trace_with_drive(self, method, kwargs):
        """Trace preserved even with a constant drive perturbing the system."""
        sigma_x = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        drive = ConstantMatrixDrive(0.1, sigma_x)
        sim = _make_eigenbasis_simulator(drive=drive)
        rho_0 = torch.tensor([[0.5, 0.2 + 0.1j], [0.2 - 0.1j, 0.5]], dtype=torch.complex128)
        times = np.linspace(0, 1.0, 15)
        result = sim.run(times, method=method, initial_state=rho_0, **kwargs)
        for i in range(len(times)):
            trace = torch.trace(result.states[i]).real.item()
            assert abs(trace - 1.0) < 1e-8, \
                f"{method} with drive: Tr[rho] = {trace} at t={times[i]}"


class TestEigenbasisHermiticity:
    """rho(t) must be Hermitian at all times."""

    @pytest.mark.parametrize("method,kwargs", [
        ('redfield', {}),
        ('secular_lindblad', {}),
        ('nonadiabatic_cgme', {'delta_tau': 0.5}),
    ])
    def test_hermitian(self, method, kwargs):
        sim = _make_eigenbasis_simulator()
        rho_0 = torch.tensor([[0.6, 0.2 + 0.1j], [0.2 - 0.1j, 0.4]], dtype=torch.complex128)
        times = np.linspace(0, 1.0, 10)
        result = sim.run(times, method=method, initial_state=rho_0, **kwargs)
        for i in range(len(times)):
            rho = result.states[i]
            assert torch.allclose(rho, rho.mH, atol=1e-10), \
                f"{method}: rho not Hermitian at t={times[i]}"


class TestPopulationDecay:
    """Excited-state population should decrease under dissipation."""

    def test_secular_decay_from_excited(self):
        sim = _make_eigenbasis_simulator(eta=0.5)
        rho_0 = torch.tensor([[0.0, 0.0], [0.0, 1.0]], dtype=torch.complex128)
        times = np.linspace(0, 50.0, 50)
        result = sim.run(times, method='secular_lindblad', initial_state=rho_0)

        pop_excited = result.get_populations(1)
        assert pop_excited[-1].item() < pop_excited[0].item(), \
            "Excited population should decay"
        assert result.get_populations(0)[-1].item() > 0, \
            "Ground state should gain population"

    def test_redfield_decay_from_excited(self):
        sim = _make_eigenbasis_simulator(eta=0.5)
        rho_0 = torch.tensor([[0.0, 0.0], [0.0, 1.0]], dtype=torch.complex128)
        times = np.linspace(0, 50.0, 50)
        result = sim.run(times, method='redfield', initial_state=rho_0)

        pop_excited = result.get_populations(1)
        assert pop_excited[-1].item() < pop_excited[0].item(), \
            "Excited population should decay"


class TestZeroDrivePopulationAgreement:
    """With V=0 and diagonal initial state, all three methods give the
    same population dynamics (the pop-to-pop block of R is identical)."""

    def test_populations_agree(self):
        omega = 5.0
        rho_0 = torch.tensor([[0.3, 0.0], [0.0, 0.7]], dtype=torch.complex128)
        times = np.linspace(0, 3.0, 30)

        results = {}
        for method, kw in [('redfield', {}),
                           ('secular_lindblad', {}),
                           ('nonadiabatic_cgme', {'delta_tau': 0.5})]:
            sim = _make_eigenbasis_simulator(omega=omega)
            results[method] = sim.run(times, method=method, initial_state=rho_0, **kw)

        for state_idx in [0, 1]:
            pop_red = results['redfield'].get_populations(state_idx)
            pop_sec = results['secular_lindblad'].get_populations(state_idx)
            pop_cgme = results['nonadiabatic_cgme'].get_populations(state_idx)

            assert torch.allclose(pop_red, pop_sec, atol=1e-8), \
                f"State {state_idx}: Redfield and secular populations differ"
            assert torch.allclose(pop_red, pop_cgme, atol=1e-8), \
                f"State {state_idx}: Redfield and CGME populations differ"


class TestNonadiabaticTracking:
    """The nonadiabatic_cgme method stores b_k population data alongside rho."""

    def test_populations_present(self):
        sim = _make_eigenbasis_simulator()
        times = np.linspace(0, 1.0, 10)
        result = sim.run(times, method='nonadiabatic_cgme', delta_tau=0.5)
        assert result.nonadiabatic_populations is not None
        assert result.nonadiabatic_populations.shape == (len(times), 2)

    def test_populations_sum_to_one(self):
        sim = _make_eigenbasis_simulator()
        times = np.linspace(0, 1.0, 10)
        result = sim.run(times, method='nonadiabatic_cgme', delta_tau=0.5)
        assert result.nonadiabatic_populations is not None
        for i in range(len(times)):
            total = result.nonadiabatic_populations[i].sum().item()
            assert abs(total - 1.0) < 1e-10, \
                f"Sum of nonadiabatic populations = {total} at step {i}"

    def test_not_present_for_other_methods(self):
        sim = _make_eigenbasis_simulator()
        times = np.linspace(0, 1.0, 5)
        result = sim.run(times, method='secular_lindblad')
        assert result.nonadiabatic_populations is None

    def test_plateau_invariance(self):
        """During a plateau (dV/dt = 0), nonadiabatic populations are constant.
        Uses SpinLockDrive with a ramp then plateau."""
        omega = 5.0
        Ix = torch.tensor([[0, 0.5], [0.5, 0]], dtype=torch.complex128)
        drive = SpinLockDrive(
            spin_lock_power_hz=100.0,
            Ix_operator=Ix,
            ramp_time_us=1.0,
            plateau_time_us=10.0,
        )

        sys = System([0.0, omega], [[0, 1], [1, 0]])
        S = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        J = DebyeSpectralDensity(eta=0.1, tau_c=1.0)
        bath = SpectralDensityBath(J, 1e-10, [S])
        energies = torch.tensor([0.0, omega], dtype=torch.float64)
        eigvecs = torch.eye(2, dtype=torch.complex128)
        bath.precompute(energies, eigvecs)

        sim = Simulator(sys, drive, bath)

        # Ramp: 0-1 us, Plateau: 1-11 us. Sample well into plateau.
        times = np.linspace(0, 6.0, 60)
        result = sim.run(times, method='nonadiabatic_cgme', delta_tau=0.5)

        # All plateau samples (t > 1.5 us) should have identical b_k populations
        assert result.nonadiabatic_populations is not None
        plateau_start = np.searchsorted(times, 1.5)
        plateau_pops = result.nonadiabatic_populations[plateau_start:]
        first = plateau_pops[0]
        for i in range(1, len(plateau_pops)):
            assert torch.allclose(plateau_pops[i], first, atol=1e-10), \
                f"Nonadiabatic populations changed during plateau at index {plateau_start + i}"


class TestConfigurableInitialState:

    def test_custom_initial_state(self):
        sim = _make_eigenbasis_simulator()
        rho_0 = torch.tensor([[0.2, 0.0], [0.0, 0.8]], dtype=torch.complex128)
        times = np.array([0.0, 0.01])
        result = sim.run(times, method='secular_lindblad', initial_state=rho_0)
        assert torch.allclose(result.states[0], rho_0, atol=1e-12)

    def test_default_ground_state(self):
        sim = _make_eigenbasis_simulator()
        times = np.array([0.0, 0.01])
        result = sim.run(times, method='secular_lindblad')
        expected = torch.zeros(2, 2, dtype=torch.complex128)
        expected[0, 0] = 1.0
        assert torch.allclose(result.states[0], expected, atol=1e-12)


class TestMethodDispatch:
    """Verify that the method parameter selects the correct code path."""

    def test_unknown_method_raises(self):
        sim = _make_eigenbasis_simulator()
        times = np.array([0.0, 1.0])
        with pytest.raises(ValueError, match="Unknown method"):
            sim.run(times, method='bogus')

    def test_unprecomputed_bath_raises(self):
        """Eigenbasis methods require bath.precompute() first."""
        sys = System([0.0, 5.0], [[0, 1], [1, 0]])
        S = torch.tensor([[0, 1], [1, 0]], dtype=torch.complex128)
        J = DebyeSpectralDensity(eta=0.01, tau_c=1.0)
        bath = SpectralDensityBath(J, 300.0, [S])
        # Do NOT call bath.precompute()
        sim = Simulator(sys, ZeroDrive(2), bath)
        times = np.array([0.0, 1.0])
        with pytest.raises(AssertionError, match="precompute"):
            sim.run(times, method='redfield')