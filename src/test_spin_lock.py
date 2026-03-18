import pytest
import torch
from signals import SpinLockDrive
from spin_system import SpinSystem

# Single spin-1/2 operators for a minimal test system
Ix = torch.tensor([[0, 0.5], [0.5, 0]], dtype=torch.complex128)
Iy = torch.tensor([[0, -0.5j], [0.5j, 0]], dtype=torch.complex128)
Iz = torch.tensor([[0.5, 0], [0, -0.5]], dtype=torch.complex128)


@pytest.fixture
def on_resonance_drive():
    """On-resonance spin-lock: 1 kHz power, 2 us ramp, 50 us plateau."""
    return SpinLockDrive(
        spin_lock_power_hz=1000.0,
        Ix_operator=Ix,
        Iz_operator=Iz,
        offset_hz=0.0,
        ramp_time_us=2.0,
        plateau_time_us=50.0,
    )


@pytest.fixture
def off_resonance_drive():
    """Off-resonance spin-lock: 1 kHz power, 1 kHz offset."""
    return SpinLockDrive(
        spin_lock_power_hz=1000.0,
        Ix_operator=Ix,
        Iz_operator=Iz,
        offset_hz=1000.0,
        ramp_time_us=2.0,
        plateau_time_us=50.0,
    )


@pytest.fixture
def no_ramp_drive():
    """Instant-on spin-lock (no ramp)."""
    return SpinLockDrive(
        spin_lock_power_hz=500.0,
        Ix_operator=Ix,
        offset_hz=0.0,
        ramp_time_us=0.0,
        plateau_time_us=50.0,
    )


def test_perturbation_matrix_on_resonance(on_resonance_drive):
    """During plateau, V(t) = omega1 * Ix."""
    t_plateau = torch.tensor(10.0)  # well within plateau (ramp ends at 2 us)
    V = on_resonance_drive(t_plateau)
    omega1 = on_resonance_drive.omega1
    expected = omega1 * Ix
    assert torch.allclose(V, expected, atol=1e-15)


def test_perturbation_matrix_off_resonance(off_resonance_drive):
    """During plateau, V(t) = omega1 * Ix + Delta * Iz."""
    t_plateau = torch.tensor(10.0)
    V = off_resonance_drive(t_plateau)
    omega1 = off_resonance_drive.omega1
    Delta = off_resonance_drive.Delta
    expected = omega1 * Ix + Delta * Iz
    assert torch.allclose(V, expected, atol=1e-15)


def test_perturbation_hermitian(on_resonance_drive):
    """V(t) must be Hermitian at all times."""
    for t in [0.0, 1.0, 10.0, 51.0, 53.0, 60.0]:
        V = on_resonance_drive(torch.tensor(t))
        assert torch.allclose(V, V.mH, atol=1e-15)


def test_zero_derivative_during_plateau(on_resonance_drive):
    """dV/dt = 0 during the plateau — the key nonadiabatic property."""
    for t in [5.0, 10.0, 25.0, 50.0]:
        dV = on_resonance_drive.derivative(torch.tensor(t))
        assert torch.allclose(dV, torch.zeros_like(dV), atol=1e-15)


def test_nonzero_derivative_during_ramp(on_resonance_drive):
    """dV/dt is nonzero during ramp-up."""
    t_ramp = torch.tensor(1.0)  # midpoint of ramp (0 to 2 us)
    dV = on_resonance_drive.derivative(t_ramp)
    assert torch.max(torch.abs(dV)).item() > 1e-10


def test_zero_before_pulse(on_resonance_drive):
    """V(t) = 0 before t = 0."""
    V = on_resonance_drive(torch.tensor(-1.0))
    assert torch.allclose(V, torch.zeros_like(V), atol=1e-15)


def test_zero_after_pulse(on_resonance_drive):
    """V(t) = 0 well after ramp-down completes."""
    t_after = torch.tensor(60.0)  # ramp-down ends at 2 + 50 + 2 = 54 us
    V = on_resonance_drive(t_after)
    assert torch.allclose(V, torch.zeros_like(V), atol=1e-15)


def test_effective_field_angle_on_resonance(on_resonance_drive):
    """On resonance: theta = pi/2 (effective field along x)."""
    theta = on_resonance_drive.effective_field_angle()
    assert torch.allclose(theta, torch.tensor(torch.pi / 2), atol=1e-10)


def test_effective_field_angle_equal_offset():
    """omega1 = Delta: theta = pi/4 (45 degree tilt)."""
    drive = SpinLockDrive(1000.0, Ix, Iz, offset_hz=1000.0)
    theta = drive.effective_field_angle()
    assert torch.allclose(theta, torch.tensor(torch.pi / 4), atol=1e-10)


def test_effective_field_magnitude(off_resonance_drive):
    """omega_eff = sqrt(omega1^2 + Delta^2)."""
    omega1 = off_resonance_drive.omega1
    Delta = off_resonance_drive.Delta
    expected = (omega1 ** 2 + Delta ** 2) ** 0.5
    assert abs(off_resonance_drive.effective_field_magnitude() - expected) < 1e-15


def test_no_ramp_instant_on(no_ramp_drive):
    """With ramp_time=0, envelope is 1 immediately at t=0."""
    V = no_ramp_drive(torch.tensor(0.0))
    expected = no_ramp_drive.omega1 * Ix
    assert torch.allclose(V, expected, atol=1e-15)


def test_envelope_ramp_midpoint(on_resonance_drive):
    """At ramp midpoint (t = ramp/2), envelope = sin^2(pi/4) = 0.5."""
    t_mid = torch.tensor(1.0)  # ramp/2 = 2/2 = 1.0 us
    V = on_resonance_drive(t_mid)
    expected = 0.5 * on_resonance_drive.omega1 * Ix
    assert torch.allclose(V, expected, atol=1e-15)


def test_derivative_numerical_vs_analytical(on_resonance_drive):
    """Central finite difference should match analytical derivative during ramp."""
    t_test = torch.tensor(1.0)
    dt = 1e-6

    V_plus = on_resonance_drive(t_test + dt)
    V_minus = on_resonance_drive(t_test - dt)
    numeric_dV = (V_plus - V_minus) / (2 * dt)

    analytic_dV = on_resonance_drive.derivative(t_test)
    assert torch.allclose(numeric_dV, analytic_dV, atol=1e-4)


def test_multispin_operator_shape():
    """V(t) must have correct shape when using 4x4 operators from a two-spin system."""
    sys = SpinSystem(
        spins=[
            {'nucleus': '15N', 'chemical_shift_ppm': 120.0, 'gamma': -27.116e6},
            {'nucleus': '1H', 'chemical_shift_ppm': 8.5, 'gamma': 267.522e6},
        ],
        B0_field=11.7,
        J_couplings={(0, 1): -92.0},
    )
    Ix_N = sys.get_spin_operator(0, 'x')  # 4x4
    Iz_N = sys.get_spin_operator(0, 'z')  # 4x4

    drive = SpinLockDrive(1000.0, Ix_N, Iz_N, offset_hz=500.0,
                          ramp_time_us=2.0, plateau_time_us=50.0)

    V = drive(torch.tensor(10.0))
    dV = drive.derivative(torch.tensor(1.0))

    assert V.shape == (4, 4)
    assert dV.shape == (4, 4)
    assert torch.allclose(V, V.mH, atol=1e-12)

    expected = drive.omega1 * Ix_N + drive.Delta * Iz_N
    assert torch.allclose(V, expected, atol=1e-15)
