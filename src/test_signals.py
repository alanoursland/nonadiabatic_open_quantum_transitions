import pytest
import torch
from signals import PlateauPulse

@pytest.fixture
def pulse():
    return PlateauPulse(
        amplitude=0.15,
        carrier_freq=5.05,
        ramp_time=10.0,
        plateau_time=100.0
    )

def test_pulse_envelope_plateau(pulse):
    # During the plateau (t=50), the envelope should be exactly 1.0
    # So the value should be exactly A0 * cos(omega * t)
    t = torch.tensor(50.0, dtype=torch.float64)
    expected_val = pulse.A0 * torch.cos(pulse.omega * t)
    actual_val = pulse(t)
    assert torch.allclose(actual_val, expected_val)

def test_pulse_envelope_boundaries(pulse):
    # Before t=0, envelope should be 0
    assert torch.abs(pulse(-1.0)) < 1e-12
    # Long after pulse, envelope should be 0
    assert torch.abs(pulse(150.0)) < 1e-12

def test_analytical_derivative_vs_numerical(pulse):
    # Test derivative during the ramp-up phase
    t_test = 5.0
    dt = 1e-6
    
    # Central finite difference approximation
    val_plus = pulse(t_test + dt)
    val_minus = pulse(t_test - dt)
    numeric_deriv = (val_plus - val_minus) / (2 * dt)
    
    # Analytical derivative
    analytic_deriv = pulse.derivative(t_test)
    
    # Should match closely
    assert torch.allclose(numeric_deriv, analytic_deriv, atol=1e-4)