import pytest
import torch
from system import System

@pytest.fixture
def qutrit():
    energy_levels = [0.0, 5.0, 9.8]
    dipole_matrix = [
        [0.0, 1.0, 0.0],
        [1.0, 0.0, 1.4],
        [0.0, 1.4, 0.0]
    ]
    return System(energy_levels, dipole_matrix)

def test_system_initialization(qutrit):
    assert qutrit.dim == 3
    assert qutrit.H0.shape == (3, 3)
    assert qutrit.dipole.shape == (3, 3)

def test_angular_frequency_conversion(qutrit):
    # Energy levels input in GHz should be converted to angular frequency (2*pi*f)
    expected_E1 = 5.0 * 2 * torch.pi
    assert torch.allclose(torch.real(qutrit.H0[1, 1]), torch.tensor(expected_E1, dtype=torch.float64))

def test_bohr_frequencies(qutrit):
    omega = qutrit.bohr_frequencies
    # omega_10 = E1 - E0
    expected_w10 = (5.0 - 0.0) * 2 * torch.pi
    assert torch.allclose(omega[1, 0], torch.tensor(expected_w10, dtype=torch.float64))
    # omega_ij = -omega_ji
    assert torch.allclose(omega[2, 1], -omega[1, 2])