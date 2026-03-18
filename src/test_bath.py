import pytest
import torch
from system import System
from bath import ThermalBath

@pytest.fixture
def qubit_system():
    return System([0.0, 5.0], [[0.0, 1.0], [1.0, 0.0]])

@pytest.fixture
def thermal_bath():
    return ThermalBath(
        temperature=0.015,
        t1_times={(1, 0): 50.0},
        t2_times={(1, 0): 30.0}
    )

@pytest.fixture
def density_matrix():
    # A valid pure state |+><+|
    rho = torch.tensor([[0.5, 0.5], [0.5, 0.5]], dtype=torch.complex128)
    return rho

def test_dissipator_trace_preservation(qubit_system, thermal_bath, density_matrix):
    D_rho = thermal_bath.get_lindblad_dissipator(qubit_system, density_matrix)
    # The trace of the derivative of the density matrix must be 0
    trace_D = torch.trace(D_rho)
    assert torch.abs(trace_D) < 1e-10

def test_dissipator_hermiticity(qubit_system, thermal_bath, density_matrix):
    D_rho = thermal_bath.get_lindblad_dissipator(qubit_system, density_matrix)
    # If rho is Hermitian, D(rho) must be Hermitian to keep rho(t) Hermitian
    assert torch.allclose(D_rho, D_rho.mH, atol=1e-10)

def test_population_decay(qubit_system, thermal_bath):
    # Start entirely in excited state
    rho = torch.tensor([[0.0, 0.0], [0.0, 1.0]], dtype=torch.complex128)
    D_rho = thermal_bath.get_lindblad_dissipator(qubit_system, rho)
    
    # Excited state population should decrease (negative derivative)
    assert torch.real(D_rho[1, 1]) < 0
    # Ground state population should increase (positive derivative)
    assert torch.real(D_rho[0, 0]) > 0