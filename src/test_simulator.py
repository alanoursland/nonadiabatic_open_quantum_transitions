import pytest
import numpy as np
import torch
from system import System
from signals import PlateauPulse
from bath import ThermalBath
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