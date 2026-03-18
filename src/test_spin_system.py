import pytest
import torch
from spin_system import SpinSystem


@pytest.fixture
def single_spin():
    """Single 15N spin in 11.7 T field."""
    return SpinSystem(
        spins=[{'nucleus': '15N', 'chemical_shift_ppm': 120.0, 'gamma': -27.116e6}],
        B0_field=11.7,
    )


@pytest.fixture
def two_spin():
    """15N-1H pair in 11.7 T field with J-coupling."""
    return SpinSystem(
        spins=[
            {'nucleus': '15N', 'chemical_shift_ppm': 120.0, 'gamma': -27.116e6},
            {'nucleus': '1H', 'chemical_shift_ppm': 8.5, 'gamma': 267.522e6},
        ],
        B0_field=11.7,
        J_couplings={(0, 1): -92.0},
        csa_tensors={0: -170.0},
    )


def test_dimensions(single_spin, two_spin):
    assert single_spin.dim == 2
    assert two_spin.dim == 4


def test_larmor_frequency(single_spin):
    gamma = -27.116e6
    B0 = 11.7
    delta = 120.0
    expected = -gamma * B0 * (1 + delta * 1e-6)
    assert abs(single_spin.larmor_frequencies[0] - expected) < 1.0


def test_hamiltonian_hermiticity(two_spin):
    H = two_spin.get_H0()
    assert torch.allclose(H, H.mH, atol=1e-10)


def test_spin_commutation_relations(two_spin):
    """[Ix, Iy] = i*Iz for each spin."""
    for spin_idx in range(two_spin.n_spins):
        Ix = two_spin.get_spin_operator(spin_idx, 'x')
        Iy = two_spin.get_spin_operator(spin_idx, 'y')
        Iz = two_spin.get_spin_operator(spin_idx, 'z')
        comm = Ix @ Iy - Iy @ Ix
        assert torch.allclose(comm, 1j * Iz, atol=1e-12)


def test_spin_operators_hermitian(two_spin):
    for spin_idx in range(two_spin.n_spins):
        for comp in ('x', 'y', 'z'):
            op = two_spin.get_spin_operator(spin_idx, comp)
            assert torch.allclose(op, op.mH, atol=1e-12)


def test_spin_operators_traceless(two_spin):
    for spin_idx in range(two_spin.n_spins):
        for comp in ('x', 'y', 'z'):
            op = two_spin.get_spin_operator(spin_idx, comp)
            assert torch.abs(torch.trace(op)) < 1e-12


def test_bohr_frequencies_single_spin(single_spin):
    """For a single spin, Bohr frequencies are +-omega_Larmor."""
    omega = single_spin.larmor_frequencies[0]
    bohr = single_spin.bohr_frequencies
    # eigvalsh returns ascending order: [-omega/2, +omega/2] when omega > 0
    # bohr[0,1] = E_0 - E_1 = -omega, bohr[1,0] = omega
    assert torch.allclose(
        bohr[1, 0], torch.tensor(omega, dtype=torch.float64), atol=1.0
    )
    assert torch.allclose(
        bohr[0, 1], torch.tensor(-omega, dtype=torch.float64), atol=1.0
    )


def test_j_coupling_splits_levels(two_spin):
    """J-coupling lifts degeneracy: the four eigenvalues should all be distinct."""
    eigs = two_spin.energy_levels
    # Sort and check pairwise differences are nonzero
    sorted_eigs, _ = torch.sort(eigs)
    diffs = sorted_eigs[1:] - sorted_eigs[:-1]
    assert torch.all(torch.abs(diffs) > 1.0)  # at least 1 rad/s apart


def test_different_spins_commute(two_spin):
    """Operators on different spins must commute: [I_{a,i}, I_{b,j}] = 0 for i != j."""
    for comp_a in ('x', 'y', 'z'):
        for comp_b in ('x', 'y', 'z'):
            op_a = two_spin.get_spin_operator(0, comp_a)
            op_b = two_spin.get_spin_operator(1, comp_b)
            comm = op_a @ op_b - op_b @ op_a
            assert torch.allclose(comm, torch.zeros_like(comm), atol=1e-12)
