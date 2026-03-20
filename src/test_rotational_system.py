"""Tests for rotational_system.py: rigid rotor energy levels and coupling operator."""
import math
import pytest
import torch

from rotational_system import (
    RotationalSystem, build_collisional_coupling_operator,
    C_CM, CM_TO_RAD_S,
)


# --- Fixtures ---

@pytest.fixture
def co_system():
    """CO molecule, B = 1.9313 cm^-1, J_max = 10."""
    return RotationalSystem(B_cm=1.9313, J_max=10)


@pytest.fixture
def co2_system():
    """CO2 molecule, B = 0.3902 cm^-1, J_max=20."""
    return RotationalSystem(B_cm=0.3902, J_max=20)


# --- Energy level tests ---

class TestEnergyLevels:

    def test_ground_state_energy_zero(self, co_system):
        """E(J=0) = B * 0 * 1 = 0."""
        E_cm = co_system.energy_levels_cm()
        assert E_cm[0].item() == pytest.approx(0.0)

    def test_energy_level_formula(self, co_system):
        """E(J) = B * J * (J+1) for several J values."""
        E_cm = co_system.energy_levels_cm()
        B = co_system.B_cm
        for J in range(co_system.J_max + 1):
            assert E_cm[J].item() == pytest.approx(B * J * (J + 1), rel=1e-10)

    def test_co_energy_levels(self, co_system):
        """Verify specific CO energy levels in cm^-1."""
        E_cm = co_system.energy_levels_cm()
        assert E_cm[1].item() == pytest.approx(3.8626, rel=1e-4)
        assert E_cm[2].item() == pytest.approx(11.5878, rel=1e-4)
        assert E_cm[3].item() == pytest.approx(23.1756, rel=1e-4)

    def test_unit_conversion(self, co_system):
        """E(J=1) in rad/s = CM_TO_RAD_S * B * 2."""
        E_rad = co_system.energy_levels[1].item()
        expected = CM_TO_RAD_S * co_system.B_cm * 2.0
        assert E_rad == pytest.approx(expected, rel=1e-10)

    def test_energy_levels_ascending(self, co_system):
        """Energy levels are strictly increasing with J."""
        E = co_system.energy_levels
        for J in range(co_system.J_max):
            assert E[J + 1].item() > E[J].item()


# --- Hamiltonian tests ---

class TestHamiltonian:

    def test_dimension(self, co_system):
        """dim = J_max + 1."""
        assert co_system.dim == 11

    def test_hamiltonian_diagonal(self, co_system):
        """H0 is diagonal (rigid rotor in the |J> eigenbasis)."""
        H0 = co_system.get_H0()
        off_diag = H0 - torch.diag(torch.diag(H0))
        assert torch.allclose(off_diag, torch.zeros_like(off_diag), atol=1e-20)

    def test_hamiltonian_hermitian(self, co_system):
        """H0 = H0†."""
        H0 = co_system.get_H0()
        assert torch.allclose(H0, H0.conj().T, atol=1e-20)

    def test_hamiltonian_diagonal_matches_energy_levels(self, co_system):
        """Diagonal of H0 matches energy_levels."""
        H0 = co_system.get_H0()
        diag = torch.diag(H0).real
        assert torch.allclose(diag, co_system.energy_levels, atol=1e-10)


# --- Bohr frequency tests ---

class TestBohrFrequencies:

    def test_antisymmetric(self, co_system):
        """bohr[J,J'] = -bohr[J',J]."""
        bf = co_system.bohr_frequencies
        assert torch.allclose(bf, -bf.T, atol=1e-10)

    def test_diagonal_zero(self, co_system):
        """bohr[J,J] = 0."""
        bf = co_system.bohr_frequencies
        assert torch.allclose(torch.diag(bf), torch.zeros(co_system.dim, dtype=torch.float64), atol=1e-10)

    def test_adjacent_line_spacing(self, co_system):
        """R-branch transitions are spaced by 2*B cm^-1 (uniform for rigid rotor)."""
        transitions = co_system.transition_frequencies_cm(delta_J=1)
        B = co_system.B_cm
        for i in range(len(transitions) - 1):
            spacing = transitions[i + 1][1] - transitions[i][1]
            assert spacing == pytest.approx(2.0 * B, rel=1e-10)

    def test_r_branch_frequency(self, co_system):
        """R-branch frequency for J->J+1 is 2*B*(J+1) cm^-1."""
        transitions = co_system.transition_frequencies_cm(delta_J=1)
        B = co_system.B_cm
        for J, freq in transitions:
            assert freq == pytest.approx(2.0 * B * (J + 1), rel=1e-10)


# --- Coupling operator tests ---

class TestCouplingOperator:

    def test_hermitian(self):
        """Coupling operator is Hermitian."""
        S = build_collisional_coupling_operator(J_max=5, coupling_strength_cm=0.065)
        assert torch.allclose(S, S.conj().T, atol=1e-20)

    def test_tridiagonal(self):
        """Only nearest-neighbor elements (DeltaJ = +-1) are nonzero."""
        S = build_collisional_coupling_operator(J_max=5, coupling_strength_cm=0.065)
        dim = 6
        for i in range(dim):
            for j in range(dim):
                if abs(i - j) > 1:
                    assert S[i, j].item() == 0.0

    def test_zero_diagonal(self):
        """Diagonal elements are zero (no self-coupling)."""
        S = build_collisional_coupling_operator(J_max=5, coupling_strength_cm=0.065)
        assert torch.allclose(torch.diag(S), torch.zeros(6, dtype=torch.complex128), atol=1e-20)

    def test_magnitude(self):
        """Off-diagonal elements match coupling_strength * CM_TO_RAD_S."""
        coupling_cm = 0.065
        S = build_collisional_coupling_operator(J_max=5, coupling_strength_cm=coupling_cm)
        expected = coupling_cm * CM_TO_RAD_S
        assert S[0, 1].real.item() == pytest.approx(expected, rel=1e-10)
        assert S[3, 4].real.item() == pytest.approx(expected, rel=1e-10)

    def test_dimension(self):
        """Shape is (J_max+1, J_max+1)."""
        S = build_collisional_coupling_operator(J_max=8, coupling_strength_cm=0.1)
        assert S.shape == (9, 9)


# --- CO2 cross-checks ---

class TestCO2:

    def test_co2_dimension(self, co2_system):
        assert co2_system.dim == 21

    def test_co2_smaller_spacing(self, co_system, co2_system):
        """CO2 (smaller B) has smaller energy level spacings than CO."""
        E_co = co_system.energy_levels_cm()
        E_co2 = co2_system.energy_levels_cm()
        # E(J=1) for CO2 < E(J=1) for CO
        assert E_co2[1].item() < E_co[1].item()
