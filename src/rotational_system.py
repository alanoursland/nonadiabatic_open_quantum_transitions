"""
Rigid rotor quantum system for molecular rotational spectroscopy.

Provides energy levels, Bohr frequencies, and collisional coupling operators
for a linear rigid rotor molecule (e.g. CO, CO2, N2O, HCl). Used in the
HITRAN line-mixing feasibility estimate (Phase 0A).

The rigid rotor has energy levels E_J = B * J * (J + 1) in cm^-1, where B
is the rotational constant and J = 0, 1, 2, ... is the angular momentum
quantum number. We work in the |J> basis without m_J degeneracy, since the
perturbation ratio V/DeltaE depends only on J.
"""
import math
import torch

# Physical constants for spectroscopic unit conversion
C_CM = 2.99792458e10       # speed of light in cm/s
CM_TO_RAD_S = 2.0 * math.pi * C_CM  # 1 cm^-1 = this many rad/s (~1.884e11)


class RotationalSystem:
    """
    Linear rigid rotor with energy levels E_J = B * J * (J + 1).

    Works in the |J> basis (J = 0, 1, ..., J_max) without m_J degeneracy.
    The Hamiltonian is diagonal: H0 = diag(omega_0, omega_1, ..., omega_Jmax)
    with omega_J in rad/s.

    Interface-compatible with SpinSystem for use with Simulator and
    SpectralDensityBath (provides get_H0(), dim, device, energy_levels,
    bohr_frequencies).

    Args:
        B_cm: rotational constant in cm^-1
        J_max: maximum angular momentum quantum number (dim = J_max + 1)
        device: torch device ('cpu' or 'cuda')
    """

    def __init__(self, B_cm, J_max, device='cpu'):
        self.B_cm = B_cm
        self.J_max = J_max
        self.device = device
        self.dim = J_max + 1

        # Energy levels in cm^-1 and rad/s
        J_values = torch.arange(J_max + 1, dtype=torch.float64)
        self._energy_levels_cm = B_cm * J_values * (J_values + 1)
        self.energy_levels = (CM_TO_RAD_S * self._energy_levels_cm).to(device)

        # Diagonal Hamiltonian in rad/s
        self.H0 = torch.diag(self.energy_levels).to(torch.complex128).to(device)

        # Bohr frequency matrix: omega_{J,J'} = omega_J - omega_{J'}
        E_i = self.energy_levels.unsqueeze(1)
        E_j = self.energy_levels.unsqueeze(0)
        self.bohr_frequencies = E_i - E_j

    def get_H0(self):
        """Return the system Hamiltonian (diagonal, in rad/s)."""
        return self.H0

    def energy_levels_cm(self):
        """Return energy levels in cm^-1 (for display and HITRAN comparison)."""
        return self._energy_levels_cm

    def bohr_frequencies_cm(self):
        """Return Bohr frequency matrix in cm^-1."""
        return self._energy_levels_cm.unsqueeze(1) - self._energy_levels_cm.unsqueeze(0)

    def transition_frequencies_cm(self, delta_J=1):
        """Return R-branch (DeltaJ = +delta_J) transition frequencies in cm^-1.

        For a rigid rotor, the R-branch frequency for J -> J + delta_J is:
            nu = E(J + delta_J) - E(J)

        Returns:
            list of (J, frequency_cm) tuples
        """
        transitions = []
        for J in range(self.J_max + 1 - delta_J):
            J_upper = J + delta_J
            freq = self.B_cm * (J_upper * (J_upper + 1) - J * (J + 1))
            transitions.append((J, freq))
        return transitions


def build_collisional_coupling_operator(J_max, coupling_strength_cm, device='cpu'):
    """Build a nearest-neighbor (DeltaJ = +-1) collisional coupling operator.

    Models the dominant dipolar contribution to buffer gas collisional coupling
    between adjacent rotational levels. The coupling strength is taken as
    J-independent for the Phase 0 estimate (proportional to the pressure-
    broadening coefficient).

    Args:
        J_max: maximum J quantum number (operator dimension = J_max + 1)
        coupling_strength_cm: off-diagonal element magnitude in cm^-1
        device: torch device

    Returns:
        Hermitian (J_max+1) x (J_max+1) torch tensor (complex128)
    """
    dim = J_max + 1
    S = torch.zeros(dim, dim, dtype=torch.complex128, device=device)
    coupling_rad = coupling_strength_cm * CM_TO_RAD_S
    for J in range(J_max):
        S[J, J + 1] = coupling_rad
        S[J + 1, J] = coupling_rad
    return S
