import torch

class System:
    def __init__(self, energy_levels, dipole_matrix, device='cpu'):
        self.device = device
        self.dim = len(energy_levels)
        
        # H0: Diagonal unperturbed Hamiltonian (in angular frequency, rad/ns)
        # Assuming input is in GHz, multiply by 2*pi for angular frequency
        self.energy_levels = torch.tensor(energy_levels, dtype=torch.float64, device=device) * 2 * torch.pi
        self.H0 = torch.diag(self.energy_levels).to(torch.complex128)
        
        # Dipole matrix (unitless or scaled to match field units)
        self.dipole = torch.tensor(dipole_matrix, dtype=torch.complex128, device=device)
        
        # Calculate Bohr frequencies (omega_ij = E_i - E_j)
        E_i = self.energy_levels.unsqueeze(1)
        E_j = self.energy_levels.unsqueeze(0)
        self.bohr_frequencies = E_i - E_j

    def get_H0(self):
        return self.H0