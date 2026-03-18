import torch


class SpinSystem:
    """
    NMR spin system: collection of coupled spin-1/2 particles in a static B0 field.

    Constructs the Zeeman + J-coupling Hamiltonian and all spin operators
    in the product basis (|αα...⟩, |αβ...⟩, ...).
    """

    # Single spin-1/2 operators
    _Ix1 = torch.tensor([[0, 0.5], [0.5, 0]], dtype=torch.complex128)
    _Iy1 = torch.tensor([[0, -0.5j], [0.5j, 0]], dtype=torch.complex128)
    _Iz1 = torch.tensor([[0.5, 0], [0, -0.5]], dtype=torch.complex128)
    _Id1 = torch.eye(2, dtype=torch.complex128)

    def __init__(self, spins, B0_field, J_couplings=None, csa_tensors=None, device='cpu'):
        """
        Args:
            spins: list of dicts, each with:
                'nucleus': str (e.g. '15N', '1H')
                'chemical_shift_ppm': float
                'gamma': float, gyromagnetic ratio in rad/s/T
            B0_field: static field strength in Tesla
            J_couplings: dict mapping (i, j) tuples to scalar coupling in Hz
            csa_tensors: dict mapping spin index to CSA anisotropy in ppm
            device: torch device
        """
        self.device = device
        self.n_spins = len(spins)
        self.dim = 2 ** self.n_spins
        self.spins = spins
        self.B0 = B0_field
        self.J_couplings = J_couplings or {}
        self.csa_tensors = csa_tensors or {}

        # Larmor frequencies: omega_i = -gamma_i * B0 * (1 + delta_i * 1e-6)
        self.larmor_frequencies = []
        for spin in spins:
            omega = -spin['gamma'] * B0_field * (1.0 + spin['chemical_shift_ppm'] * 1e-6)
            self.larmor_frequencies.append(omega)

        # Build spin operators in product basis
        self._spin_ops = {}
        for i in range(self.n_spins):
            for comp in ('x', 'y', 'z'):
                self._spin_ops[(i, comp)] = self._build_product_operator(i, comp).to(device)

        # Build Hamiltonian: H0 = H_Zeeman + H_J
        self.H0 = self._build_hamiltonian().to(device)

        # Eigenvalues and Bohr frequencies
        self.energy_levels = torch.linalg.eigvalsh(self.H0)
        E_i = self.energy_levels.unsqueeze(1)
        E_j = self.energy_levels.unsqueeze(0)
        self.bohr_frequencies = E_i - E_j

    def _build_product_operator(self, spin_index, component):
        """Build I_{component, spin_index} in the full product basis via Kronecker products."""
        single = {'x': self._Ix1, 'y': self._Iy1, 'z': self._Iz1}[component]
        op = single if spin_index == 0 else self._Id1
        for i in range(1, self.n_spins):
            op = torch.kron(op, single if i == spin_index else self._Id1)
        return op

    def _build_hamiltonian(self):
        """H0 = sum_i omega_i * Iz_i  +  sum_{i<j} 2*pi*J_ij * (Ii . Ij)"""
        H = torch.zeros((self.dim, self.dim), dtype=torch.complex128)

        # Zeeman
        for i in range(self.n_spins):
            H = H + self.larmor_frequencies[i] * self._spin_ops[(i, 'z')]

        # Scalar (J) couplings
        for (i, j), J_hz in self.J_couplings.items():
            J_rad = 2 * torch.pi * J_hz
            for comp in ('x', 'y', 'z'):
                H = H + J_rad * (self._spin_ops[(i, comp)] @ self._spin_ops[(j, comp)])

        return H

    def get_H0(self):
        return self.H0

    def get_spin_operator(self, spin_index, component):
        """Return I_{component} for the given spin in the product basis."""
        return self._spin_ops[(spin_index, component)]
