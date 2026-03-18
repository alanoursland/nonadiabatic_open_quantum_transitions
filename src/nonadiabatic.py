import torch


class NonadiabaticDecomposition:
    """
    Landau-Lifshitz decomposition of first-order time-dependent perturbation
    theory coefficients into adiabatic and nonadiabatic parts.

    Given a system with Hamiltonian H0 and time-dependent perturbation V(t),
    the first-order transition amplitude is:

        c_k(t) = -i * integral_0^t <k|V(t')|0> exp(i*omega_{k0}*t') dt'

    Integration by parts decomposes this as c_k = a_k + b_k:

        a_k(t) = -<k|V(t)|0> / omega_{k0} * exp(i*omega_{k0}*t)
            (adiabatic: responds instantaneously to V, oscillates with Bohr freq)

        b_k(t) = (1/omega_{k0}) * integral_0^t <k|dV/dt'|0> exp(i*omega_{k0}*t') dt'
            (nonadiabatic: accumulates history through dV/dt)

    Key property: during a plateau (constant V, dV/dt = 0), b_k does not change.

    The constructor takes energy eigenvalues and eigenvectors of H0. The
    perturbation V(t) and its derivative dV/dt can be in any basis — they
    are transformed to the eigenbasis internally.
    """

    def __init__(self, energies, eigenvectors, initial_state=0, device='cpu'):
        """
        Args:
            energies: 1D tensor of energy eigenvalues (ascending order from eigh)
            eigenvectors: unitary matrix whose columns are the energy eigenstates
            initial_state: index of the initially occupied eigenstate (default: 0 = ground)
            device: torch device
        """
        self.dim = len(energies)
        self.device = device
        self.initial_state = initial_state

        self.energies = energies.to(device)
        self.U = eigenvectors.to(device)

        # Bohr frequencies: omega_k = E_k - E_initial
        E0 = self.energies[initial_state]
        self.bohr = (self.energies - E0)

        # Mask for states != initial (avoid division by zero in omega_{00}=0)
        self.active = torch.arange(self.dim, device=device) != initial_state

        # b-coefficients: start in pure initial state
        self.b_coeffs = torch.zeros(self.dim, dtype=torch.complex128, device=device)
        self.b_coeffs[initial_state] = 1.0

    def _to_eigenbasis(self, matrix):
        """Transform an operator from the computational basis to the energy eigenbasis."""
        return self.U.mH @ matrix @ self.U

    def step(self, t, dt, dV_dt_matrix):
        """
        Advance b_k over the interval [t, t+dt] using the trapezoidal rule.

        The integrand for b_k is:
            (1/omega_{k0}) * <k|dV/dt|initial> * exp(i * omega_{k0} * t')

        The matrix elements <k|dV/dt|0> are treated as constant over [t, t+dt],
        while the phase factor exp(i*omega*t') is evaluated at both endpoints.
        This gives second-order accuracy with a single eigenbasis transform.

        During a plateau (dV/dt = 0), this is exactly zero and b_k doesn't change.

        Args:
            t: current time
            dt: step size
            dV_dt_matrix: dV/dt evaluated at time t (in computational basis)
        """
        dV_eig = self._to_eigenbasis(dV_dt_matrix)
        dV_col = dV_eig[:, self.initial_state]  # <k|dV/dt|initial>

        # Trapezoidal: average phase at t and t+dt
        phases_left = torch.exp(1j * self.bohr * t)
        phases_right = torch.exp(1j * self.bohr * (t + dt))
        avg_phase = 0.5 * (phases_left + phases_right)

        increment = torch.zeros(self.dim, dtype=torch.complex128, device=self.device)
        increment[self.active] = (
            dV_col[self.active] * avg_phase[self.active] / self.bohr[self.active]
        )
        self.b_coeffs = self.b_coeffs + dt * increment

    def get_adiabatic_coefficients(self, t, V_matrix):
        """
        Return a_k(t) = -<k|V(t)|initial> / omega_{k0} * exp(i*omega_{k0}*t).

        Args:
            t: current time
            V_matrix: V(t) in the computational basis
        """
        V_eig = self._to_eigenbasis(V_matrix)
        V_col = V_eig[:, self.initial_state]

        phases = torch.exp(1j * self.bohr * t)

        a = torch.zeros(self.dim, dtype=torch.complex128, device=self.device)
        a[self.active] = -V_col[self.active] / self.bohr[self.active] * phases[self.active]
        return a

    def get_nonadiabatic_coefficients(self):
        """Return the current raw b_k coefficients (b_0 uncorrected)."""
        return self.b_coeffs.clone()

    def _normalized_b(self):
        """Return b with b_0 corrected so that sum |b_k|^2 = 1.

        The b_k for k != initial are first-order accurate from the integral.
        b_0 is set to sqrt(1 - sum_{k>0} |b_k|^2) to conserve probability.
        This is a partial second-order correction that ensures Tr(sigma_nad) = 1.
        """
        b = self.b_coeffs.clone()
        sum_excited = torch.sum(torch.abs(b[self.active]) ** 2)
        b0_sq = torch.clamp(1.0 - sum_excited, min=0.0)
        b[self.initial_state] = torch.sqrt(b0_sq)
        return b

    def get_nonadiabatic_density_matrix(self):
        """Return sigma_nad = |b><b| with trace = 1 (b_0 normalized)."""
        b = self._normalized_b()
        return b.unsqueeze(1) * b.conj().unsqueeze(0)

    def get_nonadiabatic_populations(self):
        """Return |b_k|^2 with sum = 1 (b_0 normalized)."""
        b = self._normalized_b()
        return torch.abs(b) ** 2
