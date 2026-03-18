import torch
import numpy as np


# --- Physical constants ---
HBAR = 1.054571817e-34   # J·s
KB = 1.380649e-23        # J/K


class ThermalBath:
    """Phenomenological bath with T1/T2 rates. Kept for backward compatibility."""

    def __init__(self, temperature, t1_times, t2_times, device='cpu'):
        self.device = device
        self.temperature = temperature
        self.t1 = t1_times
        self.t2 = t2_times

    def get_lindblad_dissipator(self, system, rho):
        D_rho = torch.zeros_like(rho, dtype=torch.complex128, device=self.device)

        for (i, j), t1 in self.t1.items():
            rate = 1.0 / t1
            L = torch.zeros((system.dim, system.dim), dtype=torch.complex128, device=self.device)
            L[j, i] = 1.0
            L_dag = L.mH
            term1 = L @ rho @ L_dag
            term2 = 0.5 * (L_dag @ L @ rho + rho @ L_dag @ L)
            D_rho += rate * (term1 - term2)

        for (i, j), t2 in self.t2.items():
            rate_t1_i = sum([1 / t for (k, l), t in self.t1.items() if k == i])
            rate_t1_j = sum([1 / t for (k, l), t in self.t1.items() if k == j])
            gamma_phi = (1.0 / t2) - 0.5 * rate_t1_i - 0.5 * rate_t1_j
            if gamma_phi > 0:
                L_z = torch.zeros((system.dim, system.dim), dtype=torch.complex128, device=self.device)
                L_z[i, i] = 1.0
                L_z[j, j] = -1.0
                term1 = L_z @ rho @ L_z.mH
                term2 = 0.5 * (L_z.mH @ L_z @ rho + rho @ L_z.mH @ L_z)
                D_rho += (gamma_phi / 2.0) * (term1 - term2)

        return D_rho


class SpectralDensityBath:
    """
    Microscopic bath defined by a spectral density function J(omega).

    Computes dissipators via the Redfield tensor R_{abcd} built element-wise
    in the energy eigenbasis of H0. The master equation is:

        d(rho_ab)/dt = -i * omega_ab * rho_ab + sum_{cd} R_{ab,cd} * rho_{cd}

    where the first term is the coherent evolution (handled by the simulator)
    and the R-tensor encodes all dissipative effects.

    Three dissipator methods are provided:
        1. Full Redfield   — uses R_{abcd} directly (may violate positivity)
        2. Secular Lindblad — zeros R_{abcd} when |omega_ab - omega_cd| > tol
        3. CGME             — damps R_{abcd} by sinc((omega_ab - omega_cd)*Tc/2)

    Usage:
        bath = SpectralDensityBath(J_fn, T, [S_operator])
        bath.precompute(energies, eigenvectors)
        D_rho = bath.get_redfield_dissipator(rho_in_eigenbasis)

    Reference: Breuer & Petruccione, "The Theory of Open Quantum Systems",
    Eqs. 3.155-3.168; Campaioli et al., PRX Quantum 5, 020202 (2024).
    """

    def __init__(self, spectral_density_fn, temperature_K,
                 coupling_operators, device='cpu'):
        """
        Args:
            spectral_density_fn: callable J(omega) -> float, the bath spectral
                density. Must accept omega in the same angular frequency units
                as the system Hamiltonian eigenvalues.
            temperature_K: bath temperature in Kelvin
            coupling_operators: list of system operators S_alpha that couple
                to the bath, in the computational basis
            device: torch device
        """
        self.J = spectral_density_fn
        self.T = temperature_K
        self.coupling_ops = [op.to(torch.complex128).to(device) for op in coupling_operators]
        self.device = device
        self._precomputed = False

    # --- Rate functions ---

    def _bose_einstein(self, omega):
        """Bose-Einstein occupation number n(omega).

        n(omega) = 1 / (exp(hbar*omega / kT) - 1)

        Handles edge cases: omega -> 0 (classical limit kT/hbar*omega),
        large positive x (n -> 0), large negative x (n -> -1).
        """
        if abs(omega) < 1e-30:
            return 0.0
        x = HBAR * omega / (KB * self.T)
        if abs(x) < 1e-10:
            # Taylor expansion: 1/(e^x - 1) ~ 1/x - 1/2 + x/12 - ...
            return 1.0 / x - 0.5
        if x > 500:
            return 0.0
        if x < -500:
            return -1.0
        return 1.0 / (np.exp(x) - 1.0)

    def _gamma(self, omega):
        """One-sided Fourier transform of the bath correlation function (real part).

        gamma(omega) = 2 * Re[Gamma(omega)]

        For omega > 0 (emission):  gamma = J(omega) * [n(omega) + 1]
        For omega < 0 (absorption): gamma = J(|omega|) * n(|omega|)
        For omega = 0:              gamma = J(0) * kT / (hbar * omega) limit

        Satisfies detailed balance: gamma(omega) / gamma(-omega) = exp(hbar*omega/kT).
        """
        if abs(omega) < 1e-30:
            if self.T > 0:
                return self.J(0.0) * KB * self.T / HBAR
            return 0.0
        j_val = self.J(abs(omega))
        n_val = self._bose_einstein(abs(omega))
        if omega > 0:
            return j_val * (n_val + 1.0)
        else:
            return j_val * n_val

    # --- Redfield tensor construction ---

    def precompute(self, energies, eigenvectors):
        """Build the Redfield tensor from the system's eigendata.

        Must be called before any dissipator method.

        Args:
            energies: 1D tensor, energy eigenvalues (ascending, from eigh)
            eigenvectors: 2D tensor, unitary matrix whose columns are eigenstates

        The Redfield tensor is built from the standard element-wise formula
        in the energy eigenbasis. For each coupling operator S_alpha:

            R_{ab,cd} = sum_alpha [
                gamma(omega_{ca}) * S_{ac} * S_{db}       (term 1: "+" half-FT)
              + gamma(omega_{db})  * S_{db} * S_{ac}       (term 2: "-" half-FT, real gamma)
              - delta_{bd} * sum_n gamma(omega_{na}) * S_{an} * S_{nc}    (term 3: trace-preserving)
              - delta_{ac} * sum_n gamma(omega_{nb}) * S_{dn} * S_{nb}    (term 4: trace-preserving)
            ]

        where S_{ac} = <a|S|c> in the eigenbasis and omega_{ca} = E_c - E_a.

        When gamma is real (no Lamb shift), terms 1 and 2 have the same operator
        content but different frequency arguments. Terms 3 and 4 enforce
        Tr[D[rho]] = 0 for any rho.
        """
        self.energies = energies.to(self.device)
        self.U = eigenvectors.to(torch.complex128).to(self.device)
        dim = len(energies)
        self.dim = dim

        # Bohr frequency matrix: bohr[a,b] = E_a - E_b
        E = energies.to(torch.float64).to(self.device)
        self.bohr = E.unsqueeze(0) - E.unsqueeze(1)

        # Transform coupling operators to eigenbasis: S_eig = U^dag S U
        self.coupling_ops_eig = []
        for S in self.coupling_ops:
            S_eig = self.U.mH @ S @ self.U
            self.coupling_ops_eig.append(S_eig)

        # Precompute gamma at every Bohr frequency
        # gamma_mat[a,b] = gamma(omega_{ab}) = gamma(E_a - E_b)
        gamma_mat = torch.zeros(dim, dim, dtype=torch.float64, device=self.device)
        for a in range(dim):
            for b in range(dim):
                gamma_mat[a, b] = self._gamma(self.bohr[a, b].item())

        # Build the full Redfield tensor R_{ab,cd}
        R = torch.zeros(dim, dim, dim, dim, dtype=torch.complex128, device=self.device)

        for S in self.coupling_ops_eig:
            for a in range(dim):
                for b in range(dim):
                    for c in range(dim):
                        for d in range(dim):
                            # Term 1: gamma(omega_{ca}) * S_{ac} * S_{db}
                            val = gamma_mat[c, a] * S[a, c] * S[d, b]

                            # Term 2: gamma(omega_{db}) * S_{ac} * S_{db}
                            # (same operator content, different freq argument)
                            val = val + gamma_mat[d, b] * S[a, c] * S[d, b]

                            # Term 3: -delta_{bd} * sum_n gamma(omega_{cn}) * S_{an} * S_{nc}
                            if b == d:
                                for n in range(dim):
                                    val = val - gamma_mat[c, n] * S[a, n] * S[n, c]

                            # Term 4: -delta_{ac} * sum_n gamma(omega_{dn}) * S_{dn} * S_{nb}
                            if a == c:
                                for n in range(dim):
                                    val = val - gamma_mat[d, n] * S[d, n] * S[n, b]

                            R[a, b, c, d] += val

        self.R_full = R

        # --- Secular mask ---
        # Keep only elements where |omega_ab - omega_cd| < tol
        bohr_diff = torch.zeros(dim, dim, dim, dim, dtype=torch.float64, device=self.device)
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        bohr_diff[a, b, c, d] = abs(self.bohr[a, b].item() - self.bohr[c, d].item())

        self._bohr_diff = bohr_diff
        secular_tol = max(1e-6, 1e-6 * torch.max(torch.abs(self.bohr)).item())
        self.secular_mask = bohr_diff < secular_tol

        self.R_secular = R.clone()
        self.R_secular[~self.secular_mask] = 0.0

        # Reshape for fast matmul: R as (dim^2, dim^2) superoperator
        self.R_full_super = R.reshape(dim * dim, dim * dim)
        self.R_secular_super = self.R_secular.reshape(dim * dim, dim * dim)

        self._precomputed = True

    # --- Dissipator methods ---

    def _apply_superoperator(self, R_super, rho):
        """Compute D[rho]_{ab} = sum_{cd} R_{ab,cd} * rho_{cd}."""
        dim = self.dim
        rho_vec = rho.reshape(dim * dim)
        D_vec = R_super @ rho_vec
        return D_vec.reshape(dim, dim)

    def get_redfield_dissipator(self, rho):
        """Full non-secular Redfield dissipator.

        rho must be in the energy eigenbasis. May violate positivity.
        """
        assert self._precomputed, "Call precompute() first"
        return self._apply_superoperator(self.R_full_super, rho)

    def get_secular_lindblad_dissipator(self, rho):
        """Secular (GKSL) Lindblad dissipator.

        rho must be in the energy eigenbasis. Guaranteed completely positive.
        Equivalent to the full Redfield tensor with all non-secular terms
        (|omega_ab - omega_cd| > 0) set to zero.
        """
        assert self._precomputed, "Call precompute() first"
        return self._apply_superoperator(self.R_secular_super, rho)

    def get_cgme_dissipator(self, rho, delta_tau):
        """Coarse-grained master equation dissipator.

        R_cgme_{abcd} = R_{abcd} * sinc((omega_ab - omega_cd) * delta_tau / 2)

        where sinc(x) = sin(x)/x (the unnormalized sinc).

        Interpolates between full Redfield (delta_tau -> 0) and
        secular Lindblad (delta_tau -> inf). Preserves complete positivity
        for appropriate choice of delta_tau (Schaller & Brandes, PRA 2008).

        rho must be in the energy eigenbasis.
        """
        assert self._precomputed, "Call precompute() first"
        dim = self.dim

        # Build sinc damping factor for each tensor element
        # sinc(x) = sin(x)/x; np.sinc(y) = sin(pi*y)/(pi*y), so use np.sinc(x/pi)
        sinc_mat = torch.zeros(dim, dim, dim, dim, dtype=torch.float64, device=self.device)
        for a in range(dim):
            for b in range(dim):
                for c in range(dim):
                    for d in range(dim):
                        x = self._bohr_diff[a, b, c, d].item() * delta_tau / 2.0
                        if abs(x) < 1e-30:
                            sinc_mat[a, b, c, d] = 1.0
                        else:
                            sinc_mat[a, b, c, d] = np.sinc(x / np.pi)

        R_cgme = self.R_full * sinc_mat.to(torch.complex128)
        R_cgme_super = R_cgme.reshape(dim * dim, dim * dim)
        return self._apply_superoperator(R_cgme_super, rho)


# --- Spectral density models ---

class DebyeSpectralDensity:
    """Debye (Lorentzian) spectral density: J(omega) = eta * omega * tau_c / (1 + omega^2 * tau_c^2).

    This is the simplest physically motivated spectral density, corresponding
    to exponentially decaying bath correlations with correlation time tau_c.
    """

    def __init__(self, eta, tau_c):
        """
        Args:
            eta: dimensionless coupling strength (overall prefactor)
            tau_c: correlation time in the same time units as 1/omega
        """
        self.eta = eta
        self.tau_c = tau_c

    def __call__(self, omega):
        return self.eta * abs(omega) * self.tau_c / (1.0 + omega ** 2 * self.tau_c ** 2)


class LipariSzaboSpectralDensity:
    """Lipari-Szabo model-free spectral density for NMR relaxation.

    J(omega) = (2/5) * [S^2 * tau_c / (1 + omega^2 * tau_c^2)
                       + (1 - S^2) * tau / (1 + omega^2 * tau^2)]

    where 1/tau = 1/tau_c + 1/tau_e.
    """

    def __init__(self, tau_c, S2, tau_e=0.0):
        """
        Args:
            tau_c: overall rotational correlation time (seconds)
            S2: generalized order parameter (0 <= S2 <= 1)
            tau_e: internal motion correlation time (seconds)
        """
        self.tau_c = tau_c
        self.S2 = S2
        self.tau_e = tau_e
        if tau_e > 0 and tau_c > 0:
            self.tau = 1.0 / (1.0 / tau_c + 1.0 / tau_e)
        else:
            self.tau = 0.0

    def __call__(self, omega):
        w2 = omega ** 2
        term1 = self.S2 * self.tau_c / (1.0 + w2 * self.tau_c ** 2)
        term2 = 0.0
        if self.tau > 0:
            term2 = (1.0 - self.S2) * self.tau / (1.0 + w2 * self.tau ** 2)
        return 0.4 * (term1 + term2)