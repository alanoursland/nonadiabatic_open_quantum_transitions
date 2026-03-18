import pytest
import torch
from nonadiabatic import NonadiabaticDecomposition


# --- Helpers for a two-level system (H0 already diagonal) ---

def make_two_level(omega):
    """Two-level system with gap omega. H0 = diag(0, omega), eigenvectors = I."""
    energies = torch.tensor([0.0, omega], dtype=torch.float64)
    eigvecs = torch.eye(2, dtype=torch.complex128)
    return energies, eigvecs


def off_diag_perturbation(v):
    """V = [[0, v], [v, 0]] — purely off-diagonal coupling."""
    return torch.tensor([[0, v], [v, 0]], dtype=torch.complex128)


# --- Tests ---


class TestPlateauInvariance:
    """During a plateau (dV/dt = 0), b_k must not change."""

    def test_b_unchanged_after_plateau_steps(self):
        energies, eigvecs = make_two_level(5.0)
        nad = NonadiabaticDecomposition(energies, eigvecs)

        # Artificially set b_1 to a nonzero value (as if a ramp just finished)
        nad.b_coeffs[1] = 0.1 + 0.05j
        b_before = nad.get_nonadiabatic_coefficients()

        # Step through 100 plateau steps with dV/dt = 0
        dV_zero = torch.zeros(2, 2, dtype=torch.complex128)
        for i in range(100):
            nad.step(t=10.0 + i * 0.01, dt=0.01, dV_dt_matrix=dV_zero)

        b_after = nad.get_nonadiabatic_coefficients()
        assert torch.allclose(b_before, b_after, atol=1e-15)


class TestLinearRamp:
    """
    Linear ramp: V(t) = v * (t/T) * sigma_x for 0 <= t <= T.
    dV/dt = (v/T) * sigma_x (constant during ramp).

    Analytical result:
        b_1(T) = v / (omega^2 * T) * (exp(i*omega*T) - 1) / i
    """

    def setup_method(self):
        self.omega = 5.0
        self.v = 0.3
        self.T = 2.0
        self.N = 5000  # steps for accuracy

    def analytical_b1(self):
        omega, v, T = self.omega, self.v, self.T
        phase = torch.exp(1j * torch.tensor(omega * T, dtype=torch.float64))
        return v * (phase - 1) / (1j * omega**2 * T)

    def run_ramp(self):
        energies, eigvecs = make_two_level(self.omega)
        nad = NonadiabaticDecomposition(energies, eigvecs)

        dt = self.T / self.N
        dV = off_diag_perturbation(self.v / self.T)  # constant dV/dt during ramp

        for i in range(self.N):
            t = i * dt
            nad.step(t, dt, dV)

        return nad

    def test_b1_matches_analytical(self):
        nad = self.run_ramp()
        b = nad.get_nonadiabatic_coefficients()
        expected = self.analytical_b1()
        assert torch.allclose(b[1], expected, atol=1e-4)

    def test_b0_unchanged(self):
        """b_0 (initial state) should remain exactly 1."""
        nad = self.run_ramp()
        b = nad.get_nonadiabatic_coefficients()
        assert torch.allclose(b[0], torch.tensor(1.0 + 0j, dtype=torch.complex128))


class TestDecompositionConsistency:
    """
    Verify a_k + b_k = c_k where c_k is computed by direct numerical
    integration of the TDPT integral.

    Uses a sin^2 ramp: V(t) = v * sin^2(pi*t / (2T)) * sigma_x for 0 <= t <= T.
    """

    def test_a_plus_b_equals_c(self):
        omega = 4.0
        v = 0.2
        T = 3.0
        N = 5000

        energies, eigvecs = make_two_level(omega)
        nad = NonadiabaticDecomposition(energies, eigvecs)

        dt = T / N

        # Also accumulate c_k by direct integration: c_k = -i * integral V_{k0} e^{iwt} dt
        c_direct = torch.zeros(2, dtype=torch.complex128)

        for i in range(N):
            t = i * dt

            # Envelope: f(t) = sin^2(pi*t / (2*T))
            f = torch.sin(torch.tensor(torch.pi * t / (2 * T))) ** 2
            df = (torch.pi / T) * torch.sin(torch.tensor(torch.pi * t / (2 * T))) * \
                 torch.cos(torch.tensor(torch.pi * t / (2 * T)))

            V = off_diag_perturbation(v * f.item())
            dV = off_diag_perturbation(v * df.item())

            nad.step(t, dt, dV)

            # Direct integration: c_1 += -i * V_{10} * exp(i*omega*t) * dt
            V10 = v * f
            c_direct[1] += -1j * V10 * torch.exp(1j * torch.tensor(omega * t)) * dt

        # Get final V(T) for adiabatic coefficients
        f_T = torch.sin(torch.tensor(torch.pi * T / (2 * T))) ** 2  # = 1.0
        V_T = off_diag_perturbation(v * f_T.item())

        a = nad.get_adiabatic_coefficients(T, V_T)
        b = nad.get_nonadiabatic_coefficients()

        # a_1 + b_1 should equal c_1
        c_decomposed = a[1] + b[1]
        assert torch.allclose(c_decomposed, c_direct[1], atol=1e-4), \
            f"a+b = {c_decomposed}, c_direct = {c_direct[1]}"


class TestDensityMatrixAndPopulations:
    """Test the derived quantities from b_k."""

    def test_density_matrix_trace_one(self):
        """Density matrix must have trace = 1 (b_0 normalized)."""
        energies, eigvecs = make_two_level(5.0)
        nad = NonadiabaticDecomposition(energies, eigvecs)
        nad.b_coeffs = torch.tensor([1.0 + 0j, 0.3 - 0.2j], dtype=torch.complex128)

        sigma = nad.get_nonadiabatic_density_matrix()
        assert torch.allclose(
            torch.trace(sigma).real,
            torch.tensor(1.0, dtype=torch.float64),
            atol=1e-14,
        )

    def test_density_matrix_is_normalized_outer_product(self):
        """sigma_nad = |b_norm><b_norm| where b_norm has corrected b_0."""
        energies, eigvecs = make_two_level(5.0)
        nad = NonadiabaticDecomposition(energies, eigvecs)
        nad.b_coeffs = torch.tensor([1.0 + 0j, 0.3 - 0.2j], dtype=torch.complex128)

        sigma = nad.get_nonadiabatic_density_matrix()
        b_norm = nad._normalized_b()
        expected = b_norm.unsqueeze(1) * b_norm.conj().unsqueeze(0)
        assert torch.allclose(sigma, expected, atol=1e-15)

    def test_density_matrix_hermitian(self):
        energies, eigvecs = make_two_level(5.0)
        nad = NonadiabaticDecomposition(energies, eigvecs)
        nad.b_coeffs = torch.tensor([1.0 + 0j, 0.3 - 0.2j], dtype=torch.complex128)

        sigma = nad.get_nonadiabatic_density_matrix()
        assert torch.allclose(sigma, sigma.mH, atol=1e-15)

    def test_populations_sum_to_one(self):
        """Populations must sum to 1 after normalization."""
        energies, eigvecs = make_two_level(5.0)
        nad = NonadiabaticDecomposition(energies, eigvecs)
        nad.b_coeffs = torch.tensor([1.0 + 0j, 0.3 - 0.2j], dtype=torch.complex128)

        pops = nad.get_nonadiabatic_populations()
        assert torch.allclose(
            torch.sum(pops),
            torch.tensor(1.0, dtype=torch.float64),
            atol=1e-14,
        )
        # Excited state population is |0.3-0.2j|^2 = 0.13
        assert torch.allclose(pops[1], torch.tensor(0.13, dtype=torch.float64), atol=1e-15)
        # Ground state population is 1 - 0.13 = 0.87
        assert torch.allclose(pops[0], torch.tensor(0.87, dtype=torch.float64), atol=1e-14)

    def test_initial_state_populations(self):
        """Initially, all population is in the ground state."""
        energies, eigvecs = make_two_level(5.0)
        nad = NonadiabaticDecomposition(energies, eigvecs)

        pops = nad.get_nonadiabatic_populations()
        assert pops[0].item() == 1.0
        assert pops[1].item() == 0.0

    def test_raw_coefficients_unchanged_by_normalization(self):
        """get_nonadiabatic_coefficients() returns raw b (b_0 still 1.0)."""
        energies, eigvecs = make_two_level(5.0)
        nad = NonadiabaticDecomposition(energies, eigvecs)
        nad.b_coeffs = torch.tensor([1.0 + 0j, 0.3 - 0.2j], dtype=torch.complex128)

        b_raw = nad.get_nonadiabatic_coefficients()
        assert b_raw[0] == 1.0 + 0j  # raw, not corrected
        assert b_raw[1] == 0.3 - 0.2j


class TestNonDiagonalH0:
    """
    Test with a Hamiltonian that requires nontrivial eigenvectors,
    confirming the eigenbasis transformation works correctly.
    """

    def test_with_coupled_system(self):
        # H0 = [[1, 0.5], [0.5, 3]] — not diagonal
        H0 = torch.tensor([[1.0, 0.5], [0.5, 3.0]], dtype=torch.complex128)
        energies, eigvecs = torch.linalg.eigh(H0)

        nad = NonadiabaticDecomposition(energies.to(torch.float64), eigvecs)

        omega = (energies[1] - energies[0]).item()

        # Perturbation in computational basis
        V_comp = torch.tensor([[0, 0.2], [0.2, 0]], dtype=torch.complex128)
        # Transform to eigenbasis to get V_{10}
        V_eig = eigvecs.mH @ V_comp @ eigvecs
        V10 = V_eig[1, 0]

        # Linear ramp dV/dt = V_comp / T (constant)
        T = 2.0
        N = 5000
        dt = T / N
        dV = V_comp / T

        for i in range(N):
            nad.step(i * dt, dt, dV)

        b = nad.get_nonadiabatic_coefficients()

        # Analytical: b_1 = V10 * (exp(i*omega*T) - 1) / (i * omega^2 * T)
        expected = V10 * (torch.exp(1j * torch.tensor(omega * T)) - 1) / (1j * omega**2 * T)
        assert torch.allclose(b[1], expected, atol=1e-4), \
            f"b[1] = {b[1]}, expected = {expected}"
