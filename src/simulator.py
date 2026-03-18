import math
import torch
import matplotlib.pyplot as plt


class Result:
    def __init__(self, times, states):
        self.times = times
        self.states = states  # Tensor of shape (len(times), dim, dim)
        self.nonadiabatic_populations = None  # type: torch.Tensor | None

    def get_populations(self, state_index):
        # Extract diagonal elements (probabilities)
        return torch.real(self.states[:, state_index, state_index])

    def plot_populations(self, states=[0, 1, 2], label_prefix=""):
        times_np = self.times.cpu().numpy()
        for idx in states:
            pop = self.get_populations(idx).cpu().numpy()
            plt.plot(times_np, pop, label=f"{label_prefix} |{idx}><{idx}|")
        plt.xlabel("Time")
        plt.ylabel("Population")
        plt.legend()
        plt.grid(True)


class Simulator:
    def __init__(self, system, drive, bath):
        self.system = system
        self.drive = drive
        self.bath = bath

    def commutator(self, A, B):
        return A @ B - B @ A

    def liouvillian(self, t, rho, method, **kwargs):
        """
        Computes d(rho)/dt for the selected mathematical framework.

        Methods:
            'lindblad'          — legacy: scalar drive * dipole, ThermalBath dissipator
            'redfield'          — full non-secular Redfield (SpectralDensityBath)
            'secular_lindblad'  — secular Lindblad / GKSL (SpectralDensityBath)
            'nonadiabatic_cgme' — coarse-grained ME (SpectralDensityBath)

        For the three new methods, rho is in the energy eigenbasis and
        self._U / self._H0_eig must be set by _run_eigenbasis() before calling.
        """
        if method == 'lindblad':
            # Legacy path: scalar drive * dipole matrix, ThermalBath
            H0 = self.system.get_H0()
            V_t = self.system.dipole * self.drive(t)
            H_t = H0 + V_t

            unitary = -1j * self.commutator(H_t, rho)
            dissipator = self.bath.get_lindblad_dissipator(self.system, rho)

            return unitary + dissipator

        elif method in ('redfield', 'secular_lindblad', 'nonadiabatic_cgme'):
            # New path: matrix-valued drive, eigenbasis propagation
            V_comp = self.drive(t)
            V_eig = self._U.mH @ V_comp @ self._U
            H_eig = self._H0_eig + V_eig

            unitary = -1j * self.commutator(H_eig, rho)

            if method == 'redfield':
                dissipator = self.bath.get_redfield_dissipator(rho)
            elif method == 'secular_lindblad':
                dissipator = self.bath.get_secular_lindblad_dissipator(rho)
            else:
                delta_tau = kwargs.get('delta_tau', 1.0)
                dissipator = self.bath.get_cgme_dissipator(rho, delta_tau)

            return unitary + dissipator

        else:
            raise ValueError(f"Unknown method: {method}")

    def run(self, times_np, method='lindblad', initial_state=None, **kwargs):
        """
        Propagate the master equation.

        Args:
            times_np: 1D array of output times
            method: 'lindblad', 'redfield', 'secular_lindblad', or 'nonadiabatic_cgme'
            initial_state: optional initial density matrix (default: |0><0|)
            **kwargs: passed to liouvillian (e.g. delta_tau for CGME)

        Returns:
            Result object with times, states, and (for CGME) nonadiabatic_populations
        """
        times = torch.tensor(times_np, dtype=torch.float64, device=self.system.device)
        dim = self.system.dim

        if initial_state is not None:
            rho_0 = initial_state.to(torch.complex128).to(self.system.device)
        else:
            rho_0 = torch.zeros((dim, dim), dtype=torch.complex128, device=self.system.device)
            rho_0[0, 0] = 1.0

        if method in ('redfield', 'secular_lindblad', 'nonadiabatic_cgme'):
            return self._run_eigenbasis(times, rho_0, method, **kwargs)
        else:
            return self._run_legacy(times, rho_0, method, **kwargs)

    def _run_eigenbasis(self, times, rho_0, method, **kwargs):
        """Propagation in the energy eigenbasis for SpectralDensityBath methods."""
        from nonadiabatic import NonadiabaticDecomposition

        assert hasattr(self.bath, '_precomputed') and self.bath._precomputed, \
            "SpectralDensityBath.precompute() must be called before running"

        # Use the bath's eigenbasis for consistency with the Redfield tensor
        energies = self.bath.energies
        self._U = self.bath.U
        self._H0_eig = torch.diag(energies.to(torch.complex128))

        # Transform initial state to eigenbasis
        rho = self._U.mH @ rho_0 @ self._U

        # Nonadiabatic decomposition (tracks b_k alongside rho propagation)
        decomp = None
        if method == 'nonadiabatic_cgme':
            decomp = NonadiabaticDecomposition(energies, self._U)

        # RK4 stability: max_dt from spectral range of H0
        spectral_range = (energies[-1] - energies[0]).item()
        max_dt = 1.5 / spectral_range if spectral_range > 0 else float('inf')

        # Storage (computational basis for states, eigenbasis for nad populations)
        states = [rho_0.clone()]
        nad_pops = [decomp.get_nonadiabatic_populations().clone()] if decomp else None

        for i in range(1, len(times)):
            t_start = times[i - 1].item()
            interval = times[i].item() - t_start

            n_substeps = max(1, math.ceil(abs(interval) / max_dt))
            dt_sub = interval / n_substeps

            for j in range(n_substeps):
                t = t_start + j * dt_sub

                k1 = dt_sub * self.liouvillian(t, rho, method, **kwargs)
                k2 = dt_sub * self.liouvillian(t + 0.5 * dt_sub, rho + 0.5 * k1, method, **kwargs)
                k3 = dt_sub * self.liouvillian(t + 0.5 * dt_sub, rho + 0.5 * k2, method, **kwargs)
                k4 = dt_sub * self.liouvillian(t + dt_sub, rho + k3, method, **kwargs)

                rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

                # Update nonadiabatic decomposition (once per substep, not per RK4 stage)
                if decomp is not None:
                    dV_dt = self.drive.derivative(t)
                    decomp.step(t, dt_sub, dV_dt)

            # Trace check — no renormalization; drift indicates a dissipator bug
            trace_val = torch.trace(rho).real.item()
            if abs(trace_val - 1.0) > 1e-6:
                raise RuntimeError(
                    f"Trace drift at t={times[i].item():.6g}: Tr(rho) = {trace_val:.10g}")

            # Transform back to computational basis for storage
            rho_comp = self._U @ rho @ self._U.mH
            states.append(rho_comp.clone())

            if decomp is not None and nad_pops is not None:
                nad_pops.append(decomp.get_nonadiabatic_populations().clone())

        result = Result(times, torch.stack(states))
        if nad_pops is not None:
            result.nonadiabatic_populations = torch.stack(nad_pops)  # type: ignore[assignment]
        return result

    def _run_legacy(self, times, rho_0, method, **kwargs):
        """Original propagation for the 'lindblad' method with ThermalBath."""
        states = [rho_0.clone()]
        rho = rho_0.clone()

        # Estimate max stable step from system frequencies
        H0 = self.system.get_H0()
        eigs = torch.linalg.eigvalsh(H0.real.to(torch.float64))
        spectral_range = (eigs[-1] - eigs[0]).item()
        max_dt = 1.5 / spectral_range if spectral_range > 0 else float('inf')

        for i in range(1, len(times)):
            t_start = times[i - 1]
            interval = (times[i] - t_start).item()

            n_substeps = max(1, math.ceil(abs(interval) / max_dt))
            dt_sub = interval / n_substeps

            for j in range(n_substeps):
                t = t_start + j * dt_sub

                k1 = dt_sub * self.liouvillian(t, rho, method, **kwargs)
                k2 = dt_sub * self.liouvillian(t + 0.5 * dt_sub, rho + 0.5 * k1, method, **kwargs)
                k3 = dt_sub * self.liouvillian(t + 0.5 * dt_sub, rho + 0.5 * k2, method, **kwargs)
                k4 = dt_sub * self.liouvillian(t + dt_sub, rho + k3, method, **kwargs)

                rho = rho + (k1 + 2 * k2 + 2 * k3 + k4) / 6.0

            # Ensure trace preservation (correcting minor numerical drift)
            rho = rho / torch.trace(rho)
            states.append(rho.clone())

        return Result(times, torch.stack(states))
