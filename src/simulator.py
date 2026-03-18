import math
import torch
import matplotlib.pyplot as plt

class Result:
    def __init__(self, times, states):
        self.times = times
        self.states = states # Tensor of shape (len(times), dim, dim)

    def get_populations(self, state_index):
        # Extract diagonal elements (probabilities)
        return torch.real(self.states[:, state_index, state_index])

    def plot_populations(self, states=[0, 1, 2], label_prefix=""):
        times_np = self.times.cpu().numpy()
        for idx in states:
            pop = self.get_populations(idx).cpu().numpy()
            plt.plot(times_np, pop, label=f"{label_prefix} |{idx}><{idx}|")
        plt.xlabel("Time (ns)")
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
        """
        H0 = self.system.get_H0()
        
        if method == 'lindblad':
            # Standard Dirac Basis: H(t) = H0 + V(t)
            V_t = self.system.dipole * self.drive(t)
            H_t = H0 + V_t
            
            # Unitary evolution: -i [H, rho]
            unitary = -1j * self.commutator(H_t, rho)
            
            # Markovian Secular Dissipation
            dissipator = self.bath.get_lindblad_dissipator(self.system, rho)
            
            return unitary + dissipator
            
        elif method == 'nonadiabatic_cgme':
            # 1. Compute V(t) and its exact time derivative dV/dt
            V_t = self.system.dipole * self.drive(t)
            dV_dt = self.system.dipole * self.drive.derivative(t)
            
            H_t = H0 + V_t
            
            # 2. Structural Hook for Nonadiabatic Subspace Projection P_nad(t)
            # In a full implementation, dV_dt is used here to compute the 
            # adiabatic polarization a_k(t) and subtract it from the memory kernel.
            # unitary_nad = -i [H_eff(t), rho_nad]
            unitary = -1j * self.commutator(H_t, rho) 
            
            # 3. Structural Hook for Coarse-Grained Master Equation (CGME)
            # P_CG operates over time window \Delta\tau to ensure positivity
            # while retaining non-secular coherence transfers.
            dt_cg = kwargs.get('coarse_grain_timescale', 0.5)
            
            # Fallback to standard dissipator for the structural skeleton
            dissipator = self.bath.get_lindblad_dissipator(self.system, rho)
            
            return unitary + dissipator
            
        else:
            raise ValueError(f"Unknown method: {method}")

    def run(self, times_np, method='lindblad', **kwargs):
        times = torch.tensor(times_np, dtype=torch.float64, device=self.system.device)

        # Initial state: Ground state |0><0|
        rho_0 = torch.zeros((self.system.dim, self.system.dim), dtype=torch.complex128, device=self.system.device)
        rho_0[0, 0] = 1.0

        states = [rho_0]
        rho = rho_0.clone()

        # Estimate max stable step from system frequencies
        # RK4 stability on imaginary axis: |lambda * dt| < 2*sqrt(2)
        H0 = self.system.get_H0()
        eigs = torch.linalg.eigvalsh(H0.real.to(torch.float64))
        spectral_range = (eigs[-1] - eigs[0]).item()
        max_dt = 1.5 / spectral_range if spectral_range > 0 else float('inf')

        # 4th-Order Runge-Kutta with automatic sub-stepping
        for i in range(1, len(times)):
            t_start = times[i-1]
            interval = (times[i] - t_start).item()

            n_substeps = max(1, math.ceil(abs(interval) / max_dt))
            dt_sub = interval / n_substeps

            for j in range(n_substeps):
                t = t_start + j * dt_sub

                k1 = dt_sub * self.liouvillian(t, rho, method, **kwargs)
                k2 = dt_sub * self.liouvillian(t + 0.5*dt_sub, rho + 0.5*k1, method, **kwargs)
                k3 = dt_sub * self.liouvillian(t + 0.5*dt_sub, rho + 0.5*k2, method, **kwargs)
                k4 = dt_sub * self.liouvillian(t + dt_sub, rho + k3, method, **kwargs)

                rho = rho + (k1 + 2*k2 + 2*k3 + k4) / 6.0

            # Ensure trace preservation (correcting minor numerical drift)
            rho = rho / torch.trace(rho)
            states.append(rho.clone())

        return Result(times, torch.stack(states))