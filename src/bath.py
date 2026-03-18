import torch

class ThermalBath:
    def __init__(self, temperature, t1_times, t2_times, device='cpu'):
        self.device = device
        self.temperature = temperature
        self.t1 = t1_times
        self.t2 = t2_times

    def get_lindblad_dissipator(self, system, rho):
        """
        Standard secular Lindblad dissipator D[rho].
        Constructs jump operators from phenomenological T1/T2 times.
        """
        D_rho = torch.zeros_like(rho, dtype=torch.complex128, device=self.device)
        
        # Population relaxation (T1)
        for (i, j), t1 in self.t1.items():
            rate = 1.0 / t1
            # Jump operator L = |j><i| (decay from i to j)
            L = torch.zeros((system.dim, system.dim), dtype=torch.complex128, device=self.device)
            L[j, i] = 1.0 
            
            L_dag = L.mH
            # D = rate * (L rho L^dagger - 0.5 * {L^dagger L, rho})
            term1 = L @ rho @ L_dag
            term2 = 0.5 * (L_dag @ L @ rho + rho @ L_dag @ L)
            D_rho += rate * (term1 - term2)
            
        # Pure dephasing (T2* contribution calculated from T1 and T2)
        for (i, j), t2 in self.t2.items():
            rate_t1_i = sum([1/t for (k, l), t in self.t1.items() if k == i])
            rate_t1_j = sum([1/t for (k, l), t in self.t1.items() if k == j])
            
            # 1/T2 = 1/(2*T1_i) + 1/(2*T1_j) + 1/Tphi
            gamma_phi = (1.0 / t2) - 0.5 * rate_t1_i - 0.5 * rate_t1_j
            
            if gamma_phi > 0:
                L_z = torch.zeros((system.dim, system.dim), dtype=torch.complex128, device=self.device)
                L_z[i, i] = 1.0
                L_z[j, j] = -1.0
                
                term1 = L_z @ rho @ L_z.mH
                term2 = 0.5 * (L_z.mH @ L_z @ rho + rho @ L_z.mH @ L_z)
                D_rho += (gamma_phi / 2.0) * (term1 - term2)
                
        return D_rho