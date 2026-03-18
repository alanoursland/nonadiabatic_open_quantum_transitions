import torch

class Signal:
    def __call__(self, t):
        raise NotImplementedError
        
    def derivative(self, t):
        raise NotImplementedError

class PlateauPulse(Signal):
    def __init__(self, amplitude, carrier_freq, ramp_time, plateau_time):
        """
        amplitude: Peak drive strength in GHz
        carrier_freq: Drive frequency in GHz
        ramp_time: Ramping duration in ns
        plateau_time: Flat constant drive duration in ns
        """
        # Convert GHz to angular frequency
        self.A0 = amplitude * 2 * torch.pi
        self.omega = carrier_freq * 2 * torch.pi
        self.ramp = ramp_time
        self.plateau = plateau_time
        
    def _envelope(self, t):
        # Smooth step up, flat plateau, smooth step down
        t_tensor = torch.as_tensor(t, dtype=torch.float64)
        env = torch.zeros_like(t_tensor)
        
        # Ramp up
        mask_up = (t_tensor >= 0) & (t_tensor < self.ramp)
        env[mask_up] = torch.sin((torch.pi / 2) * (t_tensor[mask_up] / self.ramp))**2
        
        # Plateau
        mask_plat = (t_tensor >= self.ramp) & (t_tensor <= self.ramp + self.plateau)
        env[mask_plat] = 1.0
        
        # Ramp down
        mask_down = t_tensor > self.ramp + self.plateau
        decay_t = t_tensor[mask_down] - (self.ramp + self.plateau)
        env[mask_down] = torch.cos((torch.pi / 2) * torch.clamp(decay_t / self.ramp, max=1.0))**2
        
        return env

    def _envelope_derivative(self, t):
        t_tensor = torch.as_tensor(t, dtype=torch.float64)
        denv = torch.zeros_like(t_tensor)
        
        # Ramp up derivative
        mask_up = t_tensor < self.ramp
        denv[mask_up] = (torch.pi / self.ramp) * torch.sin((torch.pi / 2) * (t_tensor[mask_up] / self.ramp)) * \
                        torch.cos((torch.pi / 2) * (t_tensor[mask_up] / self.ramp))
        
        # Plateau derivative is 0
        
        # Ramp down derivative
        mask_down = (t_tensor > self.ramp + self.plateau) & (t_tensor < 2*self.ramp + self.plateau)
        decay_t = t_tensor[mask_down] - (self.ramp + self.plateau)
        denv[mask_down] = -(torch.pi / self.ramp) * torch.cos((torch.pi / 2) * (decay_t / self.ramp)) * \
                          torch.sin((torch.pi / 2) * (decay_t / self.ramp))
        
        return denv

    def __call__(self, t):
        t_tensor = torch.as_tensor(t, dtype=torch.float64)
        env = self._envelope(t_tensor)
        oscillation = torch.cos(self.omega * t_tensor)
        return self.A0 * env * oscillation

    def derivative(self, t):
        t_tensor = torch.as_tensor(t, dtype=torch.float64)
        env = self._envelope(t_tensor)
        denv = self._envelope_derivative(t_tensor)
        oscillation = torch.cos(self.omega * t_tensor)
        d_oscillation = -self.omega * torch.sin(self.omega * t_tensor)
        
        # Product rule: d(env * osc)/dt = denv * osc + env * d_osc
        return self.A0 * (denv * oscillation + env * d_oscillation)