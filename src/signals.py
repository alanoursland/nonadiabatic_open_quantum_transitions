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


class SpinLockDrive(Signal):
    """
    NMR spin-lock drive in the rotating frame.

    In the rotating frame at the carrier frequency, the perturbation is:
        V(t) = omega_1 * f(t) * Ix  +  Delta_omega * Iz
    where f(t) is the sin^2 ramp-up / plateau / cos^2 ramp-down envelope,
    omega_1 is the spin-lock field strength, and Delta_omega is the offset.

    Time unit for all calls: microseconds (us).
    """

    def __init__(self, spin_lock_power_hz, Ix_operator, Iz_operator=None,
                 offset_hz=0.0, ramp_time_us=0.0, plateau_time_us=50.0):
        """
        Args:
            spin_lock_power_hz: omega_1 / (2*pi) in Hz
            Ix_operator: Ix spin operator matrix for the target spin
            Iz_operator: Iz spin operator matrix (needed for off-resonance)
            offset_hz: off-resonance offset Delta_omega / (2*pi) in Hz
            ramp_time_us: ramp duration in microseconds
            plateau_time_us: plateau duration in microseconds
        """
        # Convert Hz to rad/us: 1 Hz = 2*pi rad/s = 2*pi * 1e-6 rad/us
        self.omega1 = 2 * torch.pi * spin_lock_power_hz * 1e-6
        self.Delta = 2 * torch.pi * offset_hz * 1e-6
        self.Ix = Ix_operator
        self.Iz = Iz_operator
        self.ramp = ramp_time_us
        self.plateau = plateau_time_us

    def _envelope(self, t):
        """Smooth envelope: sin^2 ramp-up, flat plateau, cos^2 ramp-down."""
        t_tensor = torch.as_tensor(t, dtype=torch.float64)
        env = torch.zeros_like(t_tensor)

        mask_up = (t_tensor >= 0) & (t_tensor < self.ramp)
        if self.ramp > 0:
            env[mask_up] = torch.sin((torch.pi / 2) * (t_tensor[mask_up] / self.ramp)) ** 2

        mask_plat = (t_tensor >= self.ramp) & (t_tensor <= self.ramp + self.plateau)
        env[mask_plat] = 1.0

        mask_down = t_tensor > self.ramp + self.plateau
        if self.ramp > 0:
            decay_t = t_tensor[mask_down] - (self.ramp + self.plateau)
            env[mask_down] = torch.cos((torch.pi / 2) * torch.clamp(decay_t / self.ramp, max=1.0)) ** 2

        return env

    def _envelope_derivative(self, t):
        """Time derivative of the envelope (rad/us units)."""
        t_tensor = torch.as_tensor(t, dtype=torch.float64)
        denv = torch.zeros_like(t_tensor)

        if self.ramp > 0:
            mask_up = (t_tensor >= 0) & (t_tensor < self.ramp)
            denv[mask_up] = (torch.pi / self.ramp) * \
                torch.sin((torch.pi / 2) * (t_tensor[mask_up] / self.ramp)) * \
                torch.cos((torch.pi / 2) * (t_tensor[mask_up] / self.ramp))

            mask_down = (t_tensor > self.ramp + self.plateau) & \
                        (t_tensor < 2 * self.ramp + self.plateau)
            decay_t = t_tensor[mask_down] - (self.ramp + self.plateau)
            denv[mask_down] = -(torch.pi / self.ramp) * \
                torch.cos((torch.pi / 2) * (decay_t / self.ramp)) * \
                torch.sin((torch.pi / 2) * (decay_t / self.ramp))

        return denv

    def __call__(self, t):
        """V(t) = omega1 * f(t) * Ix + Delta * Iz.  Time t in us."""
        f = self._envelope(t)
        V = self.omega1 * f * self.Ix
        if self.Iz is not None and abs(self.Delta) > 0:
            V = V + self.Delta * self.Iz
        return V

    def derivative(self, t):
        """dV/dt = omega1 * df/dt * Ix.  Time t in us."""
        df = self._envelope_derivative(t)
        return self.omega1 * df * self.Ix

    def effective_field_angle(self):
        """Tilt angle theta of the effective field from z-axis (radians).
        theta = pi/2 on resonance (field along x), 0 far off resonance."""
        if abs(self.Delta) < 1e-30:
            return torch.tensor(torch.pi / 2)
        return torch.atan2(torch.tensor(self.omega1), torch.tensor(self.Delta))

    def effective_field_magnitude(self):
        """Effective field magnitude omega_eff in rad/us."""
        return (self.omega1 ** 2 + self.Delta ** 2) ** 0.5