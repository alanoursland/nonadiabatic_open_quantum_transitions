import torch
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


class R1rhoExtractor:
    """
    Extracts R1rho from a simulation Result by fitting the decay of
    <O(t)> = Tr[rho(t) @ O] to M0 * exp(-R1rho * t) over a specified
    time window (typically the spin-lock plateau).
    """

    def __init__(self, result, observable, method_name=''):
        """
        Args:
            result: Result object from Simulator.run()
            observable: Hermitian operator O (torch tensor, same dim as rho)
            method_name: label for plotting / display
        """
        self.result = result
        self.observable = observable
        self.method_name = method_name

    def get_expectation_values(self):
        """Compute <O(t)> = Tr[rho(t) @ O] for all stored time steps.

        Returns:
            torch.Tensor of shape (n_times,), real-valued
        """
        # states: (n_times, dim, dim), observable: (dim, dim)
        # Tr[rho @ O] = sum_{ij} rho_{ij} O_{ji} = einsum('tij,ji->t')
        vals = torch.einsum('tij,ji->t', self.result.states, self.observable)
        return vals.real

    def fit(self, t_start=None, t_end=None):
        """Fit <O(t)> to M0 * exp(-R1rho * t) in [t_start, t_end].

        Args:
            t_start: start of fitting window (default: first time point)
            t_end: end of fitting window (default: last time point)

        Returns:
            dict with keys 'R1rho', 'M0', 'residual_std'
        """
        times_np = self.result.times.cpu().numpy()
        expectation = self.get_expectation_values().cpu().numpy()

        # Select fitting window
        mask = np.ones(len(times_np), dtype=bool)
        if t_start is not None:
            mask &= times_np >= t_start
        if t_end is not None:
            mask &= times_np <= t_end

        t_fit = times_np[mask]
        y_fit = expectation[mask]

        if len(t_fit) < 3:
            raise ValueError(
                f"Need at least 3 points in fitting window, got {len(t_fit)}")

        # Shift time origin to t_fit[0] for numerical stability
        t0 = t_fit[0]
        t_shifted = t_fit - t0

        def exp_decay(t, M0, R):
            return M0 * np.exp(-R * t)

        # Initial guesses
        M0_guess = abs(y_fit[0]) if abs(y_fit[0]) > 1e-30 else 1.0
        R_guess = 0.1

        try:
            popt, _ = curve_fit(
                exp_decay, t_shifted, y_fit,
                p0=[M0_guess, R_guess],
                bounds=([0.0, 0.0], [np.inf, np.inf]),
                maxfev=10000,
            )
        except RuntimeError as e:
            raise RuntimeError(f"Exponential fit failed: {e}")

        M0_fit, R1rho_fit = popt
        residuals = y_fit - exp_decay(t_shifted, *popt)

        return {
            'R1rho': float(R1rho_fit),
            'M0': float(M0_fit),
            'residual_std': float(np.std(residuals)),
        }


class R1rhoDispersion:
    """Container for R1rho vs spin-lock power dispersion data.

    Stores (omega1, R1rho) pairs for one or more methods and
    provides plotting.
    """

    def __init__(self):
        # method_name -> {'omega1': list, 'R1rho': list}
        self._data = {}

    def add_point(self, method_name, omega1_hz, R1rho):
        """Add a single (omega1, R1rho) measurement.

        Args:
            method_name: str label for the method
            omega1_hz: spin-lock power in Hz
            R1rho: extracted relaxation rate
        """
        if method_name not in self._data:
            self._data[method_name] = {'omega1': [], 'R1rho': []}
        self._data[method_name]['omega1'].append(float(omega1_hz))
        self._data[method_name]['R1rho'].append(float(R1rho))

    def get_methods(self):
        """Return list of method names with stored data."""
        return list(self._data.keys())

    def get_dispersion(self, method_name):
        """Return (omega1_array, R1rho_array) for the given method.

        Returns:
            tuple of numpy arrays
        """
        d = self._data[method_name]
        return np.array(d['omega1']), np.array(d['R1rho'])

    def plot(self, ax=None, experimental=None):
        """Plot R1rho dispersion curves for all methods.

        Args:
            ax: matplotlib Axes (created if None)
            experimental: optional dict {'omega1': array, 'R1rho': array}
                for overlay of experimental data
        """
        if ax is None:
            _, ax = plt.subplots()

        for method_name, d in self._data.items():
            omega1 = np.array(d['omega1'])
            r1rho = np.array(d['R1rho'])
            sort_idx = np.argsort(omega1)
            ax.plot(omega1[sort_idx], r1rho[sort_idx], 'o-', label=method_name)

        if experimental is not None:
            ax.plot(experimental['omega1'], experimental['R1rho'],
                    's', color='black', label='Experimental', markersize=8)

        ax.set_xlabel(r'$\omega_1 / 2\pi$ (Hz)')
        ax.set_ylabel(r'$R_{1\rho}$ (s$^{-1}$)')
        ax.legend()
        ax.grid(True)
        return ax
