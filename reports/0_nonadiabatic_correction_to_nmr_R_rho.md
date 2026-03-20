# Phase 0 Feasibility Report: Nonadiabatic Correction to NMR $R_{1\rho}$

## Abstract

We estimated the magnitude of the nonadiabatic correction to rotating-frame relaxation rates ($R_{1\rho}$) for $^{15}$N backbone amide nuclei in proteins, using analytical perturbation theory and numerical simulation. The correction arises from the Landau-Lifshitz decomposition of Dirac's transition coefficients into adiabatic (polarization) and nonadiabatic (genuine transition) components, as developed by Mandal and Hunt (2012–2023). The key question was whether this correction produces a measurable difference in $R_{1\rho}$ that could be tested against published NMR data.

The answer is no. The maximum relative correction is $6.3 \times 10^{-11}$, eleven orders of magnitude below the experimental precision threshold of $\sim 10^{-2}$. NMR $R_{1\rho}$ cannot distinguish the nonadiabatic framework from standard Redfield/Lindblad theory. This negative result is a consequence of the extremely small ratio of the spin-lock field strength ($\sim 10^2$–$10^3$ Hz) to the Larmor frequency ($\sim 5 \times 10^7$ Hz), which makes the adiabatic polarization component negligible.

---

## Background

### The Nonadiabatic Framework

The Landau-Lifshitz separation decomposes the time-dependent perturbation theory coefficient $c_k(t)$ for an excited state into an adiabatic term $a_k(t)$ (reversible polarization response to the instantaneous field) and a nonadiabatic term $b_k(t)$ (genuine transitions driven by the field's time variation). Mandal and Hunt demonstrated through energy balance, statistical consistency, and gauge invariance arguments that $|b_k(t)|^2$, not $|c_k(t)|^2$, is the physically correct transition probability.

When a quantum system is simultaneously driven by an external field and coupled to a thermal bath, the choice between $|c_k|^2$ and $|b_k|^2$ has physical consequences: using Dirac's coefficients causes the bath to act on the adiabatic polarization as though it were a real excitation, producing fictitious heat flow and incorrect thermal equilibrium. The nonadiabatic framework avoids these pathologies by ensuring the bath couples only to genuine transitions.

### The NMR $R_{1\rho}$ Experiment

In a rotating-frame relaxation experiment, a nuclear spin is subjected to a continuous RF spin-lock field of amplitude $\omega_1$ while simultaneously relaxing through coupling to the molecular tumbling (the thermal bath). This is a direct realization of the driven-dissipative scenario: the spin-lock field is the perturbation $V$, it is constant during the measurement plateau ($\partial V/\partial t = 0$), and the bath acts simultaneously.

The standard analytical result for the on-resonance $R_{1\rho}$ is:

$$R_{1\rho} = R_1 \cos^2\theta + R_2 \sin^2\theta$$

where $\theta = \arctan(\omega_1/\Delta\omega)$ is the effective field tilt angle, and $R_1$, $R_2$ are computed from the Solomon-Bloembergen equations using the spectral density function $J(\omega)$.

The nonadiabatic framework predicts a correction to this formula: the adiabatic polarization induced by the spin-lock field should be excluded from the bath-coupled dynamics. The magnitude of this correction depends on the perturbation ratio $\omega_1 / \omega_0$, where $\omega_0$ is the Larmor frequency.

### Purpose of This Estimate

Before investing in a detailed numerical comparison against published experimental data (the full Phase 3 of the experiment plan), we needed to determine whether the predicted correction is large enough to be measurable. This Phase 0 estimate provides that determination.

---

## Methods

### System Parameters

We used standard parameters for a $^{15}$N backbone amide in a small globular protein at 11.7 T (500 MHz $^1$H frequency):

| Parameter | Value | Source |
|:---|:---|:---|
| Static field $B_0$ | 11.7 T | Standard NMR |
| Correlation time $\tau_c$ | 5.0 ns | Typical small protein |
| Order parameter $S^2$ | 0.85 | Lipari-Szabo model-free |
| Internal correlation time $\tau_e$ | 50 ps | Lipari-Szabo |
| $^{15}$N gyromagnetic ratio $\gamma_N$ | $-2.7116 \times 10^7$ rad/s/T | CODATA |
| $^1$H gyromagnetic ratio $\gamma_H$ | $2.6752 \times 10^8$ rad/s/T | CODATA |
| N-H bond length $r_{\text{NH}}$ | 1.02 Å | Standard |
| $^{15}$N CSA $\Delta\sigma$ | $-170$ ppm | Standard backbone |
| Spin-lock powers $\omega_1/2\pi$ | 25–1000 Hz | Korzhnev et al. (2005) |

### Derived Quantities

From these parameters, the following were computed:

| Quantity | Value |
|:---|:---|
| $^{15}$N Larmor frequency $\omega_N$ | $3.173 \times 10^8$ rad/s (50.49 MHz) |
| $^1$H Larmor frequency $\omega_H$ | $3.130 \times 10^9$ rad/s (498.15 MHz) |
| Dipolar coupling constant $d$ | 72,087 rad/s |
| CSA coupling constant $c$ | 31,139 rad/s |

### Relaxation Rates (Standard Redfield)

The Solomon-Bloembergen equations with the Lipari-Szabo spectral density give:

| Rate | Value |
|:---|:---|
| $R_1$ | 2.449 s$^{-1}$ ($T_1 = 0.408$ s) |
| $R_2$ | 6.789 s$^{-1}$ ($T_2 = 0.147$ s) |
| NOE | 0.790 |

These values are consistent with published experimental data for ubiquitin-class proteins at 500 MHz.

### Correction Mechanism 1: Perturbative Adiabatic Mixing

In the Landau-Lifshitz decomposition, the adiabatic coefficient for a two-level system under a transverse perturbation of strength $\omega_1$ is:

$$|a_1|^2 = \left(\frac{\omega_1}{2\omega_N}\right)^2$$

This represents the fraction of the density matrix that is adiabatically polarized by the spin-lock field. In the nonadiabatic framework, this fraction is excluded from bath-induced relaxation. The resulting correction to $R_{1\rho}$ is bounded by:

$$\delta R_{1\rho} \leq |a_1|^2 \times |R_2 - R_1|$$

### Correction Mechanism 2: CGME vs. Secular Lindblad

The coarse-grained master equation (CGME) damps each non-secular Redfield tensor element by a factor:

$$\text{sinc}\left(\frac{(\omega_{ab} - \omega_{cd})\Delta\tau}{2}\right)$$

For NMR Bohr frequencies ($\sim 10^8$–$10^9$ rad/s) and a coarse-graining timescale $\Delta\tau = 1$ μs, the argument of the sinc is large ($\sim 100$–$1000$), and the damping factor is negligible. The CGME is effectively identical to the secular Lindblad equation for these parameters.

### Numerical Verification

In addition to the analytical estimates, we ran numerical simulations using the full Redfield tensor machinery (Work Packages 1–6 of the codebase) on a toy two-level system ($\omega = 5$ rad/s) to verify that the secular Lindblad and CGME dissipators produce identical population dynamics when the non-secular sinc factors are small, and differ in the coherence sector when they are not.

---

## Results

### Correction 1: Perturbative Adiabatic Mixing

| $\omega_1/2\pi$ (Hz) | $\|a_1\|^2$ | $\delta R_{1\rho}$ (s$^{-1}$) | $\delta R_{1\rho} / R_{1\rho}$ |
|:---|:---|:---|:---|
| 25 | $6.13 \times 10^{-14}$ | $2.66 \times 10^{-13}$ | $3.92 \times 10^{-14}$ |
| 50 | $2.45 \times 10^{-13}$ | $1.06 \times 10^{-12}$ | $1.57 \times 10^{-13}$ |
| 100 | $9.81 \times 10^{-13}$ | $4.25 \times 10^{-12}$ | $6.27 \times 10^{-13}$ |
| 250 | $6.13 \times 10^{-12}$ | $2.66 \times 10^{-11}$ | $3.92 \times 10^{-12}$ |
| 500 | $2.45 \times 10^{-11}$ | $1.06 \times 10^{-10}$ | $1.57 \times 10^{-11}$ |
| 1000 | $9.81 \times 10^{-11}$ | $4.25 \times 10^{-10}$ | $6.27 \times 10^{-11}$ |

The maximum relative correction, at the strongest accessible spin-lock power of 1000 Hz, is $6.27 \times 10^{-11}$.

### Correction 2: CGME Sinc Retention Factors

| Frequency pair | $\|\Delta\omega\|$ (rad/s) | $\|\text{sinc}\|$ |
|:---|:---|:---|
| $\omega_N$ | $3.173 \times 10^8$ | $6.30 \times 10^{-3}$ |
| $\omega_H$ | $3.130 \times 10^9$ | $2.94 \times 10^{-4}$ |
| $\omega_H - \omega_N$ | $2.813 \times 10^9$ | $6.24 \times 10^{-4}$ |
| $\omega_H + \omega_N$ | $3.447 \times 10^9$ | $5.21 \times 10^{-4}$ |
| $2\omega_N$ | $6.345 \times 10^8$ | $1.38 \times 10^{-4}$ |
| $2\omega_H$ | $6.260 \times 10^9$ | $2.61 \times 10^{-4}$ |

All sinc retention factors are below 1%, confirming that the CGME and secular Lindblad dissipators are functionally identical for NMR parameters.

### Numerical Simulator Comparison

For the toy two-level system ($\omega = 5$ rad/s, $\Delta\tau = 0.5$), the secular and CGME methods produced identical population dynamics (agreement to $10^{-8}$) but differed in coherence decay rates, confirming that the CGME code correctly retains non-secular terms when the sinc factors are $O(1)$. Both methods decayed from the excited state and reached the same thermal equilibrium.

---

## Discussion

### Why the Correction Is So Small

The nonadiabatic correction to $R_{1\rho}$ is proportional to $(\omega_1/\omega_N)^2$. For NMR:

$$\frac{\omega_1}{\omega_N} = \frac{2\pi \times 1000 \text{ Hz}}{2\pi \times 50.49 \times 10^6 \text{ Hz}} \approx 2 \times 10^{-5}$$

This ratio squared is $\sim 4 \times 10^{-10}$. The spin-lock field is not a small perturbation compared to the internal dynamics — it is a *vanishingly* small perturbation. The entire edifice of the nonadiabatic correction (adiabatic polarization, gauge invariance, thermodynamic consistency) is physically correct, but the effect it corrects is negligible for NMR because the perturbation is negligible.

This is a property of NMR specifically, not of the nonadiabatic framework in general. The Mandal-Hunt papers demonstrate 35%+ differences between $|c_k|^2$ and $|b_k|^2$ for overlapping electromagnetic pulses with perturbation ratios of order $10^{-1}$. The framework produces large corrections when $V/\Delta E \sim 0.1$–$1$. NMR spin-locks have $V/\Delta E \sim 10^{-5}$.

### Why NMR Was a Reasonable Candidate

NMR was selected because it satisfies the experimental requirements identified in the experiment plan: (1) simultaneously driven and dissipating, (2) independently known Hamiltonian parameters, (3) the perturbative regime applies, and (4) precision measurements exist. It fails on an additional criterion that was not sufficiently weighted: the perturbation must be large enough relative to the energy gap for the correction to be detectable.

The spin-lock in NMR is a plateau pulse — exactly the scenario where the nonadiabatic framework predicts constant populations versus the Dirac framework's oscillating populations. But the oscillation amplitude is $(\omega_1/\omega_N)^2 \sim 10^{-10}$, so the difference between "constant" and "oscillating with amplitude $10^{-10}$" is unmeasurable.

### Implications for the Experimental Search

The perturbation ratio $V/\Delta E$ is the key parameter for identifying systems where the nonadiabatic correction is detectable. The remaining experimental candidates have the following approximate perturbation ratios:

| System | Perturbation | Energy gap | $V/\Delta E$ |
|:---|:---|:---|:---|
| NMR $R_{1\rho}$ (this work) | Spin-lock ($\sim$kHz) | Larmor ($\sim$50 MHz) | $\sim 10^{-5}$ |
| FTMW spectroscopy | Microwave pulse ($\sim$MHz) | Rotational splitting ($\sim$GHz) | $\sim 10^{-3}$–$10^{-1}$ |
| HITRAN line mixing | Collisional coupling ($\sim$0.1 cm$^{-1}$) | Rotational spacing ($\sim$1 cm$^{-1}$) | $\sim 10^{-1}$ |
| 2D electronic spectroscopy | Laser pulse ($\sim$100 cm$^{-1}$) | Exciton splitting ($\sim$200 cm$^{-1}$) | $\sim 10^{-1}$–$1$ |

The correction scales as $(V/\Delta E)^2$, so moving from NMR ($10^{-5}$) to HITRAN ($10^{-1}$) increases the correction by a factor of $10^8$, from $10^{-10}$ to $10^{-2}$ — potentially measurable.

---

## Conclusion

The nonadiabatic correction to NMR $R_{1\rho}$ is $6.3 \times 10^{-11}$ relative, far below the experimental precision of $\sim 10^{-2}$. This definitively rules out NMR rotating-frame relaxation as a test bed for the nonadiabatic framework.

The result does not invalidate the nonadiabatic framework itself — it identifies a regime where the correction is negligible due to the extremely small perturbation ratio. Systems with larger perturbation ratios ($V/\Delta E \gtrsim 10^{-2}$) remain viable candidates, with HITRAN pressure-broadened line mixing and FTMW coherent transients as the most promising alternatives.

---

## Code and Reproducibility

All calculations were performed using the `nonadiabatic_estimate.py` module in the `nonadiabatic_open_quantum_transitions` codebase. The full result is reproduced by running:

```
python nonadiabatic_estimate.py
```

The codebase includes 200+ automated tests verifying the spin system construction, spectral density models, Redfield tensor properties (trace preservation, Hermiticity, detailed balance, secular vs. non-secular structure), nonadiabatic decomposition (plateau invariance, analytical benchmarks, decomposition consistency), and the Solomon-Bloembergen analytical relaxation rates.

---

## References

- Mandal, A., & Hunt, K. L. C. (2012). Adiabatic and nonadiabatic contributions to the energy of a system subject to a time-dependent perturbation. *J. Chem. Phys.*, 137, 164109.
- Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
- Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields. *J. Chem. Phys.*, 158, 164107.
- Korzhnev, D. M., Orekhov, V. Y., & Kay, L. E. (2005). Off-resonance $R_{1\rho}$ NMR studies of exchange dynamics in proteins with low spin-lock fields. *JACS*, 127, 713–721.
- Lipari, G., & Szabo, A. (1982). Model-free approach to the interpretation of nuclear magnetic resonance relaxation in macromolecules. *JACS*, 104, 4546–4559.
- Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Physical Review A*, 78, 022106.
- Campaioli, F., Cole, J. H., & Hapuarachchi, H. (2024). A tutorial on quantum master equations. *PRX Quantum*, 5, 020202.