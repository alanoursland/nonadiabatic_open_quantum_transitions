# Nonadiabatic Correction to Pressure-Broadened Spectral Line Shapes

## Overview

This document derives how the Mandal-Hunt nonadiabatic framework modifies the relaxation matrix $\mathbf{W}$ that governs pressure-broadened line mixing in gas-phase molecular spectroscopy. The result connects the time-domain nonadiabatic decomposition (Landau-Lifshitz / Mandal-Hunt) to the frequency-domain line-shape formalism (Anderson-Tsao-Curnutte / Ben-Reuven / Fano).

The central result is:

$$\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \, \mathcal{F}$$

where $\mathcal{F} = (\mathbf{I} - \mathbf{M}) \otimes (\mathbf{I}^* - \mathbf{M}^*)$ is a superoperator that filters out the adiabatic polarization created by the mean collisional field before the bath can act. The filter $\mathbf{M}$ is built from the same pressure-shift physics already present in the standard theory.

This produces two testable predictions:

1. **First-order inter-branch coupling.** The nonadiabatic filter mixes R-branch and P-branch line spaces at order $V/\Delta E$, creating off-diagonal relaxation pathways that the standard secular approximation forbids.

2. **Second-order modification of intra-branch line mixing.** The standard off-diagonal $W_{JJ'}$ within a branch are modified at order $(V/\Delta E)^2$, with a $J$-dependent correction that is largest for low-$J$ lines.

Both predictions are quantitatively significant for molecules with small rotational constants (CO₂, N₂O) at pressures above ~5 atm, and are testable against published HITRAN line-mixing data.

---

## The Problem: Two Formalisms, One Physics

### The Mandal-Hunt framework (time domain)

A quantum system with Hamiltonian $H_0$ is driven by a time-dependent perturbation $V(t)$ and coupled to a thermal bath. The Landau-Lifshitz integration by parts decomposes each transition amplitude into an adiabatic component $a_k(t)$ (reversible polarization response to $V$) and a nonadiabatic component $b_k(t)$ (genuine irreversible transition). The bath should relax only the nonadiabatic part: $\mathcal{D}[\sigma_{\text{nad}}]$, not $\mathcal{D}[\sigma]$.

### The ATC line-shape formalism (frequency domain)

The absorption coefficient of a pressure-broadened molecular gas is:

$$I(\omega) \propto \frac{1}{\pi} \operatorname{Im} \left[ \mathbf{d}^T \cdot \frac{1}{\omega \mathbf{1} - \boldsymbol{\omega}_0 - i\mathbf{W}} \cdot \boldsymbol{\rho}_{\text{eq}} \right]$$

where $\boldsymbol{\omega}_0$ is the diagonal matrix of line positions, $\mathbf{d}$ is the transition dipole vector, $\boldsymbol{\rho}_{\text{eq}}$ is the thermal population vector, and $\mathbf{W}$ is the relaxation matrix computed from the Redfield tensor projected onto the space of optical coherences. The relaxation matrix $\mathbf{W}$ encodes all collisional effects: pressure broadening (diagonal), line mixing (off-diagonal), and pressure shifts (imaginary parts).

### The question

How does the nonadiabatic rule — "the bath relaxes only $\sigma_{\text{nad}}$" — modify $\mathbf{W}$?

---

## The Linear Response Constraint

The most natural first attempt is to identify the probe radiation field as the perturbation $V(t)$ and apply the Landau-Lifshitz decomposition to the probe-induced coherences. This fails.

In linear response, the probe field creates an optical coherence with both a dispersive (in-phase, adiabatic) component and an absorptive (out-of-phase, nonadiabatic) component. The standard line-shape formula already isolates the absorptive component by taking $\operatorname{Im}[\ldots]$. The dispersive component contributes to the refractive index, not the absorption.

Therefore, applying the nonadiabatic decomposition to the probe-matter interaction produces zero correction to the absorption line shape. The standard formula already computes the nonadiabatic (resonant) response of the molecule to the probe. The correction must come from elsewhere.

---

## The Mean-Field / Fluctuation Decomposition

The relevant perturbation is not the probe field but the collisional interaction itself. In a gas cell, buffer-gas collisions simultaneously mix rotational states (perturbation) and provide thermal relaxation (bath). The Mandal-Hunt framework requires these roles to be separated.

The separation is achieved by decomposing the collisional interaction:

$$V_{\text{int}}(t) = \langle V_{\text{int}} \rangle + \delta V_{\text{int}}(t)$$

**The mean field $\langle V_{\text{int}} \rangle$** is the time-averaged collisional coupling, proportional to the buffer gas pressure $P$. It is a static perturbation that creates a steady-state adiabatic mixing of the rotational states. In the standard ATC formalism, this same mean field produces the pressure-induced line shifts $\delta_J$. Its strength is $V = \gamma_L P$ where $\gamma_L$ is the pressure-broadening coefficient.

**The fluctuations $\delta V_{\text{int}}(t)$** are the rapid, stochastic, collision-to-collision variations around the mean. Their correlation function generates the Redfield tensor and therefore the standard relaxation matrix $\mathbf{W}^{\text{standard}}$. They play the role of the thermal bath.

**The nonadiabatic rule** then states: the fluctuation bath $\delta V_{\text{int}}(t)$ must relax only the nonadiabatic components of the density matrix — the part that represents genuine rotational transitions — and must not act on the adiabatic polarization created by the mean field $\langle V_{\text{int}} \rangle$.

This decomposition is not ad hoc. The mean-field and fluctuation components already appear separately in the standard ATC theory: the mean field produces pressure shifts (real part of the forward-scattering amplitude), while the fluctuations produce pressure broadening and line mixing (the relaxation matrix). The nonadiabatic framework adds a physical constraint on how these two components interact: the bath (fluctuations) must not relax the state created by the perturbation (mean field).

---

## Construction of the Nonadiabatic Filter

### The adiabatic mixing matrix

The mean-field perturbation $\langle V_{\text{int}} \rangle$ mixes the rotational eigenstates of $H_0$. In first-order perturbation theory, the adiabatic coefficient for state $|k\rangle$ due to coupling to state $|j\rangle$ is:

$$a_k^{(j)} = \frac{\langle k | \langle V_{\text{int}} \rangle | j \rangle}{E_j - E_k}$$

Define the mixing matrix $\mathbf{M}$ with elements:

$$M_{kj} = \frac{\langle k | \langle V_{\text{int}} \rangle | j \rangle}{E_j - E_k} \quad (k \neq j), \qquad M_{kk} = 0$$

The nonadiabatic state amplitudes are:

$$\mathbf{b} = (\mathbf{I} - \mathbf{M}) \, \mathbf{c}$$

where $\mathbf{c}$ are the full (Dirac) amplitudes.

### The filter superoperator

The density matrix transforms as an outer product: $\rho = \mathbf{c} \otimes \mathbf{c}^*$. The nonadiabatic density matrix is:

$$\sigma_{\text{nad}} = \mathbf{b} \otimes \mathbf{b}^* = [(\mathbf{I} - \mathbf{M}) \otimes (\mathbf{I}^* - \mathbf{M}^*)] \, \rho$$

Define the **nonadiabatic filter superoperator**:

$$\mathcal{F} = (\mathbf{I} - \mathbf{M}) \otimes (\mathbf{I}^* - \mathbf{M}^*)$$

so that $\sigma_{\text{nad}} = \mathcal{F} \rho$.

Note: $\mathcal{F}$ is not idempotent ($\mathcal{F}^2 \neq \mathcal{F}$), so it is a filter, not a projector. This is consistent with the discussion in the gaps analysis (Gap 1): the nonadiabatic operation is affine at the density matrix level. At the amplitude level, $(\mathbf{I} - \mathbf{M})$ is idempotent only if $\mathbf{M}^2 = \mathbf{M}$, which is not generally true. For the perturbative regime where $\|\mathbf{M}\| \ll 1$, the distinction between $\mathcal{F}$ and a true projector is second-order in $\mathbf{M}$.

### The modified master equation

The standard master equation for optical coherences is:

$$\frac{d\rho}{dt} = -i\boldsymbol{\omega}_0 \rho - \mathbf{W}^{\text{standard}} \rho$$

The nonadiabatic modification replaces the argument of the dissipator:

$$\frac{d\rho}{dt} = -i\boldsymbol{\omega}_0 \rho - \mathbf{W}^{\text{standard}} (\mathcal{F} \rho)$$

The coherent evolution ($-i\boldsymbol{\omega}_0 \rho$) acts on the full density matrix — the probe interacts with the whole molecule, adiabatic polarization included. The dissipator ($\mathbf{W}$) acts only on the nonadiabatic part ($\mathcal{F}\rho$). The adiabatic polarization evolves coherently but is shielded from collisional relaxation.

By direct identification, the effective nonadiabatic relaxation matrix is:

$$\boxed{\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \, \mathcal{F}}$$

The standard relaxation matrix $\mathbf{W}^{\text{standard}}$ is a property of the bath (collision fluctuations) and is structurally unchanged. The filter $\mathcal{F}$ modifies what the bath acts on, producing an effective relaxation matrix with qualitatively new structure.

---

## Structure of the Correction

### Expansion in orders of $\mathbf{M}$

Expanding $\mathcal{F}$ to second order:

$$\mathcal{F} = \mathbf{I} \otimes \mathbf{I} - (\mathbf{M} \otimes \mathbf{I} + \mathbf{I} \otimes \mathbf{M}^*) + \mathbf{M} \otimes \mathbf{M}^* + \ldots$$

The correction to $W$ is:

$$\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} - \mathbf{W}^{\text{standard}}(\mathbf{M} \otimes \mathbf{I} + \mathbf{I} \otimes \mathbf{M}^*) + \mathbf{W}^{\text{standard}}(\mathbf{M} \otimes \mathbf{M}^*) + \ldots$$

### First-order correction: inter-branch coupling

The first-order term $-\mathbf{W}^{\text{standard}}(\mathbf{M} \otimes \mathbf{I} + \mathbf{I} \otimes \mathbf{M}^*)$ generates new off-diagonal elements in the relaxation matrix.

Consider an R-branch optical coherence $\rho_J = |J+1\rangle\langle J|$ (a line with $\Delta J = +1$). The first-order correction generates:

$$\boldsymbol{\Delta}\rho_J = \mathbf{M}|J+1\rangle\langle J| + |J+1\rangle\langle J|\mathbf{M}^{\dagger}$$

For a collisional potential dominated by $P_2(\cos\theta)$ anisotropy (the leading term for most linear molecule–atom interactions), $\mathbf{M}$ connects states with $\Delta J = \pm 2$. Applying $\mathbf{M}$ to the R-branch coherence generates four new coherences:

| Generated coherence | $\Delta J$ | Branch |
|:---|:---|:---|
| $\|J+3\rangle\langle J\|$ | +3 | Neither R nor P |
| $\|J-1\rangle\langle J\|$ | $-1$ | **P-branch** |
| $\|J+1\rangle\langle J+2\|$ | $-1$ | **P-branch** |
| $\|J+1\rangle\langle J-2\|$ | +3 | Neither R nor P |

The first-order nonadiabatic correction creates coupling between R-branch and P-branch lines. Standard ATC theory treats these branches independently because their Bohr frequencies are well-separated and the secular approximation decouples them. The nonadiabatic framework breaks this isolation through the adiabatic polarization created by the mean collisional field.

The magnitude of this inter-branch coupling is of order $|M| \sim V/\Delta E = \gamma_L P / [2B(J+1)]$. For CO₂ at 5 atm: $|M| \sim 0.45$ at $J = 0$, decreasing as $1/(J+1)$ for higher $J$.

### Second-order correction: modified intra-branch mixing

The second-order term $\mathbf{W}^{\text{standard}}(\mathbf{M} \otimes \mathbf{M}^*)$ generates corrections within the R-branch line space. Applying $\mathbf{M}$ twice to an R-branch coherence returns to the R-branch:

$$\mathbf{M}|J+1\rangle\langle J|\mathbf{M}^{\dagger} \to |J+3\rangle\langle J+2|, \quad |J-1\rangle\langle J-2|, \quad \ldots$$

These are standard R-branch lines. The correction to the intra-branch relaxation matrix element $W_{JJ'}$ is of order $|M|^2 \sim (V/\Delta E)^2$. For CO₂ at 5 atm, this is $\sim 20\%$ at $J = 0$ and decreases as $1/(J+1)^2$ for higher $J$.

This $J$-dependent modification of the standard line-mixing parameters is a distinctive prediction: low-$J$ lines are more affected than high-$J$ lines, with the correction falling off quadratically.

---

## Consistency Checks

### Sum rule (intensity conservation)

The standard sum rule states that collisions conserve total band intensity:

$$\sum_J d_J \rho_J^{\text{eq}} \, W_{JJ'} = 0 \quad \text{for all } J'$$

In vector notation: $\mathbf{v}^T \mathbf{W} = 0$ where $v_J = d_J \rho_J^{\text{eq}}$.

For the nonadiabatic relaxation matrix:

$$\mathbf{v}^T \mathbf{W}^{\text{nad}} = \mathbf{v}^T \mathbf{W}^{\text{standard}} \mathcal{F} = (\mathbf{v}^T \mathbf{W}^{\text{standard}}) \mathcal{F} = \mathbf{0} \cdot \mathcal{F} = \mathbf{0}$$

**The sum rule is satisfied.** The proof is that the sum rule is a property of the left null space of $\mathbf{W}^{\text{standard}}$, and right-multiplication by $\mathcal{F}$ cannot affect the left null space. Intensity conservation is guaranteed by the structure of the bath, independent of the nonadiabatic filter.

### Detailed balance

The standard detailed balance condition is $W_{JJ'} d_{J'} \rho_{J'} = W_{J'J} d_J \rho_J$, ensuring relaxation toward the Boltzmann distribution at the bare rotational energies.

$\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$ violates this condition at order $|\mathbf{M}|^2$. This violation is expected and physically correct: the mean collisional field $\langle V_{\text{int}} \rangle$ shifts the effective rotational energies (pressure shifts), so the appropriate thermal equilibrium is with respect to the dressed energies $E_J^{\text{dressed}} = E_J + \langle J | \langle V_{\text{int}} \rangle | J \rangle$, not the bare energies.

The detailed balance condition for the nonadiabatic relaxation matrix should be:

$$W_{JJ'}^{\text{nad}} \, d_{J'}^{\text{eff}} \, \rho_{J'}^{\text{dressed}} = W_{J'J}^{\text{nad}} \, d_J^{\text{eff}} \, \rho_J^{\text{dressed}}$$

where $\rho_J^{\text{dressed}} \propto (2J+1) e^{-E_J^{\text{dressed}}/k_BT}$. The violation of bare-energy detailed balance is second-order in $V/\Delta E$ — the same order as the pressure shift correction to the thermal populations — and does not affect the first-order inter-branch prediction.

---

## Connection to Known Spectroscopic Phenomena

### Sub-Lorentzian far wings

Pressure-broadened line shapes are experimentally observed to fall off faster than Lorentzian in the far wings. The standard theory requires empirical $\chi$-factors to correct this. Various physical mechanisms have been proposed (duration-of-collision effects, breakdown of the impact approximation).

The nonadiabatic framework provides an additional mechanism. In the far wings (large detuning from line center), the optical coherence created by the probe is predominantly off-resonant — it is an adiabatic polarization response rather than a genuine transition. The nonadiabatic filter $\mathcal{F}$ reduces the fraction of this coherence that is subject to collisional relaxation, causing the effective broadening to decrease with detuning. This produces sub-Lorentzian wing behavior without empirical correction factors.

### Known discrepancies in CO₂ line mixing

Hartmann and collaborators have documented systematic discrepancies between standard Redfield predictions and measured HITRAN line shapes for CO₂ at elevated pressures, particularly in Q-branches. These discrepancies have been attributed to various beyond-Redfield effects. The nonadiabatic correction provides a candidate explanation: the standard Redfield treatment over-relaxes the density matrix by including adiabatic polarization in the bath coupling, producing incorrect line-mixing coefficients at low $J$.

---

## Testable Predictions

### Prediction 1: Inter-branch mixing signature

The nonadiabatic framework predicts that R-branch and P-branch relaxation matrices are not independent at elevated pressures. The coupling is first-order in $V/\Delta E$ and mediated by the $P_2(\cos\theta)$ anisotropy of the intermolecular potential.

**Observable:** At pressures where $V/\Delta E \gtrsim 0.1$ (e.g., CO₂ in N₂ at $\geq 5$ atm), the R-branch and P-branch line shapes should show correlated deviations from independent-branch predictions. Specifically, intensity "borrowed" from R-branch lines should appear in P-branch lines, and vice versa, with the magnitude proportional to $\gamma_L P / [2B(J+1)]$.

**How to test:** Compute the absorption spectrum of the CO₂ $\nu_3$ band at 5–10 atm using $\mathbf{W}^{\text{nad}}$ (with inter-branch coupling) and $\mathbf{W}^{\text{standard}}$ (without), and compare both against measured spectra. The inter-branch coupling should improve the residuals in the overlap region between R and P branches.

### Prediction 2: $J$-dependent line-mixing correction

Within the R-branch (or P-branch), the nonadiabatic correction to $W_{JJ'}$ scales as $1/(J+1)^2$. Low-$J$ lines have the largest correction; high-$J$ lines are essentially uncorrected.

**Observable:** The first-order line-mixing $Y$-coefficients tabulated in HITRAN should show systematic $J$-dependent deviations from standard Redfield predictions, with the largest deviations at low $J$. The deviations should scale as $(\gamma_L P)^2 / [2B(J+1)]^2$.

**How to test:** Compare measured $Y_J$ against standard Redfield predictions and nonadiabatic predictions for CO₂, N₂O, and CO at multiple pressures. The nonadiabatic prediction has no free parameters — $\gamma_L$, $B$, and the intermolecular potential are all independently known.

### Prediction 3: Sub-Lorentzian wings without $\chi$-factors

The nonadiabatic framework predicts sub-Lorentzian far-wing behavior as a consequence of the adiabatic/nonadiabatic decomposition, without requiring empirical correction factors.

**How to test:** Compute the far-wing absorption (5–50 cm⁻¹ from line center) for CO₂ at elevated pressures using $\mathbf{W}^{\text{nad}}$ and compare against measured absorption coefficients. If the nonadiabatic prediction matches the measured sub-Lorentzian wings without $\chi$-factors, it provides strong evidence for the framework.

---

## Implementation

The nonadiabatic line-shape calculation requires:

1. **The standard relaxation matrix $\mathbf{W}^{\text{standard}}$** — computed from the Redfield tensor using published intermolecular potentials and the existing `SpectralDensityBath` infrastructure.

2. **The mixing matrix $\mathbf{M}$** — computed from the mean collisional field $\langle V_{\text{int}} \rangle$ and the rotational energy levels. The matrix elements $M_{kj} = \langle k | \langle V_{\text{int}} \rangle | j \rangle / (E_j - E_k)$ require the same intermolecular potential used for $\mathbf{W}$.

3. **The filter superoperator $\mathcal{F} = (\mathbf{I} - \mathbf{M}) \otimes (\mathbf{I}^* - \mathbf{M}^*)$** — constructed from $\mathbf{M}$.

4. **The product $\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$** — a matrix multiplication in Liouville space.

5. **The absorption line shape** from the standard formula with $\mathbf{W}^{\text{nad}}$ replacing $\mathbf{W}^{\text{standard}}$.

Steps 1 and 2 are the computationally intensive parts. Steps 3–5 are algebraic operations on the matrices produced by steps 1–2.

---

## Regime of Validity

The construction requires:

- **Perturbative mixing:** $\|\mathbf{M}\| < 1$, i.e., $\gamma_L P / [2B(J+1)] < 1$ for the relevant $J$ values. For CO₂, this limits pressures to roughly $\leq 10$ atm for the lowest $J$ levels. Above this, the first-order Landau-Lifshitz decomposition breaks down and higher-order or non-perturbative methods are needed.

- **Markovian bath:** The collision fluctuations must be fast compared to the system dynamics. This is the standard Born-Markov condition, satisfied when the collision duration ($\sim 1$ ps) is short compared to the inverse Bohr frequencies ($\sim 10$ ps for CO₂).

- **Mean-field validity:** The gas must be dense enough that many simultaneous long-range interactions create a well-defined average field. This is satisfied above $\sim 1$ atm for most buffer gases.

---

## References

- Mandal, A. & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
- Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields. *J. Chem. Phys.*, 158, 164107.
- Hartmann, J.-M., Boulet, C., & Robert, D. (2021). *Collisional Effects on Molecular Spectra* (2nd ed.). Elsevier.
- Ben-Reuven, A. (1966). Impact broadening of microwave spectra. *Phys. Rev.*, 145, 7.
- Fano, U. (1963). Pressure broadening as a prototype of relaxation. *Phys. Rev.*, 131, 259.
- Schaller, G. & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Phys. Rev. A*, 78, 022106.
- Gordon, I. E. et al. (2022). The HITRAN2020 molecular spectroscopic database. *JQSRT*, 277, 107949.