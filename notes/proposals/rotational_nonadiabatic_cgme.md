# The Nonadiabatic CGME for a Rotational Line-Mixing System

## Purpose

The three-level nonadiabatic CGME document reduced the general nonadiabatic Nakajima-Zwanzig framework to a 3-level system driven by a classical field and weakly coupled to a thermal bath. That reduction led to the NMR $R_{1\rho}$ experiment plan, which returned NEGATIVE at Phase 0 because $V/\Delta E \sim 10^{-5}$.

This document performs the analogous reduction for a different physical system: pressure-broadened rotational line mixing in a gas-phase linear molecule (CO$_2$, N$_2$O, CO) colliding with a buffer gas (N$_2$). The Phase 0A feasibility estimate returned POSITIVE for CO$_2$ at $\geq 5$ atm, with a 20% correction at $J = 0$. The goal here is to write out the master equation explicitly for this system, with all matrix elements hardcoded for a small truncation ($J_{\max} = 2$, giving a 3-level system with $J = 0, 1, 2$), so that the equation can be implemented, tested against published HITRAN data, and compared with the standard (Dirac-basis) treatment.

The structure follows the three-level CGME document step by step, substituting the rotational physics for the generic three-level system at each stage.

---

## The Physical Scenario

In pressure-broadened molecular spectroscopy, a gas-phase molecule absorbs infrared radiation while simultaneously undergoing collisions with buffer gas molecules. These two processes happen concurrently:

**The radiation field** drives transitions between rotational levels. For a linear molecule like CO$_2$, the allowed transitions are $\Delta J = \pm 1$ (R-branch and P-branch). The radiation field is the perturbation $V$ — it couples the molecular dipole to the electromagnetic field.

**Buffer gas collisions** provide the thermal bath. Collisions interrupt the molecule's coherent interaction with the radiation field, broaden the spectral lines, and — crucially — couple different rotational transitions to each other. This coupling between transitions is **line mixing**: the absorption profile of a group of closely spaced lines deviates from the sum of independent Lorentzians because collisions transfer coherence between transitions.

The standard theoretical treatment of line mixing uses the Redfield relaxation matrix (called the $W$-matrix in the spectroscopy literature). The off-diagonal elements $W_{JJ'}$ describe how collisions couple the $J \to J+1$ and $J' \to J'+1$ transitions. The standard treatment computes these elements using Dirac's transition coefficients. The nonadiabatic framework computes them using the Landau-Lifshitz decomposition, which removes the adiabatic polarization and retains only the genuine collisional transitions.

The observable is the absorption line shape, which depends on the $W$-matrix through a resolvent:

$$F(\nu) \propto \text{Im}\left[\mathbf{d}^T \left( i(\nu - \boldsymbol{\nu}_0) + \mathbf{W} \right)^{-1} \boldsymbol{\rho}_{\text{eq}} \right]$$

where $\boldsymbol{\nu}_0$ is the vector of line-center frequencies, $\mathbf{d}$ is the vector of transition dipole strengths, and $\boldsymbol{\rho}_{\text{eq}}$ is the thermal population vector. The nonadiabatic correction enters through the matrix $\mathbf{W}$.

---

## The System: Truncated Rigid Rotor ($J = 0, 1, 2$)

### Energy Levels

A linear rigid rotor has energy levels:

$$E_J = B \, J(J+1) \quad \text{(in cm}^{-1}\text{)}$$

where $B$ is the rotational constant. For CO$_2$, $B = 0.3902$ cm$^{-1}$.

Truncating at $J_{\max} = 2$ gives a 3-level system:

| State | $J$ | $E_J$ (cm$^{-1}$) | $E_J$ (rad/s) |
|:---|:---|:---|:---|
| $\|0\rangle$ | 0 | 0 | 0 |
| $\|1\rangle$ | 1 | $2B = 0.7804$ | $1.47 \times 10^{11}$ |
| $\|2\rangle$ | 2 | $6B = 2.3412$ | $4.41 \times 10^{11}$ |

### Bohr Frequencies

$$\omega_{10} = 2B \cdot C_{\text{conv}} = 0.7804 \text{ cm}^{-1} = 1.47 \times 10^{11} \text{ rad/s}$$

$$\omega_{20} = 6B \cdot C_{\text{conv}} = 2.3412 \text{ cm}^{-1} = 4.41 \times 10^{11} \text{ rad/s}$$

$$\omega_{21} = 4B \cdot C_{\text{conv}} = 1.5608 \text{ cm}^{-1} = 2.94 \times 10^{11} \text{ rad/s}$$

where $C_{\text{conv}} = 2\pi c = 1.884 \times 10^{11}$ rad/s/cm$^{-1}$.

Note the specific structure of the rigid rotor: $\omega_{20} = 3\omega_{10}$ and $\omega_{21} = 2\omega_{10}$. The level spacing increases linearly with $J$: the gap between $J = 0$ and $J = 1$ is $2B$, between $J = 1$ and $J = 2$ is $4B$. This means the nonadiabatic correction is largest for the lowest transition, where the gap is smallest.

### The Unperturbed Hamiltonian

$$H_0 = \begin{pmatrix} 0 & 0 & 0 \\ 0 & \hbar\omega_{10} & 0 \\ 0 & 0 & \hbar\omega_{20} \end{pmatrix}$$

This is already diagonal in the $|J\rangle$ basis. The eigenvectors are trivial: $|0\rangle = (1,0,0)^T$, $|1\rangle = (0,1,0)^T$, $|2\rangle = (0,0,1)^T$. No basis transformation is needed.

---

## Step 1: The Perturbation and the Nonadiabatic Coefficients

### Identifying the Perturbation

In the three-level CGME document, the perturbation $V(t)$ was a generic time-dependent external field. For the HITRAN system, the perturbation has two components that must be carefully distinguished:

**The radiation field** couples the molecule to the electromagnetic field through the transition dipole moment. For a linear molecule, this drives $\Delta J = \pm 1$ transitions. In the 3-level system, the radiation field couples $|0\rangle \leftrightarrow |1\rangle$ and $|1\rangle \leftrightarrow |2\rangle$, but not $|0\rangle \leftrightarrow |2\rangle$ (since $|\Delta J| = 2$ is forbidden for dipole transitions).

**The collisional coupling** is the perturbation that causes line mixing. Buffer gas collisions couple adjacent rotational states with a strength proportional to the pressure-broadening coefficient $\gamma_L$ times the gas pressure $P$. This is the perturbation whose nonadiabatic decomposition we need.

The key realization: in the line-mixing problem, the roles are inverted relative to the NMR problem. In NMR, the spin-lock field was the perturbation $V$ and molecular tumbling was the bath. Here, **collisions play both roles simultaneously** — they are both the perturbation that couples rotational states (the off-diagonal $W$-matrix elements) and the bath that causes relaxation and dephasing (the diagonal elements, i.e., the pressure-broadened linewidths).

More precisely: the collisional interaction during a single collision event can be decomposed into a coherent part (which mixes rotational states — this is the perturbation $V$) and a stochastic part (which causes dephasing and population relaxation — this is the bath). The coherent part is what the Landau-Lifshitz decomposition acts on. The stochastic part provides the spectral density for the Redfield tensor.

In the Anderson-Tsao-Curnutte formalism, this decomposition is standard: the $S$-matrix for a single collision is expanded in terms of the intermolecular potential, and the off-diagonal elements give the line-coupling coefficients while the diagonal elements give the broadening and shifting.

### The Collisional Coupling Operator

The coherent collisional coupling is modeled as a nearest-neighbor ($\Delta J = \pm 1$) interaction:

$$V_{\text{coll}} = V_c \begin{pmatrix} 0 & 1 & 0 \\ 1 & 0 & 1 \\ 0 & 1 & 0 \end{pmatrix}$$

where $V_c = \gamma_L \times P$ is the coupling strength in cm$^{-1}$, converted to energy units. For CO$_2$ with N$_2$ buffer gas: $\gamma_L = 0.070$ cm$^{-1}$/atm.

| Pressure (atm) | $V_c$ (cm$^{-1}$) | $V_c$ (rad/s) |
|:---|:---|:---|
| 1 | 0.070 | $1.32 \times 10^{10}$ |
| 5 | 0.350 | $6.59 \times 10^{10}$ |
| 10 | 0.700 | $1.32 \times 10^{11}$ |

Note: this $J$-independent coupling is the Phase 0 approximation. A realistic treatment would use $J$-dependent coupling elements from the intermolecular potential, but for the reduced master equation the structure is the same — only the numerical values of the matrix elements change.

### The Landau-Lifshitz Decomposition

Applying the integration by parts to the first-order perturbation theory expression for the collisional perturbation $V_{\text{coll}}$:

**Adiabatic coefficients** (instantaneous polarization response to collisions):

$$a_1^{(1)}(t) = \frac{\langle 1 | V_{\text{coll}}(t) | 0 \rangle}{\hbar\omega_{10}} \, e^{i\omega_{10} t} = \frac{V_c}{\hbar\omega_{10}} \, e^{i\omega_{10} t}$$

$$a_2^{(1)}(t) = \frac{\langle 2 | V_{\text{coll}}(t) | 0 \rangle}{\hbar\omega_{20}} \, e^{i\omega_{20} t} = 0$$

The second equation vanishes because $\langle 2 | V_{\text{coll}} | 0 \rangle = 0$ — the collisional coupling is nearest-neighbor only, so there is no direct coupling between $J = 0$ and $J = 2$.

Similarly, for transitions from $|1\rangle$:

$$a_2^{(1)}(t) \big|_{\text{from }|1\rangle} = \frac{\langle 2 | V_{\text{coll}}(t) | 1 \rangle}{\hbar\omega_{21}} \, e^{i\omega_{21} t} = \frac{V_c}{\hbar\omega_{21}} \, e^{i\omega_{21} t}$$

**Adiabatic mixing fractions** (the quantities that matter for the correction):

$$|a_1|^2 = \left(\frac{V_c}{\hbar\omega_{10}}\right)^2 = \left(\frac{V_c}{2B}\right)^2$$

$$|a_2|^2\big|_{\text{from }|1\rangle} = \left(\frac{V_c}{\hbar\omega_{21}}\right)^2 = \left(\frac{V_c}{4B}\right)^2$$

For CO$_2$ at 5 atm:

| Transition pair | $V_c$ (cm$^{-1}$) | $\Delta E$ (cm$^{-1}$) | $V_c/\Delta E$ | $\|a\|^2$ |
|:---|:---|:---|:---|:---|
| $J=0 \leftrightarrow J=1$ | 0.350 | 0.780 | 0.449 | 0.201 |
| $J=1 \leftrightarrow J=2$ | 0.350 | 1.561 | 0.224 | 0.050 |

These are the fractions of the density matrix that are adiabatic polarization — the part the bath should *not* act on in the nonadiabatic framework.

**Nonadiabatic coefficients** $b_k(t)$: these are the residual integrals that depend on $dV_{\text{coll}}/dt$. During steady-state absorption (where the collisional perturbation is statistically stationary), the nonadiabatic populations reach a steady state determined by the balance between the coherent collisional driving and the dissipative dephasing. The exact values require propagating the master equation; what matters for constructing the equation is that $\sigma_{\text{nad}}$ is built from $b_k b_l^*$, not $c_k c_l^*$.

---

## Step 2: The Nonadiabatic Density Matrix

The nonadiabatic density matrix for the 3-level rotational system is:

$$\sigma_{\text{nad}} = \begin{pmatrix} 1 - |b_1|^2 - |b_2|^2 & b_0^* b_1 & b_0^* b_2 \\ b_1^* b_0 & |b_1|^2 & b_1^* b_2 \\ b_2^* b_0 & b_2^* b_1 & |b_2|^2 \end{pmatrix}$$

where $b_0 = \sqrt{1 - |b_1|^2 - |b_2|^2}$ (normalization), and $b_1$, $b_2$ are the nonadiabatic amplitudes for states $|J=1\rangle$ and $|J=2\rangle$.

The populations $[\sigma_{\text{nad}}]_{JJ} = |b_J|^2$ are the genuine collisional-transfer populations, with the adiabatic polarization removed. The coherences $[\sigma_{\text{nad}}]_{JJ'}$ describe quantum coherence between states that have been genuinely populated by collisions.

The critical property: these populations do not include the $|a_J|^2$ contribution. When the nonadiabatic Redfield tensor acts on this density matrix, it does not attempt to relax the adiabatic polarization — only the genuine excitations are subject to bath-induced relaxation.

---

## Step 3: The Nonadiabatic Redfield Tensor

### The Bath: Collisional Dephasing and Relaxation

The bath for the rotational system is the stochastic component of buffer gas collisions. It is characterized by a spectral density $J_{\text{coll}}(\omega)$ that reflects the collision dynamics.

For gas-phase collisions at room temperature, the bath correlation time is the collision duration $\tau_c \sim 1$ ps. The collisional spectral density is modeled as a Lorentzian:

$$J_{\text{coll}}(\omega) = \frac{2\gamma_0 \tau_c}{1 + \omega^2 \tau_c^2}$$

where $\gamma_0$ is the zero-frequency collision rate, proportional to pressure. For CO$_2$ in N$_2$ at 296 K: $\gamma_0 \approx \gamma_L \times P \times C_{\text{conv}} \approx 0.070 \times P \times 1.884 \times 10^{11}$ rad/s.

### The System-Bath Coupling Operator

The system-bath coupling operator $S$ describes how collisions couple the rotational states. For the dominant dipolar contribution ($\Delta J = \pm 1$):

$$S = \begin{pmatrix} 0 & s_{01} & 0 \\ s_{01} & 0 & s_{12} \\ 0 & s_{12} & 0 \end{pmatrix}$$

where $s_{01} = \langle 0 | S | 1 \rangle$ and $s_{12} = \langle 1 | S | 2 \rangle$ are the coupling matrix elements. In the simplest (Phase 0) approximation these are $J$-independent and equal to a dimensionless coupling constant. A realistic treatment uses $J$-dependent matrix elements derived from the intermolecular potential.

### Transition Rates

For each pair of states $(k, l)$ with $E_k > E_l$, the downward and upward rates are:

$$\gamma_{k \to l} = \frac{2\pi}{\hbar^2} |s_{lk}|^2 \, J_{\text{coll}}(\omega_{kl}) \, [\bar{n}(\omega_{kl}) + 1]$$

$$\gamma_{l \to k} = \frac{2\pi}{\hbar^2} |s_{lk}|^2 \, J_{\text{coll}}(\omega_{kl}) \, \bar{n}(\omega_{kl})$$

where $\bar{n}(\omega) = (e^{\hbar\omega/k_BT} - 1)^{-1}$ is the Bose-Einstein occupation number.

For the 3-level rotational system, the three transition pairs and their Bohr frequencies are:

| Pair | Bohr frequency | $\Delta E$ (cm$^{-1}$) | Coupling element |
|:---|:---|:---|:---|
| $(1,0)$: $J=1 \to J=0$ | $\omega_{10} = 2B \cdot C_{\text{conv}}$ | 0.780 | $s_{01}$ |
| $(2,1)$: $J=2 \to J=1$ | $\omega_{21} = 4B \cdot C_{\text{conv}}$ | 1.561 | $s_{12}$ |
| $(2,0)$: $J=2 \to J=0$ | $\omega_{20} = 6B \cdot C_{\text{conv}}$ | 2.341 | 0 (forbidden) |

The $(2,0)$ pair has zero coupling because $|\Delta J| = 2$ transitions are not directly coupled by the dominant dipolar collisional mechanism. This simplifies the Redfield tensor: there are only four nonzero rates (two downward, two upward) instead of six.

### The High-Temperature Regime

At room temperature ($T = 296$ K) and rotational energy scales ($\hbar\omega_{10}/k_BT \approx 0.0038$ for CO$_2$), the thermal occupation numbers are large:

$$\bar{n}(\omega_{10}) = \frac{1}{e^{0.0038} - 1} \approx 263$$

$$\bar{n}(\omega_{21}) = \frac{1}{e^{0.0076} - 1} \approx 131$$

In this limit, $\bar{n} \gg 1$, and the upward and downward rates are nearly equal: $\gamma_{l \to k} / \gamma_{k \to l} = \bar{n}/(\bar{n}+1) \approx 1 - \hbar\omega/k_BT$. Detailed balance is satisfied but the system is near the classical (equipartition) regime.

This has an important consequence: the thermal equilibrium populations are nearly uniform across the three levels (weighted by degeneracy $2J+1$, which we are not including in this simplified treatment — see the note below on degeneracy).

---

## Step 4: The Coarse-Graining Procedure

### The CGME Sinc Factors for the Rotational System

Each element of the Redfield tensor $R_{abcd}$ is multiplied by:

$$R_{abcd}^{\text{CG}} = R_{abcd} \cdot \text{sinc}\left(\frac{(\omega_{ab} - \omega_{cd})\Delta\tau}{2}\right)$$

For the 3-level rotational system, the relevant frequency differences between transitions are:

$$|\omega_{10} - \omega_{21}| = |2B - 4B| \cdot C_{\text{conv}} = 2B \cdot C_{\text{conv}} = \omega_{10}$$

$$|\omega_{10} - \omega_{20}| = |2B - 6B| \cdot C_{\text{conv}} = 4B \cdot C_{\text{conv}} = \omega_{21}$$

$$|\omega_{21} - \omega_{20}| = |4B - 6B| \cdot C_{\text{conv}} = 2B \cdot C_{\text{conv}} = \omega_{10}$$

### Choosing $\Delta\tau$

The coarse-graining timescale must satisfy:

$$\tau_c \ll \Delta\tau \ll T_{\text{relax}}$$

where $\tau_c \sim 1$ ps is the collision duration and $T_{\text{relax}}$ is the population relaxation time.

For CO$_2$ at 5 atm, the pressure-broadened linewidth is $\gamma_L P \approx 0.35$ cm$^{-1}$, giving a dephasing time $T_2 \sim 1/(\gamma_L P \cdot C_{\text{conv}}) \sim 15$ ps. The population relaxation time is longer: $T_1 \sim 100$ ps (estimated from the relaxation rates). So the window is:

$$1 \text{ ps} \ll \Delta\tau \ll 100 \text{ ps}$$

A choice of $\Delta\tau \sim 5$–$10$ ps is appropriate.

### The Sinc Factors Are Near Unity

For CO$_2$ with $\Delta\tau = 5$ ps:

| Frequency pair | $\|\Delta\omega\|$ (cm$^{-1}$) | $\|\Delta\omega\|$ (rad/s) | sinc argument | $\|\text{sinc}\|$ |
|:---|:---|:---|:---|:---|
| $\omega_{10}$ vs $\omega_{21}$ | 0.780 | $1.47 \times 10^{11}$ | 0.37 | 0.955 |
| $\omega_{10}$ vs $\omega_{20}$ | 1.561 | $2.94 \times 10^{11}$ | 0.74 | 0.831 |
| $\omega_{21}$ vs $\omega_{20}$ | 0.780 | $1.47 \times 10^{11}$ | 0.37 | 0.955 |

All sinc factors are of order 1. **The CGME retains essentially all non-secular Redfield elements.** This is the opposite of the NMR case (where all sinc factors were $< 10^{-3}$) and is the reason the HITRAN system is physically interesting: the secular approximation fails here, and the three-way comparison (full Redfield / CGME / secular Lindblad) produces genuinely different predictions.

The secular Lindblad approximation — setting all non-secular elements to zero — corresponds to **ignoring line mixing entirely**. Each rotational transition is treated as independent. This is known to fail for closely spaced lines at elevated pressures, which is precisely where the HITRAN data shows discrepancies with simple models.

---

## Step 5: The Complete Nonadiabatic CGME for the Rotational System

Assembling the pieces, the master equation is:

$$\frac{d\sigma_{\text{nad}}}{dt} = -\frac{i}{\hbar}[H_0 + H_{LS} + V_{\text{nad}}(t), \, \sigma_{\text{nad}}] + \mathcal{L}_{\text{CG}}[\sigma_{\text{nad}}]$$

### The Coherent Terms

**$H_0$** is the diagonal rigid-rotor Hamiltonian defined above.

**$H_{LS}$** is the collisional Lamb shift — the pressure-induced shift of the rotational line positions. In spectroscopic terms, this is the pressure shift coefficient $\delta_L$ (typically $\sim 10^{-3}$ cm$^{-1}$/atm for CO$_2$, much smaller than the broadening coefficient $\gamma_L$). It enters as a diagonal correction to the energy levels:

$$H_{LS} = \text{diag}(0, \, \delta_{10} P, \, \delta_{20} P) \times \hbar C_{\text{conv}}$$

This is a small correction but is included because HITRAN tabulates pressure shifts alongside broadening coefficients, and the full comparison must account for them.

**$V_{\text{nad}}(t)$** is the coherent nonadiabatic driving — the part of the collisional coupling that drives genuine transitions (not adiabatic polarization). For steady-state absorption, this term describes how the radiation field, filtered through the collisional dynamics, drives population transfer between rotational levels. Its matrix elements are related to the nonadiabatic coefficients $b_k(t)$ and the time derivative of the collisional coupling.

### The Coarse-Grained Dissipator

The dissipator $\mathcal{L}_{\text{CG}}$ has Lindblad form:

$$\mathcal{L}_{\text{CG}}[\sigma_{\text{nad}}] = \sum_{\alpha} \gamma_\alpha^{\text{CG}} \left( L_\alpha \, \sigma_{\text{nad}} \, L_\alpha^\dagger - \frac{1}{2}\{L_\alpha^\dagger L_\alpha, \, \sigma_{\text{nad}}\} \right)$$

For the 3-level rotational system with nearest-neighbor coupling, the Lindblad operators and their coarse-grained rates are:

**Population transfer (4 operators):**

$$L_1 = |0\rangle\langle 1|, \quad \gamma_1^{\text{CG}} = \gamma_{1 \to 0} \cdot f_{10}$$

$$L_2 = |1\rangle\langle 0|, \quad \gamma_2^{\text{CG}} = \gamma_{0 \to 1} \cdot f_{10}$$

$$L_3 = |1\rangle\langle 2|, \quad \gamma_3^{\text{CG}} = \gamma_{2 \to 1} \cdot f_{21}$$

$$L_4 = |2\rangle\langle 1|, \quad \gamma_4^{\text{CG}} = \gamma_{1 \to 2} \cdot f_{21}$$

where $f_{kl} = \text{sinc}(0) = 1$ for the secular (diagonal) elements. The secular population-transfer rates are unmodified by coarse-graining.

**Coherence transfer (non-secular terms, retained by the CGME):**

The non-secular elements of the Redfield tensor couple different transition pairs. For the 3-level system, the most important non-secular coupling is between the $0 \leftrightarrow 1$ and $1 \leftrightarrow 2$ transitions. These are the terms that produce line mixing — they are the off-diagonal elements of the $W$-matrix.

In the standard (Dirac-basis) treatment, the non-secular coherence-transfer rate between the $0 \leftrightarrow 1$ and $1 \leftrightarrow 2$ transitions is:

$$W_{01,12}^{\text{Dirac}} = \text{(full Redfield tensor element)}$$

In the nonadiabatic treatment, this rate is modified. The bath acts on $\sigma_{\text{nad}}$ instead of $\sigma$. Because $\sigma_{\text{nad}}$ has the adiabatic polarization removed, the effective off-diagonal relaxation matrix element is reduced:

$$W_{01,12}^{\text{nad}} = W_{01,12}^{\text{Dirac}} \times \left(1 - |a_{01}|^2\right) \times \left(1 - |a_{12}|^2\right)$$

where $|a_{01}|^2 = (V_c/2B)^2$ and $|a_{12}|^2 = (V_c/4B)^2$ are the adiabatic mixing fractions for each transition pair.

This is the central result. The nonadiabatic correction reduces the off-diagonal $W$-matrix elements by a factor that depends on the perturbation ratio $V_c/\Delta E$ for each pair of coupled transitions. The reduction is:

For CO$_2$ at 5 atm:

| $W$-matrix element | $\|a\|^2$ factor 1 | $\|a\|^2$ factor 2 | Reduction | Surviving fraction |
|:---|:---|:---|:---|:---|
| $W_{01,12}$ ($J=0\text{-}1$ coupled to $J=1\text{-}2$) | 0.201 | 0.050 | $0.799 \times 0.950$ | 0.759 (24% reduction) |

At 1 atm:

| $W$-matrix element | $\|a\|^2$ factor 1 | $\|a\|^2$ factor 2 | Reduction | Surviving fraction |
|:---|:---|:---|:---|:---|
| $W_{01,12}$ | 0.008 | 0.002 | $0.992 \times 0.998$ | 0.990 (1% reduction) |

The correction is pressure-dependent (because $V_c = \gamma_L P$) and $J$-dependent (because $\Delta E = 2B(J+1)$ increases with $J$). Both dependencies are distinctive predictions that can be tested against HITRAN data.

---

## The Nonadiabatic Correction to the $W$-Matrix

### Connection to the Line-Shape Formalism

The absorption line shape for a manifold of overlapping transitions is given by the resolvent formula:

$$F(\nu) = \frac{1}{\pi} \text{Im}\left[\sum_{J,J'} d_J \left[ (\nu - \nu_0) \mathbf{I} - i\mathbf{W} \right]^{-1}_{JJ'} \rho_{J'} \, d_{J'} \right]$$

where:

- $\nu$ is the observation frequency
- $\nu_0$ is the band center
- $d_J$ is the transition dipole matrix element for the $J \to J+1$ transition
- $\rho_J$ is the thermal population of level $J$
- $\mathbf{W}$ is the relaxation matrix

The diagonal elements of $\mathbf{W}$ are the pressure-broadened linewidths (plus pressure shifts):

$$W_{JJ} = \gamma_J P - i\delta_J P$$

where $\gamma_J$ and $\delta_J$ are the broadening and shift coefficients for the $J \to J+1$ transition.

The off-diagonal elements of $\mathbf{W}$ are the line-mixing coefficients. These are the Redfield tensor elements that couple different transitions through collisions:

$$W_{JJ'} = \text{(collisional coupling between transitions } J \text{ and } J' \text{)}$$

### The Standard vs. Nonadiabatic $W$-Matrix

**Standard (Dirac-basis) $W$-matrix:** The off-diagonal elements are computed from the bath correlation functions using Dirac's transition coefficients. The bath acts on the full density matrix, including the adiabatic polarization. This is what the Anderson-Tsao-Curnutte formalism and its modern extensions compute.

**Nonadiabatic $W$-matrix:** The off-diagonal elements are computed from the same bath correlation functions, but the bath acts on $\sigma_{\text{nad}}$ — the density matrix with adiabatic polarization removed. For each off-diagonal element coupling transitions $J$ and $J'$:

$$W_{JJ'}^{\text{nad}} = W_{JJ'}^{\text{Dirac}} \times \mathcal{C}_{JJ'}$$

where $\mathcal{C}_{JJ'}$ is the nonadiabatic correction factor. At first order in perturbation theory:

$$\mathcal{C}_{JJ'} = \left(1 - \frac{V_c^2}{[2B(J+1)]^2}\right)^{1/2} \left(1 - \frac{V_c^2}{[2B(J'+1)]^2}\right)^{1/2}$$

This factor is less than 1 for all $J, J'$, and is closest to 1 for large $J$ (where $\Delta E \gg V_c$). For small $J$ and high pressure, the correction is significant.

**The diagonal elements are unchanged.** The pressure-broadened linewidths $\gamma_J P$ are properties of the individual transitions and do not involve the adiabatic polarization of other levels. (More precisely: the diagonal elements involve $|c_J|^2$ for the same state, and the adiabatic correction to the diagonal is a higher-order effect.) The nonadiabatic correction acts specifically on the off-diagonal (line-mixing) elements.

### Explicit $W$-Matrix for the 3-Level System

For CO$_2$ at pressure $P$ in N$_2$ at 296 K, the $2 \times 2$ $W$-matrix for the $R(0)$ and $R(1)$ transitions ($J = 0 \to 1$ and $J = 1 \to 2$) is:

**Standard:**

$$\mathbf{W}^{\text{Dirac}} = \begin{pmatrix} \gamma_0 P - i\delta_0 P & W_{01}^{\text{Dirac}} \\ W_{10}^{\text{Dirac}} & \gamma_1 P - i\delta_1 P \end{pmatrix}$$

**Nonadiabatic:**

$$\mathbf{W}^{\text{nad}} = \begin{pmatrix} \gamma_0 P - i\delta_0 P & W_{01}^{\text{Dirac}} \times \mathcal{C}_{01} \\ W_{10}^{\text{Dirac}} \times \mathcal{C}_{10} & \gamma_1 P - i\delta_1 P \end{pmatrix}$$

where:

$$\mathcal{C}_{01} = \mathcal{C}_{10} = \left(1 - \frac{(\gamma_L P)^2}{(2B)^2}\right)^{1/2} \left(1 - \frac{(\gamma_L P)^2}{(4B)^2}\right)^{1/2}$$

The two frameworks predict the same diagonal elements (linewidths and shifts) but different off-diagonal elements (line mixing). The difference is measurable when $\gamma_L P / 2B$ is not negligible — which is the regime CO$_2$ at $\geq 5$ atm occupies.

---

## What This Equation Satisfies

### Complete Positivity

Guaranteed by the Lindblad form of the coarse-grained dissipator, exactly as in the three-level CGME document. The sinc factors are all positive and less than 1, so the coarse-grained rates are non-negative.

### Gauge Invariance

Inherited from the nonadiabatic coefficients $b_k(t)$. The collisional perturbation does not have a gauge ambiguity in the same sense as electromagnetic coupling (it is a short-range molecular interaction), but the principle is the same: the nonadiabatic decomposition removes the unphysical adiabatic component that would contaminate the relaxation dynamics.

### Thermodynamic Consistency

When collisions are the only process (no radiation field), the system thermalizes to the Boltzmann distribution. The nonadiabatic CGME ensures that the bath drives the system toward the correct thermal equilibrium — the Boltzmann distribution of the field-dressed system — without the fictitious heat flow that arises when the bath acts on the adiabatic polarization.

### Non-Secular Accuracy

The CGME retains the non-secular terms (the off-diagonal $W$-matrix elements) that describe line mixing. The secular Lindblad approximation, which discards these terms, corresponds to treating each line as independent — a known failure for closely spaced transitions at elevated pressures. The nonadiabatic CGME provides the correct line-mixing description with the nonadiabatic correction applied.

---

## What Needs to Be Computed

For a comparison against HITRAN data, the inputs are:

1. **The energy levels** $E_J = B \cdot J(J+1)$ for CO$_2$ ($B = 0.3902$ cm$^{-1}$), truncated at $J_{\max}$ large enough to capture the transitions of interest. For the $\nu_3$ R-branch, $J = 0$ to $J = 20$ covers the important lines; the correction is largest for $J = 0$–$3$ and falls below 1% for $J \geq 4$ at 5 atm.

2. **The standard $W$-matrix elements** $W_{JJ'}^{\text{Dirac}}$ from published calculations. Sources: Hartmann et al. (*Collisional Effects on Molecular Spectra*, 2021), Gamache et al. (ab initio broadening and line-mixing parameters), or the HITRAN database itself (first-order line-mixing Y-parameters).

3. **The nonadiabatic correction factors** $\mathcal{C}_{JJ'}$ computed from $V_c = \gamma_L P$ and $\Delta E_J = 2B(J+1)$.

4. **The line-shape parameters** $d_J$ (transition dipole strengths), $\nu_J$ (line positions), and $\rho_J$ (thermal populations) from HITRAN.

5. **Comparison:** Compute the line shape $F(\nu)$ using $\mathbf{W}^{\text{Dirac}}$ and $\mathbf{W}^{\text{nad}}$, and compare both against the measured HITRAN absorption cross-sections at elevated pressures.

The computation is a matrix inversion at each frequency point — a $J_{\max} \times J_{\max}$ complex matrix — which is trivial for $J_{\max} \leq 30$. The entire calculation runs on a laptop.

---

## Note on $m_J$ Degeneracy

The treatment above works in the $|J\rangle$ basis without $m_J$ degeneracy. This is appropriate because the line-mixing $W$-matrix is conventionally defined in the transition basis (each transition $J \to J+1$ is one index), not the state basis, and the $m_J$-degeneracy is handled through the thermal population weights $\rho_J \propto (2J+1) \exp(-E_J/k_BT)$ and the transition dipole strengths $d_J \propto (J+1)/(2J+1)$ (for R-branch). The master equation acts in the transition space, and the $2J+1$ degeneracy enters through the initial conditions and the observable, not through the dimension of the density matrix.

---

## Note on CO$_2$ Symmetry

CO$_2$ is a symmetric linear molecule with no permanent dipole moment. The $\nu_3$ antisymmetric stretch has a transition dipole and produces R- and P-branch absorption ($\Delta J = \pm 1$), but no Q-branch ($\Delta J = 0$). Additionally, the nuclear spin statistics of $^{16}$O ($I = 0$, bosonic) require that only even-$J$ or only odd-$J$ levels are populated for a given vibrational state. This means the $W$-matrix for the $\nu_3$ R-branch couples only even-$J$ transitions to even-$J$ transitions, and odd-$J$ to odd-$J$. The nonadiabatic correction applies within each symmetry block independently.

For the minimal 3-level system ($J = 0, 1, 2$), this symmetry constraint means that $J = 0$ and $J = 2$ are in the same symmetry block (both even), while $J = 1$ is in a separate block. The physically meaningful line-mixing coupling is between even-$J$ R-branch transitions: $R(0)$ coupled to $R(2)$, $R(2)$ coupled to $R(4)$, etc. The minimal system for a genuine line-mixing comparison is therefore $J = 0, 2, 4$ (three even-$J$ levels), not $J = 0, 1, 2$.

This does not change the structure of the master equation — only the specific energy levels and Bohr frequencies that enter. The master equation for $J = 0, 2, 4$ has the same form as above, with:

$$E_0 = 0, \quad E_2 = 6B, \quad E_4 = 20B$$

$$\omega_{20} = 6B \cdot C_{\text{conv}}, \quad \omega_{40} = 20B \cdot C_{\text{conv}}, \quad \omega_{42} = 14B \cdot C_{\text{conv}}$$

and the coupling is between $R(0)$ ($J=0 \to J=1$) and $R(2)$ ($J=2 \to J=3$), with a frequency separation of $4B = 1.561$ cm$^{-1}$.

---

## Comparison with the Three-Level CGME

The structure is identical to the generic three-level nonadiabatic CGME. The differences are:

| Feature | Generic 3-level | Rotational (HITRAN) |
|:---|:---|:---|
| Energy levels | Arbitrary $E_0, E_1, E_2$ | $E_J = BJ(J+1)$, specific spacing |
| Perturbation | Generic $V(t)$ | Collisional coupling $V_c = \gamma_L P$ |
| Bath | Generic spectral density $J(\omega)$ | Collisional Lorentzian, $\tau_c \sim 1$ ps |
| Coupling structure | Arbitrary $S_{kl}$ | Nearest-neighbor $\Delta J = \pm 1$ |
| Sinc factors | System-dependent | Near unity for CO$_2$ (non-secular terms survive) |
| Observable | Time-domain populations | Frequency-domain line shape $F(\nu)$ |
| $V/\Delta E$ | System-dependent | $\gamma_L P / 2B(J+1)$, tunable via pressure |

The observable is the key difference. The NMR experiment extracted $R_{1\rho}$ from exponential decay in the time domain. The HITRAN experiment extracts line-mixing parameters from the absorption line shape in the frequency domain. The master equation is the same; the connection to experiment goes through the resolvent formula rather than through exponential fitting.

---

## What This Document Does Not Cover

**Realistic $J$-dependent coupling elements.** The Phase 0 estimate uses $J$-independent $V_c = \gamma_L P$. A full comparison requires $J$-dependent off-diagonal $W$-matrix elements from published intermolecular potential calculations (Gamache, Hartmann). The structure of the master equation is unchanged; only the numerical values of the matrix elements change.

**Centrifugal distortion.** The rigid-rotor energy levels $E_J = BJ(J+1)$ do not include the centrifugal distortion correction $-DJ^2(J+1)^2$. For CO$_2$, $D \approx 1.3 \times 10^{-7}$ cm$^{-1}$, which is negligible for $J \leq 20$.

**Vibrational dependence.** The rotational constant $B$ depends slightly on the vibrational state. For the $\nu_3$ band, $B'' = 0.3902$ cm$^{-1}$ (ground state) and $B' \approx 0.3866$ cm$^{-1}$ (excited state). This shifts the line positions but does not affect the structure of the $W$-matrix.

**Higher-order nonadiabatic corrections.** At $V_c/\Delta E = 0.45$ (CO$_2$ at 5 atm, $J = 0$), first-order perturbation theory is being stretched. A second-order treatment would introduce additional correction terms. The first-order result provides the leading correction and establishes whether the effect is detectable; refinement to higher order is warranted only if the first-order comparison with data is promising.

**Non-Markovian effects.** The CGME operates within the Born-Markov regime. For very high pressures where the collisional perturbation is comparable to or larger than the energy spacing (the "impact limit" breaking down), non-Markovian effects become important. This corresponds to $V_c > 2B$, or $P > 2B/\gamma_L \approx 11$ atm for CO$_2$. The useful operating regime is below this threshold.

---

## References

* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields: Dephasing and population relaxation due to contact with a bath. *J. Chem. Phys.*, 158, 164107.
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
* Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Physical Review A*, 78, 022106.
* Hartmann, J.-M., Boulet, C., & Robert, D. (2021). *Collisional Effects on Molecular Spectra* (2nd ed.). Elsevier.
* Gordon, I. E. et al. (2022). The HITRAN2020 molecular spectroscopic database. *JQSRT*, 277, 107949.
* Anderson, P. W. (1949). Pressure broadening in the microwave and infra-red regions. *Physical Review*, 76, 647.
* Rosenkranz, P. W. (1975). Shape of the 5 mm oxygen band in the atmosphere. *IEEE Trans. Antennas Propag.*, 23, 498–506.
* Strow, L. L., & Reuter, D. (1988). Effect of line mixing on atmospheric brightness temperatures near 15 µm. *Applied Optics*, 27, 872–878.