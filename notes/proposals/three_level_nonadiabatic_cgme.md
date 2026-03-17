# The Nonadiabatic CGME for a 3-Level System

## Purpose

The nonadiabatic Nakajima-Zwanzig proposal describes a general framework for driven, dissipative quantum systems. This document reduces that framework to its simplest nontrivial case: a 3-level system driven by a classical electromagnetic field and weakly coupled to a thermal bath. This is the regime Hunt's group works in, and the goal is a concrete, implementable master equation that improves on what Jovanovski, Mandal, and Hunt (2023) achieved with nonadiabatic Redfield theory — specifically, by adding positivity preservation without sacrificing the nonadiabatic decomposition's gauge invariance and thermodynamic consistency.

At this scale, the BQP structure of the general framework is irrelevant. The density matrix is 3×3. The memory kernel can be computed analytically or by direct numerical integration. The entire calculation runs on a laptop. What matters here is not computational complexity but *physical correctness*: does the master equation satisfy all four criteria (complete positivity, non-secular accuracy, gauge invariance, thermodynamic consistency) for a driven, dissipative 3-level system?

---

## The System

Consider three energy eigenstates of the unperturbed Hamiltonian $H_0$:

$$H_0 |n\rangle = E_n |n\rangle, \quad n \in \{0, 1, 2\}$$

with $E_0 < E_1 < E_2$. The Bohr frequencies are:

$$\omega_{10} = \frac{E_1 - E_0}{\hbar}, \quad \omega_{20} = \frac{E_2 - E_0}{\hbar}, \quad \omega_{21} = \frac{E_2 - E_1}{\hbar}$$

Note that $\omega_{20} = \omega_{10} + \omega_{21}$. The system is driven by a time-dependent perturbation $V(t)$ (a laser pulse, a microwave field, or any classical electromagnetic interaction) and coupled to a thermal bath at temperature $T$ through a system-bath interaction $H_{SB} = S \otimes B$, where $S$ is a system operator and $B$ is a bath operator.

---

## Step 1: The Nonadiabatic Coefficients

The time-evolving state is expanded in the unperturbed basis:

$$|\Psi(t)\rangle = \sum_{n=0}^{2} c_n(t) \, e^{-iE_n t/\hbar} |n\rangle$$

Following the Landau-Lifshitz integration by parts (described in the Landau-Lifshitz separation notes), each excited-state coefficient is decomposed:

$$c_n^{(1)}(t) = a_n^{(1)}(t) + b_n^{(1)}(t), \quad n \in \{1, 2\}$$

**The adiabatic term** $a_n^{(1)}(t)$ is the boundary term from the integration by parts:

$$a_n^{(1)}(t) = \frac{\langle n | V(t) | 0 \rangle}{E_n - E_0} \, e^{i\omega_{n0} t}$$

This is identical to what static perturbation theory gives for the instantaneous perturbation $V(t)$. It represents the polarization of the ground state — the admixture of $|n\rangle$ into the instantaneous ground state due to the field — without any real transition.

**The nonadiabatic term** $b_n^{(1)}(t)$ is the residual integral:

$$b_n^{(1)}(t) = -\frac{1}{i\hbar} \frac{1}{i\omega_{n0}} \int_{-\infty}^{t} \frac{d}{dt'}\langle n | V(t') | 0 \rangle \, e^{i\omega_{n0} t'} \, dt'$$

This depends on the time derivative of the perturbation matrix element. It represents genuine transitions driven by the field's time variation. It vanishes when the field is constant, is gauge-invariant, and gives the correct power absorption.

For the 3-level system, we have two nonadiabatic amplitudes: $b_1(t)$ (transition from ground to first excited state) and $b_2(t)$ (transition from ground to second excited state). At first order, transitions between the two excited states are mediated through the ground state and appear at second order.

**These coefficients are what Hunt's group already computes.** Nothing new is required here.

---

## Step 2: The Nonadiabatic Density Matrix

The nonadiabatic reduced density matrix is the 3×3 matrix whose elements are constructed from the nonadiabatic coefficients:

$$[\sigma_{\text{nad}}]_{00} = 1 - |b_1(t)|^2 - |b_2(t)|^2$$
$$[\sigma_{\text{nad}}]_{11} = |b_1(t)|^2$$
$$[\sigma_{\text{nad}}]_{22} = |b_2(t)|^2$$
$$[\sigma_{\text{nad}}]_{10} = b_1(t) \cdot f_0^*(t), \quad [\sigma_{\text{nad}}]_{20} = b_2(t) \cdot f_0^*(t)$$
$$[\sigma_{\text{nad}}]_{21} = b_2(t) \, b_1^*(t)$$

where $f_0(t)$ represents the ground-state amplitude corrected for normalization. The off-diagonal elements $[\sigma_{\text{nad}}]_{n0}$ are the coherences between the genuinely excited states and the ground state. The element $[\sigma_{\text{nad}}]_{21}$ is the coherence between the two excited states.

The critical difference from the standard density matrix: the populations $[\sigma_{\text{nad}}]_{nn}$ are $|b_n(t)|^2$, not $|c_n(t)|^2$. The adiabatic polarization has been removed. When the field is constant, these populations are constant — there is nothing for the bath to oscillate against.

---

## Step 3: The Nonadiabatic Redfield Tensor

The Redfield relaxation tensor describes how the bath drives transitions between the system's states. In the nonadiabatic framework, the tensor elements are computed in the standard way — from the bath correlation functions evaluated at the system's Bohr frequencies — but they act on the nonadiabatic density matrix $\sigma_{\text{nad}}$ rather than the Dirac density matrix.

The bath is characterized by its spectral density $J(\omega)$ and the thermal occupation number $\bar{n}(\omega) = (e^{\hbar\omega/k_BT} - 1)^{-1}$. The Redfield tensor elements relevant for the 3-level system are built from the one-sided Fourier transforms of the bath correlation functions:

$$\Gamma(\omega) = \int_0^\infty d\tau \; e^{i\omega\tau} \; \langle B(\tau) B(0) \rangle_{\text{eq}}$$

The real part of $\Gamma(\omega)$ gives the relaxation rates; the imaginary part gives the Lamb shift corrections.

For each pair of states $(k, l)$ with $k > l$, the transition rates are:

$$\gamma_{k \to l} = \frac{2}{\hbar^2} |S_{lk}|^2 \, \text{Re}[\Gamma(\omega_{kl})] = \frac{2\pi}{\hbar^2} |S_{lk}|^2 \, J(\omega_{kl}) \, [\bar{n}(\omega_{kl}) + 1]$$

$$\gamma_{l \to k} = \frac{2\pi}{\hbar^2} |S_{lk}|^2 \, J(\omega_{kl}) \, \bar{n}(\omega_{kl})$$

where $S_{lk} = \langle l | S | k \rangle$ are the matrix elements of the system coupling operator. These rates satisfy detailed balance:

$$\frac{\gamma_{l \to k}}{\gamma_{k \to l}} = e^{-\hbar\omega_{kl}/k_BT}$$

which ensures the correct Boltzmann thermal equilibrium.

For the 3-level system, there are three transition pairs: $(1,0)$, $(2,0)$, and $(2,1)$, giving six rates total (three downward, three upward).

**This is also essentially what the 2023 Jovanovski-Mandal-Hunt paper did.** The new step is next.

---

## Step 4: The Coarse-Graining Procedure

The nonadiabatic Redfield tensor, like any Redfield tensor, can violate the positivity of the density matrix. The coarse-grained master equation (CGME) fixes this by averaging the Redfield tensor over a finite time window $\Delta\tau$.

### The Mechanism

In the Redfield equation, different density matrix elements are coupled by terms that oscillate at frequency differences $|\omega_{ab} - \omega_{cd}|$. In the standard (secular) Lindblad equation, all such oscillating terms are dropped — equivalent to averaging over $\Delta\tau \to \infty$. In the CGME, the averaging is over a finite $\Delta\tau$, which retains terms whose oscillation period is longer than $\Delta\tau$ and suppresses those that oscillate faster.

Concretely, each element of the Redfield tensor $R_{abcd}$ is multiplied by a damping factor:

$$R_{abcd}^{\text{CG}} = R_{abcd} \cdot \text{sinc}\left(\frac{(\omega_{ab} - \omega_{cd})\Delta\tau}{2}\right)$$

where $\text{sinc}(x) = \sin(x)/x$. When $|\omega_{ab} - \omega_{cd}| \ll 1/\Delta\tau$, the sinc factor is approximately 1 and the term is retained. When $|\omega_{ab} - \omega_{cd}| \gg 1/\Delta\tau$, the sinc factor is approximately 0 and the term is suppressed.

The coarse-grained tensor $R_{abcd}^{\text{CG}}$ generates a master equation in Lindblad form, guaranteeing complete positivity.

### Application to the 3-Level System

For the 3-level system, the relevant frequency differences between transitions are:

$$|\omega_{10} - \omega_{20}| = \omega_{21}, \quad |\omega_{10} - \omega_{21}| = \omega_{20} - 2\omega_{21}, \quad |\omega_{20} - \omega_{21}| = \omega_{10}$$

**Case 1: Well-separated levels** ($\omega_{10}$, $\omega_{20}$, $\omega_{21}$ all much larger than the relaxation rates). Choose $\Delta\tau$ such that $1/\Delta\tau$ is smaller than all three frequency differences but larger than the relaxation rates. All non-secular terms are suppressed, and the CGME reduces to the nonadiabatic secular Lindblad equation:

$$\frac{d\sigma_{\text{nad}}}{dt} = -\frac{i}{\hbar}[H_{\text{eff}}(t), \sigma_{\text{nad}}] + \sum_{(k,l)} \gamma_{k \to l} \, \mathcal{D}[|l\rangle\langle k|] \, \sigma_{\text{nad}} + \gamma_{l \to k} \, \mathcal{D}[|k\rangle\langle l|] \, \sigma_{\text{nad}}$$

where the sum runs over all three pairs $(1,0)$, $(2,0)$, $(2,1)$, and $\mathcal{D}[L]\rho = L\rho L^\dagger - \frac{1}{2}\{L^\dagger L, \rho\}$ is the standard Lindblad dissipator. This is the simplest version — six Lindblad operators with rates determined by the spectral density and temperature, acting on the nonadiabatic density matrix.

**Case 2: Two levels nearly degenerate** ($\omega_{21} \ll \omega_{10}, \omega_{20}$, meaning states $|1\rangle$ and $|2\rangle$ are close in energy). Now $|\omega_{10} - \omega_{20}| = \omega_{21}$ is small, and the non-secular terms coupling the $0 \leftrightarrow 1$ and $0 \leftrightarrow 2$ transitions are retained by the coarse-graining. The CGME includes coherence transfer between these two transitions — the bath-mediated coupling that the secular approximation would miss. This is where the CGME provides genuine improvement over the secular Lindblad equation, and it matters for molecular systems with closely spaced electronic or vibrational states.

### Choosing $\Delta\tau$

The coarse-graining timescale must satisfy:

$$\tau_B \ll \Delta\tau \ll T_1, T_2$$

where $\tau_B$ is the bath correlation time (how quickly the bath "forgets"), and $T_1$, $T_2$ are the system's relaxation and dephasing times. This window exists whenever the bath is fast compared to the system's dissipative dynamics — the same regime where the Born-Markov approximation is valid. The CGME does not extend the validity beyond Born-Markov; it fixes the positivity problem *within* that regime.

For Hunt's systems (atomic/molecular transitions coupled to thermal radiation or a solvent bath), this hierarchy of timescales is typically well satisfied.

---

## Step 5: The Complete Nonadiabatic CGME

Assembling the pieces, the master equation for the 3-level system is:

$$\frac{d\sigma_{\text{nad}}}{dt} = -\frac{i}{\hbar}[H_S + H_{LS}(t) + V_{\text{nad}}(t), \, \sigma_{\text{nad}}] + \mathcal{L}_{\text{CG}}[\sigma_{\text{nad}}]$$

where:

**$H_S = H_0$** is the unperturbed system Hamiltonian.

**$H_{LS}(t)$** is the Lamb shift Hamiltonian — the small energy-level shifts induced by the bath, extracted from the imaginary parts of $\Gamma(\omega)$. This is typically a minor correction but is included for completeness.

**$V_{\text{nad}}(t)$** represents the coherent nonadiabatic dynamics — the driving of genuine transitions by the time-varying field. This is not the full perturbation $V(t)$; it is the part of the dynamics associated with the nonadiabatic coefficients. The adiabatic part (the polarization response) has been absorbed into the definition of the nonadiabatic subspace and does not appear in the master equation.

**$\mathcal{L}_{\text{CG}}$** is the coarse-grained dissipator — the Lindblad-form dissipative term constructed from the coarse-grained Redfield tensor. It acts on the nonadiabatic density matrix elements, relaxing the nonadiabatic populations toward the thermal equilibrium of the driven system and dephasing the nonadiabatic coherences.

---

## What This Equation Satisfies

### Complete Positivity

Guaranteed by the Lindblad form of the coarse-grained dissipator. The density matrix $\sigma_{\text{nad}}(t)$ will never develop negative eigenvalues, regardless of the driving protocol, the temperature, or the coupling strength (within the Born-Markov validity regime).

### Gauge Invariance

Inherited from the nonadiabatic coefficients $b_k(t)$. The density matrix elements, and therefore all physical predictions, are independent of whether the electromagnetic interaction is written in the length gauge or the velocity gauge. This was proven by Mandal and Hunt (2016) for the nonadiabatic coefficients; the CGME inherits it because the dissipator acts on these gauge-invariant quantities.

### Thermodynamic Consistency

When the field is constant ($\partial V / \partial t = 0$), the nonadiabatic populations $|b_k|^2$ are constant and $V_{\text{nad}} = 0$. The dissipator drives the system toward the Boltzmann distribution at the bath temperature. The heat flow is zero. The first law is satisfied.

When the field is varying, the nonadiabatic populations change (driven by $V_{\text{nad}}$), and the dissipator relaxes them toward the instantaneous thermal target. The power balance — work done by the field equals the rate of change of system energy plus heat dissipated — holds because the energy separation (adiabatic energy + nonadiabatic energy, with no cross-terms) is built into the density matrix.

### Non-Secular Accuracy

Controlled by the coarse-graining timescale $\Delta\tau$. For well-separated levels, the CGME reduces to the secular Lindblad equation (no non-secular terms needed). For nearly degenerate levels, the CGME retains the coherence-transfer terms that the secular approximation would discard. The level of non-secular accuracy is adjustable through $\Delta\tau$ and is always at least as good as the secular Lindblad equation.

---

## What Needs to Be Computed

For a specific 3-level system (e.g., the HCl vibration-rotation system that Mandal and Hunt used for numerical illustrations), the inputs are:

1. **The energy levels** $E_0$, $E_1$, $E_2$ and the transition matrix elements $\langle n | V(t) | m \rangle$ — from the molecular Hamiltonian.

2. **The nonadiabatic coefficients** $b_1(t)$, $b_2(t)$ — computed by the Landau-Lifshitz integration by parts from the perturbation $V(t)$. Already established in the Mandal-Hunt program.

3. **The bath spectral density** $J(\omega)$ and temperature $T$ — model inputs that characterize the environment (Drude-Lorentz for a solvent, Ohmic for a generic thermal bath, etc.).

4. **The system-bath coupling matrix elements** $S_{kl} = \langle k | S | l \rangle$ — from the model of how the system couples to the bath.

5. **The coarse-graining timescale** $\Delta\tau$ — chosen to satisfy $\tau_B \ll \Delta\tau \ll T_1$.

From these inputs, the six relaxation rates, the Lamb shifts, and the coarse-grained dissipator are fully determined. The master equation can be integrated numerically using standard ODE methods (it's a system of 9 coupled ODEs for the real and imaginary parts of the density matrix elements).

---

## Comparison with Existing Approaches

For the same 3-level system, the nonadiabatic CGME can be compared against:

**Standard Redfield** (what most textbooks use): Uses $|c_k(t)|^2$ as populations, includes non-secular terms, may violate positivity. Predicts oscillating populations and nonzero heat flow during a plateau pulse. Gauge-dependent for driven systems.

**Standard secular Lindblad** (what most quantum optics codes use): Uses $|c_k(t)|^2$ as populations, discards non-secular terms, guarantees positivity. Predicts oscillating populations during a plateau pulse (from the Dirac coefficients), no coherence transfer between nearly degenerate transitions. Gauge-dependent for driven systems.

**Nonadiabatic Redfield** (Jovanovski, Mandal, Hunt 2023): Uses $|b_k(t)|^2$ as populations, includes non-secular terms, may violate positivity. Predicts constant populations during a plateau pulse, zero heat flow, correct thermal equilibrium. Gauge-invariant. This is the current state of the art from Hunt's group.

**Nonadiabatic CGME** (this proposal): Uses $|b_k(t)|^2$ as populations, retains non-secular terms for nearly degenerate transitions, guarantees positivity. Predicts constant populations during a plateau pulse, zero heat flow, correct thermal equilibrium. Gauge-invariant. This is the proposed next step.

The nonadiabatic CGME adds positivity preservation to everything the 2023 paper already achieved, at the cost of introducing the coarse-graining timescale $\Delta\tau$ — a single, physically motivated parameter with a well-defined validity regime.

---

## What This Document Does Not Cover

**Higher-order nonadiabatic corrections.** The treatment here uses the first-order Landau-Lifshitz decomposition. Mandal and Hunt have extended the energy separation to third order (2012) and the variance analysis to second order (2020). Incorporating higher-order nonadiabatic coefficients into the CGME is straightforward in principle but adds algebraic complexity.

**Non-Markovian effects.** The CGME operates within the Born-Markov regime. For baths with long memory (structured spectral densities, strong coupling), the full nonadiabatic GQME with a numerically extracted memory kernel (or transfer tensors) would be needed. That is a more involved calculation described in the general proposal.

**Systems larger than 3 levels.** The framework generalizes to $n$ levels without conceptual change — more Bohr frequencies, more transition pairs, a larger density matrix — but the computational cost grows as $n^4$ (the size of the Redfield tensor). For large $n$, the BQP constraints discussed in the general proposal become relevant.

---

## References

* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields: Dephasing and population relaxation due to contact with a bath. *J. Chem. Phys.*, 158, 164107.
* Mandal, A., & Hunt, K. L. C. (2012). Adiabatic and nonadiabatic contributions to the energy of a system subject to a time-dependent perturbation: Complete separation and physical interpretation. *J. Chem. Phys.*, 137, 164109.
* Mandal, A., & Hunt, K. L. C. (2016). Gauge-invariant expectation values of the energy of a molecule in an electromagnetic field. *J. Chem. Phys.*, 144, 044109.
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
* Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Physical Review A*, 78, 022106.
* Redfield, A. G. (1957). On the theory of relaxation processes. *IBM J. Res. Dev.*, 1, 19–31.