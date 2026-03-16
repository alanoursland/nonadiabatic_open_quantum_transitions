In quantum mechanics, **Dirac's standard approach** to transition probabilities is the foundation of **Time-Dependent Perturbation Theory (TDPT)**. It provides a framework for calculating how a system transitions between the eigenstates of an unperturbed Hamiltonian $H_0$ when subjected to a time-varying potential $V(t)$.

### **The Standard Approach**

Dirac's method, developed in the late 1920s, treats the total Hamiltonian as $H = H_0 + V(t)$. The goal is to solve the time-dependent Schrödinger equation:
$$i\hbar \frac{\partial}{\partial t} |\Psi(t)\rangle = (H_0 + V(t)) |\Psi(t)\rangle$$

1.  **Interaction Picture:** The state is expanded in the basis of the unperturbed eigenstates $|n\rangle$:
    $|\Psi(t)\rangle = \sum_n c_n(t) e^{-i E_n t / \hbar} |n\rangle$
2.  **Probability Interpretation:** Dirac defined the transition probability from an initial state $|i\rangle$ to a final state $|f\rangle$ as the squared modulus of the expansion coefficient $|c_f(t)|^2$ (Haque & Zhang, 2014; Quantum transition probabilities..., 2021).
3.  **Fermi's Golden Rule:** When the perturbation is periodic or constant and the final states form a continuum, this approach leads to **Fermi's Golden Rule**, which calculates the transition *rate* rather than just the probability (Jang & Rhee, 2023).

-----

### **Core Assumptions**

Dirac’s approach relies on several critical simplifications to remain mathematically tractable:

  * **Weak Perturbation:** The most fundamental assumption is that the perturbation $V(t)$ is small compared to $H_0$. This allows the use of a power-series expansion where higher-order terms are neglected.
  * **Initial State Purity:** The system is usually assumed to start in a single, well-defined eigenstate $|i\rangle$ of $H_0$.
  * **Adiabatic/Sudden Switching:** The perturbation is often assumed to be "turned on" either infinitely slowly (adiabatically) or instantly at $t=0$ to avoid complex transient effects.
  * **Coherence:** The standard derivation assumes the system remains perfectly coherent. It does not account for interaction with a "bath" or environment, which would cause dephasing (Nonadiabatic transition..., 2023).

-----

### **Where the Approach Breaks Down**

While highly successful (e.g., for light-matter interaction in weak fields), Dirac's method fails in several physical regimes:

#### **1. Strong Field Regime**

When the applied field is strong (e.g., high-intensity lasers), the perturbation $V(t)$ is no longer "small." In these cases, the transition probability can exceed 1 if calculated using only first-order theory, which is physically impossible. This requires non-perturbative methods like the **Strong-Field Approximation (SFA)** or direct numerical solutions of the Dirac equation (Reiss, 1990; Telnov et al., 2018).

#### **2. Long Time Scales & Dephasing**

Dirac’s coefficients $|c_k(t)|^2$ often oscillate indefinitely while a field is constant, which can be unphysical in real systems. Furthermore, Dirac's approach cannot account for **dephasing**—the loss of quantum interference due to environmental coupling. In these cases, the probability calculated by Dirac can differ from experimental results by over 35% (Quantum transition probabilities..., 2021).

#### **3. Non-Hermitian and Dissipative Systems**

The standard theory is built for Hermitian Hamiltonians where probability is conserved. In systems with gain or loss (non-Hermitian), the standard approach fails to describe the inherent asymmetry in forward versus reverse transition probabilities (Choi, 2020; Non-Hermitian Dynamics..., 2026).

#### **4. Level Crossings (Dirac Points)**

In materials like graphene, where energy bands cross (Dirac points), the adiabatic assumptions break down. Transitions at these points are better described by the **Landau-Zener** formula rather than standard first-order perturbation theory (Faraj & Jin, 2015).

-----

### **References**

  * Choi, J. R. (2020). Perturbation Theory for Time-Dependent Quantum Systems Involving Complex Potentials. *Frontiers in Physics*, *8*(189). [https://doi.org/10.3389/fphy.2020.00189](https://doi.org/10.3389/fphy.2020.00189)
    Cited by: 14
  * Faraj, A., & Jin, S. (2015). *The Landau-Zener transition and the surface hopping method for the 2D Dirac equation for graphene*. arXiv. [https://doi.org/10.48550/arxiv.1505.05988](https://www.google.com/search?q=https://doi.org/10.48550/arxiv.1505.05988)
    Cited by: 4
  * Haque, M., & Zhang, J. M. (2014). Nonsmooth and level-resolved dynamics illustrated with a periodically driven tight binding model. *ScienceOpen Research*. [https://doi.org/10.14293/s2199-1006.1.sor-phys.a2cem4.v1](https://www.google.com/search?q=https://doi.org/10.14293/s2199-1006.1.sor-phys.a2cem4.v1)
    Cited by: 5
  * Jang, S. J., & Rhee, Y. M. (2023). *Modified Fermi's golden rule rate expressions*. arXiv. [https://doi.org/10.1063/5.0152804](https://doi.org/10.1063/5.0152804)
    Cited by: 8
  * Mendoza, M., Araújo, N. A. M., Succi, S., & Herrmann, H. J. (2012). Transition in the Equilibrium Distribution Function of Relativistic Particles. *Scientific Reports*, *2*. [https://doi.org/10.1038/srep00611](https://doi.org/10.1038/srep00611)
    Cited by: 23
  * Non-Hermitian Dynamics in Three-Level Systems: A Perturbative Approach for Time-Dependent Hamiltonians. (2026). *Entropy*, *28*(3), 268. [https://doi.org/10.3390/e28030268](https://www.google.com/search?q=https://doi.org/10.3390/e28030268)
  * Nonadiabatic transition probabilities for quantum systems in electromagnetic fields. (2023). *NSF PAR*.
  * Quantum transition probabilities due to overlapping electromagnetic pulses: Persistent differences between Dirac's form and no. (2021). *AIP Publishing*. [https://doi.org/10.1063/5.0020169](https://www.google.com/search?q=https://doi.org/10.1063/5.0020169)
  * Reiss, H. R. (1990). Relativistic strong-field photoionization. *Journal of the Optical Society of America B*, *7*(4), 574. [https://doi.org/10.1364/josab.7.000574](https://www.google.com/search?q=https://doi.org/10.1364/josab.7.000574)
    Cited by: 225
  * Telnov, D. A., Krapivin, D. A., Heslar, J., & Chu, S.-I. (2018). Multiphoton Ionization of One-Electron Relativistic Diatomic Quasimolecules in Strong Laser Fields. *The Journal of Physical Chemistry A*, *122*, 8026-8036. [https://doi.org/10.1021/acs.jpca.8b07463](https://doi.org/10.1021/acs.jpca.8b07463)
    Cited by: 17

Would you like me to explain how the Landau-Zener formula specifically addresses the breakdown at level crossings?