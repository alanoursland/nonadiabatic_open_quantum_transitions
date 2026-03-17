# The Landau-Lifshitz Separation of Excited-State Coefficients

## The Key Idea

In Dirac's time-dependent perturbation theory, the coefficient $c_k(t)$ of an excited state $|k\rangle$ in the evolving wavefunction measures the overlap between the time-evolved state and the unperturbed eigenstate. As discussed in the Dirac transition probabilities notes, the standard interpretation treats $|c_k(t)|^2$ as the probability that the system has transitioned to state $|k\rangle$. But this interpretation conflates two physically distinct processes: the reversible adjustment of the ground state to the perturbation (without any real excitation), and the irreversible excitation of the system into a genuinely different quantum state.

Landau and Lifshitz, in their *Quantum Mechanics: Non-Relativistic Theory* textbook, introduced a method for separating these two contributions. By performing an integration by parts on Dirac's first-order expression for $c_k^{(1)}(t)$, they decomposed it into an **adiabatic term** $a_k(t)$ and a **nonadiabatic term** $b_k(t)$:

$$c_k^{(1)}(t) = a_k^{(1)}(t) + b_k^{(1)}(t)$$

Landau and Lifshitz stated that the true transition probability into the excited state $|k\rangle$ is given by $|b_k(t)|^2$, not by $|c_k(t)|^2$. This is the **Landau-Lifshitz separation**, and it is the mathematical starting point for a line of research — developed extensively by Mandal, Hunt, and collaborators — that has reexamined and extended the foundations of quantum transition probability theory.

---

## The Integration by Parts

To see how the separation works, start from Dirac's first-order result. For a system initially in the ground state $|0\rangle$, the first-order coefficient of an excited state $|k\rangle$ is:

$$c_k^{(1)}(t) = \frac{1}{i\hbar} \int_{-\infty}^{t} \langle k | H'(t') | 0 \rangle \, e^{i\omega_{k0} t'} \, dt'$$

where $H'(t)$ is the time-dependent perturbation, $\omega_{k0} = (E_k - E_0)/\hbar$ is the Bohr frequency between the excited state and the ground state, and the lower limit reflects the assumption that the perturbation was absent in the distant past.

The key step is to integrate by parts, differentiating the matrix element $\langle k | H'(t') | 0 \rangle$ and integrating the oscillatory exponential $e^{i\omega_{k0} t'}$. This yields two terms:

The **boundary term** evaluates the integrand at the current time $t$. It is proportional to $\langle k | H'(t) | 0 \rangle / (E_k - E_0)$ — the matrix element of the perturbation at the present instant, divided by the energy gap. This is precisely what you would get from static (time-independent) perturbation theory applied to the instantaneous Hamiltonian $H_0 + H'(t)$. It represents the ground state's adjustment to the perturbation at time $t$: the system's charge distribution has distorted to partially accommodate the field, incorporating some amplitude of the excited state $|k\rangle$ into the instantaneous ground state, but no real transition has occurred. This is the adiabatic term $a_k^{(1)}(t)$.

The **integral term** depends on the time derivative $\partial H'(t') / \partial t'$ evaluated at all earlier times $t' \leq t$. It vanishes if the perturbation is constant (because $\partial H'/\partial t = 0$), and it is sensitive to how rapidly the perturbation changes. This is the nonadiabatic term $b_k^{(1)}(t)$. It represents the actual excitation of the system — the component that would survive even if you could perfectly subtract out the instantaneous polarization response.

The physical logic is clean: the adiabatic term tracks the system's reversible response to whatever the field is doing right now; the nonadiabatic term accumulates the irreversible effects of the field's history of change. A perturbation that varies slowly produces large adiabatic coefficients (strong polarization) but small nonadiabatic coefficients (few transitions). A perturbation that varies rapidly produces the opposite.

---

## What the Adiabatic Term Contains

The adiabatic coefficient $a_k(t)$ follows directly from the **adiabatic theorem of Born and Fock**. The theorem states that a system initially in an eigenstate of a slowly varying Hamiltonian will remain in the corresponding instantaneous eigenstate at all later times, acquiring only a phase. The "corresponding instantaneous eigenstate" of $H_0 + H'(t)$ is not the same as the original eigenstate $|0\rangle$ of $H_0$ — it is the ground state of the full instantaneous Hamiltonian, which includes admixtures of the excited states $|k\rangle$ of $H_0$. Those admixtures are exactly what $a_k(t)$ measures.

In other words, when you expand the instantaneous ground state $|0'(t)\rangle$ of $H_0 + H'(t)$ in the unperturbed basis $\{|k\rangle\}$, the coefficients of the excited states are the $a_k(t)$. These are nonzero — the instantaneous ground state is a superposition in the old basis — but they do not represent transitions. They represent the fact that the ground state has changed character in response to the field. If the field were suddenly removed, these admixtures would vanish as the system relaxes back to the unperturbed ground state, without any energy having been deposited into the excited states.

This is the quantum mechanical version of polarization. An atom in an external electric field develops an induced dipole moment because the field mixes states of different parity into the ground state. The atom has not been excited; it has been polarized. Dirac's $|c_k(t)|^2$ counts this polarization as part of the transition probability. The Landau-Lifshitz separation removes it.

---

## What the Nonadiabatic Term Contains

The nonadiabatic coefficient $b_k(t)$ captures genuine transitions — changes in the system's state that are not reversed when the perturbation is removed. Its structure reveals when and why transitions happen:

The dependence on $\partial H'/\partial t'$ means that transitions are driven by the *rate of change* of the perturbation, not by its instantaneous magnitude. A strong but constant field polarizes the system without causing transitions; it is the turning on, turning off, or modulation of the field that produces real excitation. This is the content of the adiabatic theorem stated in reverse: departures from adiabatic following occur precisely when the Hamiltonian changes faster than the system can track.

The nonadiabatic term is obtained as a time integral over the entire history of $\partial H'/\partial t'$, weighted by the oscillatory factor $e^{i\omega_{k0} t'}$. The oscillatory weighting means that the dominant contributions come from times when the rate of change of the perturbation has Fourier components near the resonance frequency $\omega_{k0}$. This connects the Landau-Lifshitz framework back to the familiar resonance condition of Fermi's Golden Rule, but without the long-time or constant-perturbation approximations that Fermi's Golden Rule requires.

Crucially, $|b_k(t)|^2$ vanishes when the perturbation is turned on and off smoothly and the system returns to its original Hamiltonian $H_0$, provided the process is adiabatic throughout. Dirac's $|c_k(t)|^2$ does not generally vanish in this scenario — it can retain a nonzero value from the transient dynamics of the pulse, even if no net energy was deposited. This is one of the clearest demonstrations that $|b_k(t)|^2$, not $|c_k(t)|^2$, is the physically meaningful transition probability.

---

## Extension Beyond First Order

Landau and Lifshitz presented the separation at first order in the perturbation. A natural question is whether it survives at higher orders, where the algebra becomes considerably more involved and the adiabatic and nonadiabatic contributions might mix.

Mandal and Hunt (2012) proved that the separation is complete to all orders: the energy of the system at any time $t$ decomposes exactly into an adiabatic part and a nonadiabatic part, with no cross-terms. The adiabatic part is identical to the result of static perturbation theory applied to the instantaneous Hamiltonian — it is the energy the system would have if it had adjusted perfectly to the current value of the field. The nonadiabatic part is a sum over excited states of the transition probability $|b_k(t)|^2$ multiplied by the transition energy $(E_k - E_0)$. They verified this explicitly through third order in the perturbation and showed that the cross-terms between $a_k(t)$ and $b_k(t)$ vanish when computing expectation values.

This is a nontrivial result. It means the Landau-Lifshitz separation is not an artifact of first-order perturbation theory — it reflects a genuine physical decomposition that persists to arbitrary accuracy. The adiabatic and nonadiabatic contributions to the energy do not interfere with each other. They are, in a meaningful sense, independent channels.

In subsequent work, Mandal and Hunt (2018, 2020) showed that the variance and higher moments of the energy distribution depend entirely on the nonadiabatic transition probabilities $|b_k(t)|^2$. A standard statistical analysis of the energy distribution gives consistent results only if the excitation probability is identified with $|b_k(t)|^2$, not with $|c_k(t)|^2$. This provides an additional, statistics-based argument that the nonadiabatic coefficients are the correct measure of transition probability.

---

## The Connection to Gauge Invariance

One of the more subtle virtues of the Landau-Lifshitz separation is its behavior under gauge transformations.

In electrodynamics, the interaction between a charged system and an electromagnetic field can be written in different gauges — the length gauge (where the interaction is $-\mathbf{d} \cdot \mathbf{E}$, with $\mathbf{d}$ the dipole operator and $\mathbf{E}$ the electric field) or the velocity gauge (where the interaction involves the vector potential $\mathbf{A}$), among others. The physical predictions of the theory must not depend on this choice. For a completed process — a pulse that has come and gone — this invariance is guaranteed: the S-matrix is gauge-invariant. But for the instantaneous coefficients during a pulse, gauge invariance is not automatic.

Dirac's coefficients $c_k(t)$ are gauge-dependent while the field is on. The instantaneous value of $|c_k(t)|^2$ can change depending on whether you work in the length gauge or the velocity gauge, which means it cannot be a physical observable during the pulse. The nonadiabatic coefficients $b_k(t)$, by contrast, are gauge-invariant at all times. Mandal and Hunt demonstrated this explicitly: the same $|b_k(t)|^2$ is obtained regardless of the gauge used to describe the electromagnetic interaction. This is a strong formal argument that $|b_k(t)|^2$ is the physically correct transition probability at every instant, not just after the perturbation has ended.

---

## Why This Matters for the Broader Program

The Landau-Lifshitz separation might appear to be a technical refinement — a cleaner way to bookkeep the same underlying physics. But its consequences propagate far beyond bookkeeping once the system is coupled to an environment.

As discussed in the nonadiabatic bath coupling notes, the choice between $|c_k(t)|^2$ and $|b_k(t)|^2$ as the "transition probability" has hard physical consequences when the system interacts with a thermal bath. A master equation built on Dirac's coefficients treats the adiabatic polarization as a real population, causing the bath to relax it, generating fictitious heat flow, and driving the system toward the wrong thermal equilibrium. A master equation built on the nonadiabatic coefficients avoids all of these pathologies.

The Landau-Lifshitz separation is what makes that correction possible. It identifies which part of $c_k(t)$ is "real" (in the sense of representing an actual change in the system's quantum state) and which part is "virtual" (in the sense of representing a reversible response to the instantaneous field). Without this decomposition, there is no principled way to couple the driven dynamics to the bath. With it, the coupling becomes thermodynamically consistent.

This is the sense in which a textbook observation by Landau and Lifshitz — an integration by parts and a remark about the correct identification of transition probabilities — has turned out to be the foundation of a substantial modern research program. The original insight was brief and, by Landau and Lifshitz's standards, almost offhand. Its full implications for open quantum systems, gauge invariance, quantum thermodynamics, and the limitations of Fermi's Golden Rule have only been developed in the past two decades, primarily through the work of Mandal, Hunt, Jovanovski, and their collaborators.

---

## References

* Landau, L. D., & Lifshitz, E. M. (1977). *Quantum Mechanics: Non-Relativistic Theory* (3rd ed.). Pergamon Press.
* Mandal, A., & Hunt, K. L. C. (2012). Adiabatic and nonadiabatic contributions to the energy of a system subject to a time-dependent perturbation: Complete separation and physical interpretation. *J. Chem. Phys.*, 137, 164109. https://doi.org/10.1063/1.4750045
* Mandal, A., & Hunt, K. L. C. (2018). Quantum transition probabilities during a perturbing pulse: Differences between the nonadiabatic results and Fermi's golden rule forms. *J. Chem. Phys.*, 148, 194107. https://doi.org/10.1063/1.5019172
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse: Toward experimental tests of the differences from Dirac's transition probabilities. *J. Chem. Phys.*, 149, 204110. https://doi.org/10.1063/1.5054313
* Mandal, A., & Hunt, K. L. C. (2020). Variance of the energy of a quantum system in a time-dependent perturbation: Determination by nonadiabatic transition probabilities. *J. Chem. Phys.*, 152, 104110. https://doi.org/10.1063/1.5140009
* Mandal, A., & Hunt, K. L. C. (2021). Quantum transition probabilities due to overlapping electromagnetic pulses: Persistent differences between Dirac's form and nonadiabatic perturbation theory. *J. Chem. Phys.*, 154, 024116. https://doi.org/10.1063/5.0020169