# Dirac's Transition Probabilities: The Standard Approach, Its Assumptions, and Its Limits

## Historical Context

In 1927, Paul Dirac posed a question that would shape the next century of quantum mechanics: given a system in a known initial state, what is the probability that a time-varying perturbation causes it to appear in a different state at a later time? His answer — time-dependent perturbation theory (TDPT) — remains the default framework for computing transition probabilities in quantum mechanics. It is the derivation behind Fermi's Golden Rule, the theoretical backbone of spectroscopy, and the first tool reached for in any textbook treatment of light-matter interaction.

But Dirac's framework is built on specific assumptions, and those assumptions have boundaries. Understanding what the framework actually does — not just the final formulas, but the mechanical steps and the physical content of each approximation — is essential for recognizing when it can be trusted and when it cannot.

---

## The Setup: Splitting the Hamiltonian

The starting point is a system whose total Hamiltonian can be written as two pieces:

$$H = H_0 + V(t)$$

$H_0$ is the **unperturbed Hamiltonian** — the part of the physics you can solve exactly. For an atom in free space, this is the Coulomb Hamiltonian. For a spin in a static magnetic field, it is the Zeeman Hamiltonian. The eigenstates of $H_0$ are the states $|n\rangle$ with energies $E_n$:

$$H_0 |n\rangle = E_n |n\rangle$$

These are the states you know. They form a complete basis, and any state of the system can be expanded in terms of them.

$V(t)$ is the **perturbation** — the time-dependent interaction you are turning on. It might be an oscillating electric field (light), a collision with another particle, or a pulse from a microwave source. The entire framework is built around the idea that $V(t)$ is an addition to a problem you have already solved.

---

## The Interaction Picture

The next step is the key formal move. The full time-dependent Schrödinger equation is:

$$i\hbar \frac{\partial}{\partial t} |\Psi(t)\rangle = (H_0 + V(t)) |\Psi(t)\rangle$$

Since we know the eigenstates of $H_0$, we expand the evolving state in that basis:

$$|\Psi(t)\rangle = \sum_n c_n(t) \, e^{-iE_n t/\hbar} |n\rangle$$

The exponential factors $e^{-iE_n t/\hbar}$ carry the "free" time evolution — the phase accumulation that would happen even without the perturbation. By factoring this out, the coefficients $c_n(t)$ track *only* the changes caused by $V(t)$. If the perturbation is absent, every $c_n(t)$ is constant.

This is the **interaction picture** (sometimes called the Dirac picture, appropriately enough). Its purpose is to isolate the effect of the perturbation from the background dynamics of $H_0$.

Substituting the expansion into the Schrödinger equation and projecting onto a particular final state $|f\rangle$ yields the coupled differential equations for the coefficients:

$$i\hbar \, \dot{c}_f(t) = \sum_n V_{fn}(t) \, e^{i\omega_{fn} t} \, c_n(t)$$

where $V_{fn}(t) = \langle f | V(t) | n \rangle$ are the matrix elements of the perturbation and $\omega_{fn} = (E_f - E_n)/\hbar$ are the Bohr transition frequencies. This is an exact set of equations — no approximations have been made yet. Every $c_n$ is coupled to every other through the matrix elements of $V(t)$, and solving the full coupled system is equivalent to solving the original Schrödinger equation.

---

## The Perturbative Expansion

The system of coupled equations is generally not solvable in closed form. Dirac's approach is to solve it iteratively under the assumption that $V(t)$ is small.

Suppose the system starts at $t = 0$ in a definite eigenstate $|i\rangle$, so that $c_n(0) = \delta_{ni}$. The perturbative expansion proceeds in orders of $V$:

**Zeroth order:** Nothing happens. The system stays in $|i\rangle$, and $c_f^{(0)}(t) = \delta_{fi}$.

**First order:** Substitute the zeroth-order solution into the right-hand side of the coupled equations. Since only $c_i^{(0)} = 1$ is nonzero, the equation for any $f \neq i$ simplifies to:

$$i\hbar \, \dot{c}_f^{(1)}(t) = V_{fi}(t) \, e^{i\omega_{fi} t}$$

This can be integrated directly:

$$c_f^{(1)}(t) = \frac{1}{i\hbar} \int_0^t V_{fi}(t') \, e^{i\omega_{fi} t'} \, dt'$$

The first-order transition probability is then:

$$P_{i \to f}^{(1)}(t) = |c_f^{(1)}(t)|^2 = \frac{1}{\hbar^2} \left| \int_0^t V_{fi}(t') \, e^{i\omega_{fi} t'} \, dt' \right|^2$$

This is the central result of Dirac's approach at first order. It says: the transition probability is determined by the Fourier component of the perturbation matrix element at the Bohr frequency $\omega_{fi}$. Perturbations that oscillate in resonance with the energy gap between $|i\rangle$ and $|f\rangle$ are effective at driving transitions; off-resonant perturbations are not.

**Higher orders** follow the same pattern — substitute the $(n-1)$th-order solution into the right-hand side to get the $n$th-order correction. Each order introduces one more time integral and one more factor of $V$. In practice, most applications stop at first order, because the algebra and the physical interpretation both become considerably more involved at second order and beyond.

---

## The Route to Fermi's Golden Rule

The first-order result takes a particularly clean form when the perturbation is monochromatic ($V(t) = V_0 e^{-i\omega t} + V_0^\dagger e^{i\omega t}$) or constant, and when the final states form a dense continuum rather than a set of discrete levels.

For a constant perturbation turned on at $t = 0$, the first-order probability becomes:

$$P_{i \to f}^{(1)}(t) = \frac{|V_{fi}|^2}{\hbar^2} \cdot \frac{\sin^2(\omega_{fi} t / 2)}{(\omega_{fi}/2)^2}$$

The function $\sin^2(x t/2) / (x/2)^2$ is sharply peaked around $x = 0$ — that is, around $\omega_{fi} = 0$, meaning $E_f = E_i$. As $t$ grows, the peak becomes narrower and taller, approaching $2\pi t \, \delta(\omega_{fi})$ in the long-time limit. This is energy conservation emerging from the time integral: only transitions that conserve energy (to within $\hbar/t$) accumulate probability over long times.

Summing over a continuum of final states with density $\rho(E_f)$ and taking the long-time limit gives the transition rate:

$$W_{i \to f} = \frac{2\pi}{\hbar} |V_{fi}|^2 \, \rho(E_f)$$

This is **Fermi's Golden Rule**. It is not a separate result from Dirac's theory — it is a specific limiting case of the first-order perturbative transition probability, valid when the perturbation has been on long enough for the energy-conserving delta function to be well resolved, and when the final-state spectrum is dense enough to be treated as continuous.

---

## The Assumptions, Made Explicit

Every step above involves a choice or an approximation. Collecting them:

**The perturbation is weak.** The entire expansion is organized in powers of $V/H_0$. At first order, the transition probability grows as $|V_{fi}|^2$. If $V$ is not small compared to the energy scales of $H_0$, the first-order result can exceed unity — a probability greater than 1, which signals that higher-order terms are not negligible and the expansion is unreliable. This is not a soft failure; it is a qualitative breakdown. The perturbation series does not simply become less accurate; it ceases to converge.

**The system starts in a single eigenstate.** The initial condition $c_n(0) = \delta_{ni}$ means the system is in a pure eigenstate of $H_0$. This excludes thermal mixtures, coherent superpositions, and entangled states as starting points — though the formalism can be extended to handle these by working with the density matrix rather than the state vector.

**The unperturbed basis is the right basis.** By expanding in the eigenstates of $H_0$, Dirac's approach defines "transition" as a change in the quantum numbers of the unperturbed system. This is natural when the perturbation is brief or weak, but it embeds a choice: the coefficients $c_n(t)$ track the system's projection onto the *field-free* states, not onto the instantaneous eigenstates of the full Hamiltonian $H_0 + V(t)$. While the perturbation is active, $|c_f(t)|^2$ includes contributions from the system's reversible polarization response to the field — the distortion of the charge distribution that would relax immediately if the field were removed. This is not a transition in the physical sense, but it contributes to the mathematical quantity $|c_f(t)|^2$. For an isolated system observed only after the perturbation has ended, the distinction is immaterial. In other contexts, it matters; this is discussed further in the nonadiabatic bath coupling notes.

**The system is isolated.** No bath, no environment, no decoherence. The Schrödinger equation preserves probability and phase relationships perfectly. If the system is coupled to a thermal environment while the perturbation is active, the coefficients $c_n(t)$ are no longer sufficient to describe the dynamics; a density matrix and a master equation are needed, and the way Dirac's coefficients enter that master equation is nontrivial.

**The switching is idealized.** The standard derivations assume the perturbation appears suddenly at $t = 0$ (sudden approximation) or rises infinitely slowly (adiabatic switching). Real perturbations have finite switch-on times, which introduce transient effects not captured by either limit. For pulses with smooth envelopes, the first-order integral can be evaluated directly, but the clean asymptotic forms (Fermi's Golden Rule, the energy-conserving delta function) apply only in the appropriate limits.

**Energy conservation is approximate at finite times.** The delta function $\delta(E_f - E_i)$ is a long-time idealization. At finite times, transitions to states that do not exactly conserve energy are allowed, with a probability that oscillates and decays as the energy mismatch grows. This is sometimes called the "energy-time uncertainty" aspect of perturbation theory. It means that Fermi's Golden Rule is not valid at short times, and that the concept of a well-defined transition "rate" only emerges after the system has been perturbed for long enough that the sinc-squared function is well approximated by a delta function.

---

## Where the Framework Breaks Down

The assumptions above define the boundaries of Dirac's approach. Each boundary corresponds to a regime where the formalism fails and different methods are required.

### Strong Fields

When the perturbation strength becomes comparable to the internal energy scales of the system — as in high-intensity laser physics, where the electric field of the laser can rival the Coulomb field binding an electron to a nucleus — the perturbative expansion fails at every order. The transition probability computed at first order can exceed unity, and summing higher orders does not rescue the series because it does not converge. This regime requires non-perturbative methods: direct numerical propagation of the time-dependent Schrödinger equation, the strong-field approximation (SFA) developed by Keldysh and extended by Reiss (1990) and others, or Floquet theory for periodically driven systems. These methods do not treat the field as a small correction; they incorporate it into the zeroth-order description of the system.

### Long Times and Persistent Oscillations

For a constant or periodic perturbation, the first-order coefficients $|c_f(t)|^2$ oscillate indefinitely — the well-known Rabi-type oscillations in the perturbative limit. In a real physical system, these oscillations are eventually damped by environmental interactions, but Dirac's formalism has no mechanism for this. The coefficients simply keep oscillating, and the transition probability never settles to a steady value. This is not wrong, strictly speaking — it is the correct behavior for a perfectly isolated system under a perfectly coherent drive. But it means the formalism cannot describe relaxation, thermalization, or the approach to a steady state without being supplemented by additional physics (master equations, stochastic methods, or explicit environmental models).

### Level Crossings and Near-Degeneracies

Dirac's perturbation theory assumes that the unperturbed energy levels are well separated compared to the perturbation matrix elements. When two levels approach each other — at an avoided crossing in a parameter-dependent Hamiltonian, or at a true crossing (as occurs at Dirac points in graphene band structures) — the perturbative expansion diverges even for weak perturbations because the energy denominator $1/\omega_{fi}$ blows up. Transitions in the vicinity of level crossings are better described by the **Landau-Zener formula**, which treats the two-level crossing problem exactly rather than perturbatively (Faraj & Jin, 2015). The Landau-Zener transition probability depends exponentially on the ratio of the coupling strength to the sweep rate through the crossing, a dependence that no finite order of Dirac's perturbation theory can reproduce.

### Probability Non-Conservation in Non-Hermitian Systems

The standard formalism assumes a Hermitian Hamiltonian, which guarantees that the total probability $\sum_n |c_n(t)|^2 = 1$ at all times. In systems with gain or loss — modeled by non-Hermitian Hamiltonians with complex potentials — this conservation law does not hold. The probability can grow (gain) or shrink (loss), and the symmetry between forward and reverse transition probabilities breaks down. Extending Dirac's framework to non-Hermitian systems requires careful modification of the inner product and the treatment of left and right eigenstates, which differ when $H$ is not self-adjoint (Choi, 2020).

---

## What the Framework Does and Does Not Tell You

Dirac's time-dependent perturbation theory is a method for computing the coefficients of expansion in a particular basis (the eigenstates of $H_0$) under a particular approximation (weakness of $V$). Its outputs — the $|c_f(t)|^2$ — are rigorously the probabilities of finding the system in the unperturbed eigenstates at time $t$, assuming a projective measurement in that basis. Within its domain of validity (weak perturbation, isolated system, well-separated levels, times long enough for Fermi's Golden Rule but not so long that higher-order effects accumulate), it is extraordinarily successful and remains the standard tool of the field.

What it does not do is provide a decomposition of the dynamics into physically distinct processes. The coefficients $c_n(t)$ mix together the reversible polarization response, transient interference effects, and genuine irreversible transitions into a single set of numbers. For an isolated system measured after the perturbation is over, this mixing is harmless — only the irreversible part survives. But for a system that is continuously interacting with both a driving field and a dissipative environment, unpacking that mixture becomes essential.

Recognizing what is inside Dirac's coefficients — and what must be separated out before coupling to a bath — is the bridge between this formalism and the open-systems framework discussed in the nonadiabatic bath coupling notes.

---

## References

* Choi, J. R. (2020). Perturbation Theory for Time-Dependent Quantum Systems Involving Complex Potentials. *Frontiers in Physics*, 8, 189. https://doi.org/10.3389/fphy.2020.00189
* Faraj, A., & Jin, S. (2015). The Landau-Zener transition and the surface hopping method for the 2D Dirac equation for graphene. arXiv:1505.05988. https://doi.org/10.48550/arxiv.1505.05988
* Jang, S. J., & Rhee, Y. M. (2023). Modified Fermi's golden rule rate expressions. *Journal of Chemical Physics*, 159, 024115. https://doi.org/10.1063/5.0152804
* Jovanovski, Mandal, & Hunt. (2021). Quantum transition probabilities due to overlapping electromagnetic pulses: Persistent differences between Dirac's form and nonadiabatic perturbation theory. *Journal of Chemical Physics*, 155, 164101. https://doi.org/10.1063/5.0020169
* Reiss, H. R. (1990). Relativistic strong-field photoionization. *Journal of the Optical Society of America B*, 7(4), 574. https://doi.org/10.1364/josab.7.000574