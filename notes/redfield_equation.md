# Redfield Theory: The Master Equation for Open Quantum Systems

## The Problem Redfield Theory Solves

A quantum system in the real world is never perfectly isolated. An atom in a solid interacts with lattice vibrations. A molecule in solution collides with solvent molecules. A qubit in a superconducting circuit couples to electromagnetic noise in its substrate. In every case, the system of interest exchanges energy and phase information with a much larger environment — the **bath**.

The full quantum mechanical description of the system plus bath is, in principle, exact: write down the total Hamiltonian, solve the Schrödinger equation for the combined system, and extract the properties of the subsystem you care about. In practice, this is impossible. The bath has an enormous (effectively infinite) number of degrees of freedom. You cannot track them all, and you do not want to — the goal is to describe the system alone, with the bath's effects captured implicitly.

**Redfield theory** provides a systematic framework for doing this. Starting from the total system-bath Hamiltonian, it derives an equation of motion for the **reduced density matrix** of the system — the density matrix obtained by tracing over all bath degrees of freedom. The result is a **master equation**: a differential equation that describes how the system's populations and coherences evolve under the combined influence of its own internal dynamics and the dissipative effects of the environment.

---

## The Total Hamiltonian

The starting point is the decomposition of the full Hamiltonian into three parts:

$$H_{\text{total}} = H_S + H_B + H_{SB}$$

$H_S$ is the **system Hamiltonian**, whose eigenstates $|a\rangle$ and eigenvalues $E_a$ define the energy levels of the isolated system. $H_B$ is the **bath Hamiltonian**, describing the free dynamics of the environment (typically modeled as a collection of harmonic oscillators or a thermal radiation field). $H_{SB}$ is the **system-bath coupling**, the interaction that allows energy and phase information to flow between the system and its environment.

The system-bath coupling is generally written in factored form:

$$H_{SB} = \sum_\alpha S_\alpha \otimes B_\alpha$$

where $S_\alpha$ are operators acting on the system and $B_\alpha$ are operators acting on the bath. This factored form is not a restriction — any bilinear coupling can be written this way — but it makes the subsequent algebra tractable.

The reduced density matrix of the system is obtained by tracing over the bath:

$$\rho(t) = \text{Tr}_B[\rho_{\text{total}}(t)]$$

The goal is to find a closed equation of motion for $\rho(t)$ alone, without needing to track the bath's state.

---

## The Born-Markov Approximations

Deriving a closed equation for $\rho(t)$ from the exact dynamics of $\rho_{\text{total}}(t)$ requires approximations. Redfield theory makes two, and understanding what each one does — and what it costs — is essential.

### The Born Approximation (Weak Coupling)

The Born approximation assumes that the system-bath coupling $H_{SB}$ is weak enough that the bath is essentially unaffected by the system. The bath remains in its thermal equilibrium state $\rho_B^{\text{eq}}$ at all times, and the total density matrix factorizes:

$$\rho_{\text{total}}(t) \approx \rho(t) \otimes \rho_B^{\text{eq}}$$

Physically, this means the bath is so large and so well thermalized that any perturbation the system imparts to it dissipates instantly into the bath's vast number of degrees of freedom. The bath acts as an infinite reservoir — it absorbs energy from or delivers energy to the system without its own temperature or state changing appreciably.

This approximation is second order in the coupling $H_{SB}$. It captures the leading-order effects of the environment (relaxation and dephasing) but neglects higher-order processes like bath-mediated correlations between system transitions and renormalization of the system's energy levels by the bath (the Lamb shift, though a partial Lamb shift correction can be included within the Redfield framework).

### The Markov Approximation (Memoryless Bath)

The Markov approximation assumes that the bath's correlation functions decay on a timescale $\tau_B$ that is much shorter than the timescale $\tau_S$ on which the system evolves. In physical terms: the bath "forgets" its interactions with the system almost immediately.

When this condition holds, the rate of change of $\rho(t)$ at time $t$ depends only on $\rho(t)$ at that same instant — not on the system's history. The exact integro-differential equation (where the future evolution depends on an integral over the past) reduces to a first-order differential equation with time-independent coefficients. This is an enormous simplification: it means the system's dynamics are described by a linear, time-local equation, and the bath's influence is entirely encoded in a set of constant rates.

The Markov approximation is good when the bath has a broad, featureless spectral density — many modes spanning a wide frequency range, none of which is strongly coupled to the system. It breaks down when the bath has sharp spectral structure (e.g., a cavity mode or a narrow phonon band), when the temperature is very low (so that few bath modes are thermally populated), or when the system-bath coupling is strong enough that the bath's response time is no longer negligible.

---

## The Redfield Equation

Applying both approximations to the equation of motion for $\rho(t)$, working in the energy eigenbasis $\{|a\rangle\}$ of $H_S$, produces the **Redfield equation**:

$$\frac{d\rho_{ab}(t)}{dt} = -i\omega_{ab}\,\rho_{ab}(t) + \sum_{cd} R_{abcd}\,\rho_{cd}(t)$$

The first term is the free (coherent) evolution of the system: each element of the density matrix oscillates at the Bohr frequency $\omega_{ab} = (E_a - E_b)/\hbar$ between the corresponding energy levels. This is what the system would do in the absence of the bath.

The second term contains all of the bath's effects, encoded in the **Redfield relaxation tensor** $R_{abcd}$. This tensor is a fourth-rank object with indices running over the system's energy eigenstates. Its elements are determined by the bath correlation functions — the Fourier transforms of the two-time correlation functions $\langle B_\alpha(t) B_{\beta}(0) \rangle_{\text{eq}}$ evaluated at the system's Bohr frequencies.

### What the Redfield Tensor Contains

The Redfield tensor is constructed from the **one-sided Fourier transforms** of the bath correlation functions, often written as:

$$\Gamma_{\alpha\beta}(\omega) = \int_0^\infty d\tau \; e^{i\omega\tau} \; \langle B_\alpha(\tau) B_\beta(0) \rangle_{\text{eq}}$$

These are complex quantities. Their real parts give the **relaxation rates** — the rates at which populations and coherences decay. Their imaginary parts give the **Lamb shift corrections** — small frequency shifts of the system's energy levels due to the bath coupling.

The connection to the **spectral density** $J(\omega)$ of the bath is direct. For a bath of harmonic oscillators linearly coupled to the system, the bath correlation functions are determined by $J(\omega)$ and the thermal occupation numbers $\bar{n}(\omega) = (e^{\hbar\omega/k_BT} - 1)^{-1}$. The Redfield rates at a given Bohr frequency $\omega_{ab}$ sample the spectral density at that frequency, weighted by the thermal factors. A bath with more modes near a system transition frequency produces faster relaxation of that transition.

Different index combinations of $R_{abcd}$ describe different physical processes:

**Population relaxation** ($R_{aabb}$ with $a \neq b$): the rate at which population transfers from state $|b\rangle$ to state $|a\rangle$, driven by bath-induced transitions. These rates satisfy **detailed balance** — the ratio of the forward and reverse rates equals the Boltzmann factor $e^{-\hbar\omega_{ab}/k_BT}$ — which ensures that the system relaxes toward the correct thermal equilibrium. The characteristic timescale for population decay is the longitudinal relaxation time $T_1$.

**Dephasing** ($R_{abab}$ for $a \neq b$): the rate at which off-diagonal coherences $\rho_{ab}$ decay. Dephasing has two contributions: "population-transfer dephasing" (each population-relaxation event also destroys the coherence involving the states whose populations changed) and "pure dephasing" (fluctuations in the energy levels caused by the bath, without any population transfer). The characteristic timescale is the transverse relaxation time $T_2$, which satisfies the inequality $T_2 \leq 2T_1$ — coherences always decay at least as fast as populations.

**Coherence transfer** ($R_{abcd}$ with $\{a,b\} \neq \{c,d\}$): the rate at which one coherence $\rho_{cd}$ feeds into another $\rho_{ab}$. These terms couple the equations for different off-diagonal elements and are responsible for effects like quantum beats that persist in the presence of the bath. They are the terms that the secular approximation discards.

---

## The Positivity Problem

Despite its systematic derivation, Redfield theory has a well-known pathology: the solutions can violate the **positivity of the density matrix**. The diagonal elements $\rho_{aa}$, which represent the populations of the energy eigenstates, can become negative during the time evolution — an unphysical result, since a negative population has no probabilistic interpretation.

### Why It Happens

The root cause is a subtle inconsistency in how the Born-Markov approximations are applied. The derivation of the Redfield equation is second-order in the system-bath coupling: terms of order $H_{SB}^2$ are retained, while higher orders are dropped. But the *dynamics* generated by the resulting equation are not restricted to second order — they resum the second-order rates into an exponential time evolution that propagates to all times. This resummation can produce transient behaviors that are artifacts of the truncation, not features of the true dynamics.

More technically, the Redfield equation generates a dynamical semigroup only approximately. A **completely positive** dynamical map — one that guarantees the density matrix remains a valid quantum state at all times — must have the mathematical structure of a **Lindblad equation**. The Redfield equation does not automatically satisfy this structure. The non-Lindblad terms are small (they are the same order as the terms that were already neglected in the Born approximation), but they can cause the density matrix to develop small negative eigenvalues, particularly during short transients, at very low temperatures, or when energy levels are nearly degenerate.

### The Physical Significance

In many practical applications, the positivity violations are numerically tiny and confined to short times — the system quickly settles into a regime where the density matrix is well-behaved. But the violation is a matter of principle: a negative probability, however small, signals that the mathematical framework has produced something unphysical. For applications requiring formal guarantees of physicality (quantum information processing, entanglement measures, entropy calculations), this is unacceptable.

---

## The Secular Approximation and the Lindblad Form

The standard resolution is the **secular approximation**: discard all terms in the Redfield tensor that couple density matrix elements oscillating at different Bohr frequencies.

In the Redfield equation, the term $R_{abcd}\,\rho_{cd}$ is multiplied by an implicit oscillating factor $e^{i(\omega_{ab} - \omega_{cd})t}$ when working in the Schrödinger picture. When $\omega_{ab} \neq \omega_{cd}$, this factor oscillates rapidly and its time-averaged effect is small. The secular approximation drops these terms, retaining only those where $\omega_{ab} = \omega_{cd}$.

The result is an equation that takes the **Lindblad form**:

$$\frac{d\rho}{dt} = -\frac{i}{\hbar}[H_S + H_{LS},\, \rho] + \sum_k \gamma_k \left(L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k,\, \rho\}\right)$$

where $H_{LS}$ is the Lamb shift Hamiltonian and the $L_k$ are **Lindblad operators** constructed from the system's eigenstate transition operators, with rates $\gamma_k$ determined by the bath spectral density evaluated at the corresponding Bohr frequencies. This form automatically guarantees complete positivity — the density matrix remains a valid quantum state at all times.

### What the Secular Approximation Costs

The guarantee of positivity comes at a price. By dropping the coherence transfer terms ($R_{abcd}$ with $\omega_{ab} \neq \omega_{cd}$), the secular approximation eliminates the coupling between different off-diagonal elements of the density matrix. This is accurate when the system's energy levels are well separated compared to the relaxation rates — when the Bohr frequencies are all distinct and much larger than the $R_{abcd}$. But when energy levels are close together or nearly degenerate, the "fast-oscillating" terms are not actually fast, and discarding them loses real physics.

This matters in systems like photosynthetic light-harvesting complexes, where closely spaced electronic states exhibit long-lived quantum coherences. In NMR, where chemical shift differences between spins can be comparable to relaxation rates, the full (non-secular) Redfield equation is often necessary. In these contexts, researchers use the full Redfield equation and monitor for positivity violations, accepting the small risk in exchange for more accurate coherence dynamics.

Several intermediate approaches exist: partial secular approximations that keep terms coupling nearly degenerate levels while dropping the rest, and alternative derivations (coarse-grained master equations, Nakajima-Zwanzig formalisms with controlled resummations) that attempt to preserve both accuracy and positivity. None is universally satisfactory, and the tension between microscopic accuracy and mathematical consistency remains an active area of research.

---

## Connection to Driven Systems

The standard Redfield equation is written in the eigenbasis of the time-independent system Hamiltonian $H_S$. When an external driving field is present — a laser, a microwave pulse, a time-varying electric field — the situation changes in a way that connects directly to the concerns discussed in the nonadiabatic bath coupling notes.

If the driving field is simply added to the system Hamiltonian ($H_S \to H_S + V(t)$) while the Redfield relaxation tensor is kept in its original form (computed from the eigenstates of $H_S$ alone), the bath will relax the system toward the thermal equilibrium of the *undriven* system, not the correct equilibrium of the driven one. The relaxation tensor was derived under the assumption that $H_S$ is time-independent; applying it unchanged to a driven system implicitly treats the driving as part of the coherent dynamics while leaving the dissipative dynamics anchored to the wrong basis.

This is the open-systems manifestation of the same issue that the Landau-Lifshitz separation addresses for isolated systems: the Dirac coefficients, defined relative to the unperturbed basis, conflate adiabatic polarization with genuine transitions. When those coefficients enter the Redfield tensor, the bath acts on the polarization as though it were a real population, producing the thermodynamic inconsistencies described in the nonadiabatic bath coupling document. The nonadiabatic framework provides a route to correcting this: by redefining the system's transition coefficients so that only genuine excitations are exposed to the bath, the modified master equation respects energy conservation and drives the system toward the correct steady state.

---

## Summary

| Feature | Full Redfield | Secular (Lindblad) |
| :--- | :--- | :--- |
| **Derivation** | Microscopic, from $H_{SB}$ and bath correlation functions | Redfield + secular approximation, or axiomatic |
| **Positivity** | Not guaranteed; may produce small negative populations | Guaranteed at all times |
| **Coherence transfer** | Retained; important for nearly degenerate levels | Dropped; accurate only for well-separated levels |
| **Detailed balance** | Satisfied; relaxes to correct thermal state | Satisfied; relaxes to correct thermal state |
| **Typical applications** | NMR, molecular spectroscopy, photosynthetic complexes | Quantum optics, quantum information, general dissipation |

---

## References

* Redfield, A. G. (1957). On the theory of relaxation processes. *IBM Journal of Research and Development*, 1(1), 19–31.
* Redfield, A. G. (1965). The theory of relaxation processes. *Advances in Magnetic and Optical Resonance*, 1, 1–32.
* Bloch, F. (1957). Generalized theory of relaxation. *Physical Review*, 105(4), 1206–1222.
* Lindblad, G. (1976). On the generators of quantum dynamical semigroups. *Communications in Mathematical Physics*, 48, 119–130.
* Gorini, V., Kossakowski, A., & Sudarshan, E. C. G. (1976). Completely positive dynamical semigroups of N-level systems. *Journal of Mathematical Physics*, 17, 821–825.
* Breuer, H.-P., & Petruccione, F. (2002). *The Theory of Open Quantum Systems*. Oxford University Press.