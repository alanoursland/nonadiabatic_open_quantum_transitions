# Nonadiabatic Transitions in Open Quantum Systems

## The Problem: Driven Systems That Dissipate

Many situations in physics involve a quantum system that is simultaneously **driven** by an external field and **coupled to a thermal environment**. An atom in a laser field that also interacts with the surrounding electromagnetic vacuum. A molecular chromophore absorbing light while embedded in a solvent. A qubit being manipulated by microwave pulses while losing coherence to substrate phonons. In every case, two things happen at once: the field pushes the system between energy levels, and the environment pulls it toward thermal equilibrium.

Describing either process alone is well-trodden ground. Time-dependent perturbation theory (or exact propagation of the Schrödinger equation) handles the driven dynamics. Master equations — Redfield, Lindblad, and their variants — handle the dissipative coupling to the bath. The difficulty emerges when you try to do both at the same time, because the way you define "transition" in the driven problem turns out to matter enormously once the bath enters the picture.

This is the domain of **nonadiabatic bath coupling**: the problem of correctly separating field-driven quantum dynamics from environment-induced relaxation when both are active. Getting this separation wrong doesn't just introduce small numerical errors — it can produce results that violate thermodynamics.

---

## Open Quantum Systems: A Brief Orientation

An **open quantum system** is one that exchanges energy, phase information, or both with its surroundings. Unlike the idealized closed systems of textbook quantum mechanics, open systems cannot generally be described by a wavefunction evolving under the Schrödinger equation. Instead, we work with the **density matrix** $\rho$ and its equation of motion — the **master equation**.

The most widely used form is the **Lindblad master equation**, which guarantees that the density matrix remains physical (positive, trace-preserving) at all times:

$$\frac{d\rho}{dt} = -\frac{i}{\hbar}[H, \rho] + \sum_k \gamma_k \left( L_k \rho L_k^\dagger - \frac{1}{2}\{L_k^\dagger L_k, \rho\} \right)$$

The first term is the usual unitary (Hamiltonian) evolution. The second term, built from **Lindblad operators** $L_k$ and rates $\gamma_k$, describes the irreversible effects of the environment: energy relaxation (decay toward lower states), dephasing (loss of off-diagonal coherence), and thermalization (approach to the Boltzmann distribution at the bath temperature).

The **Redfield equation** is a related but microscopically derived alternative. It starts from a specific model of the system-bath coupling and derives the relaxation rates from the bath's correlation functions, rather than postulating them. It can capture more detailed physics — memory effects, frequency-dependent relaxation — but it does not automatically guarantee positivity of the density matrix, which can cause problems in some regimes.

For the present topic, the critical point is this: both formalisms require you to specify *what the system states are* — the basis in which you write $\rho$ and the operators $L_k$. When no external driving field is present, the natural choice is the energy eigenbasis of the system Hamiltonian $H_0$. The bath then relaxes populations toward the Boltzmann distribution and damps coherences, and everything is thermodynamically consistent.

The trouble starts when the driving field is turned on.

---

## Dirac's Transition Coefficients and Their Hidden Assumption

The standard textbook treatment of a driven quantum system expands the time-evolving state in the **unperturbed eigenbasis** $\{|n\rangle\}$ of the field-free Hamiltonian $H_0$:

$$|\Psi(t)\rangle = \sum_n c_n(t) \, e^{-iE_n t/\hbar} |n\rangle$$

The coefficients $c_n(t)$ evolve according to the time-dependent perturbation equations, and $|c_n(t)|^2$ is interpreted as the probability of finding the system in state $|n\rangle$ at time $t$. These are **Dirac's transition coefficients**, and they are the starting point for Fermi's Golden Rule and most of perturbation theory.

For an isolated system observed at asymptotically long times after a completed interaction (a scattering event, a pulse that has come and gone), these coefficients give physically correct transition probabilities. But during the interaction — while the field is on — they contain something extra.

When an external field acts on a quantum system, the electron cloud (or, more generally, the charge and current distribution) **polarizes**: it distorts to partially follow the field. This is not a transition. No energy has been irreversibly absorbed. The system's charge distribution has simply adjusted to the instantaneous field, the way a classical dielectric polarizes in a capacitor. If the field is turned off slowly, this distortion reverses and the system returns to its original state.

Dirac's coefficients do not distinguish between this reversible polarization and a genuine nonadiabatic transition. Both contribute to $|c_n(t)|^2$. In an isolated system observed only after the field is gone, this doesn't matter — the polarization contribution vanishes when the field is off, and what remains is the true transition probability. But if a bath is present *while the field is on*, the distinction becomes critical.

---

## Why the Standard Approach Fails Thermodynamically

Consider what happens when you plug Dirac's coefficients into a master equation. The bath "sees" the populations $|c_n(t)|^2$ and attempts to relax them toward the thermal equilibrium distribution. But those populations include the virtual polarization — the reversible, field-following part of the wavefunction. The bath treats this polarization as though it were a real excitation and tries to relax it. The result is a system that:

**Relaxes toward the wrong equilibrium.** The bath drives the populations toward the Boltzmann distribution of the *unperturbed* Hamiltonian $H_0$, not toward the thermal state of the *actual* (field-dressed) system. Since the field modifies the effective energy levels, the correct thermal state is different from the unperturbed one. Using Dirac's coefficients, the system thermalizes to the wrong target.

**Generates fictitious heat flow.** Dirac's coefficients oscillate in time even for a constant driving field — these are the well-known Rabi oscillations in the bare-state basis. The bath continually tries to damp these oscillations, which means it is perpetually extracting "heat" from the system. But much of what it's damping is the virtual polarization, not real excitation energy. The result is a steady-state heat current that has no physical source: energy appears to flow from the system to the bath without being replenished by the field. This violates the first law of thermodynamics.

**Invents transitions that don't exist.** In certain parameter regimes, the dephasing terms of the master equation, acting on Dirac's coefficients, can transfer population between states in a way that looks like a bath-induced transition. But these apparent transitions are artifacts: they arise because the coherent polarization has been misidentified as a population that the bath can act on.

These are not edge cases or small corrections. They represent a structural incompatibility between Dirac's decomposition and the requirements of thermodynamic consistency.

---

## The Nonadiabatic Decomposition: Separating Real Transitions from Polarization

The resolution is to decompose the system's time evolution differently. Instead of tracking the coefficients relative to the unperturbed states $|n\rangle$, one separates the dynamics into an **adiabatic** (field-following) part and a **nonadiabatic** (genuinely transitional) part.

The adiabatic part captures the instantaneous polarization response — the piece of the wavefunction that continuously adjusts to the field without any irreversible energy exchange. In the limit of a very slowly varying field, this is the only part that matters, and no transitions occur. The nonadiabatic part captures everything else: the component of the evolution that corresponds to real, irreversible changes in the system's quantum state. These are the true transitions.

Concretely, this means constructing a set of **nonadiabatic transition coefficients** that have the adiabatic response subtracted out. When these corrected coefficients are used in the master equation instead of Dirac's, the bath acts only on genuine excitations. The polarization is invisible to the dissipative dynamics, as it should be — you don't thermalize a polarization any more than you thermalize the displacement of a spring that is currently being held in a compressed position.

This decomposition is not new in the sense that the adiabatic theorem and the concept of nonadiabatic transitions have been part of quantum mechanics since its earliest days. What is relatively recent is the rigorous application of this idea to the open-systems problem — the recognition that the choice of decomposition is not merely a matter of convenience or interpretation, but has hard physical consequences when dissipation is present.

---

## The Work of Mandal, Hunt, and Collaborators

The program of applying the nonadiabatic decomposition to open quantum systems has been developed in a series of papers by Anirban Mandal, Katharine Hunt, and collaborators, with a particularly important synthesis by Jovanovski, Mandal, and Hunt (2021). Their contribution is both conceptual and quantitative.

**Conceptually**, they demonstrated that Dirac's transition probability is not merely an approximation that becomes exact in certain limits — it is *physically incompatible* with thermodynamics whenever dissipation is present during the driving. The nonadiabatic decomposition is not an optional refinement; it is required for the theory to respect the first law and to produce the correct thermal steady state. This is a stronger claim than simply saying one method is more accurate than another. It says the standard method is structurally wrong for this class of problems.

**Formally**, they showed that the nonadiabatic framework respects two properties that Dirac's does not in the dissipative context. **Gauge invariance**: the physical predictions do not depend on arbitrary choices in how the electromagnetic field is mathematically represented (e.g., length gauge vs. velocity gauge). In Dirac's formulation, the instantaneous transition probability during a pulse can change depending on the gauge — a clearly unphysical result. The nonadiabatic coefficients are gauge-invariant by construction. **Asymptotic consistency**: if a field is turned on and then turned off without depositing net energy in the system, the nonadiabatic transition probability returns exactly to zero. Dirac's coefficients generally do not satisfy this unless the field vanishes everywhere, because the polarization contribution can persist as a transient even after the pulse peak has passed.

**Quantitatively**, for systems driven by pulses with finite switch-on times or by overlapping pulse sequences — situations common in ultrafast spectroscopy and coherent control — they showed that the transition probabilities computed from Dirac's coefficients can differ from the nonadiabatic results by more than 35%. This is not a small correction buried in high-order terms. It is a leading-order discrepancy that propagates into the predicted heat flow, steady-state populations, and spectral response of the system.

**Physically**, their work clarifies the role of dephasing. In the standard Dirac-based master equation, dephasing can appear to cause or enhance transitions — an effect that is difficult to interpret and sometimes treated as a real phenomenon ("dephasing-assisted transport" claims must be evaluated carefully in this light). In the nonadiabatic framework, dephasing acts only on the coherence of real transitions. It destroys phase relationships between states that have genuinely been populated by the field, but it does not create population transfer on its own. The separation between what the field does and what the bath does becomes clean.

---

## Implications and Open Questions

The nonadiabatic bath coupling framework has direct consequences for any setting where a quantum system is both driven and dissipative.

In **ultrafast spectroscopy**, pump-probe experiments and multidimensional spectroscopy techniques interrogate molecular systems with sequences of short, intense pulses while the molecules are embedded in solvents or protein environments that act as thermal baths. The interpretation of these experiments — extracting relaxation rates, energy transfer pathways, and coupling constants from the measured signals — relies on theoretical models that combine driving and dissipation. If those models use Dirac's coefficients, the extracted parameters may be systematically biased.

In **quantum control and quantum computing**, gate operations are precisely shaped pulses applied to qubits that are simultaneously decohering. Error models that do not correctly separate the coherent field response from genuine nonadiabatic transitions will misattribute error sources — blaming the bath for errors that are actually artifacts of the theoretical description, or failing to account for real bath-induced effects that are masked by the larger polarization contribution.

In **quantum thermodynamics**, the emerging field that studies heat, work, and entropy at the quantum scale, the correct identification of energy transfer between the field, the system, and the bath is foundational. The nonadiabatic framework provides a decomposition in which the first law holds at every instant, not just on average or in the long-time limit. This makes it a natural starting point for consistent definitions of quantum work and quantum heat.

Several questions remain open. The framework as developed applies most naturally to few-level systems driven by classical fields. Extending it to quantized radiation fields, to many-body systems with internal interactions, and to non-Markovian baths (where the environment has memory) are active areas of development. The relationship between the nonadiabatic decomposition and other approaches to the same problem — Floquet theory for periodically driven systems, polaron transformations for strong system-bath coupling — is not yet fully mapped out.

What is clear is that the textbook approach of simply inserting Dirac's coefficients into a Lindblad or Redfield equation is not adequate when the driving and dissipation overlap in time. The nonadiabatic decomposition, as developed by Mandal, Hunt, and collaborators, provides a physically grounded alternative that respects thermodynamics and removes gauge artifacts. For any quantitative treatment of driven, dissipative quantum systems, this distinction matters.