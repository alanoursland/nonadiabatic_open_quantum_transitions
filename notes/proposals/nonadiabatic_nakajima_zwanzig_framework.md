# The Nonadiabatic Nakajima-Zwanzig Framework: A Proposed General Theory

## The Proposal in Brief

The theory gap identified in this collection — the absence of a master equation for driven, dissipative quantum systems that simultaneously satisfies complete positivity, non-secular accuracy, gauge invariance, and thermodynamic consistency — can potentially be closed by combining two existing formalisms that have not previously been integrated: the **nonadiabatic decomposition** of Landau-Lifshitz / Mandal-Hunt, and the **Nakajima-Zwanzig generalized quantum master equation** with a modified, time-dependent projection operator.

The idea is to change what the projection operator projects *onto*. Instead of projecting the total system-bath state onto the bare system eigenstates of $H_0$ (which is what standard Nakajima-Zwanzig does, and which inherits the Dirac-coefficient problem), project onto the **nonadiabatic subspace** — the space of genuinely excited states, with the adiabatic polarization already subtracted. The resulting equation of motion — the **nonadiabatic GQME** — is formally exact, automatically gauge-invariant, automatically thermodynamically consistent, and automatically completely positive (because the exact reduced dynamics of a subsystem of a unitarily evolving whole is always completely positive). The standard approximate theories emerge as controlled limits.

---

## Two Layers of the Problem

The driven open quantum systems problem has two distinct layers, and conflating them is what has made the theory gap so persistent.

**Layer 1: What the bath sees.** This is the problem Hunt and Mandal have solved. When a driving field is present, Dirac's transition coefficients $|c_k(t)|^2$ conflate the reversible adiabatic polarization with genuine nonadiabatic transitions. If these coefficients enter the master equation, the bath acts on the polarization, producing the three thermodynamic pathologies described in the nonadiabatic bath coupling notes (wrong equilibrium, fictitious heat flow, spurious dephasing-induced transitions). The solution is to use the nonadiabatic coefficients $|b_k(t)|^2$ instead.

**Layer 2: What mathematical structure couples the bath to the system.** This is the "better treatment" Hunt is seeking. Redfield theory provides a microscopic derivation but breaks positivity. Lindblad theory guarantees positivity but requires the secular approximation. The CGME and GAME provide intermediate options. All of these operate at layer 2 — they specify the *structure* of the dissipative coupling — and all of them, as currently formulated, use the standard (Dirac-basis) projection at layer 1.

The existing work has addressed layer 1 within the simplest layer-2 structure (Redfield). The proposal here is to address both layers simultaneously by modifying the projection operator in the most general layer-2 framework available: the Nakajima-Zwanzig equation.

---

## The Standard Nakajima-Zwanzig Formalism

The Nakajima-Zwanzig GQME starts from the exact Liouville-von Neumann equation for the total system-bath state $\rho_{\text{total}}(t)$ and introduces a projection operator $\mathcal{P}$ that defines what we mean by "the system's state." The complement $\mathcal{Q} = 1 - \mathcal{P}$ captures everything else.

The standard choice is:

$$\mathcal{P}\rho_{\text{total}} = \text{Tr}_B[\rho_{\text{total}}] \otimes \rho_B^{\text{eq}}$$

This traces over the bath and replaces the bath state with its equilibrium. The "system state" is the reduced density matrix $\sigma(t) = \text{Tr}_B[\rho_{\text{total}}(t)]$, expressed in the eigenbasis of the system Hamiltonian $H_S$ (or, when a driving field is present, in the eigenbasis of the unperturbed $H_0$).

Applying the projection to the Liouville-von Neumann equation yields the exact GQME:

$$\frac{d\sigma(t)}{dt} = -\frac{i}{\hbar}\langle\mathcal{L}\rangle\sigma(t) - \int_0^t d\tau\;\mathcal{K}(\tau)\,\sigma(t-\tau) + \mathcal{I}(t)$$

where $\mathcal{K}(\tau)$ is the memory kernel encoding the bath's influence and $\mathcal{I}(t)$ is the inhomogeneous term encoding initial system-bath correlations.

This equation is exact for any system, any bath, any coupling strength, any driving protocol. The problem is that the "system state" $\sigma(t)$ is defined via the standard projection — which means its matrix elements are the Dirac coefficients. When a driving field is present, $\sigma(t)$ includes the adiabatic polarization, and any subsequent approximation of the memory kernel (Born, Markov, secular) will produce a master equation in which the bath acts on that polarization.

---

## The Nonadiabatic Projection Operator

The modification is to replace the standard projection with one that projects onto the nonadiabatic subspace.

At each instant $t$, the system's state can be decomposed (following Landau-Lifshitz and Mandal-Hunt) into an adiabatic part (the instantaneous ground state of $H_0 + V(t)$, expanded in the unperturbed basis) and a nonadiabatic part (everything that represents genuine excitation). The adiabatic part is what static perturbation theory gives for the current value of the field. The nonadiabatic part is what the field's time variation has irreversibly produced.

Define the **nonadiabatic projection operator** $\mathcal{P}_{\text{nad}}(t)$ as:

$$\mathcal{P}_{\text{nad}}(t)\rho_{\text{total}}(t) = \sigma_{\text{nad}}(t) \otimes \rho_B^{\text{eq}}$$

where $\sigma_{\text{nad}}(t)$ is the reduced density matrix constructed from the nonadiabatic coefficients $b_k(t)$ rather than the full Dirac coefficients $c_k(t)$. The matrix elements of $\sigma_{\text{nad}}$ are:

$$[\sigma_{\text{nad}}]_{kl}(t) = b_k(t)\,b_l^*(t)$$

rather than $c_k(t)\,c_l^*(t)$.

This projection operator is **time-dependent** — it changes as the driving field changes, because the decomposition into adiabatic and nonadiabatic parts depends on the instantaneous value of $V(t)$. This is a departure from the standard Nakajima-Zwanzig formalism, which typically uses a time-independent projector, but the formalism has been extended to time-dependent projectors in the literature. The key consequence is an additional term in the GQME.

**A caveat on near-degeneracies.** The Landau-Lifshitz integration by parts divides by the Bohr frequency difference $\omega_{k0} = (E_k - E_0)/\hbar$. When two energy levels are nearly degenerate ($\omega_{k0} \to 0$), the adiabatic term $a_k(t)$ becomes large and the decomposition can become numerically unstable — the "adiabatic polarization" of a nearly degenerate level is not a small correction but a dominant contribution, and separating it from the nonadiabatic part requires high precision. Hunt's current systems have well-separated energy levels where this issue does not arise. For the general framework, however, near-degenerate levels may need to be treated as a coupled subspace — analogous to the treatment in degenerate perturbation theory, where the degenerate states are first diagonalized within their subspace before perturbative corrections are applied. The nonadiabatic projection would then separate the full Hilbert space into three parts: the adiabatic ground-state adjustment, the near-degenerate subspace (treated exactly, without integration by parts), and the well-separated excited states (where the Landau-Lifshitz decomposition applies cleanly). This partitioned approach is conceptually similar to the selective application of the secular approximation in the CGME, and it would extend the framework's validity to molecular systems with closely spaced electronic states.

---

## The Nonadiabatic GQME

Applying the time-dependent nonadiabatic projection to the Liouville-von Neumann equation yields the **nonadiabatic GQME**:

$$\frac{d\sigma_{\text{nad}}(t)}{dt} = -\frac{i}{\hbar}\langle\mathcal{L}\rangle_{\text{nad}}\,\sigma_{\text{nad}}(t) - \int_0^t d\tau\;\mathcal{K}_{\text{nad}}(t, \tau)\,\sigma_{\text{nad}}(t-\tau) + \mathcal{I}_{\text{nad}}(t) + \mathcal{A}(t)$$

The terms are:

**The projected Liouvillian** $\langle\mathcal{L}\rangle_{\text{nad}}$ describes the coherent dynamics of the nonadiabatic density matrix. It includes the system's Hamiltonian evolution minus the adiabatic part — that is, it drives only the genuine transitions.

**The nonadiabatic memory kernel** $\mathcal{K}_{\text{nad}}(t, \tau)$ encodes the bath's influence on the nonadiabatic populations and coherences. Because the projection has already removed the adiabatic polarization, this kernel describes the bath's coupling to genuine excitations only. Note that the kernel is now a function of two times $(t, \tau)$ rather than just the lag $\tau$, because the nonadiabatic projection is time-dependent. This two-time dependence is the mathematical reflection of the fact that the meaning of "nonadiabatic transition" changes as the field changes.

**The inhomogeneous term** $\mathcal{I}_{\text{nad}}(t)$ accounts for initial system-bath correlations, projected onto the nonadiabatic subspace.

**The adiabatic evolution term** $\mathcal{A}(t)$ arises from the time dependence of the projection operator itself: $\mathcal{A}(t) = \dot{\mathcal{P}}_{\text{nad}}(t)\rho_{\text{total}}(t)$. This term captures the rate at which the adiabatic subspace is changing — the rate at which the field is reshaping the instantaneous ground state. It is the formal expression of the adiabatic dynamics that the standard GQME conflates with the nonadiabatic transitions. In the nonadiabatic GQME, this term is separated out explicitly: it modifies the coherent evolution but does not enter the dissipative dynamics.

---

## Properties of the Nonadiabatic GQME

### Exactness

The equation is exact. No Born, Markov, or secular approximation has been made. The Nakajima-Zwanzig formalism produces an exact equation for whatever subspace the projection operator selects, and the nonadiabatic projection is a valid projection operator. The exactness holds for arbitrary system size, arbitrary coupling strength, and arbitrary driving protocols.

### Gauge Invariance

The nonadiabatic coefficients $b_k(t)$ are gauge-invariant, as established by Mandal and Hunt. The projection operator $\mathcal{P}_{\text{nad}}(t)$ is defined in terms of these gauge-invariant quantities. Therefore the nonadiabatic density matrix $\sigma_{\text{nad}}(t)$, the memory kernel $\mathcal{K}_{\text{nad}}$, and all physical predictions derived from the GQME are gauge-invariant. This is inherited, not imposed.

### Thermodynamic Consistency

When the driving field is constant ($\partial V / \partial t = 0$), the nonadiabatic populations $|b_k(t)|^2$ are constant (this is the plateau-pulse result from the Mandal-Hunt program). The adiabatic evolution term $\mathcal{A}(t)$ vanishes (because the projection is no longer changing). The memory kernel acts on a static nonadiabatic density matrix and drives it toward the thermal equilibrium of the field-dressed system. The heat flow is zero. The first law is satisfied.

When the field is varying, the adiabatic term $\mathcal{A}(t)$ is nonzero and accounts for the changing polarization, but this term enters the coherent dynamics, not the dissipative dynamics. The bath continues to act only on the nonadiabatic part. The power balance — field work equals system energy change plus heat dissipated — holds at every instant, because the energy separation (proved by Mandal and Hunt to be complete, with no cross-terms) is built into the projection.

### Complete Positivity

The exact reduced dynamics of a subsystem of a unitarily evolving composite system is always completely positive. This is a theorem — the Stinespring dilation theorem guarantees it. The nonadiabatic GQME generates the exact reduced dynamics (projected onto the nonadiabatic subspace), so it is automatically completely positive.

The positivity problem arises only when the memory kernel is *approximated*. The exact kernel generates a CPTP map; a truncated kernel may not. This is the same situation as in the standard GQME — the Redfield equation is a second-order truncation of the standard memory kernel, and it breaks positivity. The nonadiabatic Redfield equation would be a second-order truncation of the nonadiabatic memory kernel, and it might also break positivity.

But the critical point is that positivity-preserving truncation schemes — coarse-graining (CGME), geometric-arithmetic interpolation (GAME) — can be applied to the nonadiabatic kernel just as they are applied to the standard kernel. The result would be a nonadiabatic CGME or nonadiabatic GAME: a master equation that satisfies all four criteria (positivity, non-secular accuracy, gauge invariance, thermodynamic consistency) within the Born-Markov regime.

### BQP Structure

The exact nonadiabatic memory kernel $\mathcal{K}_{\text{nad}}(t, \tau)$ encodes the full quantum dynamics of the bath and its correlations with the system's nonadiabatic transitions. For a system of $n$ levels coupled to a bath, evaluating this kernel exactly requires propagating the dynamics in the $\mathcal{Q}_{\text{nad}}$-subspace — the complement of the nonadiabatic projection — which involves the full bath Hilbert space. This is BQP-hard for general $n$.

A quantum computer can evaluate the kernel efficiently: prepare the initial state, evolve under the full Hamiltonian, and project onto the nonadiabatic subspace. The memory kernel is extracted from the correlation functions of the projected dynamics. For small $n$, this can also be done classically (using HEOM, process tensors, or direct integration). For large $n$, only the quantum implementation scales polynomially.

The classical approximation hierarchy reduces the BQP-hard kernel to tractable forms:
- Born approximation → second-order nonadiabatic kernel (classically polynomial in $n$, but may break positivity)
- Born + Markov → nonadiabatic Redfield (classically polynomial, local in time, may break positivity)
- Born + Markov + coarse-graining → nonadiabatic CGME (classically polynomial, positive, non-secular, gauge-invariant, thermodynamically consistent)
- Born + Markov + secular → nonadiabatic Lindblad (classically polynomial, positive, loses non-secular accuracy, but retains gauge invariance and thermodynamic consistency)

Each step discards quantum correlations to reduce the complexity, exactly as the BQP constraints document predicts.

---

## The Hierarchy of Approximations

One of the strengths of this framework is that all existing master equations emerge as controlled limits, and the *specific physics each approximation discards* is identifiable.

### Exact Nonadiabatic GQME (BQP)
- All criteria satisfied
- Memory kernel encodes full system-bath correlations in the nonadiabatic subspace
- Efficiently evaluable on a quantum computer; exponentially costly classically for large $n$

### Nonadiabatic Redfield (Born + Markov)
- Gauge-invariant and thermodynamically consistent (from the nonadiabatic projection)
- Non-secular accuracy retained
- May violate positivity (from the second-order truncation of the kernel)
- This is approximately what Jovanovski, Mandal, and Hunt (2023) implemented

### Nonadiabatic CGME (Born + Markov + Coarse-Graining)
- All four criteria satisfied within the weak-coupling, Markovian regime
- Positivity guaranteed by the coarse-graining construction
- Non-secular terms retained for frequencies within the coarse-graining window
- This is the most immediate "better treatment" — tractable now for few-level systems

### Nonadiabatic Lindblad (Born + Markov + Secular)
- Gauge-invariant and thermodynamically consistent
- Positive by construction
- Loses non-secular accuracy
- Appropriate for systems with well-separated energy levels under driving

### Standard Redfield (Born + Markov + Dirac Basis)
- Loses gauge invariance and thermodynamic consistency for driven systems
- Retains non-secular accuracy
- May violate positivity
- The current standard approach, and the one Hunt seeks to improve upon

### Standard Lindblad (Born + Markov + Secular + Dirac Basis)
- Loses gauge invariance and thermodynamic consistency for driven systems
- Positive by construction
- Loses non-secular accuracy
- The most commonly used approach in quantum computing and quantum optics

The ordering reveals the structure: the standard theories discard physics at *both* layers (wrong basis *and* aggressive structural approximations). The nonadiabatic framework fixes layer 1 first, then applies structural approximations only as needed. Each approximation has a clear physical meaning and a clear cost.

---

## The Near-Term Path: Nonadiabatic CGME

For Hunt's immediate research program — few-level systems, classical driving fields, weak-to-moderate bath coupling — the most actionable step is the **nonadiabatic CGME**.

This requires:

1. **Compute the nonadiabatic coefficients** $b_k(t)$ for the driven system, using the Landau-Lifshitz integration-by-parts method (already established in the Mandal-Hunt program).

2. **Construct the Redfield tensor in the nonadiabatic basis.** Evaluate the bath correlation functions at the Bohr frequencies of the system, using the nonadiabatic populations and coherences as the density matrix elements that the tensor acts on. This is the step Jovanovski, Mandal, and Hunt (2023) performed.

3. **Apply the coarse-graining procedure** to the nonadiabatic Redfield tensor. Average over a time window $\Delta\tau$ chosen to be longer than the inverse of the relevant Bohr frequency differences but shorter than the relaxation timescale. This brings the equation into Lindblad form, guaranteeing positivity, while retaining the non-secular terms for nearly degenerate transitions.

4. **Verify the four criteria.** Check gauge invariance (compare length and velocity gauge results), thermodynamic consistency (verify zero heat flow for constant field, correct thermal steady state), complete positivity (verify non-negative eigenvalues of $\sigma_{\text{nad}}$ at all times), and non-secular accuracy (compare against full nonadiabatic Redfield for systems with closely spaced levels).

This is a concrete, publishable calculation. It uses tools that already exist (the nonadiabatic coefficients from the Mandal-Hunt program, the coarse-graining procedure from Schaller and Brandes), combines them in a way that has not been done before, and produces a master equation that satisfies all four criteria within a well-defined validity regime.

**A note on kernel extraction: transfer tensors vs. continuous kernels.** If the nonadiabatic memory kernel is to be extracted numerically from short-time exact simulations (rather than derived analytically via the Born approximation), the extraction method matters. The standard approach — inverting the Volterra integral equation that defines the continuous kernel $\mathcal{K}_{\text{nad}}(t, \tau)$ — is numerically ill-conditioned. Small errors in the short-time data are amplified by the inversion, and the resulting kernel can produce positivity violations when fed into the GQME, even if the underlying dynamics are exact. An alternative is the **transfer tensor method** (TTM): instead of extracting a continuous kernel, one extracts discrete dynamical maps $\mathcal{T}_k$ from the short-time simulation, where $\mathcal{T}_k$ represents the bath's influence at lag $k\Delta t$. Because each transfer tensor is derived from exact unitary evolution followed by a partial trace, it is inherently a completely positive map. The reduced dynamics are then propagated by summing the transfer tensors: $\rho(t_n) = \sum_{k=1}^{n} \mathcal{T}_k \rho(t_{n-k})$. This discrete approach avoids the Volterra inversion entirely and preserves complete positivity at each step, making it a more robust numerical strategy for the nonadiabatic GQME pipeline — particularly when the kernel is extracted from noisy simulations or from quantum hardware.

---

## The Long-Term Framework: Exact Nonadiabatic GQME

The general framework — the exact nonadiabatic Nakajima-Zwanzig equation with the time-dependent nonadiabatic projection — provides the formal umbrella under which all of this sits. Its development requires addressing several open questions:

**The time-dependent projector formalism.** The Nakajima-Zwanzig equation with time-dependent projectors has been studied in the literature (Shibata, Hashitsume; Breuer and Petruccione discuss it), but applying it with the specific nonadiabatic projection defined by the Landau-Lifshitz decomposition is new. The additional term $\mathcal{A}(t) = \dot{\mathcal{P}}_{\text{nad}}(t)\rho_{\text{total}}$ needs to be evaluated explicitly and shown to correspond to the adiabatic dynamics.

**The nonadiabatic projection beyond perturbation theory.** The Landau-Lifshitz decomposition is perturbative — it relies on integrating by parts in Dirac's first-order equation. For strong driving, a non-perturbative definition of the nonadiabatic subspace is needed. Candidates include the instantaneous eigenstates of $H_0 + V(t)$ (the adiabatic states), Floquet states for periodic driving, or a variational definition that minimizes the nonadiabatic coupling. Each choice defines a different projection operator and a different GQME.

**Quantum evaluation of the memory kernel.** For large systems, the nonadiabatic memory kernel must be evaluated on a quantum computer. This requires quantum algorithms for preparing the initial state, evolving under the projected dynamics in the $\mathcal{Q}_{\text{nad}}$-subspace, and extracting the kernel from the resulting correlation functions. The existing quantum algorithms for GQME simulation (Wang et al., 2022; Ding et al., 2024) provide a starting point, but they use the standard projection and would need to be adapted to the nonadiabatic one. The transfer tensor method may be preferable to continuous kernel extraction in this setting, since discrete dynamical maps extracted from quantum simulations naturally preserve complete positivity and avoid the numerical instabilities of Volterra inversion.

**Benchmarking against numerically exact methods.** For small systems where HEOM or process tensor calculations are feasible, the nonadiabatic GQME (at various levels of truncation) should be compared against the exact dynamics. This would establish the accuracy of each level in the approximation hierarchy and determine whether the nonadiabatic projection provides advantages even when the full dynamics are available.

---

## What This Framework Does Not Do

This proposal does not solve the theory gap by producing a new closed-form equation. It solves it by *reframing the problem*: the correct general theory is an exact equation (the nonadiabatic Nakajima-Zwanzig GQME) that is BQP-class and from which all existing approximate theories can be derived as controlled limits. The "gap" is not that we lack the equation — the Nakajima-Zwanzig equation has existed since 1958. The gap was that the equation was being applied with the wrong projection operator, which corrupted the physical content of every approximation derived from it.

The nonadiabatic projection fixes the physical content. The existing approximation machinery (Born, Markov, secular, coarse-graining) provides the complexity reductions needed for classical computation. The BQP constraint tells us that no further reduction is possible without loss of physics. The framework is complete — what remains is implementation, verification, and extension.

---

## References

* Nakajima, S. (1958). On quantum theory of transport phenomena. *Prog. Theor. Phys.*, 20, 948–959.
* Zwanzig, R. (1960). Ensemble method in the theory of irreversibility. *J. Chem. Phys.*, 33, 1338–1341.
* Shibata, F., Takahashi, Y., & Hashitsume, N. (1977). A generalized stochastic Liouville equation. *J. Stat. Phys.*, 17, 171–187.
* Breuer, H.-P., & Petruccione, F. (2002). *The Theory of Open Quantum Systems*. Oxford University Press.
* Mandal, A., & Hunt, K. L. C. (2012). Adiabatic and nonadiabatic contributions to the energy of a system subject to a time-dependent perturbation: Complete separation and physical interpretation. *J. Chem. Phys.*, 137, 164109.
* Mandal, A., & Hunt, K. L. C. (2016). Gauge-invariant expectation values of the energy of a molecule in an electromagnetic field. *J. Chem. Phys.*, 144, 044109.
* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields: Dephasing and population relaxation due to contact with a bath. *J. Chem. Phys.*, 158, 164107.
* Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Physical Review A*, 78, 022106.
* Campaioli, F., Cole, J. H., & Hapuarachchi, H. (2024). A tutorial on quantum master equations. *PRX Quantum*, 5, 020202.
* Wang, Y., et al. (2022). Simulating open quantum system dynamics on NISQ computers with generalized quantum master equations. arXiv:2209.04956.
* Ding, Z., et al. (2024). Simulating open quantum systems using Hamiltonian simulations. *PRX Quantum*, 5, 020332.