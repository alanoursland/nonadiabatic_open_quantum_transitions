# Open Technical Issues in the Nonadiabatic Nakajima-Zwanzig Framework

## Purpose

This document identifies the specific mathematical and physical gaps in the proposed nonadiabatic Nakajima-Zwanzig (NNZ) framework for driven, dissipative quantum systems. The framework proposes to close the theory gap — the absence of a master equation simultaneously satisfying complete positivity, non-secular accuracy, gauge invariance, and thermodynamic consistency — by replacing the standard projection operator in the Nakajima-Zwanzig GQME with a time-dependent nonadiabatic projection defined by the Mandal-Hunt decomposition. The proposal then derives the nonadiabatic CGME for few-level systems as a concrete, near-term deliverable, and frames the general theory within BQP complexity constraints.

The framework's overall architecture is sound. The identification of the projection operator as the source of the Dirac-coefficient pathologies is a genuine insight, and the approximation hierarchy is well-structured. But several claims require stronger justification, and at least one — the complete positivity of the exact nonadiabatic GQME — may be incorrect as stated. Each gap is described below with its scope, its consequences if unresolved, and suggested paths to resolution.

---

## Gap 1: The Nonadiabatic Projection Operator Is Not Constructively Defined

### The issue

The standard Nakajima-Zwanzig projection operator has a simple, constructive definition: trace over the bath, tensor with the bath equilibrium state. Written as a map on the total density matrix:

$$\mathcal{P}\rho_{\text{total}} = \text{Tr}_B[\rho_{\text{total}}] \otimes \rho_B^{\text{eq}}$$

This map is linear in $\rho_{\text{total}}$ and manifestly idempotent: applying it twice gives the same result as applying it once, because $\text{Tr}_B[\sigma \otimes \rho_B^{\text{eq}}] = \sigma$.

The proposed nonadiabatic projection is defined by its output:

$$\mathcal{P}_{\text{nad}}(t)\rho_{\text{total}}(t) = \sigma_{\text{nad}}(t) \otimes \rho_B^{\text{eq}}$$

where $\sigma_{\text{nad}}$ is the reduced density matrix constructed from the nonadiabatic coefficients $b_k(t)$ rather than the Dirac coefficients $c_k(t)$. But the document does not provide a constructive definition of the map $\mathcal{P}_{\text{nad}}(t)$ — the explicit procedure that takes an arbitrary $\rho_{\text{total}}$ and produces $\sigma_{\text{nad}} \otimes \rho_B^{\text{eq}}$.

### Why it matters

For the Nakajima-Zwanzig formalism to apply, $\mathcal{P}_{\text{nad}}$ must satisfy:

1. **Linearity.** $\mathcal{P}_{\text{nad}}(\alpha \rho_1 + \beta \rho_2) = \alpha \mathcal{P}_{\text{nad}}(\rho_1) + \beta \mathcal{P}_{\text{nad}}(\rho_2)$.
2. **Idempotency.** $\mathcal{P}_{\text{nad}}^2 = \mathcal{P}_{\text{nad}}$.

The first condition is likely satisfied: the adiabatic coefficients $a_k(t)$ are determined by the instantaneous Hamiltonian $V(t)$, not by the state, so subtracting them from $c_k(t)$ is a state-independent operation. The map "trace over bath, then subtract the adiabatic part" is linear.

The second condition is more delicate. The first application of $\mathcal{P}_{\text{nad}}$ produces $\sigma_{\text{nad}} \otimes \rho_B^{\text{eq}}$, which already has the adiabatic part removed. The second application traces over the bath (yielding $\sigma_{\text{nad}}$), then must "subtract the adiabatic part" from something that has no adiabatic part to subtract. The operation must be defined so that this second subtraction is a no-op — i.e., the projection recognizes that its input is already in the nonadiabatic subspace.

This can probably be arranged (for instance, by defining the projection in terms of a decomposition of the system Hilbert space into adiabatic and nonadiabatic subspaces, with $\mathcal{P}_{\text{nad}}$ projecting onto the latter), but the construction must be made explicit. Until it is, the formal exactness claim rests on an incompletely specified mathematical object.

### Path to resolution

Define $\mathcal{P}_{\text{nad}}(t)$ constructively as a two-step operation: (a) the standard partial trace $\text{Tr}_B$, followed by (b) a linear map $\mathcal{S}(t)$ on the system's density matrix space that projects onto the nonadiabatic subspace. The map $\mathcal{S}(t)$ must satisfy $\mathcal{S}^2 = \mathcal{S}$, which requires identifying the adiabatic and nonadiabatic subspaces as complementary image and kernel of a well-defined projector on the system's operator space. The Landau-Lifshitz decomposition at each order of perturbation theory provides a candidate for this projector, but it needs to be formulated as an operator-space projection rather than a coefficient-level subtraction.

---

## Gap 2: Complete Positivity of the Exact Nonadiabatic GQME

### The issue

The proposal claims that the exact nonadiabatic GQME is automatically completely positive, citing the Stinespring dilation theorem: the exact reduced dynamics of a subsystem of a unitarily evolving composite system is always completely positive.

This argument is correct for the *standard* reduced density matrix $\sigma(t) = \text{Tr}_B[\rho_{\text{total}}(t)]$. The map $\rho_{\text{total}}(0) \mapsto \text{Tr}_B[U(t)\rho_{\text{total}}(0)U^\dagger(t)]$ is CPTP by Stinespring.

But the nonadiabatic density matrix $\sigma_{\text{nad}}(t)$ is not the output of a partial trace. It is the output of a partial trace *followed by the subtraction of the adiabatic part*:

$$\sigma_{\text{nad}}(t) = \sigma(t) - \sigma_{\text{ad}}(t)$$

where $\sigma_{\text{ad}}(t)$ is the component of the reduced density matrix corresponding to the adiabatic polarization. The map from $\rho_{\text{total}}(0)$ to $\sigma_{\text{nad}}(t)$ is therefore:

$$\rho_{\text{total}}(0) \mapsto \text{Tr}_B[U(t)\rho_{\text{total}}(0)U^\dagger(t)] - \sigma_{\text{ad}}(t)$$

If $\sigma_{\text{ad}}(t)$ is independent of $\rho_{\text{total}}(0)$ — which it is, since the adiabatic coefficients depend on $V(t)$ and $H_0$, not on the initial state — then this map is *affine*, not linear. Specifically, it is a CPTP map plus a constant shift. Affine maps do not generally preserve complete positivity. A constant subtraction from a valid density matrix can produce a matrix with negative eigenvalues.

### Why it matters

This is not a marginal issue. The formal claim that the exact NNZ framework satisfies all four criteria simultaneously — and therefore closes the theory gap at the exact level — depends on complete positivity holding without approximation. If complete positivity holds only after the CGME coarse-graining is applied, then the framework's structure changes: three criteria (gauge invariance, thermodynamic consistency, non-secular accuracy) are satisfied exactly, and the fourth (complete positivity) is restored by approximation. This is still a major improvement over the status quo, but it is a different claim.

### Possible resolutions

There are several ways this could go:

**The claim might be salvageable.** If $\sigma_{\text{nad}}(t)$ can be shown to be positive semidefinite for all physical initial states (not all conceivable initial states), the positivity issue may not arise in practice. The nonadiabatic coefficients are squared moduli ($|b_k|^2 \geq 0$) and sum to at most 1 (since $|b_0|^2 + \sum_{k>0} |b_k|^2 \leq 1$ by the perturbative normalization). But positive semidefiniteness of the full matrix (including off-diagonal elements) requires more: the coherences $b_k b_l^*$ must be consistent with a valid density matrix. This needs to be checked.

**The claim might need to be weakened.** The exact NNZ equation could be described as generating dynamics that are gauge-invariant, thermodynamically consistent, and non-secular, with complete positivity guaranteed at the approximate (CGME) level. The theory gap would then be closed for the Born-Markov regime (where the CGME applies) rather than in full generality.

**The projection might need to be redefined.** Rather than subtracting the adiabatic part (an affine operation), the projection could be defined as a genuine linear projection onto a subspace — for instance, by working in a rotated basis where the adiabatic and nonadiabatic parts are orthogonal subspaces. This would require the adiabatic and nonadiabatic components to be orthogonal in the operator-space inner product, which connects to the Mandal-Hunt result that the cross-terms vanish in expectation values. If this orthogonality extends to the operator-space level, a true linear projection exists and the Stinespring argument applies.

This third option is the most promising and connects directly to Gap 1.

---

## Gap 3: Detailed Balance and the Thermal Target State

### The issue

The framework claims thermodynamic consistency: when the driving field is constant, the system relaxes to the correct thermal equilibrium, and the heat flow is zero. The zero-heat-flow claim is well-supported — it follows directly from $|b_k(t)|^2$ being constant during a plateau pulse. But the claim about the *correct thermal equilibrium* requires scrutiny.

In the three-level nonadiabatic CGME, the transition rates satisfy detailed balance:

$$\frac{\gamma_{l \to k}}{\gamma_{k \to l}} = e^{-\hbar\omega_{kl}/k_BT}$$

where $\omega_{kl} = (E_k - E_l)/\hbar$ are the Bohr frequencies of the *bare* Hamiltonian $H_0$. This means the system relaxes toward the Boltzmann distribution of $H_0$:

$$P_k^{\text{eq}} \propto e^{-E_k / k_BT}$$

But when a constant field $V$ is present, the correct thermal equilibrium is the Boltzmann distribution of the *dressed* Hamiltonian $H_0 + V$, whose eigenvalues $E_k'$ differ from the bare values $E_k$ by field-induced shifts. The system should relax to:

$$P_k^{\text{eq}} \propto e^{-E_k' / k_BT}$$

### Why it matters

For weak fields (where $V$ is a small perturbation), the dressed and bare energies are close, and the error in the thermal target is small — second order in the coupling to the field. This is consistent with the Born-Markov regime, where the Redfield tensor itself is second-order in $H_{SB}$, so the error in the thermal target is the same order as the errors already present.

But the thermodynamic consistency claim is stated as a structural property of the framework, not as an approximation. If the claim is only valid to the same order as the Born-Markov approximation, it should be stated as such.

### Path to resolution

The cleanest resolution would be to evaluate the bath correlation functions at the *dressed* Bohr frequencies $\omega_{kl}' = (E_k' - E_l')/\hbar$ rather than the bare ones. This would produce detailed balance with respect to the dressed Hamiltonian. The practical challenge is that the dressed energies are time-dependent when the field varies, so the rates become time-dependent. For a constant field this is straightforward; for a time-varying field, it introduces additional complexity.

An alternative is to accept the bare-frequency detailed balance as a controlled approximation, consistent with the Born-Markov regime, and state the thermodynamic consistency claim accordingly: the framework eliminates the fictitious heat flow and the spurious oscillations that plague the Dirac-based approach, and it produces the correct thermal equilibrium to the same order of accuracy as the underlying Redfield theory.

---

## Gap 4: The Near-Degeneracy Problem

### The issue

The Landau-Lifshitz integration by parts separates the adiabatic and nonadiabatic contributions by dividing by the Bohr frequency $\omega_{k0} = (E_k - E_0)/\hbar$:

$$a_k^{(1)}(t) = \frac{\langle k | V(t) | 0 \rangle}{E_k - E_0} \, e^{i\omega_{k0} t}$$

When two energy levels are nearly degenerate ($E_k \approx E_0$, so $\omega_{k0} \to 0$), the adiabatic term diverges. The decomposition becomes numerically unstable and physically questionable — a nearly degenerate level is not well described as either "adiabatic polarization" or "nonadiabatic transition" because the distinction loses its meaning when the energy gap is smaller than the perturbation matrix element.

### Why it matters

Hunt's current systems have well-separated energy levels, so this is not an immediate practical obstacle. But the framework is proposed as general, and molecular systems with closely spaced electronic states (J-aggregates, photosynthetic complexes, conical intersection regions) are precisely the systems where non-secular accuracy matters most — and where the framework would be most valuable if it worked.

The proposal sketches a partitioned treatment: near-degenerate states would be treated as a coupled subspace (diagonalized exactly, without the integration by parts), while well-separated states would use the standard Landau-Lifshitz decomposition. This is conceptually sensible but not developed.

### What needs to be worked out

1. A precise criterion for "near-degenerate" — presumably $|E_k - E_l| \lesssim |\langle k | V | l \rangle|$ or $|E_k - E_l| \lesssim \hbar/\Delta\tau$ — that determines which states are treated as a coupled subspace.

2. The form of the nonadiabatic projection within the coupled subspace. If the integration by parts cannot be applied, the adiabatic/nonadiabatic decomposition must come from a different source — perhaps from the instantaneous eigenstates of $H_0 + V(t)$ restricted to the degenerate subspace, or from a variational criterion.

3. Verification that the partitioned treatment preserves the four criteria, particularly gauge invariance (which relies on the specific structure of the Landau-Lifshitz decomposition) and the completeness of the energy separation (which Mandal and Hunt proved for the standard decomposition, not the partitioned version).

---

## Gap 5: BQP-Completeness of the Specific Problem

### The issue

The complexity argument claims that computing transition probabilities for general driven, dissipative quantum systems is BQP-complete. This is used to establish that no classical master equation can be both exact and general, framing the existing approximations as controlled complexity reductions from BQP toward P.

The argument relies on the fact that simulating quantum dynamics under general Hamiltonians (including open systems under general Lindblad operators) is BQP-complete. The cited references (Ding et al. 2024, Trivedi et al. 2025) establish quantum advantage for open-system simulation, but the exact scope of the hardness results matters.

### What needs verification

1. **Is the hardness for general or physically structured Lindblad operators?** If BQP-completeness requires Lindblad operators that are computationally universal (capable of encoding arbitrary quantum circuits), but physical system-bath couplings are restricted to local or low-rank operators, then the physical problem might be in a lower complexity class. The argument would still hold for sufficiently general systems, but its applicability to the specific systems studied in chemical physics (few-body, local couplings, structured baths) would be weaker.

2. **Is the hardness for the full density matrix or for specific observables?** Computing the full $2^n \times 2^n$ density matrix is clearly BQP-hard. Computing a single transition probability $|b_k(t)|^2$ might be easier — it is a specific expectation value, and there are cases where individual expectation values can be estimated more efficiently than the full state. If the transition probability problem (as opposed to the full simulation problem) is not BQP-hard, the complexity constraint would need to be restated.

3. **The cited references.** The Ding et al. (2024) paper addresses Hamiltonian simulation of open systems via Stinespring dilation, establishing quantum speedups. The Trivedi et al. (2025) paper addresses analog simulation with noise, showing superpolynomial advantage assuming BQP $\neq$ BPP. Neither may state BQP-completeness of Lindblad simulation in exactly the form the argument requires. The claim should be checked against the precise theorem statements.

### Consequences

If the BQP-completeness claim holds as stated, the complexity framing is correct and powerful. If it holds only for general (non-physical) Lindblad operators, the framing is still useful as a limiting argument but cannot be cited as a proof that the physical transition probability problem is BQP-hard. The practical consequences for the near-term program (few-level nonadiabatic CGME) are nil — the complexity argument is about scaling, not about small instances — but the claim's role as the conceptual foundation of the framework means it should be stated precisely.

---

## Gap 6: The Coarse-Grained Nonadiabatic Redfield Tensor

### The issue

The three-level document applies the Schaller-Brandes coarse-graining procedure to the nonadiabatic Redfield tensor, producing a Lindblad-form master equation. The sinc-damping of each tensor element,

$$R_{abcd}^{\text{CG}} = R_{abcd} \cdot \text{sinc}\left(\frac{(\omega_{ab} - \omega_{cd})\Delta\tau}{2}\right)$$

is claimed to produce a completely positive dynamical map.

### What needs verification

The Schaller-Brandes proof that coarse-graining yields Lindblad form relies on the positive semidefiniteness of a specific rate matrix constructed from the Redfield tensor elements and the sinc factors. This positive semidefiniteness follows from the structure of the bath correlation functions (specifically, the positivity of the bath spectral density $J(\omega) \geq 0$ for $\omega > 0$).

The proof was developed for the standard Redfield tensor acting on the standard density matrix. The nonadiabatic Redfield tensor uses the same bath correlation functions and the same system-bath coupling matrix elements, but it acts on $\sigma_{\text{nad}}$ rather than $\sigma$. The question is: does the coarse-graining proof depend only on the tensor's structure (which is the same) or also on properties of the density matrix it acts on (which is different)?

The proof should depend only on the tensor structure, since the Lindblad form is a property of the generator, not of the state. But this should be verified explicitly, particularly for the 3-level case where the verification is a finite computation.

### Path to resolution

For the 3-level system, the coarse-grained rate matrix is a finite matrix whose positive semidefiniteness can be checked numerically for any specific choice of bath parameters and coarse-graining timescale. A general proof (for arbitrary $n$) would require showing that the nonadiabatic Redfield tensor inherits the structural properties used in the Schaller-Brandes proof, which reduces to confirming that the bath correlation functions enter the tensor in the same way regardless of whether the tensor acts on $\sigma$ or $\sigma_{\text{nad}}$.

---

## Summary Table

| Gap | Severity | Near-term impact | Resolution difficulty |
| :--- | :--- | :--- | :--- |
| 1. Projection operator definition | High (formal) | Blocks rigorous NNZ derivation | Moderate — requires operator-space formulation |
| 2. Complete positivity of exact NNZ | High (conceptual) | Changes the nature of the claim | Hard — may require redefining the projection |
| 3. Detailed balance / thermal target | Moderate | Affects precision of thermodynamic consistency claim | Moderate — dressed frequencies or explicit order statement |
| 4. Near-degeneracy | Low (near-term) | Does not affect Hunt's systems | Hard — requires partitioned decomposition theory |
| 5. BQP-completeness scope | Low-Moderate | Affects framing, not results | Moderate — requires careful reference verification |
| 6. Coarse-grained tensor validity | Low | Likely fine, needs explicit check | Easy — finite computation for 3-level case |

---

## Relationship Between the Gaps

Gaps 1 and 2 are entangled. If the nonadiabatic projection can be defined as a genuine linear projector on the system's operator space (Gap 1), then the question of whether the resulting dynamical map is CPTP (Gap 2) becomes a standard question about projected dynamics, and the Stinespring argument may apply directly. The key is whether the adiabatic and nonadiabatic parts of the density matrix are orthogonal in the operator-space inner product (the Hilbert-Schmidt inner product $\text{Tr}[A^\dagger B]$). The Mandal-Hunt result — that cross-terms between adiabatic and nonadiabatic contributions vanish in energy expectation values — is suggestive but not identical to operator-space orthogonality. Establishing this orthogonality (or its failure) would resolve both gaps simultaneously.

Gap 3 is independent and can be addressed by a straightforward calculation: compare the thermal steady states of the nonadiabatic CGME with bare versus dressed Bohr frequencies, for a concrete 3-level system with a constant field, and quantify the difference.

Gap 4 is a longer-term issue that becomes relevant only when the framework is extended beyond Hunt's current systems. It can be deferred without affecting the near-term deliverables.

Gap 5 is a framing issue. The complexity argument's role is to motivate and constrain, not to prove. If the precise BQP-completeness statement needs qualification, the overall logic of the proposal survives — it just needs to be stated more carefully.

Gap 6 is the easiest to resolve and should be done first, since it directly validates the concrete 3-level calculation.