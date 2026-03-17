# Positivity vs. Physics: The Central Tension in Open Quantum Systems

## The Tradeoff

The theory of open quantum systems is built on a tension that runs through every file in this collection. On one side is the demand for **physical accuracy**: the master equation should capture the real dynamics of the system, including coherence between nearly degenerate states, the correct thermal equilibrium, and gauge-invariant transition rates. On the other side is the demand for **mathematical consistency**: the density matrix must remain a valid quantum state at all times — Hermitian, trace-preserving, and completely positive. No negative probabilities, ever.

The standard tools of the field each satisfy one demand at the expense of the other. The Redfield equation, derived microscopically from the system-bath coupling, retains the coherence-transfer terms that describe interference between transitions at different frequencies. But it can produce density matrices with negative eigenvalues — a mathematical pathology that renders the results unphysical in principle, even if the violations are numerically small in practice. The Lindblad equation, obtained by applying the secular approximation (or postulated axiomatically from the requirement of complete positivity), guarantees a valid density matrix at all times. But it discards the non-secular terms, losing the coherence physics that drives the dynamics in systems with closely spaced energy levels.

This tension is not merely academic. It determines which phenomena a simulation can predict, which experimental signatures a model can reproduce, and which theoretical claims about driven, dissipative quantum systems can be trusted. The question that animates the research discussed in this collection is whether the tension is fundamental or whether it can be resolved — and if so, what the resolution requires.

---

## Four Criteria for a Correct Theory

The work of Mandal, Hunt, Jovanovski, and others, together with the broader developments in open quantum systems theory, points toward four criteria that a fully satisfactory master equation must satisfy simultaneously. Each criterion is established in detail in the other files of this collection; here they are assembled as a unified standard.

### I. Complete Positivity

The density matrix $\rho(t)$ must satisfy $\rho(t) \geq 0$ and $\text{Tr}[\rho(t)] = 1$ at all times. More precisely, the dynamical map must be **completely positive** — it must preserve the positivity of the density matrix even when the system is entangled with an ancilla that is not being evolved. This is the mathematical content of the Lindblad theorem (discussed in the Lindblad/GKSL file): the GKSL form is the unique Markovian master equation structure that guarantees this property.

**Where the standard approaches stand.** The Lindblad equation satisfies this by construction. The Redfield equation does not — its non-secular terms can produce small negative eigenvalues during transients, at low temperatures, or for nearly degenerate systems. The positivity problem of Redfield theory (discussed in the Redfield file) is the mathematical half of the central tension.

**What a resolution requires.** Any candidate theory must either take the Lindblad form or provide an alternative guarantee of complete positivity. The coarse-grained master equation (CGME) and the geometric-arithmetic master equation (GAME), discussed in the beyond-Redfield file, achieve this by different mathematical routes — time-averaging and rate interpolation, respectively — while retaining more of the Redfield tensor's structure than the secular approximation allows.

### II. Gauge Invariance

The physical predictions of the theory — transition probabilities, steady-state populations, heat flow — must not depend on the mathematical representation of the electromagnetic field. The same physics described in the length gauge ($-\mathbf{d} \cdot \mathbf{E}$) and the velocity gauge ($-\mathbf{p} \cdot \mathbf{A}/m$) must give the same results.

**Where the standard approaches stand.** For an undriven system, gauge invariance is not an issue — there is no external field to represent. The problem arises when a driving field is added to either a Redfield or Lindblad master equation whose dissipative terms are built from the eigenstates of the field-free Hamiltonian $H_0$. In that setting, Dirac's transition coefficients $|c_k(t)|^2$ are gauge-dependent while the field is active (as discussed in the Dirac file and the Landau-Lifshitz file). Since the dissipative terms implicitly use these coefficients, the resulting master equation inherits the gauge dependence. This is not a failure of Redfield or Lindblad theory per se — it is a failure of the basis choice used to couple the driven dynamics to the dissipation.

**What a resolution requires.** The nonadiabatic coefficients $|b_k(t)|^2$ are gauge-invariant at all times, as demonstrated by Mandal and Hunt (discussed in the Landau-Lifshitz and Hunt corrections files). Using these as the basis for the master equation removes the gauge artifact without modifying the master equation's mathematical structure.

### III. Thermodynamic Consistency

If the driving field is constant (or absent), the system must relax to the correct Boltzmann distribution. The first law of thermodynamics — the power absorbed from the field equals the rate of change of the system's energy plus the heat dissipated to the bath — must hold at every instant, not just on average or in the long-time limit. There must be no fictitious steady-state heat flow.

**Where the standard approaches stand.** Both Redfield and Lindblad satisfy detailed balance and produce the correct thermal equilibrium *for an undriven system*. The thermodynamic failures arise specifically when a driving field is active and Dirac's coefficients are used to define the system populations. The three specific pathologies — wrong equilibrium, fictitious heat flow, and dephasing-induced spurious transitions — are the central topic of the nonadiabatic bath coupling file, and the power absorption argument that diagnoses them is developed in the Hunt corrections file.

**What a resolution requires.** The nonadiabatic decomposition, applied to the master equation as proposed by Jovanovski, Mandal, and Hunt (2023), eliminates all three pathologies. The bath acts on genuine excitations rather than on the adiabatic polarization, and the resulting energy balance is exact.

### IV. Non-Secular Accuracy

The theory must retain the coherence-transfer terms that the secular approximation discards — the terms that couple density matrix elements oscillating at different Bohr frequencies. These terms describe how the relaxation of one transition is influenced by another when the two transitions are nearly resonant, and they are essential for systems with closely spaced energy levels.

**Where the standard approaches stand.** The Redfield equation retains these terms; the secular Lindblad equation discards them. This is the physical half of the central tension — the accuracy that is lost when positivity is enforced by the secular approximation. The specific physics at stake (coherence transfer, transport efficiency in molecular aggregates, NMR line shapes) is discussed in both the Redfield and Lindblad files.

**What a resolution requires.** Any method that achieves complete positivity without the secular approximation — CGME, GAME, or a nonadiabatic master equation that inherits Lindblad structure from coarse-graining — satisfies this criterion. The non-secular terms are retained because the averaging timescale is finite rather than infinite.

---

## Where the Field Stands

Stating the four criteria is easier than satisfying all of them simultaneously. Here is an honest assessment of the current landscape.

### What Has Been Demonstrated

The nonadiabatic decomposition (criteria II and III) is on firm ground. The Mandal-Hunt program, spanning a decade of publications from 2012 to 2023, has established through formal proofs and quantitative calculations that:

The nonadiabatic coefficients $|b_k(t)|^2$ are gauge-invariant, produce the correct power absorption, give a complete and cross-term-free energy separation, are statistically consistent with the quantum variance, and — when used in a master equation — yield thermodynamically consistent dynamics including the correct thermal steady state. These results are demonstrated for few-level systems driven by classical fields within time-dependent perturbation theory.

The CGME and GAME (criterion I + IV) independently address the positivity-versus-accuracy tradeoff. They produce completely positive master equations that retain non-secular coherence-transfer terms, and they have been benchmarked against exact solutions (Jaynes-Cummings models, spin chains) with favorable results.

### What Has Not Yet Been Demonstrated

The full synthesis — a single master equation that satisfies all four criteria simultaneously — is a goal rather than an accomplished fact. Specifically:

**The nonadiabatic decomposition has not been formally integrated into the CGME or GAME frameworks.** The Mandal-Hunt work modifies the basis (which coefficients enter the master equation), while the CGME/GAME modify the mathematical structure (how the dissipative terms are constructed from the Redfield tensor). Combining both modifications — using nonadiabatic coefficients as input to a coarse-grained or geometric-arithmetic master equation — is a natural next step, but it has not been carried out and published in that form.

**The extension to non-perturbative regimes is open.** The nonadiabatic decomposition is developed within the perturbative framework of Dirac's time-dependent theory (extended to higher orders by Mandal and Hunt). In the strong-coupling, non-Markovian regimes handled by HEOM and process tensor methods, the question of how to separate adiabatic and nonadiabatic contributions has not been addressed. It is possible that these numerically exact methods automatically avoid the thermodynamic pathologies by treating the system-bath correlations exactly, making the decomposition unnecessary. It is also possible that the issue persists in a different form. This is not yet known.

**The extension to many-body and quantized-field settings is open.** The current framework applies to few-level systems driven by classical electromagnetic fields. Extending it to many-body systems with internal interactions, to quantized radiation fields, and to structured (non-Markovian) baths are all active areas where the framework's applicability has not been established.

### The Philosophical Point

Despite these open questions, the Mandal-Hunt program has established something important: the tension between positivity and physics is not the only tension, and it may not be the most fundamental one. Before worrying about whether the master equation is in Lindblad form, one must first ask whether the *quantities entering the master equation* are physically meaningful. If the transition probabilities are wrong — if they include the adiabatic polarization that should not be subject to dissipation — then no amount of mathematical sophistication in the master equation structure can produce physically correct results. Fixing the math (positivity) requires first fixing the physics (the definition of a transition).

This is the central lesson of the collection: the four criteria are not independent. Gauge invariance and thermodynamic consistency (criteria II and III) are prerequisites that constrain what can serve as input to a master equation. Complete positivity and non-secular accuracy (criteria I and IV) are properties of the master equation structure itself. A correct theory must address both levels — the input and the structure — and the standard approach of inserting Dirac's coefficients into a Lindblad or Redfield equation fails at the first level before the second level is even relevant.

---

## The Corrected Comparison

| Criterion | Redfield | Secular Lindblad | CGME / GAME | Nonadiabatic + CGME (proposed) |
| :--- | :--- | :--- | :--- | :--- |
| **Complete positivity** | Violated in some regimes | Guaranteed | Guaranteed | Expected (not yet demonstrated) |
| **Gauge invariance** | Fails for driven systems (Dirac basis) | Fails for driven systems (Dirac basis) | Fails for driven systems (Dirac basis) | Satisfied (nonadiabatic basis) |
| **Thermodynamic consistency** | Violated for driven systems | Violated for driven systems | Violated for driven systems | Satisfied |
| **Non-secular accuracy** | Retained | Lost | Retained | Expected to be retained |

The critical observation is that the first three columns all share the same failure in rows 2 and 3: they use Dirac's coefficients. The CGME and GAME fix the positivity and coherence problems but inherit the gauge and thermodynamic problems from the Dirac basis. Only the combination of the nonadiabatic decomposition (fixing the input) with a positivity-preserving structure like the CGME (fixing the math) is expected to satisfy all four criteria. This combination is the frontier.

---

## References

* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields: Dephasing and population relaxation due to contact with a bath. *J. Chem. Phys.*, 158, 164107.
* Mandal, A., & Hunt, K. L. C. (2012–2021). Series of papers on nonadiabatic transition probabilities. *J. Chem. Phys.* [see Landau-Lifshitz and Hunt corrections files for complete list].
* Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Physical Review A*, 78, 022106.
* Campaioli, F., Cole, J. H., & Hapuarachchi, H. (2024). A tutorial on quantum master equations. *PRX Quantum*, 5, 020202.
* Lindblad, G. (1976). On the generators of quantum dynamical semigroups. *Commun. Math. Phys.*, 48, 119–130.
* Redfield, A. G. (1957). On the theory of relaxation processes. *IBM J. Res. Dev.*, 1, 19–31.