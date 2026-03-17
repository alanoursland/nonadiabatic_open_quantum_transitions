# The Mandal-Hunt Results: Energy, Power, and Variance in the Nonadiabatic Framework

## Scope of This Document

The Landau-Lifshitz separation splits Dirac's excited-state coefficients $c_k(t)$ into an adiabatic term $a_k(t)$ (reversible polarization) and a nonadiabatic term $b_k(t)$ (genuine transitions). The nonadiabatic bath coupling notes explain why this distinction matters when the system is coupled to a thermal environment. This document covers the middle ground: the specific quantitative results, derived by Anirban Mandal and Katharine L. C. Hunt in a series of papers beginning in 2012, that demonstrate *why* the nonadiabatic transition probability $|b_k(t)|^2$ is the physically correct quantity — and where, precisely, Dirac's $|c_k(t)|^2$ goes wrong.

The arguments fall into three categories: energy balance and power absorption, the complete separation of the energy into adiabatic and nonadiabatic parts, and the statistics of energy fluctuations. Each provides an independent line of evidence that $|b_k(t)|^2$, not $|c_k(t)|^2$, is the physical transition probability.

---

## Power Absorption: The Energy-Work Theorem Test

The most direct way to test whether a proposed transition probability is physically correct is to check whether it gives a consistent energy balance. Mandal and Hunt (2015) used this approach to expose a concrete failure of Dirac's coefficients.

### The Exact Power

The average energy of the system at time $t$ is the expectation value of the full Hamiltonian:

$$\langle E(t) \rangle = \langle \Psi(t) | H(t) | \Psi(t) \rangle$$

The rate of change of this energy — the instantaneous power absorbed from the field — follows from the time-dependent Schrödinger equation. Using $i\hbar \, \partial_t |\Psi\rangle = H|\Psi\rangle$, the commutator terms cancel and the result is:

$$P(t) = \frac{d}{dt}\langle E(t) \rangle = \left\langle \Psi(t) \left| \frac{\partial H(t)}{\partial t} \right| \Psi(t) \right\rangle$$

Since the unperturbed Hamiltonian $H_0$ is time-independent, this reduces to:

$$P(t) = \left\langle \Psi(t) \left| \frac{\partial V(t)}{\partial t} \right| \Psi(t) \right\rangle$$

This is exact. It says something physically obvious but mathematically constraining: the system gains energy only when the external field is changing. A constant field, no matter how strong, does zero work on the system. The power depends on $\partial V / \partial t$ and nothing else.

### The Dirac-Based Power

Now consider what happens if you try to reconstruct the power from Dirac's transition probabilities. The natural approach is to say: the system occupies unperturbed state $|n\rangle$ with probability $|c_n(t)|^2$, and each such state has energy $E_n$, so the average energy is $\sum_n E_n |c_n(t)|^2$ and the power is:

$$P_{\text{Dirac}}(t) = \sum_n E_n \frac{d}{dt}|c_n(t)|^2$$

Mandal and Hunt showed that this expression does *not* equal the exact power $P(t)$. The discrepancy takes the form:

$$P_{\text{Dirac}}(t) = P(t) + \text{oscillatory terms}$$

The oscillatory terms arise from the interference between the adiabatic and nonadiabatic components within Dirac's coefficients. They do not vanish when the field becomes constant. This means that $P_{\text{Dirac}}(t)$ predicts the system is absorbing energy from a static field — a clear violation of energy conservation. A constant perturbation cannot do work, yet the Dirac-based accounting says it does.

### The Nonadiabatic Power

When the same calculation is performed using the nonadiabatic coefficients $b_k(t)$ instead of $c_k(t)$, the oscillatory terms cancel exactly. The power computed from the nonadiabatic transition probabilities matches the exact expression:

$$\sum_n E_n \frac{d}{dt}|b_n(t)|^2 = \left\langle \Psi(t) \left| \frac{\partial V(t)}{\partial t} \right| \Psi(t) \right\rangle$$

More precisely, Mandal and Hunt proved that the power absorbed from the electromagnetic field equals the time derivative of the nonadiabatic contribution to the energy. When the field is constant, $|b_k(t)|^2$ is constant, and the power is zero — exactly as required.

This result was first established in the Coulomb gauge (Mandal & Hunt, 2015) and later shown to be gauge-invariant (Mandal & Hunt, 2016), holding identically in both the length gauge and the velocity gauge.

---

## Complete Energy Separation

The power absorption result is a statement about time derivatives. Mandal and Hunt (2012) proved a stronger result about the energy itself: the total energy of the system separates completely into adiabatic and nonadiabatic components at every instant, with no cross-terms.

### The Structure of the Result

For a system initially in the ground state $|0\rangle$ of $H_0$, perturbed by $V(t)$, the energy at time $t$ can be written as:

$$\langle E(t) \rangle = E_{\text{ad}}(t) + E_{\text{nad}}(t)$$

The **adiabatic energy** $E_{\text{ad}}(t)$ is identical to the energy that would be obtained from time-independent perturbation theory applied to the instantaneous Hamiltonian $H_0 + V(t)$. It is the energy of a system that has perfectly tracked the field — adjusted its ground state to the current perturbation without any real excitations. This is the energy of polarization.

The **nonadiabatic energy** $E_{\text{nad}}(t)$ is a sum over excited states:

$$E_{\text{nad}}(t) = \sum_{k \neq 0} |b_k(t)|^2 \, (E_k - E_0)$$

(with perturbative corrections at higher orders). It is the energy stored in genuine transitions — states that have been truly populated by the time-varying field.

### Why the Cross-Terms Vanish

The nontrivial part of the proof is showing that the cross-terms between $a_k(t)$ and $b_k(t)$ vanish when the energy is computed as an expectation value. This is not obvious, because the wavefunction contains both adiabatic and nonadiabatic amplitudes, and in general the expectation value of an operator applied to a sum of amplitudes produces cross-terms. Mandal and Hunt demonstrated — explicitly through third order in the perturbation parameter — that these cross-terms cancel identically when the full structure of the perturbation expansion is respected.

The physical significance is that the adiabatic and nonadiabatic contributions to the energy are independent. Polarization energy and excitation energy do not mix. This is what allows the nonadiabatic framework to provide clean energy accounting: the work done by the field goes entirely into the nonadiabatic energy (real transitions), while the adiabatic energy tracks the instantaneous polarization and returns to zero when the field is removed.

---

## Variance and Higher Moments: A Statistical Test

The energy separation tells us about the mean energy. Mandal and Hunt (2020) extended the analysis to the variance and higher moments of the energy distribution, producing what is perhaps the most elegant argument for $|b_k(t)|^2$ as the correct transition probability.

### The Result

They proved that the variance of the energy depends entirely on the nonadiabatic transition probabilities:

$$\sigma^2(E) = \sum_{k \neq 0} |b_k(t)|^2 \, (E_k - E_0)^2 - \left[\sum_{k \neq 0} |b_k(t)|^2 \, (E_k - E_0)\right]^2$$

If the system's response is entirely adiabatic — if it perfectly tracks the field without any real excitations — the variance is zero. The system has a definite energy (the adiabatic energy), and measurements of the energy always return the same value. Fluctuations in the energy arise only from nonadiabatic transitions.

### The Statistical Consistency Argument

This result enables a decisive test. Consider the standard statistical framework: if the system occupies excited state $|k\rangle$ with probability $P_k$, and the excitation energy is $\Delta E_k = E_k - E_0$, then the variance should be:

$$\sigma^2 = \sum_k P_k \, (\Delta E_k)^2 - \left(\sum_k P_k \, \Delta E_k\right)^2$$

Mandal and Hunt showed that the quantum mechanical variance (computed exactly from the wavefunction) matches this standard statistical expression if and only if $P_k = |b_k(t)|^2$. Setting $P_k = |c_k(t)|^2$ produces a different, inconsistent result — the statistical formula gives the wrong variance.

The same holds for all higher moments of the energy distribution. The nonadiabatic transition probabilities are the only ones that make the quantum mechanical results consistent with a standard statistical interpretation of "the system is in state $|k\rangle$ with probability $P_k$."

### Numerical Illustration

The 2020 paper illustrates the discrepancy numerically for vibration-rotation transitions of HCl in the gas phase, showing that the variances computed from $|c_k(t)|^2$ and $|b_k(t)|^2$ differ substantially during and after a perturbing pulse. The "fictitious" fluctuations seen in the Dirac-based variance are artifacts of counting the adiabatic polarization as a real excitation — exactly the error that the Landau-Lifshitz separation is designed to remove.

---

## Behavior During Pulses: Quantitative Differences

The results above are formal proofs of consistency. But how large are the actual numerical differences between $|b_k(t)|^2$ and $|c_k(t)|^2$ in practice? Mandal and Hunt investigated this in two papers (2018) by computing both quantities for specific, experimentally relevant pulse shapes.

### Gaussian-Envelope Pulses

For a cosine wave of frequency $\omega$ in a Gaussian envelope, they found that the nonadiabatic and Dirac transition probabilities differ systematically depending on the detuning. When the carrier frequency exceeds the transition frequency ($\omega > \omega_{k0}$), the nonadiabatic probability $|b_k(t)|^2$ is larger than $|c_k(t)|^2$. When $\omega < \omega_{k0}$, the opposite is true. These differences are independent of the phase of the carrier wave relative to the envelope peak — they are structural, not accidental.

### Plateau Pulses

The differences are even starker for pulses that rise to a constant plateau and later fall off. During the plateau — while $V(t)$ is constant — the nonadiabatic transition probability $|b_k(t)|^2$ is constant, as it must be: a static perturbation cannot cause transitions. But Dirac's $|c_k(t)|^2$ continues to oscillate throughout the plateau. This is the Rabi-oscillation artifact visible in the bare-state basis, and it is the same behavior that leads to the fictitious power absorption described above.

### Overlapping Pulses and Persistent Differences

In 2021, Mandal and Hunt proved a result that goes beyond pulse-shape dependence: for overlapping pulse sequences — where a second pulse acts while the first is still on — the differences between $|b_k(t)|^2$ and $|c_k(t)|^2$ can persist even after both pulses have ended. This is a stronger statement than anything Landau and Lifshitz considered, since they assumed a single pulse that eventually returns to zero.

The mechanism involves dephasing. The first pulse populates an excited state and creates coherences. If dephasing occurs while the system is still being driven by the second pulse, the nonadiabatic framework handles the loss of coherence correctly — it affects only the phase relationships between genuinely populated states. In Dirac's framework, dephasing of the perturbed part of the wavefunction cannot be properly described, because the adiabatic and nonadiabatic components are lumped together.

For overlapping pulses with on-resonance frequencies, the lasting differences in calculated transition probabilities exceed 35%. For off-resonant perturbations, the differences are larger still.

---

## A Proposed Experimental Test

The theoretical differences between $|b_k(t)|^2$ and $|c_k(t)|^2$ are large enough to be measurable, and Mandal and Hunt (2018) proposed a specific experimental protocol to distinguish them.

The idea is a **ratio test**. Rather than measuring absolute transition probabilities — which would require precise knowledge of matrix elements and field intensities — they proposed measuring the ratio of transition probabilities at two different detunings from resonance. The Dirac and nonadiabatic frameworks predict different ratios, and the measurement depends only on relative populations, which can be extracted from spectroscopic signals without absolute calibration.

This proposal is significant because it converts a theoretical disagreement about the correct definition of transition probability into a concrete, falsifiable experimental prediction. As of their most recent publications, this experiment has not yet been performed, but the predicted differences are well within the resolution of modern ultrafast spectroscopy.

---

## Summary: Three Independent Consistency Tests

The Mandal-Hunt results provide three independent arguments that $|b_k(t)|^2$ is the physical transition probability and $|c_k(t)|^2$ is not:

**Energy balance.** The power absorbed from the field, computed from the exact Schrödinger dynamics, matches the nonadiabatic energy rate. The Dirac-based power contains spurious oscillatory terms that predict energy absorption from a static field.

**Energy separation.** The total energy decomposes into adiabatic and nonadiabatic parts with no cross-terms. The nonadiabatic energy has the structure $\sum_k |b_k|^2 (E_k - E_0)$, confirming that $|b_k|^2$ plays the role of the transition probability in the energy accounting.

**Statistical consistency.** The variance and higher moments of the energy distribution are consistent with standard statistical analysis only if $P_k = |b_k(t)|^2$. Using $P_k = |c_k(t)|^2$ gives the wrong fluctuations.

Each of these tests is independent: the power argument is about time derivatives, the energy separation is about expectation values, and the variance argument is about fluctuations. The fact that all three point to the same conclusion — and that the nonadiabatic framework additionally respects gauge invariance and the adiabatic limit — constitutes a strong case that Dirac's transition probability, while correct as a projection onto the unperturbed basis, is not the right quantity for describing the physical process of excitation.

---

## References

* Mandal, A., & Hunt, K. L. C. (2012). Adiabatic and nonadiabatic contributions to the energy of a system subject to a time-dependent perturbation: Complete separation and physical interpretation. *J. Chem. Phys.*, 137, 164109. https://doi.org/10.1063/1.4750045
* Mandal, A., & Hunt, K. L. C. (2015). Non-adiabatic current densities, transitions, and power absorbed by a molecule in a time-dependent electromagnetic field. *J. Chem. Phys.*, 143, 034102. https://doi.org/10.1063/1.4923181
* Mandal, A., & Hunt, K. L. C. (2016). Gauge-invariant expectation values of the energy of a molecule in an electromagnetic field. *J. Chem. Phys.*, 144, 044109. https://doi.org/10.1063/1.4938564
* Mandal, A., & Hunt, K. L. C. (2018). Quantum transition probabilities during a perturbing pulse: Differences between the nonadiabatic results and Fermi's golden rule forms. *J. Chem. Phys.*, 148, 194107. https://doi.org/10.1063/1.5019172
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse: Toward experimental tests of the differences from Dirac's transition probabilities. *J. Chem. Phys.*, 149, 204110. https://doi.org/10.1063/1.5054313
* Mandal, A., & Hunt, K. L. C. (2020). Variance of the energy of a quantum system in a time-dependent perturbation: Determination by nonadiabatic transition probabilities. *J. Chem. Phys.*, 152, 104110. https://doi.org/10.1063/1.5140009
* Mandal, A., & Hunt, K. L. C. (2021). Quantum transition probabilities due to overlapping electromagnetic pulses: Persistent differences between Dirac's form and nonadiabatic perturbation theory. *J. Chem. Phys.*, 154, 024116. https://doi.org/10.1063/5.0020169