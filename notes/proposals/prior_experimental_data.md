# Prior Experimental Data for Testing the Nonadiabatic Framework

## Purpose

This document identifies published experimental data sets that could be used to test the nonadiabatic Nakajima-Zwanzig framework and its CGME approximation against real physical measurements, without requiring new laboratory experiments. The strategy is to find systems where (1) the Hamiltonian parameters are known independently of the data being compared against, (2) the system is simultaneously driven and dissipating, (3) the perturbative regime applies so that the Landau-Lifshitz decomposition is valid, and (4) the measured observable depends sensitively enough on transition probabilities to distinguish $|b_k(t)|^2$ from $|c_k(t)|^2$.

Each candidate is assessed for feasibility, theoretical match to the nonadiabatic framework, data availability, and the specific researchers who could be contacted for collaboration or data access.

---

## Why Published Data, Not New Experiments

The Mandal-Hunt program has established through internal consistency tests (power absorption, energy separation, variance) that $|b_k(t)|^2$ is the physically correct transition probability. But all of these tests are theory-against-theory comparisons: they show that the nonadiabatic coefficients are consistent with exact quantum mechanics in ways that Dirac's coefficients are not. What has not been done is a comparison against independent experimental measurements.

The proposed experimental ratio test (Mandal & Hunt, 2018) has not been performed. Designing and executing a purpose-built experiment requires laboratory resources and collaborators. However, decades of precision spectroscopy and relaxation measurements have produced published data sets with the right characteristics — known Hamiltonians, characterized driving fields, measured transition-dependent observables — that can be compared against nonadiabatic predictions computationally. If the nonadiabatic framework matches published data better than the standard Dirac-based framework using the same independently known parameters, that constitutes evidence from an independent source, arguably stronger than a purpose-built test because it rules out experimenter bias.

The key constraint: the comparison must use independently known parameters, not fitted ones. If the nonadiabatic framework introduces even one additional adjustable parameter, an improved fit proves nothing. The test is meaningful only if both frameworks use identical inputs and differ only in whether $|c_k(t)|^2$ or $|b_k(t)|^2$ enters the master equation.

---

## Candidate 1: NMR $R_{1\rho}$ Relaxation Dispersion

### Why This Is the Strongest Candidate

In a rotating-frame relaxation experiment ($R_{1\rho}$ or $T_{1\rho}$), a spin system is prepared in a transverse state and then subjected to a continuous RF spin-lock field while relaxing due to coupling to the molecular tumbling (the "bath"). This is a literal realization of the driven-dissipative scenario: the spin-lock field is the perturbation $V$, it is constant during the measurement (a plateau pulse), and the bath acts simultaneously.

The Dirac framework predicts that the spin-lock field mixes the spin states in a way that includes the adiabatic polarization response, and the bath relaxes all of it. The nonadiabatic framework predicts that the polarization response is not subject to relaxation — only genuine transitions contribute. The measured quantity, $R_{1\rho}$, depends on which picture is correct.

The critical advantage: $R_{1\rho}$ experiments can be performed with spin-lock field strengths as low as 25 Hz (demonstrated by Korzhnev, Orekhov, and Kay, JACS 2005), which is genuinely weak compared to the Larmor frequency (hundreds of MHz). This puts the system solidly in the perturbative regime where the Landau-Lifshitz decomposition applies. The relaxation rate depends on the spin-lock field strength, with the dependence being more pronounced at weak fields — exactly the regime where the nonadiabatic correction would be largest.

### What Makes NMR Ideal

The spin Hamiltonian parameters (chemical shifts, scalar couplings, dipolar couplings, chemical shift anisotropy tensors) are independently measurable to high precision from separate experiments. The bath is characterized by the rotational correlation time $\tau_c$ of the molecule, measurable from $T_1/T_2$ ratios. The spectral density function $J(\omega)$ is given by well-established models (Lipari-Szabo model-free formalism for proteins, or the simple $\tau_c / (1 + \omega^2 \tau_c^2)$ for rigid molecules). None of these parameters need to be fitted from the $R_{1\rho}$ data itself.

The standard theoretical prediction for $R_{1\rho}$ is derived from Redfield theory (with or without the secular approximation). For systems with nearly degenerate transitions (e.g., two spins with similar chemical shifts), the non-secular Redfield terms contribute to $R_{1\rho}$, and the CGME vs. secular Lindblad distinction matters. This is the regime where the nonadiabatic CGME would provide the most significant correction.

### Available Data and Contacts

**Lewis Kay's group (University of Toronto):** The leading NMR relaxation methodology group. They have published $R_{1\rho}$ dispersion data (relaxation rate as a function of spin-lock power) for numerous protein systems with fully characterized parameters. Key papers: Korzhnev, Orekhov, Kay (JACS 2005) on off-resonance $R_{1\rho}$ with low spin-lock fields; Korzhnev et al. (J. Biomol. NMR 2002) on accurate measurement of heteronuclear spin-lock relaxation rates.

**Arthur Palmer's group (Columbia University):** Extensive published $R_{1\rho}$ data on protein backbone and side-chain dynamics. Palmer's group has developed much of the theoretical framework for analyzing $R_{1\rho}$ dispersion in terms of chemical exchange, and they publish raw dispersion curves with full experimental parameters. Key papers: Massi, Johnson, Wang, Palmer (JACS 2004); Miloushev and Palmer (J. Magn. Reson. 2005).

**Dmitry Korzhnev (UConn Health):** Has published $R_{1\rho}$ dispersion data on well-characterized protein folding systems (Fyn SH3 domain, aB-crystallin). The G48M Fyn SH3 domain is particularly well characterized — the thermodynamic and kinetic parameters of the folding transition are known independently from other methods.

### What to Compute

For a two-spin or three-spin subsystem within a protein (e.g., an amide $^{15}$N-$^1$H pair with a nearby $^{13}$C), with known chemical shifts, J-couplings, and dipolar couplings:

1. Compute $R_{1\rho}$ as a function of spin-lock power using standard Redfield theory (Dirac coefficients).
2. Compute $R_{1\rho}$ as a function of spin-lock power using the nonadiabatic CGME (nonadiabatic coefficients).
3. Compare both predictions against the published experimental dispersion curve, using the same independently known parameters.

The prediction is: at low spin-lock powers (deep perturbative regime), the nonadiabatic correction will be largest, and the nonadiabatic CGME should match the experimental data more closely than standard Redfield theory. At high spin-lock powers (strong driving, outside the perturbative regime), the two predictions should converge or the perturbative framework may break down.

### Feasibility Assessment

**Physics match:** Excellent. Simultaneous driving and dissipation, perturbative regime accessible, constant drive (plateau).

**Parameter independence:** Excellent. All Hamiltonian and bath parameters are independently measured.

**Data availability:** Good. Published in journal supplementary materials, and the research groups are accustomed to sharing raw data with theorists.

**Computational difficulty:** Moderate. The spin Hamiltonian for a few-spin subsystem is small (4×4 or 8×8), and the Redfield tensor can be constructed analytically. The nonadiabatic decomposition for a spin system under a spin-lock field is a new calculation but straightforward in principle.

**Risk:** The main risk is that the difference between Dirac and nonadiabatic predictions for $R_{1\rho}$ may be small for the available data — comparable to or smaller than the experimental error bars. The magnitude of the correction needs to be estimated before committing to a full comparison.

---

## Candidate 2: Coherent Transients in Pulsed Fourier-Transform Microwave Spectroscopy

### The Opportunity

Fourier-transform microwave (FTMW) spectroscopy works by applying a short microwave pulse (typically a few microseconds) to a molecular beam, then observing the free induction decay (FID). The pulse has a well-characterized shape and frequency content, the molecular rotational Hamiltonian is known to extraordinary precision, and the FID amplitude depends directly on the transition amplitudes produced by the pulse.

For molecules with closely spaced rotational transitions — due to hyperfine structure (nuclear quadrupole coupling), internal rotation, or tunneling splittings — the non-secular terms in the relaxation theory matter, and the standard vs. nonadiabatic treatments should differ in their predictions for the FID envelope and the coherence dynamics.

The modern chirped-pulse FTMW (CP-FTMW) technique, pioneered by Brooks Pate at the University of Virginia, sweeps the microwave frequency across a broad band (up to 12 GHz) in less than a microsecond. The chirped pulse has a precisely known time-frequency profile, making it possible to compute the excitation amplitudes for every transition in the molecular spectrum. The data sets from CP-FTMW experiments are rich, publicly documented, and contain information about coherences between closely spaced transitions.

### Available Data and Contacts

**Brooks Pate's group (University of Virginia):** Pioneered CP-FTMW spectroscopy and has published extensive data sets on molecules with complex hyperfine structure. The pulse shapes, frequencies, and molecular parameters are all reported in detail. The group has a culture of open data sharing.

**Melanie Schnell's group (DESY Hamburg / University of Kiel):** High-resolution rotational spectroscopy with published data on molecules with internal dynamics (internal rotation, tunneling), where closely spaced transitions create non-secular physics.

### What to Compute

For a molecule with nuclear quadrupole hyperfine structure (e.g., a nitrogen-containing molecule like pyridine or acetonitrile), with known rotational constants, quadrupole coupling constants, and transition dipole moments:

1. Compute the excitation amplitudes produced by the chirped pulse using Dirac's time-dependent perturbation theory.
2. Compute the excitation amplitudes using the nonadiabatic decomposition.
3. Predict the FID envelope and compare against the measured signal.

The prediction is: for transitions with small frequency splittings (comparable to the inverse pulse duration), the Dirac and nonadiabatic amplitudes will differ, and the FID envelope shape will reflect this difference.

### Feasibility Assessment

**Physics match:** Good. The chirped pulse provides a well-characterized driving field, and collisional relaxation in the molecular beam provides (weak) dissipation.

**Parameter independence:** Excellent. Rotational constants and quadrupole coupling constants are known to kHz precision from independent microwave spectroscopy.

**Data availability:** Good. Published in journal supplementary materials and sometimes available in online databases.

**Computational difficulty:** Moderate to hard. Computing the nonadiabatic coefficients for a chirped pulse (whose frequency sweeps in time) requires numerical evaluation of the Landau-Lifshitz integral, which is more involved than for a monochromatic pulse.

**Risk:** The dissipation in a molecular beam is very weak (collisional broadening is small at the low pressures used), which means the driven-dissipative regime may not be fully realized. The nonadiabatic correction may manifest primarily in the coherent dynamics rather than in the relaxation, which is a different test than what the framework was designed for.

---

## Candidate 3: HITRAN Line Mixing Parameters

### The Opportunity

The HITRAN database contains measured spectral line parameters (positions, intensities, pressure-broadening coefficients, line-mixing parameters) for gas-phase molecules at various pressures and temperatures. For closely spaced rotational lines, collisions with buffer gas molecules cause "line mixing" — the observed spectral profile deviates from the sum of independent Lorentzian lines because the bath (collisions) couples different transitions.

The standard theoretical treatment of line mixing uses Redfield-type relaxation theory (the Anderson-Tsao-Curnutte formalism, extended by many groups). The relaxation matrix that governs line mixing is precisely the Redfield tensor, and the non-secular terms produce the coupling between lines. The secular approximation corresponds to ignoring line mixing entirely (treating each line as independent), which fails for closely spaced lines.

The nonadiabatic framework would modify the Redfield tensor by using nonadiabatic transition probabilities rather than Dirac coefficients. For molecules observed during absorption (where the radiation field is the "driving" and the buffer gas collisions are the "bath"), this is exactly the driven-dissipative scenario.

### Available Data and Contacts

**HITRAN database (Harvard-Smithsonian Center for Astrophysics):** Publicly available at hitran.org. Contains measured line-mixing parameters (first-order and higher-order) for CO, CO₂, N₂O, CH₄, and other molecules at various temperatures and pressures.

**Laurence Rothman (Harvard-Smithsonian):** Long-time curator of HITRAN, knowledgeable about the quality and provenance of all data in the database.

**Robert Gamache (University of Massachusetts Lowell):** Has computed theoretical broadening and line-mixing parameters for many HITRAN molecules, using ab initio intermolecular potential surfaces. His calculations provide the "standard Redfield" predictions that the nonadiabatic framework would modify.

**Jean-Michel Hartmann (Laboratoire de Météorologie Dynamique, Paris):** A leading theorist on line mixing, non-Markovian effects, and beyond-Redfield corrections in molecular spectroscopy. Has published extensively on cases where standard Redfield theory fails to match measured line shapes, particularly for CO₂ Q-branches at high pressures. Author of the monograph *Collisional Effects on Molecular Spectra* (Elsevier, 2008/2021) which comprehensively covers the theory.

### What to Compute

For a molecule like CO in N₂ buffer gas at room temperature:

1. Compute the line-mixing parameters (off-diagonal relaxation matrix elements) using standard Redfield theory with known intermolecular potentials.
2. Compute the same parameters using the nonadiabatic Redfield tensor.
3. Compare both predictions against the measured HITRAN line-mixing parameters.

The prediction is: for closely spaced lines in the Q-branch (where $\Delta J = 0$ transitions have small frequency separations), the nonadiabatic correction to the relaxation matrix will shift the predicted line shape, potentially resolving known discrepancies between standard Redfield theory and experiment.

### Feasibility Assessment

**Physics match:** Good, but indirect. The "driving" here is the radiation field during absorption, and the connection between the nonadiabatic decomposition (formulated for a system driven by an external pulse) and the steady-state absorption line shape requires theoretical work to establish.

**Parameter independence:** Excellent. Molecular parameters are known spectroscopically; intermolecular potentials are computed ab initio.

**Data availability:** Excellent. HITRAN is public.

**Computational difficulty:** Hard. The relaxation matrix for a manifold of rotational states can be large, and computing it from ab initio potentials is already a substantial calculation. Adding the nonadiabatic correction on top of this is feasible but not trivial.

**Risk:** The connection between Hunt's formulation (time-domain, pulsed excitation) and the line-shape formulation (frequency-domain, steady-state absorption) needs to be established carefully. They are related by Fourier transform, but the nonadiabatic decomposition has not been developed in the frequency domain.

---

## Candidate 4: 2D Electronic Spectroscopy of Molecular Dimers

### The Opportunity

Molecular dimers (two coupled chromophores in solution) are the simplest electronic systems where non-secular effects, ultrafast driving, and solvent-mediated dissipation all matter simultaneously. The dimer Hamiltonian has two closely spaced electronic excited states (the symmetric and antisymmetric exciton states), the laser pulses provide the driving, and the solvent provides the thermal bath. The 2D electronic spectrum encodes the coherence dynamics directly through cross-peak amplitudes and oscillation frequencies.

Multiple groups have noted discrepancies between standard Redfield/Lindblad predictions and measured coherence lifetimes in molecular dimers and small aggregates. In some cases, the predicted coherence decays faster than observed; in others, the line shapes in the 2D spectrum do not match. These discrepancies have been attributed variously to vibronic coupling, correlated bath fluctuations, or non-Markovian effects. The nonadiabatic framework offers a different explanation: the standard treatment uses Dirac coefficients, which conflate adiabatic polarization with transitions, leading to incorrect relaxation dynamics.

### Available Data and Contacts

**Graham Fleming's group (UC Berkeley):** Pioneered 2D electronic spectroscopy and has published extensive data on photosynthetic complexes and model dimers. Their data on the FMO complex includes temperature-dependent coherence lifetimes that could be compared against nonadiabatic predictions.

**Gregory Engel's group (University of Chicago):** Published the original long-lived coherence observations in FMO (Science 2007) and subsequent work refining the measurements. Has developed 3D electronic spectroscopy that resolves individual exciton transition energies.

**Tobias Brixner's group (Universität Würzburg):** Extensive 2D spectroscopy data on molecular aggregates with well-characterized electronic Hamiltonians.

### Feasibility Assessment

**Physics match:** Good for the non-secular CGME. The closely spaced exciton states and the solvent bath create exactly the conditions where coherence transfer matters. But the laser pulses are intense (not perturbative), and the system-bath coupling is not always weak, so the Born-Markov and perturbative assumptions may be strained.

**Parameter independence:** Moderate. The electronic Hamiltonian can be extracted from the spectra themselves (as in the Hayes-Engel 3D spectroscopy work), but the spectral density parameters often require fitting.

**Data availability:** Good. Published in journals, sometimes with supplementary data.

**Computational difficulty:** Hard. The standard calculation (Redfield or HEOM for a dimer coupled to a Drude-Lorentz bath) is well-established but not trivial. Adding the nonadiabatic correction requires extending the framework to ultrafast pulse sequences, which is new territory.

**Risk:** The laser pulses in 2D spectroscopy are broadband and intense — far from the perturbative regime. The Landau-Lifshitz decomposition may not apply directly. This candidate is better suited as a second-generation test, after the framework has been validated in the perturbative regime (Candidate 1).

---

## Recommended Priority Order

| Priority | Candidate | Rationale |
| :--- | :--- | :--- |
| 1 | NMR $R_{1\rho}$ dispersion | Best physics match (plateau drive + dissipation), perturbative regime accessible, independently known parameters, established collaborator base |
| 2 | FTMW coherent transients | Excellent parameter knowledge, well-characterized pulses, but weaker dissipation |
| 3 | HITRAN line mixing | Vast public data, but requires frequency-domain reformulation of the nonadiabatic framework |
| 4 | 2D electronic spectroscopy | High impact if successful, but outside the perturbative regime; better as a second-generation test |

---

## Immediate Next Steps

1. **Estimate the magnitude of the nonadiabatic correction for $R_{1\rho}$.** Before contacting experimentalists, compute the expected difference between Dirac and nonadiabatic predictions for a simple two-spin system under a weak spin-lock. If the correction is smaller than typical experimental error bars (roughly 1–5% for $R_{1\rho}$ measurements), the comparison will not be informative, and a different candidate should be prioritized.

2. **Identify specific published $R_{1\rho}$ dispersion data sets.** Look for systems with: (a) well-separated chemical exchange contributions (so that the fundamental relaxation can be isolated from exchange broadening), (b) low spin-lock powers (deep perturbative regime), and (c) independently characterized Hamiltonian parameters.

3. **Contact Lewis Kay or Arthur Palmer.** Frame the request as: "We have a modified relaxation theory (nonadiabatic CGME) that predicts specific corrections to $R_{1\rho}$ at low spin-lock powers. We would like to test it against your published dispersion data using your reported parameters. Can you point us to your most precisely characterized data set?"

4. **In parallel, assess HITRAN feasibility.** Determine whether the nonadiabatic decomposition can be reformulated in the frequency domain (as a modification to the relaxation matrix in the line-shape calculation). If yes, the HITRAN comparison becomes straightforward because the data is public and the standard Redfield predictions are published alongside the experimental values.

---

## References

### NMR $R_{1\rho}$
* Korzhnev, D. M., Orekhov, V. Y., & Kay, L. E. (2005). Off-resonance $R_{1\rho}$ NMR studies of exchange dynamics in proteins with low spin-lock fields. *JACS*, 127, 713–721.
* Korzhnev, D. M., Skrynnikov, N. R., Millet, O., Torchia, D. A., & Kay, L. E. (2002). An NMR experiment for the accurate measurement of heteronuclear spin-lock relaxation rates. *JACS*, 124, 10743–10753.
* Palmer, A. G. (2004). NMR characterization of the dynamics of biomacromolecules. *Chem. Rev.*, 104, 3623–3640.
* Massi, F., Johnson, E., Wang, C., & Palmer, A. G. (2004). Microsecond timescale backbone conformational dynamics in ubiquitin studied with NMR $R_{1\rho}$ relaxation experiments. *Protein Sci.*, 13, 735–746.

### FTMW Spectroscopy
* Brown, G. G., Dian, B. C., Douglass, K. O., Geyer, S. M., Shipman, S. T., & Pate, B. H. (2008). A broadband Fourier transform microwave spectrometer based on chirped pulse excitation. *Rev. Sci. Instrum.*, 79, 053103.

### HITRAN and Line Mixing
* Gordon, I. E., et al. (2022). The HITRAN2020 molecular spectroscopic database. *JQSRT*, 277, 107949.
* Hartmann, J.-M., Boulet, C., & Robert, D. (2021). *Collisional Effects on Molecular Spectra* (2nd ed.). Elsevier.
* Gamache, R. R., et al. (2017). Total internal partition sums for 166 isotopologues of 51 molecules important in planetary atmospheres. *JQSRT*, 203, 70–87.

### 2D Electronic Spectroscopy
* Hayes, D., & Engel, G. S. (2011). Extracting the excitonic Hamiltonian of the Fenna-Matthews-Olson complex using three-dimensional third-order electronic spectroscopy. *Biophys. J.*, 100, 2043–2052.
* Duan, H.-G., et al. (2022). Quantum coherent energy transport in the Fenna-Matthews-Olson complex at low temperature. *PNAS*, 119, e2212630119.

### Nonadiabatic Framework
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse: Toward experimental tests. *J. Chem. Phys.*, 149, 204110.
* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields: Dephasing and population relaxation due to contact with a bath. *J. Chem. Phys.*, 158, 164107.