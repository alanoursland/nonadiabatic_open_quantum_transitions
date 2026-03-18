# Experiment Plan: Testing the Nonadiabatic CGME Against NMR $R_{1\rho}$ Relaxation Data

## Overview

This plan describes a computational experiment to test the nonadiabatic coarse-grained master equation (CGME) against published NMR rotating-frame relaxation ($R_{1\rho}$) data. No new laboratory measurements are required. The experiment compares two theoretical predictions — one from standard Redfield theory using Dirac's coefficients, and one from the nonadiabatic CGME using the Landau-Lifshitz decomposition — against the same published experimental data, using the same independently known molecular parameters. If the nonadiabatic prediction matches the data more closely, that constitutes evidence for the physical correctness of the framework.

---

## Phase 0: Feasibility Check (Do First — Before Everything Else)

### Goal

Estimate the magnitude of the nonadiabatic correction to $R_{1\rho}$ for a minimal spin system. If the correction is smaller than typical experimental uncertainty, this experiment cannot distinguish the two frameworks, and a different candidate should be pursued.

### Calculation

Consider a single $^{15}$N spin in a static magnetic field $B_0$, subject to a transverse spin-lock field of amplitude $\omega_1 = \gamma B_1$ (in angular frequency units). The spin relaxes through dipolar coupling to a bonded $^1$H and through $^{15}$N chemical shift anisotropy (CSA), with a bath characterized by the rotational correlation time $\tau_c$ of the molecule.

**Standard $R_{1\rho}$ (Redfield, Dirac basis):**

In the standard treatment, the rotating-frame relaxation rate for an on-resonance spin-lock is:

$$R_{1\rho} = R_1 \cos^2\theta + R_2 \sin^2\theta$$

where $\theta$ is the angle between the effective field and $B_0$ (for on-resonance, $\theta = 90°$, so $R_{1\rho} = R_2$; for off-resonance, $\theta$ depends on the offset and spin-lock power). The rates $R_1$ and $R_2$ are computed from the Redfield tensor using spectral densities evaluated at the relevant Bohr frequencies.

**Nonadiabatic $R_{1\rho}$ (CGME, nonadiabatic basis):**

In the nonadiabatic framework, the spin-lock field is the perturbation $V$. The adiabatic response is the adjustment of the spin state to the instantaneous effective field — essentially, the spin aligning with the effective field direction. The nonadiabatic contribution is the component that represents genuine transitions away from this aligned state. During a constant spin-lock ($\partial V/\partial t = 0$), the nonadiabatic populations are constant, and only the nonadiabatic part is subject to relaxation.

The correction enters through how the Redfield tensor couples to the system states. In the Dirac basis, the bath acts on the full density matrix including the adiabatic polarization from the spin-lock field. In the nonadiabatic basis, the bath acts only on the deviation from the adiabatic state. The difference is proportional to the ratio $(\omega_1 / \omega_0)^2$ where $\omega_0$ is the Larmor frequency, modulated by the spectral density at relevant frequencies.

**What to compute:**

1. For a $^{15}$N-$^1$H pair with $\tau_c = 5$ ns (a small protein at room temperature), $B_0 = 11.7$ T (500 MHz $^1$H), and spin-lock powers $\omega_1/2\pi$ = 25, 50, 100, 250, 500, 1000 Hz:
   - Compute $R_{1\rho}$ from standard Redfield theory.
   - Compute $R_{1\rho}$ from the nonadiabatic CGME.
   - Compute the difference $\Delta R_{1\rho} = R_{1\rho}^{\text{Dirac}} - R_{1\rho}^{\text{nad}}$.

2. Compare $\Delta R_{1\rho}$ against typical experimental uncertainty for $R_{1\rho}$ measurements (approximately 1–3% for well-characterized protein systems).

**Decision criterion:**

- If $|\Delta R_{1\rho}| / R_{1\rho} > 0.01$ (1%) at the lowest accessible spin-lock powers: proceed to Phase 1.
- If $|\Delta R_{1\rho}| / R_{1\rho} < 0.01$ at all accessible spin-lock powers: the single-spin system is not sensitive enough. Consider extending to a two-spin system with nearly degenerate transitions (where non-secular effects amplify the correction), or switch to a different experimental candidate (HITRAN line mixing or FTMW).

### Estimated Time

1–2 weeks for a careful analytical or numerical calculation.

---

## Phase 1: Target Data Selection

### Goal

Identify 2–3 specific published $R_{1\rho}$ data sets with the best characteristics for comparison.

### Selection Criteria

1. **Low spin-lock powers tested.** The nonadiabatic correction is largest when the spin-lock field is weakest (most perturbative). Data sets that include measurements at $\omega_1/2\pi \leq 100$ Hz are ideal.

2. **Independently characterized spin Hamiltonian.** The chemical shifts, scalar couplings, dipolar couplings, and CSA tensors of the observed spins must be known from independent measurements (not fitted from the $R_{1\rho}$ data itself).

3. **Characterized bath parameters.** The rotational correlation time $\tau_c$ and internal motion parameters (order parameters $S^2$ from the Lipari-Szabo model-free analysis) must be available from independent $R_1$/$R_2$/NOE measurements.

4. **Minimal chemical exchange.** The target residues should not be undergoing conformational exchange on the $\mu$s–ms timescale, because exchange contributions to $R_{1\rho}$ ($R_{\text{ex}}$) would confound the comparison. Select residues with flat CPMG dispersion profiles (no exchange broadening).

5. **Multiple magnetic field strengths.** Data at two or more $B_0$ fields allows additional cross-validation, since the nonadiabatic correction has a specific field-dependence that differs from standard Redfield predictions.

### Primary Data Sources

**Source A: Korzhnev, Orekhov, & Kay (JACS 2005)**

System: G48M Fyn SH3 domain. $^{15}$N $R_{1\rho}$ measured with spin-lock fields as low as 25 Hz. On- and off-resonance profiles recorded residue by residue. The protein's folding dynamics are well-characterized from CPMG experiments, so residues without exchange can be identified. Hamiltonian parameters (chemical shifts, $^{15}$N CSA, $^{15}$N-$^1$H dipolar coupling) are standard for backbone amide nitrogens. The correlation time $\tau_c$ for Fyn SH3 is independently known.

The specific value of this data set: spin-lock fields as low as 25 Hz have been tested, which is the deepest perturbative regime available in the published literature.

**Source B: Massi, Johnson, Wang, Rance, & Palmer (JACS 2004)**

System: Ubiquitin. $^{15}$N $R_{1\rho}$ measured with weak RF fields between 150 and 1000 Hz. Ubiquitin is one of the most thoroughly characterized proteins in NMR — its dynamics have been studied by every relaxation method, and all Hamiltonian and bath parameters are independently established. Multiple residues without exchange are available.

The specific value: ubiquitin is the NMR relaxation benchmark protein. Any new theoretical prediction can be compared against an extensive existing literature.

**Source C: Palmer group, various publications**

System: Various proteins with $^{15}$N and $^{13}$C $R_{1\rho}$ data at multiple field strengths. Arthur Palmer's group has published decades of relaxation data with meticulous attention to systematic errors, cross-correlated relaxation suppression, and parameter characterization. Their methodological papers often include raw relaxation rates in supplementary tables.

### Actions

1. Obtain the published $R_{1\rho}$ values (from journal supplementary materials or by contacting the authors).
2. Identify 5–10 residues per protein that have: (a) flat CPMG profiles (no exchange), (b) well-determined $R_1$, $R_2$, and NOE values from independent experiments, and (c) $R_{1\rho}$ measured at multiple spin-lock powers.
3. Extract the Hamiltonian parameters for each selected residue: $^{15}$N chemical shift (from the HSQC spectrum), $^{15}$N CSA (use the standard backbone amide value of approximately $-170$ ppm, or site-specific values if available from cross-correlated relaxation measurements), $^{15}$N-$^1$H dipolar coupling (from the known bond length of 1.02 Å), and the Lipari-Szabo parameters ($S^2$, $\tau_e$, $\tau_c$) from the published model-free analysis.

### Estimated Time

2–3 weeks to locate and extract all needed parameters from the literature.

---

## Phase 2: Theoretical Framework for Spin-Lock $R_{1\rho}$

### Goal

Derive the nonadiabatic $R_{1\rho}$ expression for a $^{15}$N-$^1$H spin pair under a spin-lock field, starting from the Landau-Lifshitz decomposition and the nonadiabatic CGME.

### Steps

**Step 2.1: Define the system.**

The system is a $^{15}$N spin (spin-1/2) coupled to a $^1$H spin (spin-1/2) by a dipolar interaction and a scalar coupling. The unperturbed Hamiltonian $H_0$ is the Zeeman interaction of both spins with $B_0$. The perturbation $V$ is the spin-lock field acting on the $^{15}$N spin: $V = -\omega_1 I_x$ in the rotating frame.

In the rotating frame (at the $^{15}$N Larmor frequency), the effective Hamiltonian during the spin-lock is:

$$H_{\text{eff}} = \Delta\omega \, I_z + \omega_1 \, I_x + H_{\text{IS}}$$

where $\Delta\omega$ is the offset from resonance and $H_{\text{IS}}$ is the $^{15}$N-$^1$H coupling.

**Step 2.2: Perform the Landau-Lifshitz decomposition.**

Treat $\omega_1 I_x$ as the perturbation $V$ acting on the eigenstates of $H_0 = \Delta\omega I_z + H_{\text{IS}}$. For on-resonance spin-lock ($\Delta\omega = 0$), the unperturbed eigenstates are the $I_z$ eigenstates. The spin-lock field mixes the $|\alpha\rangle$ and $|\beta\rangle$ states of $^{15}$N.

Apply the Landau-Lifshitz integration by parts to separate the adiabatic response (the spin aligning with the effective field direction) from the nonadiabatic transitions (genuine population changes between the dressed states).

For a constant spin-lock ($\partial V/\partial t = 0$ during the plateau), the nonadiabatic transition probability $|b_k|^2$ is constant, and the adiabatic term $a_k$ describes the equilibrium tilt of the spin toward the effective field axis. This tilt is not a transition — it is the spin's polarization response to the RF field.

**Step 2.3: Construct the nonadiabatic Redfield tensor.**

The system-bath coupling is the dipolar interaction between $^{15}$N and $^1$H, modulated by molecular tumbling. The bath correlation function is determined by the spectral density $J(\omega)$, which for isotropic tumbling is:

$$J(\omega) = \frac{2}{5} \frac{\tau_c}{1 + \omega^2 \tau_c^2}$$

(or the Lipari-Szabo extended model-free form if internal motion is present).

Construct the Redfield tensor elements from the spectral density evaluated at the Bohr frequencies of the spin system. In the nonadiabatic framework, these tensor elements act on the nonadiabatic density matrix $\sigma_{\text{nad}}$ rather than the full density matrix $\sigma$.

**Step 2.4: Apply coarse-graining.**

Multiply each Redfield tensor element by the sinc damping factor:

$$R_{abcd}^{\text{CG}} = R_{abcd} \cdot \text{sinc}\left(\frac{(\omega_{ab} - \omega_{cd})\Delta\tau}{2}\right)$$

Choose $\Delta\tau$ to satisfy $\tau_B \ll \Delta\tau \ll T_1, T_2$. For a protein in solution, $\tau_B \sim \tau_c \sim 5$ ns and $T_1 \sim 1$ s, so the window is wide.

**Step 2.5: Extract $R_{1\rho}$ from the master equation.**

Solve the nonadiabatic CGME for the decay of transverse magnetization under the spin-lock. The decay rate is $R_{1\rho}^{\text{nad}}$.

### Deliverable

An analytical or semi-analytical expression for $R_{1\rho}^{\text{nad}}$ as a function of spin-lock power $\omega_1$, offset $\Delta\omega$, magnetic field $B_0$, correlation time $\tau_c$, and the Lipari-Szabo order parameter $S^2$. This expression must reduce to the standard $R_{1\rho}$ in the limit $\omega_1 \to 0$ (no spin-lock, no perturbation, no adiabatic correction needed).

### Estimated Time

4–8 weeks for a careful derivation and numerical verification.

---

## Phase 3: Numerical Comparison

### Goal

Compare the standard and nonadiabatic predictions against the experimental $R_{1\rho}$ data, residue by residue, spin-lock power by spin-lock power.

### Procedure

For each selected residue in each target protein:

1. **Compute $R_{1\rho}^{\text{Dirac}}(\omega_1)$** using the standard Redfield expression with the independently known parameters ($\delta_N$, $\Delta\sigma$, $r_{\text{NH}}$, $S^2$, $\tau_c$). No adjustable parameters.

2. **Compute $R_{1\rho}^{\text{nad}}(\omega_1)$** using the nonadiabatic CGME expression from Phase 2, with the same parameters. No adjustable parameters.

3. **Plot both predictions alongside the experimental data** as a function of $\omega_1$ (or $\omega_{\text{eff}}$ for off-resonance data).

4. **Compute the residual** (experimental minus predicted) for both models. Assess which model has smaller systematic deviations, particularly at low spin-lock powers where the nonadiabatic correction is predicted to be largest.

5. **Compute a $\chi^2$ statistic** for each model, using the reported experimental uncertainties as weights.

### Controls

- **High spin-lock powers:** At large $\omega_1$, the Dirac and nonadiabatic predictions should converge (the perturbative correction becomes negligible). If they don't converge, there is an error in the derivation.

- **Zero spin-lock limit:** At $\omega_1 = 0$, both predictions should equal $R_2$ (the standard transverse relaxation rate). If they don't, the framework is inconsistent.

- **Residues with exchange:** For residues known to undergo conformational exchange, both models should underestimate $R_{1\rho}$ by the exchange contribution $R_{\text{ex}}$. If the nonadiabatic model matches better for non-exchanging residues but both fail equally for exchanging residues, that's a positive consistency check.

### Deliverable

A table and figure set showing, for each residue and spin-lock power:
- Experimental $R_{1\rho} \pm \sigma$
- Standard Redfield prediction $R_{1\rho}^{\text{Dirac}}$
- Nonadiabatic CGME prediction $R_{1\rho}^{\text{nad}}$
- Residuals for both models
- $\chi^2$ comparison

### Estimated Time

2–4 weeks for the numerical calculations and analysis.

---

## Phase 4: Interpretation and Publication

### Possible Outcomes

**Outcome A: The nonadiabatic prediction matches the data significantly better than standard Redfield.**

This is the positive result. It provides the first experimental evidence (from independent data) that the nonadiabatic decomposition produces physically correct relaxation rates in a driven, dissipative system. The paper would present the theoretical framework, the comparison, and the implications for NMR relaxation theory and for the broader open-quantum-systems program.

**Outcome B: Both predictions are indistinguishable within experimental error.**

This means the nonadiabatic correction is too small to detect with current $R_{1\rho}$ precision. The paper can still report the theoretical framework and the predicted magnitude of the correction, along with a proposal for experiments that would be sensitive enough (e.g., lower spin-lock powers, different nuclei, different field strengths). This also informs whether a different experimental domain (HITRAN, FTMW) should be pursued instead.

**Outcome C: The nonadiabatic prediction is worse than standard Redfield.**

This would indicate either an error in the derivation, a breakdown of the perturbative assumption at the spin-lock powers used, or a genuine failure of the nonadiabatic framework for this type of system. The paper would need to diagnose which explanation applies. If the perturbative assumption is the issue (spin-lock power too high relative to the Larmor frequency splitting in the rotating frame), this constrains the framework's regime of validity without invalidating it.

**Outcome D: The Phase 0 feasibility check shows the correction is too small to be relevant for NMR $R_{1\rho}$.**

Pivot to a different experimental candidate. The HITRAN line-mixing comparison or the FTMW coherent transient comparison become the next priorities. Document the NMR result as a negative feasibility finding and proceed.

### Publication Target

If Outcome A: *Journal of Chemical Physics* (where the Mandal-Hunt papers are published) or *Journal of Magnetic Resonance* (the NMR theory journal). Co-authorship with Hunt and with the NMR experimentalists who provided the data.

If Outcome B: Same journals, framed as "predicted magnitude of nonadiabatic corrections to NMR rotating-frame relaxation" with a clear statement of what precision is needed to test it.

---

## Timeline Summary

| Phase | Activity | Duration | Dependency |
| :--- | :--- | :--- | :--- |
| 0 | Feasibility estimate | 1–2 weeks | None |
| 1 | Data selection and parameter extraction | 2–3 weeks | Phase 0 positive |
| 2 | Theoretical derivation | 4–8 weeks | Phase 0 positive |
| 3 | Numerical comparison | 2–4 weeks | Phases 1 and 2 |
| 4 | Interpretation and writing | 4–6 weeks | Phase 3 |
| **Total** | | **13–23 weeks** | |

Phases 1 and 2 can run in parallel. The critical path is Phase 0 → Phase 2 → Phase 3.

---

## Resources Required

- **Computational:** Standard laptop or desktop. The density matrix is at most 4×4 (two spin-1/2 particles). All calculations are analytically tractable or require trivial numerical integration.

- **Software:** Python with NumPy/SciPy for numerical work. No specialized quantum chemistry or NMR simulation software is required, though packages like SpinDynamica (Mathematica) or SpinEvolution could be used for cross-validation.

- **Data access:** Published journal articles and their supplementary materials. If raw data is not in the supplements, direct contact with Kay, Palmer, or Korzhnev will be needed. NMR relaxation groups routinely share data with theorists — this is an established practice in the field.

- **Collaboration:** Contact with Hunt (for the nonadiabatic framework) and with at least one NMR experimentalist (for data access, parameter verification, and domain expertise on systematic errors in $R_{1\rho}$ measurements).

---

## References

* Korzhnev, D. M., Orekhov, V. Y., & Kay, L. E. (2005). Off-resonance $R_{1\rho}$ NMR studies of exchange dynamics in proteins with low spin-lock fields: An application to a Fyn SH3 domain. *JACS*, 127, 713–721.
* Massi, F., Johnson, E., Wang, C., Rance, M., & Palmer, A. G. (2004). NMR $R_{1\rho}$ rotating-frame relaxation with weak radio frequency fields. *JACS*, 126, 2247–2256.
* Korzhnev, D. M., Skrynnikov, N. R., Millet, O., Torchia, D. A., & Kay, L. E. (2002). An NMR experiment for the accurate measurement of heteronuclear spin-lock relaxation rates. *JACS*, 124, 10743–10753.
* Palmer, A. G. (2004). NMR characterization of the dynamics of biomacromolecules. *Chem. Rev.*, 104, 3623–3640.
* Lipari, G., & Szabo, A. (1982). Model-free approach to the interpretation of nuclear magnetic resonance relaxation in macromolecules. *JACS*, 104, 4546–4559.
* Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Phys. Rev. A*, 78, 022106.
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields: Dephasing and population relaxation due to contact with a bath. *J. Chem. Phys.*, 158, 164107.