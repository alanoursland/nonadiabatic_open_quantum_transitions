# Experiment Plan: Testing the Nonadiabatic CGME Against HITRAN Line-Mixing Data

## Overview

This plan describes a computational experiment to test the nonadiabatic coarse-grained master equation (CGME) against published HITRAN pressure-broadened line-mixing data for CO$_2$ in N$_2$ buffer gas. No new laboratory measurements are required. The experiment compares two theoretical predictions — one from the standard Redfield relaxation matrix using Dirac's coefficients, and one from the nonadiabatic CGME using the Landau-Lifshitz decomposition — against the same published spectroscopic data, using the same independently known molecular and collisional parameters. If the nonadiabatic prediction matches the data more closely, that constitutes evidence for the physical correctness of the framework.

The theoretical basis for the comparison is the rotational nonadiabatic CGME derived in `rotational_nonadiabatic_cgme.md`, which reduces the general nonadiabatic Nakajima-Zwanzig framework to a truncated rigid rotor coupled to a collisional bath. The key prediction is that the off-diagonal elements of the relaxation matrix (the line-mixing coefficients) are reduced by a $J$-dependent factor $\mathcal{C}_{JJ'}$ that scales as $(1 - (\gamma_L P / 2B(J+1))^2)^{1/2}$.

---

## Phase 0A: Feasibility Check

**Status: COMPLETE — POSITIVE**

The Phase 0A estimate (reported in `0a_nonadiabatic_correction_to_dissipative_dynamics.md`) established:

- CO$_2$ at 5 atm in N$_2$: $V/\Delta E = 0.449$ at $J = 0$, giving $|a|^2 = 0.201$ (20% correction). **POSITIVE.**
- The correction is $J$-dependent: $|a|^2 \propto 1/(J+1)^2$, falling below 1% for $J \geq 4$ at 5 atm.
- The CGME sinc factors are near 1 (0.955–0.999 for CO$_2$), meaning non-secular Redfield elements are fully retained. The secular approximation (which would discard line mixing entirely) is invalid for this system.
- First-order perturbation theory is marginally valid at 5 atm ($|a|^2 = 0.20$) and stretched at 10 atm ($|a|^2 = 0.80$). The useful operating regime is 1–7 atm.

The feasibility gate is passed. Proceed to data selection and numerical comparison.

---

## Phase 1: Target Data Selection

### Goal

Identify 2–3 specific published data sets containing CO$_2$ line-mixing parameters or absorption cross-sections at elevated pressures, with independently known molecular and collisional parameters.

### Selection Criteria

1. **Elevated pressure (2–10 atm).** The nonadiabatic correction scales as $(\gamma_L P)^2$. At 1 atm the correction is borderline (~1%); at 5 atm it is large (~20% at $J = 0$). Data at multiple pressures allows testing the predicted pressure-squared scaling.

2. **Low-$J$ transitions resolved.** The correction is largest for $J = 0$–$3$ and falls off as $1/(J+1)^2$. Data sets that resolve individual R-branch or P-branch lines at low $J$ are ideal. Band-averaged or low-resolution spectra wash out the $J$-dependent signature.

3. **Independently known $W$-matrix elements.** The comparison must use the same off-diagonal relaxation matrix elements $W_{JJ'}^{\text{Dirac}}$ for both frameworks. These should come from published calculations using ab initio intermolecular potentials, not fitted from the same line-shape data. Published $W$-matrices for CO$_2$-N$_2$ exist from several groups.

4. **CO$_2$ $\nu_3$ band (4.3 $\mu$m) preferred.** This is the most thoroughly studied CO$_2$ band for line mixing. The R-branch transitions are well-resolved at moderate pressures, and the line-mixing effects have been documented extensively. The $\nu_3$ band is also important for atmospheric remote sensing, so discrepancies between theory and measurement have practical consequences.

5. **Room temperature (296 K).** HITRAN's reference temperature. Most published $W$-matrices are computed at 296 K.

### Primary Data Sources

**Source A: HITRAN2020 database (Gordon et al., JQSRT 2022)**

The HITRAN database tabulates first-order line-mixing parameters (the Rosenkranz $Y$-parameters) for CO$_2$ transitions. These $Y$-parameters are the most direct connection to the off-diagonal $W$-matrix:

$$Y_J = \frac{2}{P \gamma_J} \sum_{J' \neq J} \frac{W_{JJ'} \, d_{J'} \, \rho_{J'}}{d_J \, \rho_J} \cdot \frac{1}{\nu_J - \nu_{J'}}$$

where the sum runs over all other transitions. The nonadiabatic correction modifies $W_{JJ'}$, which propagates into $Y_J$. HITRAN tabulates measured $Y_J$ values with uncertainties.

The specific value: HITRAN is public, freely accessible at hitran.org, and the $Y$-parameters are tabulated for every CO$_2$ transition in the $\nu_3$ band.

**Source B: Hartmann et al. — Relaxation matrices from ab initio potentials**

Jean-Michel Hartmann's group has computed full relaxation matrices $W_{JJ'}$ for CO$_2$-N$_2$ and CO$_2$-air using the energy-corrected sudden (ECS) and modified exponential gap (MEG) models, fitted to ab initio intermolecular potential surfaces. These calculations provide the standard Redfield $W_{JJ'}^{\text{Dirac}}$ values needed for both frameworks.

Key publications:
- Hartmann, Boulet, & Robert, *Collisional Effects on Molecular Spectra* (Elsevier, 2021) — comprehensive tables
- Niro, Boulet, & Hartmann, JQSRT (2004) — CO$_2$ R-branch relaxation matrix
- Lamouroux et al., JQSRT (2015) — updated line-mixing parameters for HITRAN

The specific value: Hartmann's $W$-matrices are the industry standard for CO$_2$ line mixing. They are what HITRAN's $Y$-parameters are derived from. Using them as $W_{JJ'}^{\text{Dirac}}$ provides the cleanest comparison — the same input, two frameworks, one set of data.

**Source C: Predoi-Cross et al. — Measured absorption cross-sections at elevated pressures**

Agnes Predoi-Cross (University of Lethbridge) has published high-resolution absorption measurements of CO$_2$ in N$_2$ at pressures from 1 to 50 atm. These are direct measurements of the absorption line shape $F(\nu)$, not derived parameters. Comparing the computed line shape (from the resolvent formula with $\mathbf{W}^{\text{Dirac}}$ vs. $\mathbf{W}^{\text{nad}}$) against these measured profiles is the most direct test.

Key publications:
- Predoi-Cross et al., JQSRT (2007) — CO$_2$ $\nu_3$ band at 5–50 atm
- Predoi-Cross et al., J. Mol. Spectrosc. (2010) — CO$_2$ broadening and line mixing

The specific value: measured line shapes at precisely the pressures where the nonadiabatic correction is predicted to be largest (5–10 atm).

**Source D: Devi et al. — Multi-spectrum fitting of CO$_2$ line parameters**

V. Malathy Devi and collaborators (College of William & Mary) have performed multi-spectrum fits of CO$_2$ absorption data, extracting broadening, shifting, and line-mixing parameters simultaneously from multiple pressures. Their fitted parameters can be compared against both the standard and nonadiabatic predictions.

### Actions

1. Download CO$_2$ $\nu_3$ R-branch $Y$-parameters from HITRAN2020 for transitions $R(0)$ through $R(20)$ at 296 K in N$_2$.
2. Obtain published $W_{JJ'}^{\text{Dirac}}$ matrices from Hartmann's calculations (from the monograph or supplementary data in the JQSRT papers).
3. If direct line-shape comparison is pursued: obtain the Predoi-Cross absorption cross-section data at 5 and 10 atm.
4. Extract the molecular parameters needed for the comparison: $B = 0.3902$ cm$^{-1}$, $\gamma_J$ and $\delta_J$ for each transition (from HITRAN), transition dipole strengths $d_J$, and thermal populations $\rho_J$ at 296 K.

### Estimated Time

1–2 weeks to locate and extract all needed parameters from the literature and databases.

---

## Phase 2: The Nonadiabatic Correction to the $W$-Matrix

### Goal

Apply the nonadiabatic CGME (from `rotational_nonadiabatic_cgme.md`) to compute the corrected relaxation matrix $\mathbf{W}^{\text{nad}}$ and the corrected $Y$-parameters, and verify the framework's internal consistency.

### Steps

**Step 2.1: Compute the standard $W$-matrix.**

For the CO$_2$ $\nu_3$ R-branch, assemble the $W$-matrix from Hartmann's published values. The matrix dimension is set by the number of R-branch transitions included (typically $R(0)$ through $R(40)$ or higher). For the initial comparison, truncate at $R(20)$ — the nonadiabatic correction is below 0.1% for $J > 10$ at 5 atm, so higher-$J$ transitions serve as built-in controls.

Account for CO$_2$ symmetry: only even-$J$ ground-state levels are populated (for $^{16}$O$^{12}$C$^{16}$O), so the $W$-matrix couples even-$J$ R-branch transitions to even-$J$ and odd-$J$ to odd-$J$. The even-$J$ block ($R(0), R(2), R(4), \ldots$) and odd-$J$ block ($R(1), R(3), R(5), \ldots$) can be treated independently.

**Step 2.2: Compute the nonadiabatic correction factors.**

For each off-diagonal element $W_{JJ'}$, compute:

$$\mathcal{C}_{JJ'} = \left(1 - \frac{(\gamma_L P)^2}{[2B(J+1)]^2}\right)^{1/2} \left(1 - \frac{(\gamma_L P)^2}{[2B(J'+1)]^2}\right)^{1/2}$$

This gives the nonadiabatic $W$-matrix:

$$W_{JJ'}^{\text{nad}} = W_{JJ'}^{\text{Dirac}} \times \mathcal{C}_{JJ'}$$

The diagonal elements (linewidths $\gamma_J P$ and shifts $\delta_J P$) are unchanged.

**Important note on the correction factor form.** The expression above uses $(1 - |a|^2)^{1/2}$ per transition, appropriate if the $W$-matrix elements are amplitude-level quantities. If they are rate-level quantities, the correct factor is $(1 - |a|^2)$ without the square root. The rotational CGME document flags this ambiguity. The resolution depends on how the density matrix transformation $\sigma \to \sigma_{\text{nad}}$ propagates through the Redfield tensor into the $W$-matrix. Both forms should be computed and compared against data — the data itself can discriminate between them.

**Step 2.3: Compute the corrected $Y$-parameters.**

From $\mathbf{W}^{\text{nad}}$, compute the Rosenkranz $Y$-parameters:

$$Y_J^{\text{nad}} = \frac{2}{P \gamma_J} \sum_{J' \neq J} \frac{W_{JJ'}^{\text{nad}} \, d_{J'} \, \rho_{J'}}{d_J \, \rho_J} \cdot \frac{1}{\nu_J - \nu_{J'}}$$

The prediction is: $Y_J^{\text{nad}} = Y_J^{\text{Dirac}} \times \mathcal{C}_J^{\text{eff}}$, where $\mathcal{C}_J^{\text{eff}}$ is a weighted average of $\mathcal{C}_{JJ'}$ over the contributing transitions. For low $J$, $\mathcal{C}_J^{\text{eff}}$ is significantly less than 1; for high $J$, it approaches 1.

**Step 2.4: Compute the corrected line shape (if direct comparison is pursued).**

Using the resolvent formula:

$$F(\nu) = \frac{1}{\pi} \text{Im}\left[\mathbf{d}^T \left( (\nu - \boldsymbol{\nu}_0) \mathbf{I} - i\mathbf{W}^{\text{nad}} \right)^{-1} \boldsymbol{\rho} \right]$$

This is a matrix inversion at each frequency point. For a $20 \times 20$ $W$-matrix, this is trivial.

**Step 2.5: Verify internal consistency.**

Check that the nonadiabatic $W$-matrix satisfies the sum rules that the standard $W$-matrix satisfies:

- **Sum rule for linewidths:** $\sum_{J'} W_{JJ'} \rho_{J'} d_{J'}^2 = 0$ (conservation of total intensity). The nonadiabatic correction modifies the off-diagonal elements but not the diagonal, so this sum rule will be slightly violated. The violation is of order $|a|^2 \times W_{JJ'}$, which is the correction itself. This is expected and physical — the nonadiabatic framework redistributes intensity differently from the standard one.

- **Detailed balance:** $W_{JJ'} \rho_{J'} d_{J'}^2 = W_{J'J} \rho_J d_J^2$. The correction factor $\mathcal{C}_{JJ'} = \mathcal{C}_{J'J}$ is symmetric, so detailed balance is preserved.

- **Limit checks:** At $P \to 0$, $\mathcal{C}_{JJ'} \to 1$ and $\mathbf{W}^{\text{nad}} \to \mathbf{W}^{\text{Dirac}}$. At large $J$, $\mathcal{C}_{JJ'} \to 1$ for all pressures.

### Deliverable

A table and code producing, for each R-branch transition $R(J)$ at pressures $P = 1, 2, 3, 5, 7, 10$ atm:

- $\mathcal{C}_J^{\text{eff}}$ (the effective nonadiabatic reduction factor)
- $Y_J^{\text{Dirac}}$ (from published $W$-matrix)
- $Y_J^{\text{nad}}$ (with nonadiabatic correction)
- The difference $\Delta Y_J = Y_J^{\text{Dirac}} - Y_J^{\text{nad}}$

### Estimated Time

2–3 weeks for the computation and verification.

---

## Phase 3: Numerical Comparison Against HITRAN Data

### Goal

Compare the standard and nonadiabatic predictions against published line-mixing measurements, transition by transition, pressure by pressure.

### Procedure

**Comparison A: $Y$-parameters against HITRAN**

For each R-branch transition $R(J)$ with $J = 0$ through $J = 20$:

1. **Compute $Y_J^{\text{Dirac}}$** from the published $W$-matrix (Hartmann et al.) using the standard Rosenkranz formula. No adjustable parameters.

2. **Compute $Y_J^{\text{nad}}$** from the nonadiabatic $W$-matrix, using the same published parameters with only the correction factor $\mathcal{C}_{JJ'}$ applied. No adjustable parameters.

3. **Compare both against the HITRAN-tabulated $Y_J$** values. Compute residuals for both models.

4. **Test the $J$-dependence.** The nonadiabatic prediction has a specific $J$-dependent signature: the correction is largest for low $J$ and falls off as $1/(J+1)^2$. If the HITRAN residuals (measured minus standard-Redfield-predicted) show this same $J$-dependence, that is evidence for the nonadiabatic framework.

5. **Test the pressure dependence.** If data at multiple pressures are available, test whether the correction scales as $P^2$ (as predicted) rather than $P$ or $P^0$.

**Comparison B: Line shapes against measured absorption (if Predoi-Cross data are used)**

For selected spectral windows containing low-$J$ R-branch transitions at 5 atm:

1. **Compute the absorption line shape** $F(\nu)$ using $\mathbf{W}^{\text{Dirac}}$ and the resolvent formula.

2. **Compute the line shape** using $\mathbf{W}^{\text{nad}}$.

3. **Overlay both with the measured absorption cross-section.** Assess which model matches the line-wing profile more closely, particularly in the regions between adjacent lines where line mixing has the largest effect.

### Controls

- **High-$J$ transitions ($J > 10$).** Both predictions should agree to better than 0.1% for these transitions, where $\mathcal{C}_{JJ'} > 0.99$. If they differ significantly at high $J$, there is an error.

- **Low-pressure limit ($P = 1$ atm).** Both predictions should nearly agree ($\mathcal{C}_{JJ'} > 0.99$). If the nonadiabatic prediction deviates significantly at 1 atm, the correction factor is wrong.

- **Isolated lines.** For transitions far from any neighbors (e.g., the highest-$J$ R-branch lines), line mixing is negligible and both models should match the data equally well.

- **Both forms of the correction factor.** Compare the data against both the amplitude-level correction ($(1 - |a|^2)^{1/2}$) and the rate-level correction ($(1 - |a|^2)$). The data should discriminate between them.

### Deliverable

A table and figure set showing, for each $R(J)$ transition and each pressure:

- Measured $Y_J$ (from HITRAN) $\pm$ uncertainty
- Standard Redfield prediction $Y_J^{\text{Dirac}}$
- Nonadiabatic prediction $Y_J^{\text{nad}}$ (both correction factor forms)
- Residuals for all models
- $\chi^2$ comparison

And if line-shape comparison is performed:

- Measured vs. predicted absorption profiles at 5 atm, showing the line-mixing contribution in the inter-line regions

### Estimated Time

2–4 weeks for the numerical calculations and analysis.

---

## Phase 4: Interpretation and Publication

### Possible Outcomes

**Outcome A: The nonadiabatic prediction matches HITRAN data significantly better than standard Redfield.**

This is the positive result. It provides the first experimental evidence that the nonadiabatic decomposition produces physically correct relaxation dynamics in a driven, dissipative system, using independently known parameters and publicly available data. The $J$-dependent and pressure-dependent signatures provide multiple independent tests. The paper presents the theoretical framework (rotational nonadiabatic CGME), the comparison, and the implications for molecular spectroscopy and open quantum systems theory.

This result would also have practical significance for atmospheric remote sensing: the line-mixing parameterization used in radiative transfer codes (which underpins satellite retrievals of CO$_2$ concentration) would need to be updated.

**Outcome B: Both predictions are indistinguishable within HITRAN's uncertainty.**

This means either (a) the correction is smaller than the Phase 0A estimate predicted (possibly because the correction factor form without the square root is correct, halving the effect), or (b) the HITRAN $Y$-parameter uncertainties are too large to resolve the difference. The paper reports the predicted correction magnitude and identifies what measurement precision is needed. Direct line-shape comparison (Comparison B) may still resolve the difference even if the $Y$-parameters cannot, because the line-shape residuals are continuous functions of frequency with correlated structure.

**Outcome C: The nonadiabatic prediction is worse than standard Redfield.**

Possible explanations:

- The correction factor form is wrong (amplitude-level vs. rate-level). Try the other form.
- The perturbative assumption is breaking down at the pressures used. Restrict to lower pressures.
- The $J$-independent coupling approximation is inadequate. Use $J$-dependent $W$-matrix elements.
- The dual role of collisions (perturbation and bath simultaneously) is not correctly handled by the Landau-Lifshitz decomposition in this context. This would be a genuine limitation of the framework for this type of system and would need to be understood.

**Outcome D: The correction has the wrong $J$-dependence or pressure dependence.**

This would strongly constrain the form of the nonadiabatic correction and possibly rule out the first-order perturbative treatment. The specific pattern of disagreement would guide the next theoretical step (higher-order corrections, non-perturbative treatment, or modified coupling model).

### Publication Target

If Outcome A: *Journal of Quantitative Spectroscopy and Radiative Transfer* (the HITRAN home journal, where Hartmann and Gamache publish) or *Journal of Chemical Physics* (where the Mandal-Hunt papers live). Co-authorship with Hunt. Consider contacting Hartmann for collaboration — his expertise on CO$_2$ line mixing and access to published $W$-matrices would strengthen the paper.

If Outcome B: Same journals, framed as "predicted nonadiabatic corrections to pressure-broadened line mixing" with quantitative predictions for experimentally testable signatures.

---

## Open Theoretical Questions

### The Correction Factor Form

The rotational CGME document derives a correction factor $\mathcal{C}_{JJ'}$ with square-root factors $(1 - |a|^2)^{1/2}$ per transition. This assumes the $W$-matrix elements are amplitude-level objects. If they are rate-level (squared), the correction should be $(1 - |a|^2)$ without the square root.

The $W$-matrix in the spectroscopy literature is defined as a rate matrix (units of cm$^{-1}$, acting on the density matrix in Liouville space). The Redfield tensor elements that compose it are bilinear in the system-bath coupling — they involve $|S_{kl}|^2$, not $S_{kl}$. This suggests the rate-level correction $(1 - |a|^2)$ is correct.

However, the nonadiabatic correction enters through the density matrix transformation $\sigma \to \sigma_{\text{nad}}$, which affects the amplitudes $c_k \to b_k$ (linear), not the rates. The correct propagation through the Redfield tensor to the $W$-matrix needs to be worked through explicitly for the rotational system. This is the most important theoretical subtlety and should be resolved before the comparison.

### The Dual Role of Collisions

In the NMR system, the perturbation (spin-lock field) and the bath (molecular tumbling) are physically distinct. In the HITRAN system, collisions serve both roles. The Landau-Lifshitz decomposition separates the coherent (off-diagonal, line-mixing) and stochastic (diagonal, broadening) contributions of collisions. This separation is standard in the Anderson-Tsao-Curnutte formalism, but the mapping onto the nonadiabatic framework needs to be verified explicitly.

The key question: does the adiabatic polarization $|a_k|^2$ arise from the coherent part of the collision (the off-diagonal $W$-matrix element) or from the full collision event (including the stochastic part)? If the former, the correction applies only to the off-diagonal elements, as assumed. If the latter, the diagonal elements (linewidths) are also modified, which changes the prediction.

### Perturbative Validity at 5 atm

At $V/\Delta E = 0.45$ ($J = 0$, 5 atm), $|a|^2 = 0.20$. First-order perturbation theory gives $(1 - |a|^2)^{1/2} = 0.89$ (amplitude) or $(1 - |a|^2) = 0.80$ (rate). The exact (non-perturbative) correction could differ. At 3 atm, $|a|^2 = 0.073$ and first-order perturbation theory is comfortable. The comparison should be performed at multiple pressures (1, 2, 3, 5, 7, 10 atm) to map out where the perturbative treatment works and where it breaks down. A pressure where the perturbative prediction starts to deviate from the data provides direct information about higher-order contributions.

---

## Timeline Summary

| Phase | Activity | Duration | Dependency |
|:---|:---|:---|:---|
| 0A | Feasibility estimate | — | **COMPLETE (POSITIVE)** |
| 1 | Data selection and parameter extraction | 1–2 weeks | Phase 0A positive |
| 2 | Nonadiabatic $W$-matrix and $Y$-parameters | 2–3 weeks | Phase 1 |
| 3 | Numerical comparison against HITRAN | 2–4 weeks | Phases 1 and 2 |
| 4 | Interpretation and writing | 4–6 weeks | Phase 3 |
| **Total remaining** | | **9–15 weeks** | |

Phases 1 and 2 can run in parallel: data extraction from HITRAN/literature does not depend on the theoretical computation. The critical path is Phase 2 (resolving the correction factor form) → Phase 3 (comparison).

---

## Resources Required

- **Computational:** Standard laptop or desktop. The $W$-matrix is at most $30 \times 30$ (for $R(0)$ through $R(40)$ with symmetry). The line-shape calculation is a complex matrix inversion at each frequency point. The existing PyTorch codebase (RotationalSystem, SpectralDensityBath, CollisionalLorentzian) provides the infrastructure; only the $W$-matrix assembly and line-shape computation need to be added.

- **Software:** Python with NumPy/SciPy. The HITRAN Application Programming Interface (HAPI, available at hitran.org) provides programmatic access to line parameters. No specialized spectroscopy codes are required, though LBLRTM (Line-By-Line Radiative Transfer Model) could be used for cross-validation of the line-shape calculation.

- **Data access:** HITRAN is public and free (hitran.org). Hartmann's $W$-matrices are published in journal supplements and in his monograph. Predoi-Cross's absorption data may require direct contact.

- **Collaboration:** Contact with Hunt (for the nonadiabatic framework). Consider contacting Hartmann (for $W$-matrix expertise and known discrepancies) and Gamache (for ab initio broadening parameters). These researchers are active and publish regularly in JQSRT.

---

## References

* Gordon, I. E. et al. (2022). The HITRAN2020 molecular spectroscopic database. *JQSRT*, 277, 107949.
* Hartmann, J.-M., Boulet, C., & Robert, D. (2021). *Collisional Effects on Molecular Spectra* (2nd ed.). Elsevier.
* Niro, F., Boulet, C., & Hartmann, J.-M. (2004). Spectra calculations in central and wing regions of CO$_2$ IR bands between 10 and 20 $\mu$m. I: Model and laboratory measurements. *JQSRT*, 88, 483–498.
* Lamouroux, J. et al. (2015). Updated database plus software for line-mixing in CO$_2$ infrared spectra and their test using laboratory spectra in the 1.5–2.3 $\mu$m region. *JQSRT*, 151, 88–96.
* Predoi-Cross, A. et al. (2007). Measurement of broadening and shift parameters in the $\nu_3$ band of CO$_2$. *JQSRT*, 107, 291–305.
* Rosenkranz, P. W. (1975). Shape of the 5 mm oxygen band in the atmosphere. *IEEE Trans. Antennas Propag.*, 23, 498–506.
* Strow, L. L., & Reuter, D. (1988). Effect of line mixing on atmospheric brightness temperatures near 15 $\mu$m. *Applied Optics*, 27, 872–878.
* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields. *J. Chem. Phys.*, 158, 164107.
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
* Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Physical Review A*, 78, 022106.