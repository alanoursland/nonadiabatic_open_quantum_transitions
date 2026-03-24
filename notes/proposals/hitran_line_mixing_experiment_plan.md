# Experiment Plan: Testing the Nonadiabatic CGME Against HITRAN Line-Mixing Data

## Overview

This plan describes a computational experiment to test the nonadiabatic framework against published HITRAN pressure-broadened line-mixing data for CO$_2$ in N$_2$ buffer gas. No new laboratory measurements are required. The experiment compares two theoretical predictions — one from the standard Redfield relaxation matrix, and one from the nonadiabatic relaxation matrix $\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$ — against the same published spectroscopic data, using the same independently known molecular and collisional parameters.

The theoretical basis is the nonadiabatic line-shape theory derived in `nonadiabatic_line_shape_theory.md`. The central result is that the effective relaxation matrix is the standard Redfield matrix right-multiplied by a nonadiabatic filter superoperator $\mathcal{F} = (\mathbf{I} - \mathbf{M}) \otimes (\mathbf{I}^* - \mathbf{M}^*)$, where $\mathbf{M}$ is the adiabatic mixing matrix built from the mean collisional field and the rotational energy denominators. This produces three testable predictions: first-order inter-branch (R-P) coupling, second-order modification of intra-branch line mixing with a $J$-dependent signature, and sub-Lorentzian far-wing behavior without empirical $\chi$-factors.

---

## Phase 0A: Feasibility Check

**Status: COMPLETE — POSITIVE**

The Phase 0A estimate established:

- CO$_2$ at 5 atm in N$_2$: $V/\Delta E = 0.449$ at $J = 0$, giving $|a|^2 = 0.201$ (20% correction). **POSITIVE.**
- The correction is $J$-dependent: $|a|^2 \propto 1/(J+1)^2$, falling below 1% for $J \geq 4$ at 5 atm.
- The CGME sinc factors are near 1 (0.955–0.999 for CO$_2$), meaning non-secular Redfield elements are fully retained. The secular approximation (which would discard line mixing entirely) is invalid for this system.
- First-order perturbation theory is marginally valid at 5 atm ($|a|^2 = 0.20$) and stretched at 10 atm ($|a|^2 = 0.80$). The useful operating regime is 1–7 atm.

The feasibility gate is passed.

---

## Theoretical Foundation (Resolved)

The connection between the time-domain nonadiabatic decomposition and the frequency-domain line-shape formalism has been established. The key results are:

### The linear response constraint

Applying the nonadiabatic decomposition to the probe radiation field produces zero correction to the absorption line shape. The standard formula already isolates the resonant (nonadiabatic) absorption by taking Im[...]. The adiabatic probe response contributes to the refractive index, not the absorption. The relevant perturbation is not the probe but the collisional interaction.

### The mean-field / fluctuation decomposition

The collisional interaction is decomposed as $V_{\text{int}}(t) = \langle V_{\text{int}} \rangle + \delta V_{\text{int}}(t)$. The mean field $\langle V_{\text{int}} \rangle$ (proportional to pressure, responsible for pressure shifts) plays the role of the perturbation in the Landau-Lifshitz decomposition. The fluctuations $\delta V_{\text{int}}(t)$ (responsible for pressure broadening) play the role of the bath. The nonadiabatic rule states: the fluctuation bath must not relax the adiabatic polarization created by the mean field.

### The nonadiabatic relaxation matrix

The mixing matrix $\mathbf{M}$ has elements $M_{kj} = \langle k | \langle V_{\text{int}} \rangle | j \rangle / (E_j - E_k)$ for $k \neq j$. The nonadiabatic filter superoperator is $\mathcal{F} = (\mathbf{I} - \mathbf{M}) \otimes (\mathbf{I}^* - \mathbf{M}^*)$. The effective relaxation matrix is:

$$\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \, \mathcal{F}$$

This is a right-multiplication, not a scalar correction to individual matrix elements. It changes the topology of the relaxation matrix.

### Consistency checks (completed)

- **Sum rule (intensity conservation): VERIFIED.** The population-weighted dipole vector is a left null vector of $\mathbf{W}^{\text{standard}}$, so $\mathbf{v}^T \mathbf{W}^{\text{nad}} = (\mathbf{v}^T \mathbf{W}^{\text{standard}}) \mathcal{F} = 0$. Intensity conservation is guaranteed regardless of $\mathcal{F}$.

- **Detailed balance: VIOLATED at $O(|\mathbf{M}|^2)$.** This is expected and physically correct. The mean collisional field shifts the effective energy levels (pressure shifts), so the correct thermal equilibrium is with respect to dressed energies $E_J^{\text{dressed}} = E_J + \langle J | \langle V_{\text{int}} \rangle | J \rangle$. The violation is second-order in $V/\Delta E$ and does not affect the first-order inter-branch prediction.

### Structure of the correction

Expanding $\mathcal{F}$ in orders of $\mathbf{M}$:

- **First order** ($O(|M|) \sim V/\Delta E$): generates inter-branch coupling. For $P_2(\cos\theta)$ collisional anisotropy, R-branch coherences couple to P-branch coherences. This is a qualitatively new prediction absent from standard ATC theory.

- **Second order** ($O(|M|^2) \sim (V/\Delta E)^2$): modifies intra-branch line mixing. The standard $W_{JJ'}$ within a branch are corrected with a $J$-dependent factor, largest at low $J$.

---

## Phase 1: Target Data Selection

### Goal

Identify 2–3 specific published data sets containing CO$_2$ line-mixing parameters or absorption cross-sections at elevated pressures, with independently known molecular and collisional parameters.

### Selection Criteria

1. **Elevated pressure (2–10 atm).** The nonadiabatic correction scales as $(\gamma_L P)^2$ for intra-branch and as $\gamma_L P$ for inter-branch. Data at multiple pressures allows testing both scalings.

2. **Low-$J$ transitions resolved.** The correction is largest for $J = 0$–$3$ and falls off as $1/(J+1)^2$ (intra-branch) or $1/(J+1)$ (inter-branch). Data sets that resolve individual R-branch and P-branch lines at low $J$ are ideal.

3. **Both R-branch and P-branch data.** The inter-branch coupling prediction requires data from both branches simultaneously. Data sets that cover both sides of the band center are essential for testing Prediction 1.

4. **Independently known $W$-matrix elements.** The comparison must use published $W$-matrices computed from ab initio intermolecular potentials, not fitted from the same line-shape data.

5. **CO$_2$ $\nu_3$ band (4.3 $\mu$m) preferred.** Most thoroughly studied CO$_2$ band for line mixing. The R-branch and P-branch transitions are well-resolved at moderate pressures.

6. **Room temperature (296 K).** HITRAN's reference temperature. Most published $W$-matrices are computed at 296 K.

### Primary Data Sources

**Source A: HITRAN2020 database (Gordon et al., JQSRT 2022)**

The HITRAN database tabulates first-order line-mixing parameters (the Rosenkranz $Y$-parameters) for CO$_2$ transitions:

$$Y_J = \frac{2}{P \gamma_J} \sum_{J' \neq J} \frac{W_{JJ'} \, d_{J'} \, \rho_{J'}}{d_J \, \rho_J} \cdot \frac{1}{\nu_J - \nu_{J'}}$$

HITRAN is public, freely accessible at hitran.org, and the $Y$-parameters are tabulated for every CO$_2$ transition in the $\nu_3$ band.

**Source B: Hartmann et al. — Relaxation matrices from ab initio potentials**

Jean-Michel Hartmann's group has computed full relaxation matrices $W_{JJ'}$ for CO$_2$-N$_2$ and CO$_2$-air using the energy-corrected sudden (ECS) and modified exponential gap (MEG) models, fitted to ab initio intermolecular potential surfaces. These provide $\mathbf{W}^{\text{standard}}$. Hartmann's matrices also contain the intermolecular potential information needed to construct $\mathbf{M}$.

Key publications:
- Hartmann, Boulet, & Robert, *Collisional Effects on Molecular Spectra* (Elsevier, 2021)
- Niro, Boulet, & Hartmann, JQSRT (2004) — CO$_2$ R-branch relaxation matrix
- Lamouroux et al., JQSRT (2015) — updated line-mixing parameters for HITRAN

**Source C: Predoi-Cross et al. — Measured absorption cross-sections at elevated pressures**

Agnes Predoi-Cross (University of Lethbridge) has published high-resolution absorption measurements of CO$_2$ in N$_2$ at pressures from 1 to 50 atm. These are direct measurements of the absorption line shape, not derived parameters.

Key publications:
- Predoi-Cross et al., JQSRT (2007) — CO$_2$ $\nu_3$ band at 5–50 atm
- Predoi-Cross et al., J. Mol. Spectrosc. (2010) — CO$_2$ broadening and line mixing

**Source D: Devi et al. — Multi-spectrum fitting of CO$_2$ line parameters**

V. Malathy Devi and collaborators (College of William & Mary) have performed multi-spectrum fits extracting broadening, shifting, and line-mixing parameters simultaneously from multiple pressures.

### Actions

1. Download CO$_2$ $\nu_3$ R-branch and P-branch $Y$-parameters from HITRAN2020 for transitions $R(0)$–$R(20)$ and $P(1)$–$P(20)$ at 296 K in N$_2$.
2. Obtain published $\mathbf{W}^{\text{standard}}$ matrices from Hartmann's calculations (monograph or JQSRT supplementary data).
3. Extract the intermolecular potential matrix elements $\langle k | V_{\text{int}} | j \rangle$ needed to construct $\mathbf{M}$. These appear in Hartmann's ECS/MEG parameterization.
4. If direct line-shape comparison is pursued: obtain the Predoi-Cross absorption cross-section data at 5 and 10 atm.
5. Extract molecular parameters: $B = 0.3902$ cm$^{-1}$, $\gamma_J$ and $\delta_J$ for each transition (from HITRAN), transition dipole strengths $d_J$, and thermal populations $\rho_J$ at 296 K.

### Estimated Time

1–2 weeks.

---

## Phase 2: Construction of $\mathbf{W}^{\text{nad}}$

### Goal

Construct the nonadiabatic relaxation matrix $\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$ and compute all observables.

### Steps

**Step 2.1: Assemble the standard $W$-matrix.**

For the CO$_2$ $\nu_3$ band, assemble the full $\mathbf{W}^{\text{standard}}$ from Hartmann's published values. Include both R-branch and P-branch transitions in a single combined matrix to capture inter-branch coupling. The matrix dimension is set by the number of transitions included. For the initial comparison, include $R(0)$–$R(20)$ and $P(1)$–$P(20)$, giving a $\sim 40 \times 40$ matrix.

Account for CO$_2$ symmetry: for $^{16}$O$^{12}$C$^{16}$O, only even-$J$ or odd-$J$ ground-state levels are populated depending on the nuclear spin statistics. The even-$J$ and odd-$J$ blocks can be treated independently within each branch.

In the standard ATC treatment, the R-branch and P-branch blocks of $\mathbf{W}^{\text{standard}}$ are independent (no off-diagonal elements connecting them). The nonadiabatic correction creates inter-branch elements.

**Step 2.2: Construct the mixing matrix $\mathbf{M}$.**

The mixing matrix operates in the space of rotational energy levels (not lines). For a rigid rotor with levels $|J\rangle$, the mean collisional field $\langle V_{\text{int}} \rangle$ couples states with $\Delta J = \pm 2$ (for $P_2(\cos\theta)$ anisotropy). The matrix elements are:

$$M_{kj} = \frac{\langle k | \langle V_{\text{int}} \rangle | j \rangle}{E_j - E_k} \quad (k \neq j), \qquad M_{kk} = 0$$

The numerator $\langle k | \langle V_{\text{int}} \rangle | j \rangle$ is extracted from the same intermolecular potential that Hartmann uses for the $W$-matrix (ECS/MEG parameters). The denominator is the rotational energy gap $E_j - E_k = B[j(j+1) - k(k+1)]$ in cm$^{-1}$.

The matrix $\mathbf{M}$ is sparse: only $\Delta J = \pm 2$ elements are nonzero for $P_2$ anisotropy. Higher-order anisotropy ($P_4$, etc.) would add $\Delta J = \pm 4$ elements but these are typically much smaller.

**Step 2.3: Construct the filter superoperator $\mathcal{F}$.**

$\mathcal{F}$ operates in line space (Liouville space), where each "basis vector" is an optical coherence $|J+1\rangle\langle J|$ (R-branch) or $|J-1\rangle\langle J|$ (P-branch). The superoperator is:

$$\mathcal{F} = (\mathbf{I} - \mathbf{M}) \otimes (\mathbf{I}^* - \mathbf{M}^*)$$

Concretely, the action of $\mathcal{F}$ on a line-space vector $\boldsymbol{\rho}$ is: for each optical coherence $\rho_{J+1,J}$, subtract the adiabatic polarization contributions from states mixed by $\mathbf{M}$. The result is a matrix in line space of the same dimension as $\mathbf{W}^{\text{standard}}$ ($\sim 40 \times 40$), but with nonzero elements connecting R-branch and P-branch lines.

**Step 2.4: Compute $\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$.**

Matrix multiplication. The result has three types of elements:

- **Intra-branch diagonal (linewidths):** Modified at $O(|M|^2)$ from the standard values.
- **Intra-branch off-diagonal (standard line mixing):** Modified at $O(|M|^2)$ with a $J$-dependent correction largest at low $J$.
- **Inter-branch off-diagonal (new):** Nonzero at $O(|M|)$. These elements are identically zero in $\mathbf{W}^{\text{standard}}$ but appear in $\mathbf{W}^{\text{nad}}$.

**Step 2.5: Compute $Y$-parameters from $\mathbf{W}^{\text{nad}}$.**

The Rosenkranz $Y$-parameters are:

$$Y_J^{\text{nad}} = \frac{2}{P \gamma_J} \sum_{J' \neq J} \frac{W_{JJ'}^{\text{nad}} \, d_{J'} \, \rho_{J'}}{d_J \, \rho_J} \cdot \frac{1}{\nu_J - \nu_{J'}}$$

The sum now includes both intra-branch terms (as in the standard calculation) and inter-branch terms (new from the nonadiabatic framework). The inter-branch contributions have large frequency denominators ($\nu_R - \nu_P$ is roughly twice the band center frequency), so they may be suppressed despite the first-order coupling. This needs to be computed numerically.

**Step 2.6: Compute the line shape (if direct comparison is pursued).**

Using the resolvent formula with the combined R+P $\mathbf{W}^{\text{nad}}$:

$$F(\nu) = \frac{1}{\pi} \operatorname{Im}\left[\mathbf{d}^T \left( (\nu - \boldsymbol{\nu}_0) \mathbf{I} - i\mathbf{W}^{\text{nad}} \right)^{-1} \boldsymbol{\rho} \right]$$

For a $40 \times 40$ matrix, this is a trivial computation at each frequency point.

**Step 2.7: Verify internal consistency.**

- **Sum rule:** Already proven analytically. Verify numerically that $\sum_J d_J \rho_J [W^{\text{nad}}]_{JJ'} = 0$ for all $J'$.
- **Limit checks:** At $P \to 0$, $\mathbf{M} \to 0$, $\mathcal{F} \to \mathbf{I}$, $\mathbf{W}^{\text{nad}} \to \mathbf{W}^{\text{standard}}$. At large $J$, $M_{kj} \propto 1/(J+1)$ so the correction vanishes. Verify both numerically.
- **Hermiticity structure:** $\mathbf{W}^{\text{nad}}$ need not be symmetric (since $\mathcal{F}$ is not symmetric), but the absorption formula does not require symmetry of $\mathbf{W}$.

### Deliverable

Code and tables producing, for each transition $R(J)$ and $P(J)$ at pressures $P = 1, 2, 3, 5, 7, 10$ atm:

- The mixing matrix $\mathbf{M}$ and filter $\mathcal{F}$
- $\mathbf{W}^{\text{nad}}$, including the inter-branch block
- $Y_J^{\text{standard}}$ and $Y_J^{\text{nad}}$ for each transition
- The inter-branch coupling matrix elements (the new, structurally distinct prediction)

### Estimated Time

2–3 weeks.

---

## Phase 3: Numerical Comparison Against HITRAN Data

### Goal

Compare the standard and nonadiabatic predictions against published line-mixing measurements, transition by transition, pressure by pressure.

### Comparison A: Intra-branch $Y$-parameters against HITRAN

For each R-branch transition $R(J)$ with $J = 0$ through $J = 20$:

1. Compute $Y_J^{\text{standard}}$ from Hartmann's published $\mathbf{W}^{\text{standard}}$. No adjustable parameters.
2. Compute $Y_J^{\text{nad}}$ from $\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$. No adjustable parameters.
3. Compare both against HITRAN-tabulated $Y_J$ values. Compute residuals for both models.
4. **Test the $J$-dependence.** The nonadiabatic correction is largest for low $J$ and falls off as $1/(J+1)^2$. If the HITRAN residuals (measured minus standard-predicted) show this $J$-dependence, that is evidence for the framework.
5. **Test the pressure dependence.** The intra-branch correction scales as $P^2$. Verify against data at multiple pressures.

### Comparison B: Inter-branch coupling signature

This is the most distinctive prediction. Standard ATC theory predicts zero coupling between R-branch and P-branch relaxation. The nonadiabatic framework predicts coupling at $O(\gamma_L P / 2B(J+1))$.

1. Compute the absorption line shape $F(\nu)$ across the entire $\nu_3$ band (both R and P branches) using $\mathbf{W}^{\text{standard}}$ (R and P independent) and $\mathbf{W}^{\text{nad}}$ (R and P coupled).
2. Compare both against measured absorption spectra (Predoi-Cross data at 5–10 atm).
3. Look for systematic residuals near the band center (the overlap region between the lowest-$J$ R-branch and P-branch lines), where inter-branch coupling would have the largest effect.
4. The inter-branch correction scales as $P$ (first order), distinguishable from the intra-branch correction scaling as $P^2$.

### Comparison C: Line shapes against measured absorption

For selected spectral windows containing low-$J$ transitions at 5 atm:

1. Compute $F(\nu)$ using $\mathbf{W}^{\text{standard}}$ and $\mathbf{W}^{\text{nad}}$.
2. Overlay both with measured absorption cross-sections.
3. Assess which model matches the inter-line regions (where line mixing is most visible) and the far wings (where sub-Lorentzian behavior is expected).

### Comparison D: Far-wing behavior

Compute the absorption coefficient 5–50 cm$^{-1}$ from line center for selected transitions at elevated pressures. Compare $\mathbf{W}^{\text{nad}}$ predictions against measured far-wing absorption. If the nonadiabatic framework produces sub-Lorentzian wings without empirical $\chi$-factors, matching the data that currently requires $\chi$-factors, that is strong independent evidence.

### Controls

- **High-$J$ transitions ($J > 10$).** Both predictions should agree to better than 0.1%, since $|M| \propto 1/(J+1)$ is small. Any significant difference at high $J$ indicates an error.
- **Low-pressure limit ($P = 1$ atm).** Both predictions should nearly agree ($|M| < 0.1$). Significant deviation at 1 atm indicates the correction is too aggressive.
- **Pressure scaling.** Inter-branch coupling should scale as $P$ (first order). Intra-branch modification should scale as $P^2$ (second order). Departures from these scalings indicate breakdown of the perturbative treatment.

### Deliverable

A table and figure set showing:

- Measured $Y_J$ (from HITRAN) $\pm$ uncertainty vs. $Y_J^{\text{standard}}$ vs. $Y_J^{\text{nad}}$ for each transition and pressure
- The inter-branch coupling matrix elements at each pressure
- Measured vs. predicted absorption profiles at 5 atm, showing inter-branch and intra-branch contributions
- Far-wing absorption comparison (with and without empirical $\chi$-factors for the standard prediction)
- $\chi^2$ comparison for all models
- Residual patterns testing $J$-dependence and pressure scaling

### Estimated Time

2–4 weeks.

---

## Phase 4: Interpretation and Publication

### Possible Outcomes

**Outcome A: The nonadiabatic prediction matches data significantly better than standard Redfield.**

The positive result. First experimental evidence that the nonadiabatic decomposition produces physically correct relaxation dynamics, using independently known parameters and publicly available data. The $J$-dependent, pressure-dependent, and inter-branch signatures provide multiple independent tests. Practical significance for atmospheric remote sensing: CO$_2$ line-mixing parameterizations in radiative transfer codes would need updating.

**Outcome B: Intra-branch correction matches but inter-branch coupling is undetectable.**

The inter-branch coupling ($O(V/\Delta E)$) may be suppressed by large frequency denominators ($\nu_R - \nu_P$) in the $Y$-parameter formula. The second-order intra-branch correction would still be testable. This outcome validates the framework at second order but cannot confirm the first-order topological prediction.

**Outcome C: Both predictions are indistinguishable within uncertainty.**

Either the correction is smaller than estimated (possible if the mean-field identification overstates $\langle V_{\text{int}} \rangle$) or HITRAN uncertainties are too large. Direct line-shape comparison (Comparison C) may still resolve the difference. The paper reports predicted correction magnitudes and identifies needed measurement precision.

**Outcome D: The nonadiabatic prediction is worse than standard Redfield.**

Possible explanations: the mean-field / fluctuation separation doesn't correctly map onto the adiabatic / nonadiabatic decomposition for collisional interactions; the perturbative treatment breaks down at the pressures used; or the $P_2$ anisotropy dominance assumption is inadequate. The specific pattern of disagreement constrains the framework.

**Outcome E: The correction has the wrong $J$-dependence or pressure dependence.**

Would strongly constrain the form of the nonadiabatic correction and possibly rule out the first-order perturbative treatment for this system. Guides next theoretical steps.

### Publication Target

Primary: *Journal of Quantitative Spectroscopy and Radiative Transfer* (HITRAN home journal) or *Journal of Chemical Physics* (Mandal-Hunt series). Consider contacting Hartmann for collaboration — his $W$-matrix expertise and knowledge of known discrepancies would strengthen the paper.

---

## Resolved Theoretical Questions

These questions appeared as open issues in earlier versions of this plan and have since been answered.

### The correction structure

**Resolved.** The correction is not a scalar factor on individual $W_{JJ'}$ elements. It is the matrix product $\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$, where $\mathcal{F}$ is a superoperator that changes the topology of the relaxation matrix. This was derived in `nonadiabatic_line_shape_theory.md` and is Outcome C (structural transformation) rather than Outcome A (multiplicative correction).

### The dual role of collisions

**Resolved.** The collisional interaction is decomposed into a mean field $\langle V_{\text{int}} \rangle$ (the perturbation in the Landau-Lifshitz decomposition, responsible for pressure shifts) and fluctuations $\delta V_{\text{int}}(t)$ (the bath, responsible for pressure broadening). This separation already exists in the standard ATC formalism. The nonadiabatic rule adds the constraint that the bath (fluctuations) must not relax the adiabatic polarization created by the perturbation (mean field).

### The linear response constraint

**Resolved.** Applying the nonadiabatic decomposition to the weak probe field produces zero correction to the absorption line shape, because the standard formula already isolates the resonant (nonadiabatic) response via Im[...]. The relevant perturbation for the Landau-Lifshitz decomposition is the collisional mean field, not the probe.

### Sum rule (intensity conservation)

**Resolved and verified.** $\mathbf{W}^{\text{nad}} = \mathbf{W}^{\text{standard}} \mathcal{F}$ satisfies the sum rule because the population-weighted dipole vector is a left null vector of $\mathbf{W}^{\text{standard}}$, and right-multiplication by $\mathcal{F}$ cannot affect the left null space.

### Detailed balance

**Resolved.** $\mathbf{W}^{\text{nad}}$ violates detailed balance with respect to bare rotational energies at $O(|M|^2)$. This is expected: the mean collisional field shifts the effective energies (pressure shifts), so the correct equilibrium is with respect to dressed energies $E_J^{\text{dressed}} = E_J + \langle J | \langle V_{\text{int}} \rangle | J \rangle$. The violation is the same order as the pressure-shift correction to thermal populations.

### Remaining open question: perturbative validity

At $V/\Delta E = 0.45$ ($J = 0$, 5 atm), $|a|^2 = 0.20$. First-order perturbation theory is marginally valid. At 3 atm, $|a|^2 = 0.07$, comfortably perturbative. The comparison should be performed at multiple pressures (1, 2, 3, 5, 7, 10 atm) to map where the perturbative treatment works and where it breaks down.

---

## Timeline Summary

| Phase | Activity | Duration | Dependency |
|:---|:---|:---|:---|
| 0A | Feasibility estimate | — | **COMPLETE (POSITIVE)** |
| Theory | Frequency-domain connection | — | **COMPLETE** |
| 1 | Data selection and parameter extraction | 1–2 weeks | Theory complete |
| 2 | Construction of $\mathbf{W}^{\text{nad}}$ | 2–3 weeks | Phase 1 |
| 3 | Numerical comparison against HITRAN | 2–4 weeks | Phase 2 |
| 4 | Interpretation and writing | 4–6 weeks | Phase 3 |
| **Total remaining** | | **9–15 weeks** | |

Phases 1 and 2 can run in parallel. The critical path is Phase 2 (constructing $\mathbf{M}$ and $\mathcal{F}$ from published intermolecular potential data) → Phase 3 (comparison).

---

## Resources Required

- **Computational:** Standard laptop or desktop. The combined R+P $\mathbf{W}$-matrix is $\sim 40 \times 40$. The line-shape calculation is a complex matrix inversion at each frequency point. The existing PyTorch codebase (RotationalSystem, SpectralDensityBath, CollisionalLorentzian) provides infrastructure; the $\mathbf{M}$-matrix construction, $\mathcal{F}$ assembly, and line-shape computation need to be added.

- **Software:** Python with NumPy/SciPy. The HITRAN Application Programming Interface (HAPI, available at hitran.org) provides programmatic access to line parameters.

- **Data access:** HITRAN is public and free (hitran.org). Hartmann's $W$-matrices are published in journal supplements and in his monograph. Predoi-Cross's absorption data may require direct contact.

- **Collaboration:** Consider contacting Hartmann for $W$-matrix expertise and access to the full ECS/MEG parameterization needed to construct $\mathbf{M}$. His knowledge of known discrepancies between standard theory and HITRAN data would focus the comparison on the most informative spectral regions.

---

## References

* Gordon, I. E. et al. (2022). The HITRAN2020 molecular spectroscopic database. *JQSRT*, 277, 107949.
* Hartmann, J.-M., Boulet, C., & Robert, D. (2021). *Collisional Effects on Molecular Spectra* (2nd ed.). Elsevier.
* Niro, F., Boulet, C., & Hartmann, J.-M. (2004). Spectra calculations in central and wing regions of CO$_2$ IR bands. *JQSRT*, 88, 483–498.
* Lamouroux, J. et al. (2015). Updated database plus software for line-mixing in CO$_2$ infrared spectra. *JQSRT*, 151, 88–96.
* Predoi-Cross, A. et al. (2007). Measurement of broadening and shift parameters in the $\nu_3$ band of CO$_2$. *JQSRT*, 107, 291–305.
* Rosenkranz, P. W. (1975). Shape of the 5 mm oxygen band in the atmosphere. *IEEE Trans. Antennas Propag.*, 23, 498–506.
* Ben-Reuven, A. (1966). Impact broadening of microwave spectra. *Phys. Rev.*, 145, 7.
* Fano, U. (1963). Pressure broadening as a prototype of relaxation. *Phys. Rev.*, 131, 259.
* Mandal, A., & Hunt, K. L. C. (2018). Nonadiabatic transition probabilities in a time-dependent Gaussian pulse or plateau pulse. *J. Chem. Phys.*, 149, 204110.
* Jovanovski, S. D., Mandal, A., & Hunt, K. L. C. (2023). Nonadiabatic transition probabilities for quantum systems in electromagnetic fields. *J. Chem. Phys.*, 158, 164107.
* Schaller, G., & Brandes, T. (2008). Preservation of positivity by dynamical coarse graining. *Physical Review A*, 78, 022106.