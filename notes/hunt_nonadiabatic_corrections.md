Hunt and Mandal’s work (specifically their 1991 analysis) provides a rigorous critique of the standard **Dirac transition probability** ($|c_k(t)|^2$) by comparing it to the **power absorption** approach in time-dependent quantum systems. 

Their results demonstrate that the standard Dirac form often leads to unphysical oscillations and incorrect energy balances, whereas a **nonadiabatic decomposition** aligns perfectly with the laws of thermodynamics and electrodynamics.

---

### **1. Energy Separation and the Hamiltonian**
Hunt and Mandal distinguish between the **unperturbed energy** (associated with $H_0$) and the **total energy** (associated with the full Hamiltonian $H(t)$). 

* **Dirac’s approach** typically interprets transition probabilities based on the eigenstates of $H_0$. 
* **Hunt and Mandal’s result** shows that in the presence of a time-varying field, the energy of the system is not simply the expectation value of $H_0$. Instead, there is a clear **energy separation**: the energy added to the system by the external field ($V(t)$) must be accounted for in the transition coefficients.

### **2. Power Absorption vs. Transition Probability**
One of their most striking findings involves the calculation of **power absorption** ($P(t) = \frac{d\langle H \rangle}{dt}$).

* **The Dirac Flaw:** If you calculate the rate of energy change using the standard Dirac transition probabilities, you often find that the system "absorbs" or "emits" power even when the field is constant or zero (due to persistent oscillations in $|c_k(t)|^2$).
* **The Hunt-Mandal Correction:** They proved that the **true power absorption** is directly proportional to the rate of change of the *nonadiabatic* part of the wavefunction. When the field is stationary, the nonadiabatic corrections vanish, and the power absorption correctly goes to zero, unlike the standard Dirac form.

### **3. Variance and Fluctuations**
Hunt and Mandal analyzed the **variance ($\sigma^2$)** of the energy measurements. 
* They found that the standard Dirac approach overestimates the variance in the system's energy. 
* By using the nonadiabatic form, the variance remains consistent with the physical constraints of the perturbation. The "fictitious" fluctuations seen in the Dirac coefficients are revealed to be artifacts of using a static basis ($H_0$) to describe a dynamic process.

---

### **4. Why the Nonadiabatic Form is Physically Correct**
Hunt and Mandal argue that the nonadiabatic form is the only one that maintains **gauge invariance** and physical consistency.

| Feature | Standard Dirac Form | Nonadiabatic Form (Hunt & Mandal) |
| :--- | :--- | :--- |
| **Field Sensitivity** | Often depends on the choice of gauge (e.g., length vs. velocity). | Is gauge-invariant; results don't change with math representation. |
| **Steady State** | Predicts oscillations even after a field is stabilized. | Predicts a steady state with zero net power absorption. |
| **Energy Balance** | Violates $d\langle H \rangle / dt = \langle \partial V / \partial t \rangle$. | Strictly obeys the energy-work theorem. |



### **5. Nonadiabatic Corrections**
The "corrections" Hunt and Mandal suggest involve redefining the transition coefficient $c_k(t)$. Instead of using the static $H_0$ basis, they employ a transformation that accounts for the **instantaneous state** of the system.

The correction essentially subtracts the "adiabatic tracking" (where the electron simply follows the field) from the "true transition" (where the electron jumps to a new state). 
> **Key Insight:** The nonadiabatic correction ensures that if you turn a field on and then off very slowly, the system returns to its original state with **zero** net transition probability—a result the standard Dirac first-order theory often fails to guarantee.

---

### **Summary of the Hunt-Mandal Result**
They concluded that $|c_k(t)|^2$ is a **mathematical convenience**, not a physical observable. The physical reality is found in the **flux of probability**, which is only correctly captured by the nonadiabatic terms. This resolved a long-standing paradox where quantum transition theory seemed to conflict with classical electrodynamics.

To derive the power absorption formula and show why the standard Dirac approach fails to maintain energy balance, we start with the fundamental definition of energy in a time-dependent quantum system.

### **1. The Physical Power Absorption**
The average energy of a system at time $t$ is the expectation value of the full Hamiltonian $H(t)$:
$$\langle E(t) \rangle = \langle \Psi(t) | H(t) | \Psi(t) \rangle$$
The physical **power absorption** $P(t)$ is defined as the time rate of change of this energy. Using the Schrödinger equation, it can be shown that:
$$P(t) = \frac{d}{dt} \langle \Psi(t) | H(t) | \Psi(t) \rangle = \left\langle \Psi(t) \left| \frac{\partial H(t)}{\partial t} \right| \Psi(t) \right\rangle$$
Since the unperturbed Hamiltonian $H_0$ is time-independent, this reduces to:
$$P(t) = \left\langle \Psi(t) \left| \frac{\partial V(t)}{\partial t} \right| \Psi(t) \right\rangle$$
This is the **exact physical power**. It depends only on the rate at which the external field $V(t)$ is changing.

---

### **2. The Dirac Transition Rate**
In the standard Dirac approach, we expand $|\Psi(t)\rangle$ in the eigenstates $|n\rangle$ of $H_0$:
$$|\Psi(t)\rangle = \sum_n c_n(t) e^{-i E_n t / \hbar} |n\rangle$$
If we attempt to calculate the "power absorbed" by summing the rates at which the system moves into higher-energy unperturbed states, we get:
$$P_{\text{Dirac}}(t) = \sum_n E_n \frac{d}{dt} |c_n(t)|^2$$
**The Hunt-Mandal Inconsistency:**
Hunt and Mandal demonstrated that $P(t) \neq P_{\text{Dirac}}(t)$. Specifically, if you use the standard first-order transition coefficients $c_n(t)$, you find:
$$P_{\text{Dirac}}(t) = \left\langle \Psi(t) \left| \frac{\partial V}{\partial t} \right| \Psi(t) \right\rangle + \text{Oscillatory Terms}$$
These "oscillatory terms" do not disappear even if the field becomes constant ($\partial V / \partial t = 0$). This implies the system is still "absorbing energy" from a static field, which violates the law of conservation of energy.



---

### **3. The Nonadiabatic Correction**
To fix this, Hunt and Mandal introduce a **nonadiabatic decomposition** of the transition amplitude. They redefine the coefficient $c_k(t)$ by splitting it into an **adiabatic** (tracking) part and a **nonadiabatic** (transition) part:
$$c_k(t) = c_k^{\text{ad}}(t) + c_k^{\text{nad}}(t)$$
* **$c_k^{\text{ad}}(t)$:** Represents the electron simply "adjusting" its cloud to the presence of the field (polarization).
* **$c_k^{\text{nad}}(t)$:** Represents the actual quantum jump.

When they recalculate the power using only the **nonadiabatic flux**, the unphysical oscillations cancel out perfectly. The energy balance is restored:
$$\frac{d \langle E \rangle}{dt} = \sum_n E_n \frac{d}{dt} |c_n^{\text{nad}}(t)|^2$$

### **4. Why this is the "Correct" Form**
The nonadiabatic form is considered physically correct because it satisfies three critical criteria that the standard Dirac form misses:
1.  **Work-Energy Theorem:** It ensures that the energy change in the system equals the work done by the external field.
2.  **The Adiabatic Limit:** If the field is turned on and off infinitely slowly, the nonadiabatic coefficients remain zero, ensuring the system returns to its ground state (as thermodynamics requires).
3.  **Gauge Invariance:** The nonadiabatic transition probability gives the same result whether you use the "length gauge" ($r \cdot E$) or the "velocity gauge" ($p \cdot A$), whereas the standard Dirac form often gives different answers.

---

### **The Resulting Physical Picture**
The standard Dirac coefficient $|c_k(t)|^2$ actually measures a mix of **real transitions** and **temporary polarization**. By separating them using the Landau-Lifshitz/Hunt-Mandal logic, we find that the "transition probability" we measure in the lab is only the nonadiabatic component.
