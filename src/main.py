import numpy as np
import matplotlib.pyplot as plt
from system import System
from signals import PlateauPulse
from bath import ThermalBath
from simulator import Simulator

def main():
    # 1. Define the unperturbed 3-Level System (Transmon Qutrit parameters)
    qutrit = System(
        energy_levels=[0.0, 5.0, 9.8], # GHz
        dipole_matrix=[
            [0.0, 1.0, 0.0],
            [1.0, 0.0, 1.4],
            [0.0, 1.4, 0.0]
        ]
    )

    # 2. Define the Time-Dependent Drive
    pulse = PlateauPulse(
        amplitude=0.15,        
        carrier_freq=5.05,     
        ramp_time=10.0,        
        plateau_time=100.0     
    )

    # 3. Define the Dissipative Environment (Calibrated transmon times)
    bath = ThermalBath(
        temperature=0.015,     
        t1_times={(1, 0): 50.0, (2, 1): 40.0},               
        t2_times={(1, 0): 30.0, (2, 1): 25.0, (2, 0): 20.0}   
    )

    # 4. Initialize and Run
    times = np.linspace(0, 150, 2000)
    simulator = Simulator(system=qutrit, drive=pulse, bath=bath)

    print("Running Standard Lindblad...")
    res_lindblad = simulator.run(times, method='lindblad')

    print("Running Nonadiabatic CGME Structure...")
    res_cgme = simulator.run(
        times, 
        method='nonadiabatic_cgme', 
        coarse_grain_timescale=0.5
    )

    # 5. Plot Results
    plt.figure(figsize=(10, 6))
    res_lindblad.plot_populations(states=[2], label_prefix="Lindblad")
    res_cgme.plot_populations(states=[2], label_prefix="CGME")
    plt.title("Qutrit |2> State Population During Plateau Pulse")
    plt.show()

if __name__ == "__main__":
    main()