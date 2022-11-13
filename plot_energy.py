import matplotlib.pyplot as plt
import pandas as pd

dati_euler = pd.read_csv("Total_energy_euler.csv")
dati_rk4 = pd.read_csv("Total_energy_rk4.csv")
dati_verlet = pd.read_csv("Total_energy_verlet.csv")

plt.plot(dati_euler.energy, label='euler')
plt.plot(dati_verlet.energy, label='verlet')
plt.plot(dati_rk4.energy, label='rk4')

plt.legend()
plt.ylabel("total energy")
plt.xlabel("Step")
plt.show()