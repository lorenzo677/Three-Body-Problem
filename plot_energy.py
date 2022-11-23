import matplotlib.pyplot as plt
import pandas as pd

dati_euler = pd.read_csv("Total_energy_euler.csv",sep=';')
dati_rk4 = pd.read_csv("Total_energy_rk4.csv",sep=';')
dati_verlet = pd.read_csv("Total_energy_verlet.csv",sep=';')

#plt.plot(dati_euler.k+dati_euler.p, label='euler')
# plt.plot(dati_verlet.k+dati_verlet.p, label='verlet')
plt.plot([i*0.008 for i in range(len(dati_rk4.k))] ,dati_rk4.k+dati_rk4.p, label='rk4')

plt.legend()
plt.ylabel("total energy")

plt.xlabel("Step")
plt.show()