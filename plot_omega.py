import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def modulo(a, b, c):
    return np.sqrt(pow(a,2)+pow(b,2)+pow(c,2))


dati_euler = pd.read_csv("omega_euler.csv", sep=';')
dati_leapfrog = pd.read_csv("omega_leapfrog.csv", sep=';')
# plt.plot(modulo(dati_rk4.omega1x, dati_rk4.omega1y, dati_rk4.omega1z)-modulo(dati_rk4.OMEGAx, dati_rk4.OMEGAy, dati_rk4.OMEGAz), label='leapfrog')
plt.plot(modulo(dati_leapfrog.omega1x, dati_leapfrog.omega1y, dati_leapfrog.omega1z)-modulo(dati_leapfrog.OMEGAx, dati_leapfrog.OMEGAy, dati_leapfrog.OMEGAz), label='leapfrog')
plt.plot(modulo(dati_euler.omega1x, dati_euler.omega1y, dati_euler.omega1z)-modulo(dati_euler.OMEGAx, dati_euler.OMEGAy, dati_euler.OMEGAz), label='euler')
# plt.plot(modulo(dati_rk4.omega2x, dati_rk4.omega2y, dati_rk4.omega2z))
plt.legend()
plt.show()