import matplotlib.pyplot as plt
import pandas as pd
from math import sqrt
import numpy as np

def modulo(a, b, c):
    return np.sqrt(pow(a,2)+pow(b,2)+pow(c,2))


dati_rk4 = pd.read_csv("omega_rk4.csv", sep=';')

plt.plot(modulo(dati_rk4.omega1x, dati_rk4.omega1y, dati_rk4.omega1z)-modulo(dati_rk4.OMEGAx, dati_rk4.OMEGAy, dati_rk4.OMEGAz))
# plt.plot(modulo(dati_rk4.omega2x, dati_rk4.omega2y, dati_rk4.omega2z))
plt.show()