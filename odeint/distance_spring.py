import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

dati = pd.read_csv("/Users/lorenzo/Downloads/final_beta1e8.csv")
d = np.sqrt((dati.x2-dati.x3)**2+(dati.y2-dati.y3)**2+(dati.z2-dati.z3)**2)
plt.plot(d)
plt.show()