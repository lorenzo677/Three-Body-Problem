import pandas as pd
import matplotlib.pyplot as plt

file = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output.csv")

plt.plot(file.t, file.energy)
plt.title('energy')
plt.show()