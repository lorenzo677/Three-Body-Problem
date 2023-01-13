import pandas as pd
import matplotlib.pyplot as plt

file = pd.read_csv("/Users/lorenzo/Desktop/output.csv")

plt.plot(file.t, file.energy)
plt.title('energy')
plt.show()