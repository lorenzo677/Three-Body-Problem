#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
#%%

file_odeint = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/odeint/output.csv", sep=',')
file_rk4 = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/positions_rk4.csv", sep=';')
file_leapfrog = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/positions_leapfrog.csv", sep=';')
file_euler = pd.read_csv("/Users/lorenzo/Desktop/Three-Body-Problem/positions_euler.csv", sep=';')
#%%
# plt.plot(np.sqrt((file_odeint.x1-file_rk4.xA)**2 + (file_odeint.x2-file_rk4.yA)**2 + (file_odeint.x3-file_rk4.zA)**2), label = "Rk4 " )
# plt.plot(np.sqrt((file_odeint.x1-file_leapfrog.xA)**2 + (file_odeint.x2-file_leapfrog.yA)**2 + (file_odeint.x3-file_leapfrog.zA)**2), label = "leap" )
# plt.plot(np.sqrt((file_odeint.x1-file_euler.xA)**2 + (file_odeint.x2-file_euler.yA)**2 + (file_odeint.x3-file_euler.zA)**2), label = "euler" )

plt.plot(file_odeint.x1, file_odeint.y1, label='ode')
plt.plot(file_euler.xA, file_euler.yA, label='euler')
plt.plot(file_rk4.xA, file_rk4.yA, label='rk4')
plt.plot(file_leapfrog.xA, file_leapfrog.yA, label='lea')

plt.legend()
plt.show()