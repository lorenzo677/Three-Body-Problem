#%%
import pandas as pd
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

#%%

file_positions_A = pd.read_csv("positions_A.csv", ';')
file_positions_B = pd.read_csv("positions_B.csv", ';')
file_positions_C = pd.read_csv("positions_C.csv", ';')

#%%

fig = plt.figure()
 
# syntax for 3-D projection
ax = plt.axes(projection ='3d')

ax.plot3D(file_positions_A.x, file_positions_A.y,file_positions_A.z)
ax.plot3D(file_positions_B.x, file_positions_B.y,file_positions_B.z)
ax.plot3D(file_positions_C.x, file_positions_C.y,file_positions_C.z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()
# %%
