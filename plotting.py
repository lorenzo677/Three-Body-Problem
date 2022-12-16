#%%
import pandas as pd
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import matplotlib.pyplot as plt
import numpy as np
import sys

#%%

file_positions = pd.read_csv("positions_"+sys.argv[1]+".csv", sep=';')

t = range(0, len(file_positions.xA))

xA, yA, zA = file_positions.xA, file_positions.yA,file_positions.zA
xB, yB, zB = file_positions.xB, file_positions.yB,file_positions.zB
xC, yC, zC = file_positions.xC, file_positions.yC,file_positions.zC

dataSetA = np.array([xA, yA, zA])  # Combining our position coordinates
dataSetB = np.array([xB, yB, zB])
dataSetC = np.array([xC, yC, zC])
numDataPoints = len(t)

def animate_func(num):
    x = 150 #x as speed multiplier
    num =  num*x
   
    if num < len(file_positions.xA):
        ax.clear()  # Clears the figure to update the line, point,   
                    # title, and axes
        # Updating Trajectory Line (num+1 due to Python indexing)
        ax.plot3D(dataSetA[0, :num+1], dataSetA[1, :num+1], 
                dataSetA[2, :num+1], c='blue', label='A')
        ax.plot3D(dataSetB[0, :num+1], dataSetB[1, :num+1], 
                dataSetB[2, :num+1], c='green', label='B')
        ax.plot3D(dataSetC[0, :num+1], dataSetC[1, :num+1], 
                dataSetC[2, :num+1], c='red', label='C')
        # Updating Point Location 
        ax.scatter(dataSetA[0, num], dataSetA[1, num], dataSetA[2, num], 
                c='blue', marker='o')
        ax.scatter(dataSetB[0, num], dataSetB[1, num], dataSetB[2, num], 
                c='green', marker='o')
        ax.scatter(dataSetC[0, num], dataSetC[1, num], dataSetC[2, num], 
                c='red', marker='o')
        # Adding Constant Origin
        ax.plot3D(dataSetA[0, 0], dataSetA[1, 0], dataSetA[2, 0],     
                c='black', marker='o')
        ax.plot3D(dataSetB[0, 0], dataSetB[1, 0], dataSetB[2, 0],     
                c='black', marker='o')
        ax.plot3D(dataSetC[0, 0], dataSetC[1, 0], dataSetC[2, 0],     
                c='black', marker='o')
    # Setting Axes Limits
#     ax.set_xlim3d([-100, 20])
#     ax.set_ylim3d([-20, 30])
#     ax.set_zlim3d([-40, 20])

    # Adding Figure Labels
    ax.set_title('Trajectory - ' + sys.argv[1]+ '\nTime = ' + str(np.round(t[num],    
                 decimals=2)) + ' sec')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

fig = plt.figure()

#plt.style.use('Solarize_Light2')
#print(plt.style.available)
ax = plt.axes(projection='3d')
line_ani = FuncAnimation(fig, animate_func, interval=1, frames=int(numDataPoints))
# plt.legend()
plt.show()
# Saving the Animation
f = r"animate_func.gif"
writergif = PillowWriter(fps=20)
# line_ani.save(f, writer=writergif, dpi=300)

# ax.plot3D(file_positions.xA, file_positions.yA,file_positions.zA)
# ax.plot3D(file_positions.xB, file_positions.yB,file_positions.zB)
# ax.plot3D(file_positions.xC, file_positions.yC,file_positions.zC)



# # %%
# plt.plot(file_positions.xA,file_positions.yA, '-o', markersize='1')
# plt.plot(file_positions.xB,file_positions.yB, '-o', markersize='1')
# plt.title(sys.argv[1])
# plt.plot(file_positions.xC,file_positions.yC, '-o', markersize='1')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# plt.plot(file_positions.zA,file_positions.xA, '-o', markersize='1')
# plt.plot(file_positions.zB,file_positions.xB, '-o', markersize='1')
# plt.plot(file_positions.zC,file_positions.xC, '-o', markersize='1')
# plt.xlabel('z')
# plt.ylabel('x')
# plt.show()

# plt.plot(file_positions.yA,file_positions.zA, '-o', markersize='1')
# plt.plot(file_positions.yB,file_positions.zB, '-o', markersize='1')
# plt.plot(file_positions.yC,file_positions.zC, '-o', markersize='1')
# plt.xlabel('y')
# plt.ylabel('z')
# plt.show()
# # # %%
# plt.plot(file_positions.xA, '-o')
# plt.plot(file_positions.xB, '-o')
# plt.plot(file_positions.xC, '-o')
# plt.xlabel('time')
# plt.ylabel('x')
# plt.show()

# %%
