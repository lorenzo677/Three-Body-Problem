#%%
import pandas as pd
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import matplotlib.pyplot as plt
import numpy as np

#%%

file_positions_A = pd.read_csv("positions_A.csv", sep=';')
file_positions_B = pd.read_csv("positions_B.csv", sep=';')
file_positions_C = pd.read_csv("positions_C.csv", sep=';')

t = range(0, 80000)

xA, yA, zA = file_positions_A.x, file_positions_A.y,file_positions_A.z
xB, yB, zB = file_positions_B.x, file_positions_B.y,file_positions_B.z
xC, yC, zC = file_positions_C.x, file_positions_C.y,file_positions_C.z

dataSetA = np.array([xA, yA, zA])  # Combining our position coordinates
dataSetB = np.array([xB, yB, zB])
dataSetC = np.array([xC, yC, zC])
numDataPoints = len(t)

def animate_func(num):
    x = 200 #x as speed multiplier
    num = num*x
    if num < 79900:
        ax.clear()  # Clears the figure to update the line, point,   
                    # title, and axes
        # Updating Trajectory Line (num+1 due to Python indexing)
        ax.plot3D(dataSetA[0, :num+1], dataSetA[1, :num+1], 
                dataSetA[2, :num+1], c='blue')
        ax.plot3D(dataSetB[0, :num+1], dataSetB[1, :num+1], 
                dataSetB[2, :num+1], c='green')
        ax.plot3D(dataSetC[0, :num+1], dataSetC[1, :num+1], 
                dataSetC[2, :num+1], c='red')
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
    ax.set_xlim3d([-100, 20])
    ax.set_ylim3d([-20, 30])
    ax.set_zlim3d([-40, 20])

    # Adding Figure Labels
    ax.set_title('Trajectory \nTime = ' + str(np.round(t[num],    
                 decimals=2)) + ' sec')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

fig = plt.figure()
#plt.style.use('Solarize_Light2')
#print(plt.style.available)
ax = plt.axes(projection='3d')
line_ani = FuncAnimation(fig, animate_func, interval=1, frames=int(numDataPoints/200))
#plt.savefig("animation.gif", dpi = 300)
plt.show()
# Saving the Animation
f = r"animate_func.gif"
writergif = PillowWriter(fps=20)
#line_ani.save(f, writer=writergif, dpi=300)

# ax.plot3D(file_positions_A.x, file_positions_A.y,file_positions_A.z)
# ax.plot3D(file_positions_B.x, file_positions_B.y,file_positions_B.z)
# ax.plot3D(file_positions_C.x, file_positions_C.y,file_positions_C.z)



# %%
plt.plot(file_positions_A.x,file_positions_A.y, '-o')
plt.plot(file_positions_B.x,file_positions_B.y, '-o')
plt.plot(file_positions_C.x,file_positions_C.y, '-o')
plt.show()

plt.plot(file_positions_A.z,file_positions_A.x, '-o')
plt.plot(file_positions_B.z,file_positions_B.x, '-o')
plt.plot(file_positions_C.z,file_positions_C.x, '-o')
plt.show()

plt.plot(file_positions_A.y,file_positions_A.z, '-o')
plt.plot(file_positions_B.y,file_positions_B.z, '-o')
plt.plot(file_positions_C.y,file_positions_C.z, '-o')
plt.show()
# %%
plt.plot(file_positions_A.x[0:300], '-o')
plt.plot(file_positions_B.x[0:300], '-o')
plt.plot(file_positions_C.x[0:300], '-o')
plt.plot()

# %%
