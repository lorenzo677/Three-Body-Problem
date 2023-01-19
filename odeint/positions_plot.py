#%%
import pandas as pd
from matplotlib.animation import FuncAnimation
from matplotlib.animation import PillowWriter
import matplotlib.pyplot as plt
import numpy as np

#%%

file_positions = pd.read_csv("/odeint/output.csv")

t = range(0, len(file_positions.x1))

x1, y1, z1 = file_positions.x1, file_positions.y1, file_positions.z1
x2, y2, z2 = file_positions.x2, file_positions.y2, file_positions.z2
x3, y3, z3 = file_positions.x3, file_positions.y3, file_positions.z3

dataSetA = np.array([x1, y1, z1])  # Combining our position coordinates
dataSetB = np.array([x2, y2, z2])
dataSetC = np.array([x3, y3, z3])
numDataPoints = len(t)
#%%
def animate_func(num):
    x = 200 #x as speed multiplier
    num =  num*x
   
    if num < len(file_positions.x1):
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
    ax.set_title('Trajectory')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

fig = plt.figure()

#plt.style.use('Solarize_Light2')
#print(plt.style.available)
ax = plt.axes(projection='3d')
line_ani = FuncAnimation(fig, animate_func, interval=1, frames=int(numDataPoints*100))
# plt.legend()
plt.show()
# Saving the Animation
f = r"animate_func.gif"
writergif = PillowWriter(fps=120)
# line_ani.save(f, writer=writergif, dpi=300)

# ax.plot3D(file_positions.x1, file_positions.y1,file_positions.z1)
# ax.plot3D(file_positions.x2, file_positions.y2,file_positions.z2)
# ax.plot3D(file_positions.x3, file_positions.y3,file_positions.z3)



# # %%
# plt.plot(file_positions.x1,file_positions.y1, '-o', markersize='1')
# plt.plot(file_positions.x2,file_positions.y2, '-o', markersize='1')
# plt.title(sys.argv[1])
# plt.plot(file_positions.x3,file_positions.y3, '-o', markersize='1')
# plt.xlabel('x')
# plt.ylabel('y')
# plt.show()

# plt.plot(file_positions.z1,file_positions.x1, '-o', markersize='1')
# plt.plot(file_positions.z2,file_positions.x2, '-o', markersize='1')
# plt.plot(file_positions.z3,file_positions.x3, '-o', markersize='1')
# plt.xlabel('z')
# plt.ylabel('x')
# plt.show()

# plt.plot(file_positions.y1,file_positions.z1, '-o', markersize='1')
# plt.plot(file_positions.y2,file_positions.z2, '-o', markersize='1')
# plt.plot(file_positions.y3,file_positions.z3, '-o', markersize='1')
# plt.xlabel('y')
# plt.ylabel('z')
# plt.show()
# # # %%
# plt.plot(file_positions.x1, '-o')
# plt.plot(file_positions.x2, '-o')
# plt.plot(file_positions.x3, '-o')
# plt.xlabel('time')
# plt.ylabel('x')
# plt.show()

# %%
### TO SAVE A SCREENSHOT OF THE TRAJECTORIES
# num = 100000
# fig = plt.figure()
# ax = plt.axes(projection='3d')

# ax.plot3D(dataSetA[0, :num+1], dataSetA[1, :num+1], dataSetA[2, :num+1], c='blue', label='A')
# ax.plot3D(dataSetB[0, :num+1], dataSetB[1, :num+1], dataSetB[2, :num+1], c='green', label='B')
# ax.plot3D(dataSetC[0, :num+1], dataSetC[1, :num+1], dataSetC[2, :num+1], c='red', label='C')
# # Updating Point Location 
# ax.plot3D(dataSetA[0, num], dataSetA[1, num], dataSetA[2, num], c='blue', marker='o', markersize=9)
# ax.plot3D(dataSetB[0, num], dataSetB[1, num], dataSetB[2, num], c='green', marker='o')
# ax.plot3D(dataSetC[0, num], dataSetC[1, num], dataSetC[2, num], c='red', marker='o')
# # Adding Constant Origin
# ax.plot3D(dataSetA[0, 0], dataSetA[1, 0], dataSetA[2, 0], c='black', marker='o')
# ax.plot3D(dataSetB[0, 0], dataSetB[1, 0], dataSetB[2, 0], c='black', marker='o')
# ax.plot3D(dataSetC[0, 0], dataSetC[1, 0], dataSetC[2, 0], c='black', marker='o')
# ax.set_title('Trajectories')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')
# plt.savefig('/Users/lorenzo/Desktop/traj.pdf', dpi=400)
# plt.show()
# %%
