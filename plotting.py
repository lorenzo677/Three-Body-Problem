#%%
import pandas as pd
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

#%%

file_positions_A = pd.read_csv("positions_A.csv", sep=';')
file_positions_B = pd.read_csv("positions_B.csv", sep=';')
file_positions_C = pd.read_csv("positions_C.csv", sep=';')

#%%
def update(j):
    a,b,c = file_positions_A.x, file_positions_A.y,file_positions_A.z
    return( a, b, c )



fig = plt.figure()
 
# syntax for 3-D projection
ax = plt.axes(projection ='3d')
line, = ax.plot3D([0,0,0], [0,0,0], [0,0,0], lw=3)
def init():
    line.set_data([], [])
    return line,
def animate(i):
    for i in range(len(file_positions_A.x)):
        x = file_positions_A.x[i]
        y = file_positions_A.y[i]
        z = file_positions_A.z[i]
    yield (x, y, z)
# ax.plot3D(file_positions_A.x, file_positions_A.y,file_positions_A.z)
# ax.plot3D(file_positions_B.x, file_positions_B.y,file_positions_B.z)
# ax.plot3D(file_positions_C.x, file_positions_C.y,file_positions_C.z)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

anim = FuncAnimation(fig, animate, init_func=init,
                               frames=200, interval=20, blit=True)


anim.save('trajectory.gif', writer='imagemagick')

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
