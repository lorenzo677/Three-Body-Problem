import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# Read in the data from the file
file = pd.read_csv("/Users/lorenzo/Downloads/final_beta1e8.csv") 

# Split the data into the x, y, and z coordinates for each body at each time
x1, y1, z1, x2, y2, z2, x3, y3, z3 =  file.x1, file.y1, file.z1, file.x2, file.y2, file.z2, file.x3,  file.y3,  file.z3   

def center_of_mass(pos1, pos2):
    total_mass = 1 + 1  # assume equal masses for simplicity
    com = (pos1 + pos2) / total_mass
    return com

# Calculate the angles for each instant of time
angles = []
for i in range(len(x1)):
    pos1 = np.array([x1[i], y1[i], z1[i]])
    pos2 = np.array([x2[i], y2[i], z2[i]])
    pos3 = np.array([x3[i], y3[i], z3[i]])
    com = center_of_mass(pos2, pos3)
    vec1 = com - pos1
    vec2 = pos3 - pos2
    angle = np.arctan2(np.cross(vec1, vec2), np.dot(vec1, vec2))
    angles.append(angle)

# Plot the angles as a function of time
sns.set_style('darkgrid')
plt.plot(file.t, np.rad2deg(angles))
plt.xlabel('Time [years]', fontsize=18)
plt.ylabel('Angle [degree]', fontsize=18)
plt.subplots_adjust(top=0.96, right=0.98)
plt.savefig('images/angle.png', dpi=400)
plt.show()