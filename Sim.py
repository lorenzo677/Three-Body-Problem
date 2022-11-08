from scipy.integrate import odeint
import numpy as np
import matplotlib.pyplot as plt


DIM = 3
G = 10
N_BODIES = 3
N_STEPS = 8000

def distance(r1, r2):
    np.sqrt(np.pow(r1[0]-r2[0],2)+np.pow(r1[1]-r2[1],2)+np.pow(r1[2]-r2[2],2))


class Planet:

    def __init__(self, mass, position, velocity):
        
        self.mass = mass
        self.positionx = position[0]
        self.positiony = position[1]
        self.positionz = position[2]
        self.velocityx = velocity[0]
        self.velocityy = velocity[1]
        self.velocityz = velocity[2]
 

def acceleration(A,  B,  C, axe):
    
    mass_A = A.mass
    mass_B = B.mass
    posx_A = A.positionx 
    posx_B = B.positionx 
    posx_C = C.positionx 
    posy_A = A.positiony 
    posy_B = B.positiony 
    posy_C = C.positiony 
    posz_A = A.positionz 
    posz_B = B.positionz 
    posz_C = C.positionz
    
    if (axe == 0):
        return (-1 * G * (mass_A * (posx_C-posx_A) / np.pow(np.sqrt(np.pow(posx_C-posx_A,2)+np.pow(posy_C-posy_A,2)+np.pow(posz_C-posz_A,2)), 3) + mass_B * (posx_C-posx_B) / np.pow(np.sqrt(np.pow(posx_C-posx_B,2)+np.pow(posy_C-posy_B,2)+np.pow(posz_C-posz_B,2)), 3)))
    elif (axe == 1):
        return (-1 * G * (mass_A * (posy_C-posy_A) / np.pow(np.sqrt(np.pow(posx_C-posx_A,2)+np.pow(posy_C-posy_A,2)+np.pow(posz_C-posz_A,2)), 3) + mass_B * (posy_C-posy_B) / np.pow(np.sqrt(np.pow(posx_C-posx_B,2)+np.pow(posy_C-posy_B,2)+np.pow(posz_C-posz_B,2)), 3)))
    elif (axe == 2):
        return (-1 * G * (mass_A * (posz_C-posz_A) / np.pow(np.sqrt(np.pow(posx_C-posx_A,2)+np.pow(posy_C-posy_A,2)+np.pow(posz_C-posz_A,2)), 3) + mass_B * (posz_C-posz_B) / np.pow(np.sqrt(np.pow(posx_C-posx_B,2)+np.pow(posy_C-posy_B,2)+np.pow(posz_C-posz_B,2)), 3)))
    

# HO SCAZZATO
def F(x, v, t, A, B, C, j ):

    mass_A = A.mass
    mass_B = B.mass
    posx_A = A.positionx 
    posx_B = B.positionx 
    posx_C = C.positionx 
    posy_A = A.positiony 
    posy_B = B.positiony 
    posy_C = C.positiony 
    posz_A = A.positionz 
    posz_B = B.positionz 
    posz_C = C.positionz

    if (j == 0):
        return (-1 * G * (mass_A * (x-posx_A) / pow(np.sqrt(pow(x-A.x,2)+pow(posy_C-A.y,2)+pow(posz_C-A.z,2)), 3) + mass_B * (x-posx_B) / pow(np.sqrt(pow(x-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)))
    elif (j == 1):
        return (-1 * G * (mass_A * (x-A.y) / pow(np.sqrt(pow(posx_C-A.x,2)+pow(x-A.y,2)+pow(posz_C-A.z,2)), 3) + mass_B * (x-posy_B) / pow(np.sqrt(pow(posx_C-posx_B,2)+pow(x-posy_B,2)+pow(posz_C-posz_B,2)), 3)))
    elif (j == 2): 
        return (-1 * G * (mass_A * (x-A.z) / pow(np.sqrt(pow(posx_C-A.x,2)+pow(posy_C-A.y,2)+pow(x-A.z,2)), 3) + mass_B * (x-posz_B) / pow(np.sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(x-posz_B,2)), 3)))
        

h = 0.001

A = Planet(10, (-10, 10, -11), (-3, 0, 0))  
B = Planet(10, (0, 0, 0), (0, 0, 0))
C = Planet(10, (10, 14, 12), (0, 0, 0))

time = np.arange(0,80000, N_STEPS)





mass_A = A.mass
mass_B = B.mass
mass_C = C.mass

for  i in range(len(time)):
    time[i]= i * h


    

def x_derivatives(x, t, A,  B,  C ):
    return [x[1], (-1 * G * (A.mass * (x[0]-A.positionx) / pow(np.sqrt(pow(x[0]-A.positionx,2)+pow(C.positiony-A.positiony,2)+pow(C.positionz-A.positionz,2)), 3) + B.mass * (x[0]-B.positionx) / pow(np.sqrt(pow(x[0]-B.positionx,2)+pow(C.positiony-B.positiony,2)+pow(C.positionz-B.positionz,2)), 3))), x[1], (-1 * G * (A.mass * (x[1]-A.positiony) / pow(np.sqrt(pow(C.positionx-A.positionx,2)+pow(x[1]-A.positiony,2)+pow(C.positionz-A.positionz,2)), 3) + B.mass * (x[1]-B.positiony) / pow(np.sqrt(pow(C.positionx-B.positionx,2)+pow(x[1]-B.positiony,2)+pow(C.positionz-B.positionz,2)), 3))), x[1], (-1 * G * (A.mass * (x[2]-A.positionz) / pow(np.sqrt(pow(C.positionx-A.positionx,2)+pow(C.positiony-A.positiony,2)+pow(x[2]-A.positionz,2)), 3) + B.mass * (x[2]-B.positionz) / pow(np.sqrt(pow(C.positionx-B.positionx,2)+pow(C.positiony-B.positiony,2)+pow(x[2]-B.positionz,2)), 3)))]


for i in range(80000):
    positionsAx, velocityAx, positionsAy, velocityAy, positionsAz, velocityAz = odeint(x_derivatives, (A.positionx, A.velocityy, A.positiony, A.velocityy, A.positionz, A.velocityz),[0, 0.01], args=((C, B, A))).T

    positionsBx, velocityBx, positionsBy, velocityBy, positionsBz, velocityBz = odeint(x_derivatives, (B.positionx, B.velocityy, B.positiony, B.velocityy, B.positionz, B.velocityz),[0, 0.01], args=((A, C, B))).T

    positionsCx, velocityCx, positionsCy, velocityCy, positionsCz, velocityCz = odeint(x_derivatives, (C.positionx, C.velocityy, C.positiony, C.velocityy, C.positionz, C.velocityz),[0, 0.01], args=((A, B, C))).T


    A.positionx = positionsAx
    A.positiony = positionsAy
    A.positionz = positionsAz

    B.positionx = positionsBx
    B.positiony = positionsBy
    B.positionz = positionsBz

    C.positionx = positionsCx
    C.positiony = positionsCy
    C.positionz = positionsCz

print(positionsAx, positionsAy, positionsAz)

output_file_A=open("positions_A.csv", "w")
output_file_B=open("positions_B.csv", "w")
output_file_C=open("positions_C.csv", "w")

output_file_A.write("x;y;z;t")
output_file_B.write("x;y;z;t")
output_file_C.write("x;y;z;t")
    
# for i in range(N_STEPS):
#     output_file_A.write( f"{x_A[0][i]};{x_A[1][i]}; {x_A[2][i]};{time[i]}")
#     output_file_A.write( f"{x_B[0][i]};{x_B[1][i]}; {x_B[2][i]};{time[i]}")
#     output_file_A.write( f"{x_C[0][i]};{x_C[1][i]}; {x_C[2][i]};{time[i]}")
  
output_file_A.close()
output_file_B.close() 
output_file_C.close()
