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
        

h = 0.11

A = Planet(10, (-10, 10, -11), (-3, 0, -0.4))  
B = Planet(10, (0, 0, 0), (0, 0, 0))
C = Planet(10, (10, 14, 12), (0, 0, 0.2))

time = np.arange(0,N_STEPS*h, N_STEPS)

mass_A = A.mass
mass_B = B.mass
mass_C = C.mass

for  i in range(len(time)):
    time[i]= i * h


def x_derivatives(x, t, A,  B,  C ):
    return [x[1], (-1 * G * (A.mass * (x[0]-A.positionx) / pow(np.sqrt(pow(x[0]-A.positionx,2)+pow(C.positiony-A.positiony,2)+pow(C.positionz-A.positionz,2)), 3) + B.mass * (x[0]-B.positionx) / pow(np.sqrt(pow(x[0]-B.positionx,2)+pow(C.positiony-B.positiony,2)+pow(C.positionz-B.positionz,2)), 3))), x[1], (-1 * G * (A.mass * (x[1]-A.positiony) / pow(np.sqrt(pow(C.positionx-A.positionx,2)+pow(x[1]-A.positiony,2)+pow(C.positionz-A.positionz,2)), 3) + B.mass * (x[1]-B.positiony) / pow(np.sqrt(pow(C.positionx-B.positionx,2)+pow(x[1]-B.positiony,2)+pow(C.positionz-B.positionz,2)), 3))), x[1], (-1 * G * (A.mass * (x[2]-A.positionz) / pow(np.sqrt(pow(C.positionx-A.positionx,2)+pow(C.positiony-A.positiony,2)+pow(x[2]-A.positionz,2)), 3) + B.mass * (x[2]-B.positionz) / pow(np.sqrt(pow(C.positionx-B.positionx,2)+pow(C.positiony-B.positiony,2)+pow(x[2]-B.positionz,2)), 3)))]

output_file_A=open("positions_A.csv", "w")
output_file_B=open("positions_B.csv", "w")
output_file_C=open("positions_C.csv", "w")

output_file_A.write("x;y;z;t\n")
output_file_B.write("x;y;z;t\n")
output_file_C.write("x;y;z;t\n")
output_file_A.write( f"{A.positionx}; {A.positiony}; {A.positionz}; {0}\n")
output_file_B.write( f"{B.positionx}; {B.positiony}; {B.positionz}; {0}\n")
output_file_C.write( f"{C.positionx}; {C.positiony}; {C.positionz}; {0}\n")

for i in range(N_STEPS):
    positionsAx, velocityAx, positionsAy, velocityAy, positionsAz, velocityAz = odeint(x_derivatives, (A.positionx, A.velocityy, A.positiony, A.velocityy, A.positionz, A.velocityz),[0, h*i], args=((C, B, A))).T

    positionsBx, velocityBx, positionsBy, velocityBy, positionsBz, velocityBz = odeint(x_derivatives, (B.positionx, B.velocityy, B.positiony, B.velocityy, B.positionz, B.velocityz),[0, h*i], args=((A, C, B))).T

    positionsCx, velocityCx, positionsCy, velocityCy, positionsCz, velocityCz = odeint(x_derivatives, (C.positionx, C.velocityy, C.positiony, C.velocityy, C.positionz, C.velocityz),[0, h*i], args=((A, B, C))).T

    A.positionx = positionsAx[1]
    A.positiony = positionsAy[1]
    A.positionz = positionsAz[1]

    B.positionx = positionsBx[1]
    B.positiony = positionsBy[1]
    B.positionz = positionsBz[1]

    C.positionx = positionsCx[1]
    C.positiony = positionsCy[1]
    C.positionz = positionsCz[1]

    A.velocityx = velocityAx[1]
    A.velocityy = velocityAy[1]
    A.velocityz = velocityAz[1]

    B.velocityx = velocityBx[1]
    B.velocityy = velocityBy[1]
    B.velocityz = velocityBz[1]

    C.velocityx = velocityCx[1]
    C.velocityy = velocityCy[1]
    C.velocityz = velocityCz[1]

    output_file_A.write( f"{A.positionx}; {A.positiony}; {A.positionz}; {h*i}\n")
    output_file_B.write( f"{B.positionx}; {B.positiony}; {B.positionz}; {h*i}\n")
    output_file_C.write( f"{C.positionx}; {C.positiony}; {C.positionz}; {h*i}\n")


output_file_A.close()
output_file_B.close() 
output_file_C.close()
