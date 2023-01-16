#include <iostream>
#include <array>
#include <fstream>

static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 300000;


std::array<double, 7> A = {100, -10, 10, -11, -3, 0, 0};
std::array<double, 7> B = {100, 0, 0, 0, 3, 0, 0};
std::array<double, 7> C = {0, 10, 14, 12, 3, 0, 0};
std::array<std::array<double, 7>,  3> bodies = {A, B, C};

std::array<double, 3> m = {bodies[0][0], bodies[1][0], bodies[2][0]};
std::array<double, 3> x = {bodies[0][1], bodies[0][2], bodies[0][3]};
std::array<double, 3> y = {bodies[1][1], bodies[1][2], bodies[1][3]};
std::array<double, 3> z = {bodies[2][1], bodies[2][2], bodies[2][3]};

std::array<double, 3> Vx = {bodies[0][4], bodies[0][5], bodies[0][6]};
std::array<double, 3> Vy = {bodies[1][4], bodies[1][5], bodies[1][6]};
std::array<double, 3> Vz = {bodies[2][4], bodies[2][5], bodies[2][6]};

std::array<std::array<double, 3>, 3> x0 = {x, y, z};
std::array<std::array<double, 3>, 3> v0 = {Vx, Vy, Vz};

std::array<std::array<double, 3>, 3> acc(std::array<std::array<double, 3>, 3> x){
    std::array<std::array<std::array<double, 3>, 3>, 3> dx;
    
}



//def acc(x):
    dx = x[None,:,:]-x[:,None,:]
    r3 = np.sum(dx**2,axis=2)**1.5
    return -G*np.sum(m[:,None,None]*dx/(1e-40+r3[:,:,None]), axis=0)
    

def sys(t,u):
    x,v = u.reshape([2,-1,3])
    return np.reshape([v,acc(x)],-1)


def kineticEnergy(m,v):
    imp = 0.5*m * np.sum(v**2,axis=1).T
    return np.sum(imp.T, axis=0)

def potentialEnergy(m1,x1,m2,x2):
    dx = x2-x1
    r = np.sum(dx**2,axis=0)**0.5
    return -G*m1*m2/r