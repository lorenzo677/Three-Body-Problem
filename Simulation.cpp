/*
3D Three-Body-Problem simulation.
Created on November 2022.

@authors: Lorenzo Barsotti & Keivan Amini 
*/

#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include "integrators.h"
#include <string>


static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 70000;

// Spring
static constexpr int K_CONST = 1.0e4;
static constexpr double L0 = 2.8284271247; 

double distance(std::array<double, 3> r1, std::array<double, 3> r2){
    return sqrt(pow(r1[0]-r2[0],2)+pow(r1[1]-r2[1],2)+pow(r1[2]-r2[2],2));
}

class Planet{
public:

    double m;
    std::array <double, 3> x;
    std::array <double, 3> v;
    std::array <double, 3> a;
    double energy;
    double omega;
    
    Planet (double mass, double x_position, double y_position, double z_position, double x_velocity, double y_velocity, double z_velocity) {
        m = mass;
		x[0] = x_position;
		x[1] = y_position;
        x[2] = z_position;
		v[0] = x_velocity;
        v[1] = y_velocity;
        v[2] = z_velocity;
    };
    void computeKineticEnergy(){
        energy = 0.5 * m * (pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
    }
    
};

std::array<double, 3> differenceOfArrays(std::array<double, 3>  v1, std::array<double, 3> v2){
    std::array<double, 3> difference;
    for(int j=0;j<3;j++){
        difference[j] = v1[j]-v2[j];
    }
    return difference;
}

double computeGravitationalEnergy(Planet planet1, Planet planet2, Planet planet3){
    return -1 * G * ( planet3.m * planet1.m / distance(planet3.x, planet1.x) +  planet3.m * planet2.m / distance(planet3.x, planet2.x) + planet1.m * planet2.m / distance(planet1.x, planet2.x));
        
}

double computeElasticEnergy(Planet planet1, Planet planet2){
    return 0.5 * K_CONST * pow(distance(planet2.x, planet1.x)-L0,2);
}

std::array<double,3> computeVcm(Planet planet1, Planet planet2){
    std::array<double, 3> Vcm; // Center mass velocity between the two spring-bodies.
    for(int j=0;j<DIM-1;j++){
        Vcm[j] = (planet1.m * planet1.v[j] + planet2.m * planet2.v[j]) / (planet1.m + planet2.m);
    }
    return Vcm;
}
/*
std::array<double,3> computeOmegaSpring(Planet planet1, Planet planet2){
    std::array<double, 3> OmegaSpring; //  Angular velocity ω of planet1 with respect to the center of mass.
    std::array<double, 3> Vcm; // Center mass velocity between the two spring-bodies.
    std::array<double, 3> v1cm; // Radial velocity of planet1 with respect to the center of mass

    double distanza = distance(planet1.x,planet2.x);
    double R = 0.5 * distanza;
    
    Vcm = computeVcm(planet1,planet2);
    v1cm = differenceOfArrays(planet1.v, Vcm);

    for(int j=0;j<DIM-1;j++){
        OmegaSpring[j] = (v1cm[j]/R);
    }
    return OmegaSpring;
}

std::array<double,3> computeOmegaRevolution(Planet planet1, Planet planet2, Planet planet3){
    // planet1 and planet2 are the body with the spring
    std::array<double, 3> OmegaRevolution; // Angular velocity ω of the center mass of the two-body-spring system wrt the massive planet.
    std::array<double, 3> VelocityRevolution; // Radial velocity of the center mass of the two-body-spring system wrt the massive planet.
    std::array<double, 3> Rcm; // Center mass position between the two spring-bodies.
    std::array<double, 3> Vcm; // Center mass velocity between the two spring-bodies.

    for(int j=0;j<DIM-1;j++){
        Rcm[j] = (planet1.m * planet1.x[j] + planet2.m * planet2.x[j]) / (planet1.m + planet2.m);
    }
    double R = distance(Rcm, planet3.x); //true only if mass of third planet >> wrt masses of the other planets

    Vcm = computeVcm(planet1, planet2);
    VelocityRevolution = differenceOfArrays(Vcm, planet3.v);

    for(int j=0;j<DIM-1;j++){
        OmegaRevolution[j] = (VelocityRevolution[j]/R);
    }

    return OmegaRevolution;
}
*/
double computeDistance(double A, double B, double C, double D, Planet planet){
    double d = std::abs(A * planet.x[0] + B * planet.x[1] + C * planet.x[2] + D) / sqrt(A*A + B*B + C*C);
}

std::array<double, 3> computeCM(Planet planet1, Planet planet2, Planet planet3){
    std::array<double, 3> cm;
    for(int j=0;j<DIM-1;j++){
        cm[j] = (planet1.m * planet1.x[j] + planet2.m * planet2.x[j] + planet3.m * planet3.x[j]) / (planet1.m + planet2.m + planet3.m);
    }
    return cm;
}

std::array<double, 3> AngularMomentum(std::array<double, 3> cm, Planet planet){
    std::array<double, 3> L;
    std::array<double, 3> difference;
    difference = differenceOfArrays(cm, planet.x);
    L[0] = planet.m * (difference[1] * planet.v[2] - difference[2] * planet.v[1]);
    L[1] = planet.m * (difference[2] * planet.v[0] - difference[0] * planet.v[2]);
    L[2] = planet.m * (difference[0] * planet.v[1] - difference[1] * planet.v[0]);
    return L;
}


double springC(double x, double v, double t, Planet A, Planet B, Planet C, int axe){
    // Compute the gravitational + spring acceleration of the body C, specifying the axis.
    return (-1 * (G * (A.m * (C.x[axe]-A.x[axe]) / pow(distance(A.x, C.x), 3) + B.m * (C.x[axe]-B.x[axe]) / pow(distance(B.x, C.x), 3))) - (K_CONST / C.m) * (std::abs(distance(B.x, C.x))-L0) * (C.x[axe]-B.x[axe]) / (distance(B.x, C.x)));
}

double springB(double x, double v, double t, Planet A, Planet B, Planet C, int axe){
    // Compute the gravitational + spring acceleration of the body C, specifying the axis.
    return (-1 * (G * (A.m * (C.x[axe]-A.x[axe]) / pow(distance(A.x, C.x), 3) + B.m * (C.x[axe]-B.x[axe]) / pow(distance(B.x, C.x), 3))) - (K_CONST / C.m) * (std::abs(distance(A.x, C.x))-L0) * (C.x[axe]-A.x[axe]) / (distance(A.x, C.x)));
}

double acceleration(double x, double v, double t, Planet A, Planet B, Planet C, int axe){
    // Compute the gravitational acceleration of the body C, specifying the axis.
    return (-1 * G * (A.m * (C.x[axe]-A.x[axe]) / pow(distance(A.x, C.x), 3) + B.m * (C.x[axe]-B.x[axe]) / pow(distance(B.x, C.x), 3)));
}

int main(int argc, char** argv){
    
    double h = 0.002;

    Planet A(20, 0, 0, 0, 0, 2, 0);   // corpi allineati sull'asse delle x
    Planet B(1, -20, 10, 10, 1, 0, 0);
    Planet C(1, -20, 12, 12, -1, 0, 0);
    // Planet A(10, 0, 0, 0, -1, 0, 0);   // corpi allineati sull'asse delle x
    // Planet B(10, 20, 0, 0, 0, 0, 1);
    // Planet C(10, 15, 15, 10, 0, 2, 0);
    // Planet A(1, 0, 0, 0, -1, 0, 0);   // corpi allineati sull'asse delle x
    // Planet B(0.01, 20, 0, 0, 0, 0, 1);
    // Planet C(0.01, 15, 15, 10, 0, 2, 0);
    // Planet A(20, 0, 0, 0, 0, 1, 0);   // corpi allineati sull'asse delle x
    // Planet B(10, 0, 0 , 5, -5, 0, 5);
    // Planet C(10, 0, 0, -5, 5, 0, -5);

    // CONFIGURAZIONI BELLE

    //Planet A(0, -10, 0, 0, 0, -5, 0); 
    //Planet B(100, 0, 0, 0, 0, 5, 0);
    //Planet C(100, 13, 14, 12, 0, -5, 0);

    double x_A[DIM][3];
    double x_B[DIM][3];
    double x_C[DIM][3];

    double mass_A = A.m;
    double mass_B = B.m;
    double mass_C = C.m;

    x_A[0][0] = A.x[0];
    x_B[0][0] = B.x[0];
    x_C[0][0] = C.x[0];

    x_A[1][0] = A.x[1];
    x_B[1][0] = B.x[1];
    x_C[1][0] = C.x[1];

    x_A[2][0] = A.x[2];
    x_B[2][0] = B.x[2];
    x_C[2][0] = C.x[2];


    std::array<std::array<double, 3>, 6> initial_conditions = {A.x, A.v, B.x, B.v, C.x, C.v};
    // std::ofstream file_energy("Total_energy_" + std::string(argv[1]) + ".csv");
    // std::ofstream output_file_A("positions_A_" + std::string(argv[1]) + ".csv");
    // std::ofstream output_file_B("positions_B_" + std::string(argv[1]) + ".csv");
    // std::ofstream output_file_C("positions_C_" + std::string(argv[1]) + ".csv");
    // std::ofstream file_angmom("Total_angular_momentum_" + std::string(argv[1]) + ".csv");
    
    std::ofstream file_energy_euler("Total_energy_euler.csv");
    std::ofstream file_energy_leapfrog("Total_energy_leapfrog.csv");
    std::ofstream file_energy_rk4("Total_energy_rk4.csv");
    std::ofstream output_file_euler("positions_euler.csv");
    std::ofstream output_file_leapfrog("positions_leapfrog.csv");
    std::ofstream output_file_rk4("positions_rk4.csv");
    std::ofstream file_angmom_euler("Total_angular_momentum_euler.csv");
    std::ofstream file_angmom_leapfrog("Total_angular_momentum_leapfrog.csv");
    std::ofstream file_angmom_rk4("Total_angular_momentum_rk4.csv");
    std::ofstream file_distance_euler("distance_euler.csv");
    std::ofstream file_distance_leapfrog("distance_leapfrog.csv");
    std::ofstream file_distance_rk4("distance_rk4.csv");
    
    output_file_euler<<"xA;yA;zA;xB;yB;zB;xC;yC;zC"<<std::endl;
    file_energy_euler<<"k;g;e"<<std::endl;
    file_angmom_euler<<"Lx;Ly;Lz"<<std::endl;
    output_file_leapfrog<<"xA;yA;zA;xB;yB;zB;xC;yC;zC"<<std::endl;
    file_energy_leapfrog<<"k;g;e"<<std::endl;
    file_angmom_leapfrog<<"Lx;Ly;Lz"<<std::endl;
    output_file_rk4<<"xA;yA;zA;xB;yB;zB;xC;yC;zC"<<std::endl;
    file_energy_rk4<<"k;g;e"<<std::endl;
    file_angmom_rk4<<"Lx;Ly;Lz"<<std::endl;
    file_distance_euler<<"d"<<std::endl;
    file_distance_leapfrog<<"d"<<std::endl;
    file_distance_rk4<<"d"<<std::endl;
    
    double m1[4][4];
    double k1[4][4];
    double m2[4][4];
    double k2[4][4];
    double m3[4][4];
    double k3[4][4];
    double m4[4][4];
    double k4[4][4];
    
    std::array<double,3> vA = A.v; //condizione iniziale velocita
    std::array<double,3> xA = A.x;
    std::array<double,3> vB = B.v; //condizione iniziale velocita
    std::array<double,3> xB = B.x;
    std::array<double,3> vC = C.v; //condizione iniziale velocita
    std::array<double,3> xC = C.x;
    double t;
    std::array<double,3> cm = computeCM(A, B, C);

// ==========================================================
//                          EULER
// ==========================================================

for (int i=0; i<N_STEPS-1; i++){
    output_file_euler << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< ";";
    output_file_euler << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< ";";
    output_file_euler << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;
    for(int j=0; j<DIM-1; j++){

        A.a[j] = acceleration(0, 0, 0, B, C, A, j);
        B.a[j] = springB(0, 0, 0, C, A, B, j);
        C.a[j] = springC(0, 0, 0, A, B, C, j);

    }
    for(int j=0; j<DIM-1; j++){
        
        A.v[j] += A.a[j] * h;
        B.v[j] += B.a[j] * h;
        C.v[j] += C.a[j] * h;
        
        A.x[j] += A.v[j] * h;
        B.x[j] += B.v[j] * h;
        C.x[j] += C.v[j] * h;
    } 
    
    A.computeKineticEnergy();
    B.computeKineticEnergy();
    C.computeKineticEnergy();
    cm = computeCM(A,B,C);
    file_energy_euler<<A.energy + B.energy + C.energy<<";"<< computeGravitationalEnergy(A, B, C)<<";"<< computeElasticEnergy(B, C) <<std::endl;
    file_angmom_euler<<AngularMomentum(cm, A)[0]+ AngularMomentum(cm, B)[0]+AngularMomentum(cm, C)[0]<<";"<< AngularMomentum(cm, A)[1]+ AngularMomentum(cm, B)[1]+AngularMomentum(cm, C)[1]<<";"<<AngularMomentum(cm, A)[2]+ AngularMomentum(cm, B)[2]+AngularMomentum(cm, C)[2]<<std::endl;
    
    // file_angmom<< <<";"<< <<";"<< <<std::endl;
    // file_energy<<A.energy + B.energy + C.energy <<std::endl;
    // file_energy<< computePotentialEnergy(A, B, C)<<std::endl;
}

// ==========================================================
//                          RUNGE KUTTA 4
// ==========================================================
A.x = initial_conditions[0];
A.v = initial_conditions[1];
B.x = initial_conditions[2];
B.v = initial_conditions[3];
C.x = initial_conditions[4];
C.v = initial_conditions[5];

h *= 4;

for(int i=0; i<N_STEPS/4-1; i++){
    vA = A.v;
    xA = A.x;
    vB = B.v;
    xB = B.x;
    vC = C.v;
    xC = C.x;
    for(int j=0; j<DIM-1; j++){
        // body A
        m1[0][j] = h * vA[j];
        k1[0][j] = h * acceleration(xA[j], vA[j], t, C, B, A, j);
        // body B
        m1[1][j] = h * vB[j];
        k1[1][j] = h * springB(xB[j], vB[j], t, C, A, B, j); 
        // body C
        m1[2][j] = h * vC[j];
        k1[2][j] = h * springC(xC[j], vC[j], t, A, B, C, j); 
    }
    for(int j=0; j<DIM-1;j++){
        A.v[j] = vA[j] + 0.5 * k1[0][j];
        B.v[j] = vB[j] + 0.5 * k1[1][j];
        C.v[j] = vC[j] + 0.5 * k1[2][j];          
        A.x[j] = xA[j] + 0.5 * m1[0][j];
        B.x[j] = xB[j] + 0.5 * m1[1][j];
        C.x[j] = xC[j] + 0.5 * m1[2][j];
    }
    for(int j=0; j<DIM-1; j++){
        //Body A
        m2[0][j] = h * A.v[j];
        k2[0][j] = h * acceleration(A.x[j], A.v[j], t+0.5*h, C, B, A, j);
        //Body B
        m2[1][j] = h * B.v[j];
        k2[1][j] = h * springB(B.x[j], B.v[j], t+0.5*h, C, A, B, j);
        // Body C
        m2[2][j] = h * C.v[j];
        k2[2][j] = h * springC(C.x[j], C.v[j], t+0.5*h, A, B, C, j);
    }
        for(int j=0; j<DIM-1;j++){
        A.v[j] = vA[j] + 0.5 * k2[0][j];
        B.v[j] = vB[j] + 0.5 * k2[1][j];
        C.v[j] = vC[j] + 0.5 * k2[2][j];          
        A.x[j] = xA[j] + 0.5 * m2[0][j];
        B.x[j] = xB[j] + 0.5 * m2[1][j];
        C.x[j] = xC[j] + 0.5 * m2[2][j];
    }
    for(int j=0; j<DIM-1; j++){
        //Body A
        m3[0][j] = h * A.v[j];
        k3[0][j] = h * acceleration(A.x[j], A.v[j], t+0.5*h, C, B, A, j);
        
        //Body B
        m3[1][j] = h * B.v[j];
        k3[1][j] = h * springB(B.x[j], B.v[j], t+0.5*h, C, A, B, j);
        // Body C
        m3[2][j] = h * C.v[j];
        k3[2][j] = h * springC(C.x[j], C.v[j], t+0.5*h, A, B, C, j);
    }
        for(int j=0; j<DIM-1;j++){
        A.v[j] = vA[j] + k3[0][j];
        B.v[j] = vB[j] + k3[1][j];
        C.v[j] = vC[j] + k3[2][j];          
        A.x[j] = xA[j] + m3[0][j];
        B.x[j] = xB[j] + m3[1][j];
        C.x[j] = xC[j] + m3[2][j];
    }
    for(int j=0; j<DIM-1; j++){
        //Body A
        m4[0][j] = h * A.v[j];
        k4[0][j] = h * acceleration(A.x[j], A.v[j], t + h, C, B, A, j);
        //Body B
        m4[1][j] = h * B.v[j];
        k4[1][j] = h * springB(B.x[j], B.v[j], t + h, C, A, B, j);
        // Body C
        m4[2][j] = h * C.v[j];
        k4[2][j] = h * springC(C.x[j], C.v[j], t + h, A, B, C, j);
    }

    for(int j=0; j<DIM-1;j++){
        A.v[j] = vA[j] + (k1[0][j] + 2*k2[0][j] + 2*k3[0][j] + k4[0][j])/6;
        B.v[j] = vB[j] + (k1[1][j] + 2*k2[1][j] + 2*k3[1][j] + k4[1][j])/6;
        C.v[j] = vC[j] + (k1[2][j] + 2*k2[2][j] + 2*k3[2][j] + k4[2][j])/6;          
        A.x[j] = xA[j] + (m1[0][j] + 2*m2[0][j] + 2*m3[0][j] + m4[0][j])/6;
        B.x[j] = xB[j] + (m1[1][j] + 2*m2[1][j] + 2*m3[1][j] + m4[1][j])/6;
        C.x[j] = xC[j] + (m1[2][j] + 2*m2[2][j] + 2*m3[2][j] + m4[2][j])/6;
    }

    output_file_rk4 << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< ";";
    output_file_rk4 << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< ";";
    output_file_rk4 << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;
    A.computeKineticEnergy();
    B.computeKineticEnergy();
    C.computeKineticEnergy();
    file_energy_rk4<<A.energy + B.energy + C.energy <<";"<< computeGravitationalEnergy(A, B, C)<<";"<<computeElasticEnergy(B, C)<<std::endl;
    cm=computeCM(A, B, C);
    file_angmom_rk4<<AngularMomentum(cm, A)[0]+ AngularMomentum(cm, B)[0]+AngularMomentum(cm, C)[0]<<";"<< AngularMomentum(cm, A)[1]+ AngularMomentum(cm, B)[1]+AngularMomentum(cm, C)[1]<<";"<<AngularMomentum(cm, A)[2]+ AngularMomentum(cm, B)[2]+AngularMomentum(cm, C)[2]<<std::endl;
}

// ==========================================================
//                          LEAPFROG
// ==========================================================
A.x = initial_conditions[0];
A.v = initial_conditions[1];
B.x = initial_conditions[2];
B.v = initial_conditions[3];
C.x = initial_conditions[4];
C.v = initial_conditions[5];
h/=4;
for (int i=0; i<N_STEPS-1; i++){
    output_file_leapfrog << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< ";";
    output_file_leapfrog << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< ";";
    output_file_leapfrog << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;

    for(int j=0; j<DIM-1; j++){
        A.x[j] += A.v[j] * h / 2;
        B.x[j] += B.v[j] * h / 2;
        C.x[j] += C.v[j] * h / 2;
    }
    
    for(int j=0; j<DIM-1; j++){
        A.a[j] = acceleration(0, 0, 0, B, C, A, j);
        B.a[j] = springB(0, 0, 0, C, A, B, j);
        C.a[j] = springC(0, 0, 0, A, B, C, j);

        A.v[j] += A.a[j] * h;
        B.v[j] += B.a[j] * h;
        C.v[j] += C.a[j] * h;

    } 
    for(int j=0; j<DIM-1; j++){

        A.x[j] += A.v[j] * h / 2;
        B.x[j] += B.v[j] * h / 2;
        C.x[j] += C.v[j] * h / 2;
    }

    A.computeKineticEnergy();
    B.computeKineticEnergy();
    C.computeKineticEnergy();
    cm = computeCM(A,B,C);
    file_energy_leapfrog<<A.energy + B.energy + C.energy<<";"<< computeGravitationalEnergy(A, B, C)<<";"<<computeElasticEnergy(B, C)<<std::endl;
    file_angmom_leapfrog<<AngularMomentum(cm, A)[0]+ AngularMomentum(cm, B)[0]+AngularMomentum(cm, C)[0]<<";"<< AngularMomentum(cm, A)[1]+ AngularMomentum(cm, B)[1]+AngularMomentum(cm, C)[1]<<";"<<AngularMomentum(cm, A)[2]+ AngularMomentum(cm, B)[2]+AngularMomentum(cm, C)[2]<<std::endl;
}
    
//----------------------------------------------------------------  

    output_file_euler.close();
    file_energy_euler.close();
    file_angmom_euler.close();
    output_file_rk4.close();
    file_energy_rk4.close();
    file_angmom_rk4.close();
    output_file_leapfrog.close();
    file_energy_leapfrog.close();
    file_angmom_leapfrog.close();
    file_distance_euler.close();
    file_distance_leapfrog.close();
    file_distance_rk4.close();
    #ifdef _WIN32
        std::string command1 ="python3.9 plotting.py rk4";
        std::string command2 ="python3.9 plotting.py euler";
        std::string command3 ="python3.9 plotting.py leapfrog";
    #elif __APPLE__
        std::string command1 ="python3.11 plotting.py rk4";
        std::string command2 ="python3.11 plotting.py euler";
        std::string command3 ="python3.11 plotting.py leapfrog";
    #endif
    FILE* pipe1 = popen(command1.c_str(), "w");
    pclose(pipe1);
    FILE* pipe2 = popen(command2.c_str(), "w");
    pclose(pipe2);
    FILE* pipe3 = popen(command3.c_str(), "w");
    pclose(pipe3);
   return 0;
}