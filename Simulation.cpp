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
#include <map>

// Value-Defintions of the different String values
enum StringValue { evNotDefined,
                          evStringValue1,
                          evStringValue2,
                          evStringValue3
                          };

// Map to associate the strings with the enum values
static std::map<std::string, StringValue> s_mapStringValues;

void Initialize(){
  s_mapStringValues["euler"] = evStringValue1;
  s_mapStringValues["rk4"] = evStringValue2;
  s_mapStringValues["leapfrog"] = evStringValue3;
}

static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 70000;

// Spring
static constexpr int K_CONST = 100000;
static constexpr double L0 = 10; 
// float l0 = sqrt(pow(L0X,2)+pow(L0Y,2)+pow(L0Z,2)); // modulo lunghezza a riposo molla

// std::array <double, 3> L0 = {L0X, L0Y, L0Z};


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
    
    void computeOmega(){
        omega = sqrt(K_CONST / m);
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


     Planet A(500, 0, 0, 0, 0, 0, 0);   // corpi allineati sull'asse delle x
     Planet B(10, 0, 0 , 5, 5, 0, 5);
     Planet C(10, 0, 0, -5, -5, 0, -5);

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

    Initialize();
   
    std::ofstream file_energy("Total_energy_" + std::string(argv[1]) + ".csv");
    std::ofstream output_file_A("positions_A_" + std::string(argv[1]) + ".csv");
    std::ofstream output_file_B("positions_B_" + std::string(argv[1]) + ".csv");
    std::ofstream output_file_C("positions_C_" + std::string(argv[1]) + ".csv");
    std::ofstream file_angmom("Total_angular_momentum_" + std::string(argv[1]) + ".csv");

    output_file_A<<"x;y;z"<<std::endl;
    output_file_B<<"x;y;z"<<std::endl;
    output_file_C<<"x;y;z"<<std::endl;
    file_energy<<"k;g;e"<<std::endl;
    file_angmom<<"Lx;Ly;Lz"<<std::endl;

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

    if (argc>=2){
        switch (s_mapStringValues[argv[1]]){
            case evStringValue1: 
                // ==========================================================
                //                          EULER
                // ==========================================================

                for (int i=0; i<N_STEPS-1; i++){
                    output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
                    output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
                    output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;
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
                    file_energy<<A.energy + B.energy + C.energy<<";"<< computeGravitationalEnergy(A, B, C)<<";"<< computeElasticEnergy(B, C) <<std::endl;
                    file_angmom<<AngularMomentum(cm, A)[0]+ AngularMomentum(cm, B)[0]+AngularMomentum(cm, C)[0]<<";"<< AngularMomentum(cm, A)[1]+ AngularMomentum(cm, B)[1]+AngularMomentum(cm, C)[1]<<";"<<AngularMomentum(cm, A)[2]+ AngularMomentum(cm, B)[2]+AngularMomentum(cm, C)[2]<<std::endl;
                    // file_angmom<< <<";"<< <<";"<< <<std::endl;
                    // file_energy<<A.energy + B.energy + C.energy <<std::endl;
                    // file_energy<< computePotentialEnergy(A, B, C)<<std::endl;
                }
                
                break;
            case evStringValue2:{
                // ==========================================================
                //                          RUNGE KUTTA 4
                // ==========================================================
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

                    output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
                    output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
                    output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;
                    A.computeKineticEnergy();
                    B.computeKineticEnergy();
                    C.computeKineticEnergy();
                    file_energy<<A.energy + B.energy + C.energy <<";"<< computeGravitationalEnergy(A, B, C)<<";"<<computeElasticEnergy(B, C)<<std::endl;
                    cm=computeCM(A, B, C);
                    file_angmom<<AngularMomentum(cm, A)[0]+ AngularMomentum(cm, B)[0]+AngularMomentum(cm, C)[0]<<";"<< AngularMomentum(cm, A)[1]+ AngularMomentum(cm, B)[1]+AngularMomentum(cm, C)[1]<<";"<<AngularMomentum(cm, A)[2]+ AngularMomentum(cm, B)[2]+AngularMomentum(cm, C)[2]<<std::endl;
                    // file_energy<<A.energy + B.energy + C.energy + computePotentialEnergy(A, B, C)<<std::endl;
                    // file_energy<<A.energy + B.energy + C.energy <<std::endl;
                    // file_energy<< computePotentialEnergy(A, B, C)<<std::endl;
                }
                }
                break;
            case evStringValue3:
                // ==========================================================
                //                          LEAPFROG
                // ==========================================================
            
                for (int i=0; i<N_STEPS-1; i++){
                    output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
                    output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
                    output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;

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
                    file_energy<<A.energy + B.energy + C.energy<<";"<< computeGravitationalEnergy(A, B, C)<<";"<<computeElasticEnergy(B, C)<<std::endl;
                    file_angmom<<AngularMomentum(cm, A)[0]+ AngularMomentum(cm, B)[0]+AngularMomentum(cm, C)[0]<<";"<< AngularMomentum(cm, A)[1]+ AngularMomentum(cm, B)[1]+AngularMomentum(cm, C)[1]<<";"<<AngularMomentum(cm, A)[2]+ AngularMomentum(cm, B)[2]+AngularMomentum(cm, C)[2]<<std::endl;
                }

            break;
            default:
                std::cout<<"Insert an argument between: euler, rk4 or leapfrog"<<std::endl;
                return 0;
        }
    }else{
        std::cout<<"Insert an argument between: euler, rk4 or leapfrog"<<std::endl;
        return 0;
    }
    
//----------------------------------------------------------------  

    output_file_A.close();
    output_file_B.close(); 
    output_file_C.close();
    file_energy.close();

    #ifdef _WIN32
        std::string command ="python3.9 plotting.py " + std::string(argv[1]);
    #elif __APPLE__
        std::string command ="python3.11 plotting.py " + std::string(argv[1]);
    #endif
    FILE* pipe = popen(command.c_str(), "w");
    pclose(pipe);
   return 0;
}