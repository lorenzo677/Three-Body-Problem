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
  s_mapStringValues["verlet"] = evStringValue3;
}

static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 300000;

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

double computePotentialEnergy(Planet planet1, Planet planet2, Planet planet3){
    return -1 * G * ( planet3.m * planet1.m / distance(planet3.x, planet1.x) +  planet3.m * planet2.m / distance(planet3.x, planet2.x) + planet1.m * planet2.m / distance(planet1.x, planet2.x));
        
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

double acceleration(Planet A, Planet B, Planet C, int axe){
    // Compute the acceleration of the body C, specifying the axis.
    double mass_A = A.m;
    double mass_B = B.m;
    double posx_A = A.x[0]; 
    double posx_B = B.x[0]; 
    double posx_C = C.x[0]; 
    double posy_A = A.x[1]; 
    double posy_B = B.x[1]; 
    double posy_C = C.x[1]; 
    double posz_A = A.x[2]; 
    double posz_B = B.x[2]; 
    double posz_C = C.x[2];
    if (axe == 0){
        return (-1 * G * (mass_A * (posx_C-posx_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posx_C-posx_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }else if (axe == 1){
        return (-1 * G * (mass_A * (posy_C-posy_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posy_C-posy_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }else if (axe == 2){
        return (-1 * G * (mass_A * (posz_C-posz_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posz_C-posz_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }
}

double F(double x, double v, double t, Planet A, Planet B, Planet C, int j ){
    // Function to integrate via Runge-Kutta.
    double mass_A = A.m;
    double mass_B = B.m;
    double posx_A = A.x[0]; 
    double posx_B = B.x[0]; 
    double posx_C = C.x[0]; 
    double posy_A = A.x[1]; 
    double posy_B = B.x[1]; 
    double posy_C = C.x[1]; 
    double posz_A = A.x[2]; 
    double posz_B = B.x[2]; 
    double posz_C = C.x[2];

    if (j == 0){
        return (-1 * G * (mass_A * (x-posx_A) / pow(sqrt(pow(x-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (x-posx_B) / pow(sqrt(pow(x-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }else if (j == 1) {
        return (-1 * G * (mass_A * (x-posy_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(x-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (x-posy_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(x-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }else if (j == 2) {
        return (-1 * G * (mass_A * (x-posz_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(x-posz_A,2)), 3) + mass_B * (x-posz_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(x-posz_B,2)), 3)));
    }
}

int main(int argc, char** argv){
    
    double h = 0.0005;

    Planet A(100, -10, 10, -11, -3, 0, 0);   // corpi allineati sull'asse delle x
    Planet B(100, 0, 0, 0, 3, 0, 0);
    Planet C(0, 10, 14, 12, 3, 0, 0);
    // Planet A(2e30, 0, 0, 0, 0, 0, 0);   // dati veri
    // Planet B(5.9e24, 150e9, 0, 0, 0, 3e4, 0);
    // Planet C(7.34e22, 60e9, 3.8e8, 0, -1022, 3e4, 0);
   
    // Planet A(10, -10, 10, -20, 0, 0, 0);   // corpi allineati sull'asse delle x
    // Planet B(10, 0, 0, 0, 10, 0, 0);
    // Planet C(10, 10, -10, -20, 0, 0, 0);
    
    double x_A[DIM][3];
    double x_B[DIM][3];
    double x_C[DIM][3];
    double vx_A;
    double vy_A;
    double vz_A;
    double vx_B;
    double vy_B;
    double vz_B;
    double vx_C;
    double vy_C;
    double vz_C;

    double mass_A = A.m;
    double mass_B = B.m;
    double mass_C = C.m;

    x_A[0][0] = A.x[0];
    x_B[0][0] = B.x[0];
    x_C[0][0] = C.x[0];
    
    vx_A = A.v[0];
    vx_B = B.v[0];
    vx_C = C.v[0];

    x_A[1][0] = A.x[1];
    x_B[1][0] = B.x[1];
    x_C[1][0] = C.x[1];

    vy_A = A.v[1];
    vy_B = B.v[1];
    vy_C = C.v[1];

    x_A[2][0] = A.x[2];
    x_B[2][0] = B.x[2];
    x_C[2][0] = C.x[2];

    vz_A = A.v[2];
    vz_B = B.v[2];
    vz_C = C.v[2];

    Initialize();
   
    std::ofstream file_energy("Total_energy_" + std::string(argv[1]) + ".csv");
    std::ofstream output_file_A("positions_A_" + std::string(argv[1]) + ".csv");
    std::ofstream output_file_B("positions_B_" + std::string(argv[1]) + ".csv");
    std::ofstream output_file_C("positions_C_" + std::string(argv[1]) + ".csv");
    std::ofstream file_angmom("Total_angular_momentum_" + std::string(argv[1]) + ".csv");

    output_file_A<<"x;y;z"<<std::endl;
    output_file_B<<"x;y;z"<<std::endl;
    output_file_C<<"x;y;z"<<std::endl;
    file_energy<<"k;p"<<std::endl;
    file_angmom<<"Lx;Ly;Lz"<<std::endl;

    double m1;
    double k1;
    double m2;
    double k2;
    double m3;
    double k3;
    double m4;
    double k4;
    
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

                        A.a[j] = acceleration(B, C, A, j);
                        B.a[j] = acceleration(A, C, B, j);
                        C.a[j] = acceleration(B, A, C, j);

                        A.v[j] += A.a[j] * h;
                        B.v[j] += B.a[j] * h;
                        C.v[j] += C.a[j] * h;

                        x_A[j][0] = x_A[j][0] + A.v[j] * h;
                        x_B[j][0] = x_B[j][0] + B.v[j] * h;
                        x_C[j][0] = x_C[j][0] + C.v[j] * h;
                        
                    } 
                    for (int j=0;j<DIM-1;j++){
                        A.x[j] = x_A[j][0];
                        B.x[j] = x_B[j][0];
                        C.x[j] = x_C[j][0];
                    }
                    
                    A.computeKineticEnergy();
                    B.computeKineticEnergy();
                    C.computeKineticEnergy();
                    cm = computeCM(A,B,C);
                    file_energy<<A.energy + B.energy + C.energy<<";"<< computePotentialEnergy(A, B, C)<<std::endl;
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

                for(int i=0; i<N_STEPS-1; i++){
                    for(int j=0; j<DIM-1; j++){
                        // body A
                        m1 = h*vA[j];
                        k1 = h*F(xA[j], vA[j], t, C, B, A, j);  

                        m2 = h*(vA[j] + 0.5*k1);
                        k2 = h*F(xA[j]+0.5*m1, vA[j]+0.5*k1, t+0.5*h, C, B, A, j);

                        m3 = h*(vA[j] + 0.5*k2);
                        k3 = h*F(xA[j]+0.5*m2, vA[j]+0.5*k2, t+0.5*h, C, B, A, j);

                        m4 = h*(vA[j] + k3);
                        k4 = h*F(xA[j]+m3, vA[j]+k3, t+h, C, B, A, j);

                        xA[j] += (m1 + 2*m2 + 2*m3 + m4)/6;
                        vA[j] += (k1 + 2*k2 + 2*k3 + k4)/6;
                    }

                    for(int j=0; j<DIM-1; j++){
                        // body B
                        m1 = h*vB[j];
                        k1 = h*F(xB[j], vB[j], t, A, C, B, j);  //(x, v, t)

                        m2 = h*(vB[j] + 0.5*k1);
                        k2 = h*F(xB[j]+0.5*m1, vB[j]+0.5*k1, t+0.5*h, A, C, B, j);

                        m3 = h*(vB[j] + 0.5*k2);
                        k3 = h*F(xB[j]+0.5*m2, vB[j]+0.5*k2, t+0.5*h, A, C, B, j);

                        m4 = h*(vB[j] + k3);
                        k4 = h*F(xB[j]+m3, vB[j]+k3, t+h, A, C, B, j);

                        xB[j] += (m1 + 2*m2 + 2*m3 + m4)/6;
                        vB[j] += (k1 + 2*k2 + 2*k3 + k4)/6;
                    }

                    for(int j=0; j<DIM-1; j++){   
                        // body C
                        m1 = h*vC[j];
                        k1 = h*F(xC[j], vC[j], t, A, B, C, j);  //(x, v, t)

                        m2 = h*(vC[j] + 0.5*k1);
                        k2 = h*F(xC[j]+0.5*m1, vC[j]+0.5*k1, t+0.5*h, A, B, C, j);

                        m3 = h*(vC[j] + 0.5*k2);
                        k3 = h*F(xC[j]+0.5*m2, vC[j]+0.5*k2, t+0.5*h, A, B, C, j);

                        m4 = h*(vC[j] + k3);
                        k4 = h*F(xC[j]+m3, vC[j]+k3, t+h, A, B, C, j);

                        xC[j] += (m1 + 2*m2 + 2*m3 + m4)/6;
                        vC[j] += (k1 + 2*k2 + 2*k3 + k4)/6;
                        
                    }

                    for(int j=0; j<DIM-1;j++){
                        A.v[j] = vA[j];
                        B.v[j] = vB[j];
                        C.v[j] = vC[j];          
                        A.x[j] = xA[j];
                        B.x[j] = xB[j];
                        C.x[j] = xC[j];
                    }
                    output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
                    output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
                    output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;
                    A.computeKineticEnergy();
                    B.computeKineticEnergy();
                    C.computeKineticEnergy();
                    file_energy<<A.energy + B.energy + C.energy <<";"<< computePotentialEnergy(A, B, C)<<std::endl;
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
                //                          VERLET
                // ==========================================================
                
                for(int j=0; j<DIM-1; j++){
                    A.a[j] = acceleration(B, C, A, j);
                    B.a[j] = acceleration(A, C, B, j);
                    C.a[j] = acceleration(B, A, C, j);
                    x_A[j][1] = x_A[j][0] + A.v[j] * h + 0.5 * A.a[j] * h * h;
                    x_B[j][1] = x_B[j][0] + B.v[j] * h + 0.5 * B.a[j] * h * h;
                    x_C[j][1] = x_C[j][0] + C.v[j] * h + 0.5 * C.a[j] * h * h;
            
                }
                // print initial condition
                output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
                output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
                output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;

                A.x[0] = x_A[0][1];
                A.x[1] = x_A[1][1];
                A.x[2] = x_A[2][1];
                B.x[0] = x_B[0][1];
                B.x[1] = x_B[1][1];
                B.x[2] = x_B[2][1];
                C.x[0] = x_C[0][1];
                C.x[1] = x_C[1][1];
                C.x[2] = x_C[2][1];
                for (int i=1; i<N_STEPS-1; i++){
                    output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
                    output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
                    output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;
                    for(int j=0; j<DIM-1; j++){
                        k1=A.a[j];
                        k2=B.a[j];
                        k3=C.a[j];

                        A.a[j] = acceleration(B, C, A, j);
                        B.a[j] = acceleration(A, C, B, j);
                        C.a[j] = acceleration(B, A, C, j);
                        
                        x_A[j][2]= x_A[j][1]; //x_A variabile temporanea
                        x_B[j][2]= x_B[j][1];
                        x_C[j][2]= x_C[j][1];

                        x_A[j][1] = 2 * x_A[j][1] - x_A[j][0] + A.a[j] * h * h;
                        x_B[j][1] = 2 * x_B[j][1] - x_B[j][0] + B.a[j] * h * h;
                        x_C[j][1] = 2 * x_C[j][1] - x_C[j][0] + C.a[j] * h * h;                
                        
                        x_A[j][0]= x_A[j][2];
                        x_B[j][0]= x_B[j][2];
                        x_C[j][0]= x_C[j][2];

                        // A.v[j] = (x_A[j][1] - x_A[j][0]) / h;
                        // B.v[j] = (x_B[j][1] - x_B[j][0]) / h;
                        // C.v[j] = (x_B[j][1] - x_B[j][0]) / h;
                        
                        A.v[j] += 0.5 * (A.a[j]+k1) * h;
                        B.v[j] += 0.5 * (B.a[j]+k2) * h;
                        C.v[j] += 0.5 * (C.a[j]+k3) * h;

                    }      
                    for (int j = 0; j < DIM-1; j++){
                            
                        A.x[j] = x_A[j][1];
                        B.x[j] = x_B[j][1];
                        C.x[j] = x_C[j][1];
                
                        A.computeKineticEnergy();
                        B.computeKineticEnergy();
                        C.computeKineticEnergy();
                    }
                    cm=computeCM(A, B, C);
                    file_energy<<A.energy + B.energy + C.energy <<";"<< computePotentialEnergy(A, B, C)<<std::endl;
                    file_angmom<<AngularMomentum(cm, A)[0]+ AngularMomentum(cm, B)[0]+AngularMomentum(cm, C)[0]<<";"<< AngularMomentum(cm, A)[1]+ AngularMomentum(cm, B)[1]+AngularMomentum(cm, C)[1]<<";"<<AngularMomentum(cm, A)[2]+ AngularMomentum(cm, B)[2]+AngularMomentum(cm, C)[2]<<std::endl;
                    // file_energy<<A.energy + B.energy + C.energy +computePotentialEnergy(A, B, C)<<std::endl;
                    // file_energy<<A.energy + B.energy + C.energy <<std::endl;
                    // file_energy<< computePotentialEnergy(A, B, C)<<std::endl;
                }
                break;
            default:
                std::cout<<"Inserire un argomento tra: euler, rk o verlet"<<std::endl;
                return 0;
        }
    }else{
        std::cout<<"Inserire argomento: euler, rk4 oppure verlet"<<std::endl;
        return 0;
    }
    
//----------------------------------------------------------------  

    output_file_A.close();
    output_file_B.close(); 
    output_file_C.close();
    file_energy.close();

    #ifdef _WIN32
        std::string command ="python3 plotting.py " + std::string(argv[1]);
    #elif __APPLE__
        std::string command ="python3.11 plotting.py " + std::string(argv[1]);
    #endif
    FILE* pipe = popen(command.c_str(), "w");
    pclose(pipe);
   return 0;
}