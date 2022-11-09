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


static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 80000;

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
		double getPositionX(void){
			return x[0];
		}
		double getPositionY(void){
			return x[1];
		}
        double getPositionZ(void){
			return x[2];
		}
        double getMass(void){
            return m;
        }
        double getVelocityX(void){
            return v[0];
        }
        double getVelocityY(void){
            return v[1];
        }
        double getVelocityZ(void){
            return v[2];
        }
        double getAccelX(void){
            return a[0];
        }
        double getAccelY(void){
            return a[1];
        }
        double getAccelZ(void){
            return a[2];
        }
        double getEnergy(void){
            return energy;
        }
        void computeEnergy(Planet planet1, Planet planet2){
            
            energy = 0.5 * m * (pow(v[0],2)+pow(v[1],2)+pow(v[2],2)) - G * m * planet1.m / distance(x, planet1.x) - G * m * planet2.m / distance(x, planet2.x);
        }
};

double acceleration(Planet A, Planet B, Planet C, int axe){
    // Compute the acceleration of the body C, specifying the axis.
    double mass_A = A.getMass();
    double mass_B = B.getMass();
    double posx_A = A.getPositionX(); 
    double posx_B = B.getPositionX(); 
    double posx_C = C.getPositionX(); 
    double posy_A = A.getPositionY(); 
    double posy_B = B.getPositionY(); 
    double posy_C = C.getPositionY(); 
    double posz_A = A.getPositionZ(); 
    double posz_B = B.getPositionZ(); 
    double posz_C = C.getPositionZ();
    if (axe == 0)
        return (-1 * G * (mass_A * (posx_C-posx_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posx_C-posx_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    else if (axe == 1){
        return (-1 * G * (mass_A * (posy_C-posy_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posy_C-posy_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }else if (axe == 2){
        return (-1 * G * (mass_A * (posz_C-posz_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posz_C-posz_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }
}

double Fi(double x, double v, double t, Planet A, Planet B, Planet C, int j ){
    return x;
}

double F(double x, double v, double t, Planet A, Planet B, Planet C, int j ){
    // Function to integrate via Runge-Kutta.
    double mass_A = A.getMass();
    double mass_B = B.getMass();
    double posx_A = A.getPositionX(); 
    double posx_B = B.getPositionX(); 
    double posx_C = C.getPositionX(); 
    double posy_A = A.getPositionY(); 
    double posy_B = B.getPositionY(); 
    double posy_C = C.getPositionY(); 
    double posz_A = A.getPositionZ(); 
    double posz_B = B.getPositionZ(); 
    double posz_C = C.getPositionZ();

    if (j == 0){
        return (-1 * G * (mass_A * (x-posx_A) / pow(sqrt(pow(x-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (x-posx_B) / pow(sqrt(pow(x-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }else if (j == 1) {
        return (-1 * G * (mass_A * (x-posy_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(x-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (x-posy_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(x-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    }else if (j == 2) {
        return (-1 * G * (mass_A * (x-posz_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(x-posz_A,2)), 3) + mass_B * (x-posz_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(x-posz_B,2)), 3)));
    }
}

int main(){
    
    double h = 0.001;

    Planet A(10, -10, 10, -11, -3, 0, 0);   // corpi allineati sull'asse delle x
    Planet B(10, 0, 0, 0, 0, 0, 0);
    Planet C(10, 10, 14, 12, 0, 0, 0);
    // Planet A(10, -20, 20, 0, 10, 10, 0.1);   // corpi allineati sulla bisettrice 2 e 3 con velocita perpendicolare
    // Planet B(10, 0, 0, 0, 0, 0, 1);
    // Planet C(10, 20, -20, 0.3, -10, -10, 0);
    // BE CAREFUL: If starting position of the body on one axis is the same, acceleration will be to inf.
    // Planet A(10, -10, 10, -20, 0, 0, 0);   // corpi allineati sull'asse delle x
    // Planet B(10, 0, 0, 0, 10, 0, 0);
    // Planet C(10, 10, -10, -20, 0, 0, 0);


    std::array<double, N_STEPS> time;
    
    double x_A[DIM][N_STEPS];
    double x_B[DIM][N_STEPS];
    double x_C[DIM][N_STEPS];
    double vx_A;
    double vy_A;
    double vz_A;
    double vx_B;
    double vy_B;
    double vz_B;
    double vx_C;
    double vy_C;
    double vz_C;
    // double v_A[DIM][N_STEPS];
    // double v_B[DIM][N_STEPS];
    // double v_C[DIM][N_STEPS];
    // double a_A[DIM][N_STEPS];
    // double a_B[DIM][N_STEPS];
    // double a_C[DIM][N_STEPS];

    double mass_A = A.getMass();
    double mass_B = B.getMass();
    double mass_C = C.getMass();
    for (int i=0;i<N_STEPS-1; i++){time[i]= i * h;}

    x_A[0][0] = A.getPositionX();
    x_B[0][0] = B.getPositionX();
    x_C[0][0] = C.getPositionX();

    vx_A = A.getVelocityX();
    vx_B = B.getVelocityX();
    vx_C = C.getVelocityX();

    x_A[1][0] = A.getPositionY();
    x_B[1][0] = B.getPositionY();
    x_C[1][0] = C.getPositionY();

    vy_A = A.getVelocityY();
    vy_B = B.getVelocityY();
    vy_C = C.getVelocityY();

    x_A[2][0] = A.getPositionZ();
    x_B[2][0] = B.getPositionZ();
    x_C[2][0] = C.getPositionZ();

    vz_A = A.getVelocityZ();
    vz_B = B.getVelocityZ();
    vz_C = C.getVelocityZ();

    std::ofstream file_energy("Total_energy_euler.csv");

    
    //Function for the Euler method
    // std::cout<<"Energia totale pianeti Euler Method"<<std::endl;
    
    // for (int i=0; i<N_STEPS-1; i++){
    //      for(int j=0; j<DIM-1; j++){


    //         A.a[j] = acceleration(B, C, A, j);
    //         B.a[j] = acceleration(A, C, B, j);
    //         C.a[j] = acceleration(B, A, C, j);

    //         A.v[j] += A.a[j] * h;
    //         B.v[j] += B.a[j] * h;
    //         C.v[j] += C.a[j] * h;

    //         x_A[j][i + 1] = x_A[j][i] + A.v[j] * h;
    //         x_B[j][i + 1] = x_B[j][i] + B.v[j] * h;
    //         x_C[j][i + 1] = x_C[j][i] + C.v[j] * h;
            
    //         A.x[j] = x_A[j][i + 1];
    //         B.x[j] = x_B[j][i + 1];
    //         C.x[j] = x_C[j][i + 1];
            
    //         A.computeEnergy(B, C);
    //         B.computeEnergy(A, C);
    //         C.computeEnergy(B, A);
            
    //     } 
        
    //     file_energy<<A.energy+B.energy+C.energy<<std::endl;
    // }


    // RUNGE KUTTA 4
    // Integro accelerazione su asse x
    // x'' = -G(...)
    // viene trasformato in
    // x' = v
    // v' = -G(...)

    // double m1;
    // double k1;
    // double m2;
    // double k2;
    // double m3;
    // double k3;
    // double m4;
    // double k4;

    // std::array<double,3> vA = A.v; //condizione iniziale velocita
    // std::array<double,3> xA = A.x;
    // std::array<double,3> vB = B.v; //condizione iniziale velocita
    // std::array<double,3> xB = B.x;
    // std::array<double,3> vC = C.v; //condizione iniziale velocita
    // std::array<double,3> xC = C.x;
    // double t;
    
    // for(int i=0; i<N_STEPS-1; i++){
    //     for(int j=0; j<DIM-1; j++){
    //         // body A
    //         m1 = h*vA[j];
    //         k1 = h*F(xA[j], vA[j], t, C, B, A, j);  

    //         m2 = h*(vA[j] + 0.5*k1);
    //         k2 = h*F(xA[j]+0.5*m1, vA[j]+0.5*k1, t+0.5*h, C, B, A, j);

    //         m3 = h*(vA[j] + 0.5*k2);
    //         k3 = h*F(xA[j]+0.5*m2, vA[j]+0.5*k2, t+0.5*h, C, B, A, j);

    //         m4 = h*(vA[j] + k3);
    //         k4 = h*F(xA[j]+m3, vA[j]+k3, t+h, C, B, A, j);

    //         xA[j] += (m1 + 2*m2 + 2*m3 + m4)/6;
    //         vA[j] += (k1 + 2*k2 + 2*k3 + k4)/6;
    //     }

    //     for(int j=0; j<DIM-1; j++){
    //         // body B
    //         m1 = h*vB[j];
    //         k1 = h*F(xB[j], vB[j], t, A, C, B, j);  //(x, v, t)

    //         m2 = h*(vB[j] + 0.5*k1);
    //         k2 = h*F(xB[j]+0.5*m1, vB[j]+0.5*k1, t+0.5*h, A, C, B, j);

    //         m3 = h*(vB[j] + 0.5*k2);
    //         k3 = h*F(xB[j]+0.5*m2, vB[j]+0.5*k2, t+0.5*h, A, C, B, j);

    //         m4 = h*(vB[j] + k3);
    //         k4 = h*F(xB[j]+m3, vB[j]+k3, t+h, A, C, B, j);

    //         xB[j] += (m1 + 2*m2 + 2*m3 + m4)/6;
    //         vB[j] += (k1 + 2*k2 + 2*k3 + k4)/6;
    //     }

    //     for(int j=0; j<DIM-1; j++){   
    //         // body C
    //         m1 = h*vC[j];
    //         k1 = h*F(xC[j], vC[j], t, A, B, C, j);  //(x, v, t)

    //         m2 = h*(vC[j] + 0.5*k1);
    //         k2 = h*F(xC[j]+0.5*m1, vC[j]+0.5*k1, t+0.5*h, A, B, C, j);

    //         m3 = h*(vC[j] + 0.5*k2);
    //         k3 = h*F(xC[j]+0.5*m2, vC[j]+0.5*k2, t+0.5*h, A, B, C, j);

    //         m4 = h*(vC[j] + k3);
    //         k4 = h*F(xC[j]+m3, vC[j]+k3, t+h, A, B, C, j);

    //         xC[j] += (m1 + 2*m2 + 2*m3 + m4)/6;
    //         vC[j] += (k1 + 2*k2 + 2*k3 + k4)/6;
            
    //     }

         
    //     //std::cout<<A.x[0]<<std::endl;
    //     for(int j=0; j<DIM-1;j++){
    //         x_A[j][i+1] = xA[j];
    //         x_B[j][i+1] = xB[j];
    //         x_C[j][i+1] = xC[j];            
    //         A.x[j] = xA[j];
    //         B.x[j] = xB[j];
    //         C.x[j] = xC[j];
    //     }
    //         A.computeEnergy(B, C);
    //          B.computeEnergy(A, C);
    //          C.computeEnergy(B, A);
    //     file_energy<<A.energy+B.energy+C.energy<<std::endl;
    // }





// VERLET
x_A[0][1] = x_A[0][0] + A.v[0]*h;
x_A[1][1] = x_A[1][0] + A.v[1]*h;
x_A[2][1] = x_A[2][0] + A.v[2]*h;
x_B[0][1] = x_B[0][0] + B.v[0]*h;
x_B[1][1] = x_B[1][0] + B.v[1]*h;
x_B[2][1] = x_B[2][0] + B.v[2]*h;
x_C[0][1] = x_C[0][0] + C.v[0]*h;
x_C[1][1] = x_C[1][0] + C.v[1]*h;
x_C[2][1] = x_C[2][0] + C.v[2]*h;

A.x[0]= x_A[0][1];
A.x[1]= x_A[1][1];
A.x[2]= x_A[2][1];

B.x[0]=x_B[0][1];
B.x[1]=x_B[1][1];
B.x[2]=x_B[2][1];

C.x[0]=x_C[0][1];
C.x[1]=x_C[1][1];
C.x[2]=x_C[2][1];

std::cout<<x_A[0][0]<<std::endl;
std::cout<<x_A[0][1]<<std::endl;
 for (int i=1; i<N_STEPS-1; i++){
         for(int j=0; j<DIM-1; j++){

            A.a[j] = acceleration(B, C, A, j);
            B.a[j] = acceleration(A, C, B, j);
            C.a[j] = acceleration(B, A, C, j);

            file_energy<<A.a[j]<<std::endl;
            
            x_A[j][i + 1] = 2 * x_A[j][i] - x_A[j][i-1] + A.a[j] * h * h;
            x_B[j][i + 1] = 2 * x_B[j][i] - x_B[j][i-1] + B.a[j] * h * h;
            x_C[j][i + 1] = 2 * x_C[j][i] - x_C[j][i-1] + C.a[j] * h * h;
            if(j==0 and i==1){
            std::cout<<x_A[j][i + 1] <<"= 2 * " <<x_A[j][i] << "-"<< x_A[j][i-1]<< "+"<< A.a[j] <<"*"<< h <<"*"<< h<<std::endl;
            }
        }      
        for (int j = 0; j < DIM-1; j++){
           
            A.x[j] = x_A[j][i + 1];
            B.x[j] = x_B[j][i + 1];
            C.x[j] = x_C[j][i + 1];
            
            A.computeEnergy(B, C);
            B.computeEnergy(A, C);
            C.computeEnergy(B, A);
        }
        
        
        //file_energy<<A.energy+B.energy+C.energy<<std::endl;
    }






    file_energy.close();
// RK4 LUTZ

//     std::array<double, 3> m1A;
//     std::array<double, 3> k1A;
//     std::array<double, 3> m2A;
//     std::array<double, 3> k2A;
//     std::array<double, 3> m3A;
//     std::array<double, 3> k3A;
//     std::array<double, 3> m4A;
//     std::array<double, 3> k4A;

//     std::array<double, 3> m1B;
//     std::array<double, 3> k1B;
//     std::array<double, 3> m2B;
//     std::array<double, 3> k2B;
//     std::array<double, 3> m3B;
//     std::array<double, 3> k3B;
//     std::array<double, 3> m4B;
//     std::array<double, 3> k4B;
    
//     std::array<double, 3> m1C;
//     std::array<double, 3> k1C;
//     std::array<double, 3> m2C;
//     std::array<double, 3> k2C;
//     std::array<double, 3> m3C;
//     std::array<double, 3> k3C;
//     std::array<double, 3> m4C;
//     std::array<double, 3> k4C;

//     std::array<double,3> vA = A.v; //condizione iniziale velocita
//     std::array<double,3> xA = A.x;
//     std::array<double,3> vB = B.v; //condizione iniziale velocita
//     std::array<double,3> xB = B.x;
//     std::array<double,3> vC = C.v; //condizione iniziale velocita
//     std::array<double,3> xC = C.x;

// double t;

// for(int i =0;i<N_STEPS-1;i++){
//     for(int j=0; j<DIM-1; j++){
//             // body A
//             m1A[j] = h*vA[j];
//             k1A[j] = h*F(xA[j], vA[j], t, C, B, A, j);  

//             m2A[j] = h*(vA[j] + 0.5*k1A[j]);
//             k2A[j] = h*F(xA[j]+0.5*m1A[j], vA[j]+0.5*k1A[j], t+0.5*h, C, B, A, j);

//             m3A[j] = h*(vA[j] + 0.5*k2A[j]);
//             k3A[j] = h*F(xA[j]+0.5*m2A[j], vA[j]+0.5*k2A[j], t+0.5*h, C, B, A, j);

//             m4A[j] = h*(vA[j] + k3A[j]);
//             k4A[j] = h*F(xA[j]+m3A[j], vA[j]+k3A[j], t+h, C, B, A, j);

//             xA[j] += (m1A[j] + 2*m2A[j] + 2*m3A[j] + m4A[j])/6;
//             vA[j] += (k1A[j] + 2*k2A[j] + 2*k3A[j] + k4A[j])/6;
        
//             // body B
//             m1B[j] = h*vB[j];
//             k1B[j] = h*F(xB[j], vB[j], t, A, C, B, j);  //(x, v, t)

//             m2B[j] = h*(vB[j] + 0.5*k1B[j]);
//             k2B[j] = h*F(xB[j]+0.5*m1B[j], vB[j]+0.5*k1B[j], t+0.5*h, A, C, B, j);

//             m3B[j] = h*(vB[j] + 0.5*k2B[j]);
//             k3B[j] = h*F(xB[j]+0.5*m2B[j], vB[j]+0.5*k2B[j], t+0.5*h, A, C, B, j);

//             m4B[j] = h*(vB[j] + k3B[j]);
//             k4B[j] = h*F(xB[j]+m3B[j], vB[j]+k3B[j], t+h, A, C, B, j);

//             xB[j] += (m1B[j] + 2*m2B[j] + 2*m3B[j] + m4B[j])/6;
//             vB[j] += (k1B[j] + 2*k2B[j] + 2*k3B[j] + k4B[j])/6;

//             // body C
//             m1C[j] = h*vC[j];
//             k1C[j] = h*F(xC[j], vC[j], t, A, B, C, j);  //(x, v, t)

//             m2C[j] = h*(vC[j] + 0.5*k1C[j]);
//             k2C[j] = h*F(xC[j]+0.5*m1C[j], vC[j]+0.5*k1C[j], t+0.5*h, A, B, C, j);

//             m3C[j] = h*(vC[j] + 0.5*k2C[j]);
//             k3C[j] = h*F(xC[j]+0.5*m2C[j], vC[j]+0.5*k2C[j], t+0.5*h, A, B, C, j);

//             m4C[j] = h*(vC[j] + k3C[j]);
//             k4C[j] = h*F(xC[j]+m3C[j], vC[j]+k3C[j], t+h, A, B, C, j);

//             xC[j] += (m1C[j] + 2*m2C[j] + 2*m3C[j] + m4C[j])/6;
//             vC[j] += (k1C[j] + 2*k2C[j] + 2*k3C[j] + k4C[j])/6;
            
//         }

//         //std::cout<<A.x[0]<<std::endl;
//         for(int j=0; j<DIM-1;j++){
//             x_A[j][i+1] = xA[j];
//             x_B[j][i+1] = xB[j];
//             x_C[j][i+1] = xC[j];            
//             A.x[j] = xA[j];
//             B.x[j] = xB[j];
//             C.x[j] = xC[j];
//         }
//     }




//----------------------------------------------------------------

// Print data on .csv
    std::ofstream output_file_A("positions_A.csv");
    std::ofstream output_file_B("positions_B.csv");
    std::ofstream output_file_C("positions_C.csv");
    output_file_A<<"x;y;z;t"<<std::endl;
    output_file_B<<"x;y;z"<<std::endl;
    output_file_C<<"x;y;z"<<std::endl;
    
    for(int i = 0; i<N_STEPS-1; i++){
        output_file_A << x_A[0][i] << ";" << x_A[1][i] << ";" << x_A[2][i]<< std::endl;
        output_file_B << x_B[0][i] << ";" << x_B[1][i] << ";" << x_B[2][i]<< std::endl;
        output_file_C << x_C[0][i] << ";" << x_C[1][i] << ";" << x_C[2][i]<< std::endl;
    }    
    // for(int i = 0; i<N_STEPS-1; i++){
    //     output_file_A << x_A[0][i] << ";" << x_A[1][i]<< ";" << x_A[2][i]<<";"<< time[i]<<std::endl;
    // }    
    output_file_A.close();
    output_file_B.close(); 
    output_file_C.close();

 
    // FILE* pipe = popen("conda activate ml\n python plotting.py", "w");
    // pclose(pipe);
   return 0;
      
}