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

std::array<double, 3> compute_vector_cm(double mass_1, double mass_2, double mass_3, double *vector_1, double *vector_2, double *vector_3){   
    std::array<double, 3> vector_cm;
    for (int i = 0; i < 3; i++) {   
        vector_cm[i] = ( mass_1 * vector_1[i] + mass_2 * vector_2[i] + mass_3 * vector_3[i] ) / (mass_1 + mass_2 + mass_3);
    }
    return vector_cm; 
}

// double distance(Planet planet1, Planet planet2){
//     return sqrt(pow(planet1.x[0]-planet2.x[0],2)+pow(planet1.x[1]-planet2.x[1],2)+pow(planet1.x[2]-planet2.x[2],2));
// }


double acceleration(Planet A, Planet B, Planet C, int axe){
    //compute the acceleration along one axis of the body C
    double mass_A = A.getMass();
    double mass_B = A.getMass();
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

    /*
    for(int j=0; j<DIM; j++){
            for(int i=0; i<N_STEPS; i++){
            v_1[j][i+1] = RK4(time[i], v_1[j][i], h, acceleration, mass_2, mass_3, x_2[j][i], x_3[j][i]); // probabilmente sbagliatto
            }
        }
    */



    //Function for the Euler method
    // std::cout<<"Posizione x primo pianeta Euler Method"<<std::endl;
    
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
        
    //     std::cout<<A.x[0]<<std::endl;
    // }


    // ---------------------------------------------------------------------------------------------
    // RUNGE KUTTA 4

    // Integro accelerazione su asse x
    // x'' = -G(...)
    // viene trasformato in
    // x' = v
    // v' = -G(...)

    // double function_to_integrate(double x0, double y0, double m1, double m2, double p1, double p2){
    //     return (-1 * G * (m1 * (x0-posx_A) / pow(sqrt(pow(x0-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + m2 * (x0-posx_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
    // }
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
    
    // x' = v
    // v' = -G(...)
    
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
        //std::cout<<A.x[0]<<std::endl;
        for(int j=0; j<DIM-1;j++){
        x_A[j][i+1] = xA[j];
        x_B[j][i+1] = xB[j];
        x_C[j][i+1] = xC[j];            
        A.x[j] = xA[j];
        B.x[j] = xB[j];
        C.x[j] = xC[j];
        }
    }

    // double RK4(double x0, double y0, double h, double (*func)(double, double, double, double, double, double, double), double m1, double m2, double p1, double p2){
    //     // Finds value of y for a given x using step size h
    //     // and initial value y0 at x0.
    //     double k1, k2, k3, k4;

    //     k1 = h * func(x0, y0, m1, m2, p1, p2, x0);
    //     k2 = h * func(x0 + 0.5 * h, y0 + 0.5 * k1, m1, m2, p1, p2, x0);
    //     k3 = h * func(x0 + 0.5 * h, y0 + 0.5 * k2, m1, m2, p1, p2, x0);
    //     k4 = h * func(x0 + h, y0 + k3, m1, m2, p1, p2, x0);
    
    //     // Update next value of y
    //     y0 = y0 + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);;
    
    //     // Update next value of x
    //     x0 = x0 + h;

//     return y0;
// }

// int prova; // velocità del corpo A 
// prova = RK4(A.a[0], A.v[0], h, double (*func)(double, double, double, double, double, double, double), A.m, B.m, B.x[0], C.x[0]){




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

 
 //   FILE* pipe = popen("conda activate ml\n python plotting.py", "w");
 //   pclose(pipe);
   return 0;
      
}