#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include "integrators.h"

static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 50000;

class Planet{

public:

    double m;
    std::array <double, 3> x;
    std::array <double, 3> v;
    std::array <double, 3> a;
    double energy;
    
    Planet () {m = 0; x = {0, 0, 0}; v = {0, 0, 0}; a = {0, 0, 0};};
		void setPlanet(double mass, double x_position, double y_position, double z_position, double x_velocity, double y_velocity, double z_velocity){
			m = mass;
			x[0] = x_position;
			x[1] = y_position;
            x[2] = z_position;
			v[0] = x_velocity;
            v[1] = y_velocity;
            v[2] = z_velocity;
		}
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
};

std::array<double, 3> compute_vector_cm(double mass_1, double mass_2, double mass_3, double *vector_1, double *vector_2, double *vector_3){   
    std::array<double, 3> vector_cm;
    for (int i = 0; i < 3; i++) {   
        vector_cm[i] = ( mass_1 * vector_1[i] + mass_2 * vector_2[i] + mass_3 * vector_3[i] ) / (mass_1 + mass_2 + mass_3);
    }
    return vector_cm; 
}

double distance(double pos_1[DIM], double pos_2[DIM]){
    return sqrt(pow(pos_1[0]-pos_2[0], 2) + pow(pos_1[1]-pos_2[1], 2) + pow(pos_1[2]-pos_2[2], 2));
}

// Funzione per l'accelerazione //

double acceleration_class(Planet A, Planet B, Planet C, int axe){
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

double acceleration(double mass_1, double mass_2, double posx_1, double posx_2, double posx_3, double posy_1, double posy_2, double posy_3, double posz_1, double posz_2, double posz_3){ 
    //compute the acceleration along one axis of the body 3
    return (-1 * G * (mass_1 * (posx_3-posx_1) / pow(sqrt(pow(posx_3-posx_1,2)+pow(posy_3-posy_1,2)+pow(posz_3-posz_1,2)), 3) + mass_2 * (posx_3-posx_2) / pow(sqrt(pow(posx_3-posx_2,2)+pow(posy_3-posy_2,2)+pow(posz_3-posz_2,2)), 3)));

}

int main(){
    
   
    double h = 0.01;

    Planet A;
    Planet B;
    Planet C;
    // BE CAREFUL: If starting position of the body on one axis is the same, acceleration will be to inf.
    
    A.setPlanet(10, -10, 10, -11, -3, 0, 0);   // corpi allineati sull'asse delle x
    B.setPlanet(10, 0, 0, 0, 0, 0, 0);
    C.setPlanet(10, 10, 14, 12, 0, 0, 0);

    A.setPlanet(10, -20, 20, 0, 10, 10, 0.1);   // corpi allineati sulla bisettrice 2 e 3 con velocita perpendicolare
    B.setPlanet(10, 0, 0, 0, 0, 0, 1);
    C.setPlanet(10, 20, -20, 0.3, -10, -10, 0);

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
    // for (int i=0;i<N_STEPS-1; i++){time[i]= i * h;}

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




//Function for the Euler method (contrasta con quella presente sull'header)

    for (int i=0; i<N_STEPS-1; i++){
        for(int j=0; j<DIM-1; j++){


            A.a[j] = acceleration_class(B, C, A, j);
            B.a[j] = acceleration_class(A, C, B, j);
            C.a[j] = acceleration_class(B, A, C, j);
            // a_A[j][i] = acceleration(mass_B, mass_C, x_B[j][i], x_C[j][i], x_A[j][i], x_B[j+1][i], x_C[j+1][i], x_A[j+1][i], x_B[j+2][i], x_C[j+2][i], x_A[j+2][i]);
            // a_B[j][i] = acceleration(mass_A, mass_C, x_A[j][i], x_C[j][i], x_B[j][i], x_A[j+1][i], x_C[j+1][i], x_B[j+1][i], x_A[j+2][i], x_C[j+2][i], x_B[j+2][i]);
            // a_C[j][i] = acceleration(mass_A, mass_B, x_A[j][i], x_B[j][i], x_C[j][i], x_A[j+1][i], x_B[j+1][i], x_C[j+1][i], x_A[j+2][i], x_B[j+2][i], x_C[j+2][i]);
            x_A[j][i + 1] = x_A[j][i] + A.v[j] * h;
            x_B[j][i + 1] = x_B[j][i] + B.v[j] * h;
            x_C[j][i + 1] = x_C[j][i] + C.v[j] * h;
            
            A.x[j] = x_A[j][i + 1];
            B.x[j] = x_B[j][i + 1];
            C.x[j] = x_C[j][i + 1];

            A.v[j] += A.a[j] * h;
            B.v[j] += B.a[j] * h;
            C.v[j] += C.a[j] * h;
            
            
        // if (j==0 and i==1){ // l'array inizia a sporcarsi per j=1 e i=1
        //     std::cout<<v_C[j][i-1]<<'+'<< a_C[j][i]<< '*'<< h;
        //     std::cout<<"I valori utilizzati per calcolare l'accelerazione valgono "<<x_B[j][i]<<", "<<x_3[j][i]<<", "<<x_1[j][i]<<std::endl;
        //     std::cout<<"L'accelerazione utilizzata per calcolare la velocità vale "<<a_3[j][i-1]<<std::endl;
        //     std::cout<<"La velocità utilizzata per calcolare la posizione vale "<<v_3[j][i]<<std::endl;
        //     std::cout<<"La posizione vale "<<x_3[j][i+1]<<std::endl<<std::endl;
        // }
        std::cout<<"i = " <<i<<"\tj="<<j<<std::endl;
    }
}

    // Alla ricerca del bug perduto
    std::cout<<"Considero la poszione del corpo 3 nei primi step:\n"; //Il problema è a_3 sull'asse y.
    std::cout<<"x, "<<"y, "<<"z"<<std::endl;
    std::cout<<x_A[0][0]<<", "<<x_A[1][0]<<", "<<x_A[2][0]<<std::endl;
    std::cout<<x_A[0][1]<<", "<<x_A[1][1]<<", "<<x_A[2][1]<<std::endl;
    std::cout<<x_A[0][2]<<", "<<x_A[1][2]<<", "<<x_A[2][2]<<std::endl;


    

    std::ofstream output_file_A("positions_A.csv");
    std::ofstream output_file_B("positions_B.csv");
    std::ofstream output_file_C("positions_C.csv");
    output_file_A<<"x;y;z"<<std::endl;
    output_file_B<<"x;y;z"<<std::endl;
    output_file_C<<"x;y;z"<<std::endl;
    
    for(int i = 0; i<N_STEPS-1; i++){
        output_file_A << x_A[0][i] << ";" << x_A[1][i] << ";" << x_A[2][i]<< std::endl;
        output_file_B << x_B[0][i] << ";" << x_B[1][i] << ";" << x_B[2][i]<< std::endl;
        output_file_C << x_C[0][i] << ";" << x_C[1][i] << ";" << x_C[2][i]<< std::endl;
    }    
    // for(int i = 0; i<N_STEPS-1; i++){
    //     output_file_A << a_1[0][i] << ";" << a_1[1][i] << ";" << a_1[2][i]<< std::endl;
    //     output_file_B << a_2[0][i] << ";" << a_2[1][i] << ";" << a_2[2][i]<< std::endl;
    //     output_file_C << a_3[0][i] << ";" << a_3[1][i] << ";" << a_3[2][i]<< std::endl;
    // }    
    output_file_A.close();
    output_file_B.close(); 
    output_file_C.close();

 
 //   FILE* pipe = popen("conda activate ml\n python plotting.py", "w");
 //   pclose(pipe);
   return 0;
      
}