#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include "integrators.h"

static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 5000;

class Planet{
private:
    double m;
    double x;
    double y;
    double z;
    double v_x;
    double v_y;
    double v_z;
    double a_x;
    double a_y;
    double a_z;
    double energy;
public:
    Planet () {m=0; x=0; y=0; z=0, v_x=0; v_y=0; v_z=0;}
		void setPlanet(double mass, double x_position, double y_position, double z_position, double x_velocity, double y_velocity, double z_velocity){
			m = mass;
			x = x_position;
			y = y_position;
            z = z_position;
			v_x = x_velocity;
            v_y = y_velocity;
            v_z = z_velocity;
		}
		double getPositionX(void){
			return x;
		}
		double getPositionY(void){
			return y;
		}
        double getPositionZ(void){
			return z;
		}
        double getMass(void){
            return m;
        }
        double getVelocityX(void){
            return v_x;
        }
        double getVelocityY(void){
            return v_y;
        }
        double getVelocityZ(void){
            return v_z;
        }
        double getAccelX(void){
            return a_x;
        }
        double getAccelY(void){
            return a_y;
        }
        double getAccelZ(void){
            return a_z;
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

double acceleration_class(Planet A, Planet B, Planet C){
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
    return (-1 * G * (mass_A * (posx_C-posx_A) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posx_C-posx_B) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
}

double acceleration(double mass_1, double mass_2, double posx_1, double posx_2, double posx_3, double posy_1, double posy_2, double posy_3, double posz_1, double posz_2, double posz_3){ 
    //compute the acceleration along one axis of the body 3
    return (-1 * G * (mass_1 * (posx_3-posx_1) / pow(sqrt(pow(posx_3-posx_1,2)+pow(posy_3-posy_1,2)+pow(posz_3-posz_1,2)), 3) + mass_2 * (posx_3-posx_2) / pow(sqrt(pow(posx_3-posx_2,2)+pow(posy_3-posy_2,2)+pow(posz_3-posz_2,2)), 3)));

}

int main(){
    
    // BE CAREFUL: If starting position of the body on one axis is the same, acceleration will be to inf.
    double mass_1 = 10;
    double mass_2 = 10;
    double mass_3 = 10;    
    double h = 0.001;

    // std::array<double, DIM-1> x0_1 = {-10, 10, -11};
    // std::array<double, DIM-1> x0_2 = {0, 0, 0};
    // std::array<double, DIM-1> x0_3 = {10, 14, 12};                                  
    // std::array<double, DIM-1> v0_1 = {-3, 0, 0};
    // std::array<double, DIM-1> v0_2 = {0, 0, 0};
    // std::array<double, DIM-1> v0_3 = {0, 0, 0};
    std::array<double, DIM-1> x0_1 = {-20, 20, 0};
    std::array<double, DIM-1> x0_2 = {0, 0, 0};
    std::array<double, DIM-1> x0_3 = {20, -20, 0};                                  
    std::array<double, DIM-1> v0_1 = {10, 10, 0};
    std::array<double, DIM-1> v0_2 = {0, 0, 0};
    std::array<double, DIM-1> v0_3 = {-10, -10, 0};
    std::array<double, N_STEPS> time;
    
    double x_1[DIM][N_STEPS];
    double x_2[DIM][N_STEPS];
    double x_3[DIM][N_STEPS];
    double v_1[DIM][N_STEPS];
    double v_2[DIM][N_STEPS];
    double v_3[DIM][N_STEPS];
    double a_1[DIM][N_STEPS];
    double a_2[DIM][N_STEPS];
    double a_3[DIM][N_STEPS];

    for (int i=0;i<N_STEPS; i++){time[i]= i * h;}

    x_1[0][0] = x0_1[0];
    x_2[0][0] = x0_2[0];
    x_3[0][0] = x0_3[0];

    v_1[0][0] = v0_1[0];
    v_2[0][0] = v0_2[0];
    v_3[0][0] = v0_3[0];

    x_1[1][0] = x0_1[1];
    x_2[1][0] = x0_2[1];
    x_3[1][0] = x0_3[1];

    v_1[1][0] = v0_1[1];
    v_2[1][0] = v0_2[1];
    v_3[1][0] = v0_3[1];

    x_1[2][0] = x0_1[2];
    x_2[2][0] = x0_2[2];
    x_3[2][0] = x0_3[2];

    v_1[2][0] = v0_1[2];
    v_2[2][0] = v0_2[2];
    v_3[2][0] = v0_3[2];

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
    a_1[j][i] = acceleration(mass_2, mass_3, x_2[j][i], x_3[j][i], x_1[j][i], x_2[j+1][i], x_3[j+1][i], x_1[j+1][i], x_2[j+2][i], x_3[j+2][i], x_1[j+2][i]);
    a_2[j][i] = acceleration(mass_1, mass_3, x_1[j][i], x_3[j][i], x_2[j][i], x_1[j+1][i], x_3[j+1][i], x_2[j+1][i], x_1[j+2][i], x_3[j+2][i], x_2[j+2][i]);
    a_3[j][i] = acceleration(mass_1, mass_2, x_1[j][i], x_2[j][i], x_3[j][i], x_1[j+1][i], x_2[j+1][i], x_3[j+1][i], x_1[j+2][i], x_2[j+2][i], x_3[j+2][i]);
    
    v_1[j][i + 1] = v_1[j][i] + a_1[j][i] * h;
	v_2[j][i + 1] = v_2[j][i] + a_2[j][i] * h;
	v_3[j][i + 1] = v_3[j][i] + a_3[j][i] * h;
    
    x_1[j][i + 1] = x_1[j][i] + v_1[j][i] * h;
    x_2[j][i + 1] = x_2[j][i] + v_2[j][i] * h;
    x_3[j][i + 1] = x_3[j][i] + v_3[j][i] * h;
    
    if (j==0 and i==1){ // l'array inizia a sporcarsi per j=1 e i=1
        std::cout<<v_3[j][i-1]<<'+'<< a_3[j][i]<< '*'<< h;
        std::cout<<"I valori utilizzati per calcolare l'accelerazione valgono "<<x_2[j][i]<<", "<<x_3[j][i]<<", "<<x_1[j][i]<<std::endl;
        std::cout<<"L'accelerazione utilizzata per calcolare la velocità vale "<<a_3[j][i-1]<<std::endl;
        std::cout<<"La velocità utilizzata per calcolare la posizione vale "<<v_3[j][i]<<std::endl;
        std::cout<<"La posizione vale "<<x_3[j][i+1]<<std::endl;
        std::cout<<"\n";

    }
    }
}

    // Alla ricerca del bug perduto
    std::cout<<"Considero la poszione del corpo 3 nei primi step:\n"; //Il problema è a_3 sull'asse y.
    std::cout<<"x, "<<"y, "<<"z"<<std::endl;
    std::cout<<x_1[0][0]<<", "<<x_1[1][0]<<", "<<x_1[2][0]<<std::endl;
    std::cout<<x_1[0][1]<<", "<<x_1[1][1]<<", "<<x_1[2][1]<<std::endl;
    std::cout<<x_1[0][2]<<", "<<x_1[1][2]<<", "<<x_1[2][2]<<std::endl;


    

    std::ofstream output_file_A("positions_A.csv");
    std::ofstream output_file_B("positions_B.csv");
    std::ofstream output_file_C("positions_C.csv");
    output_file_A<<"x;y;z"<<std::endl;
    output_file_B<<"x;y;z"<<std::endl;
    output_file_C<<"x;y;z"<<std::endl;
    
    for(int i = 0; i<N_STEPS-1; i++){
        output_file_A << x_1[0][i] << ";" << x_1[1][i] << ";" << x_1[2][i]<< std::endl;
        output_file_B << x_2[0][i] << ";" << x_2[1][i] << ";" << x_2[2][i]<< std::endl;
        output_file_C << x_3[0][i] << ";" << x_3[1][i] << ";" << x_3[2][i]<< std::endl;
    }    
    // for(int i = 0; i<N_STEPS-1; i++){
    //     output_file_A << a_1[0][i] << ";" << a_1[1][i] << ";" << a_1[2][i]<< std::endl;
    //     output_file_B << a_2[0][i] << ";" << a_2[1][i] << ";" << a_2[2][i]<< std::endl;
    //     output_file_C << a_3[0][i] << ";" << a_3[1][i] << ";" << a_3[2][i]<< std::endl;
    // }    
    output_file_A.close();
    output_file_B.close(); 
    output_file_C.close();

 /*
   FILE* pipe = popen("conda activate ml\n python plotting.py", "w");
    pclose(pipe);
    return 0;
*/       
}