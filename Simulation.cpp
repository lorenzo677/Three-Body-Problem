#include <iostream>
#include <cmath>
#include <array>
#include <fstream>
#include "integrators.h"

static constexpr int DIM = 4;
static constexpr double G = 10;
static constexpr int N_BODIES = 3;
static constexpr int N_STEPS = 1000;

std::array<float, 3> compute_vector_cm(float mass_1, float mass_2, float mass_3, float *vector_1, float *vector_2, float *vector_3){   
    std::array<float, 3> vector_cm;
    for (int i = 0; i < 3; i++) {   
        vector_cm[i] = ( mass_1 * vector_1[i] + mass_2 * vector_2[i] + mass_3 * vector_3[i] ) / (mass_1 + mass_2 + mass_3);
    }
    return vector_cm; 
}

double distance(double pos_1[DIM], double pos_2[DIM]){
    return sqrt(pow(pos_1[0]-pos_2[0], 2) + pow(pos_1[1]-pos_2[1], 2) + pow(pos_1[2]-pos_2[2], 2));
}

// Funzione per l'accelerazione //

double acceleration(float mass_1, float mass_2, float pos_1, float pos_2, float pos_3){ 
    //compute the acceleration along one axis of the body 3
    return -1 * G * (mass_1 * (pos_3-pos_1) / pow(std::abs(pos_3-pos_1), 3) + mass_2 * (pos_3-pos_2) / pow(std::abs(pos_3-pos_2), 3));
}

int main(){
    // BE CAREFUL: If starting position of the body on one axis is the same, acceleration will be to inf.
    double mass_1 = 10;
    double mass_2 = 20;
    double mass_3 = 30;                       
    double x0_1[DIM-1] = {-10, 10, -11};
    double x0_2[DIM-1] = {0, 0, 0};
    double x0_3[DIM-1] = {10, 14, 12};                                  
    double v0_1[DIM-1] = {-3, 0, 0};
    double v0_2[DIM-1] = {0, 0, 0};
    double v0_3[DIM-1] = {3, 0, 0};                                 
    double h = 0.0005;

    double time[N_STEPS];
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
    a_1[j][i] = acceleration(mass_2, mass_3, x_2[j][i], x_3[j][i], x_1[j][i]);
    a_2[j][i] = acceleration(mass_1, mass_3, x_1[j][i], x_3[j][i], x_2[j][i]);
    a_3[j][i] = acceleration(mass_1, mass_2, x_1[j][i], x_2[j][i], x_3[j][i]);
    
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
    std::cout<<"Considero la velocita del corpo 3 nei primi step:\n"; //Il problema è a_3 sull'asse y.
    std::cout<<"x, "<<"y, "<<"z"<<std::endl;
    std::cout<<v_3[0][0]<<", "<<v_3[1][0]<<", "<<v_3[2][0]<<std::endl;
    std::cout<<v_3[0][1]<<", "<<v_3[1][1]<<", "<<v_3[2][1]<<std::endl;
    std::cout<<v_3[0][2]<<", "<<v_3[1][2]<<", "<<v_3[2][2]<<std::endl;
    std::cout<<v_3[0][3]<<", "<<v_3[1][3]<<", "<<v_3[2][3]<<std::endl;
    std::cout<<v_3[0][4]<<", "<<v_3[1][4]<<", "<<v_3[2][4]<<std::endl;
    std::cout<<v_3[0][5]<<", "<<v_3[1][5]<<", "<<v_3[2][5]<<std::endl;
    std::cout<<v_3[0][25]<<", "<<v_3[1][25]<<", "<<v_3[2][25]<<std::endl;

    std::ofstream output_file_A("positions_A.csv");
    std::ofstream output_file_B("positions_B.csv");
    std::ofstream output_file_C("positions_C.csv");
    output_file_A<<"x;y;z"<<std::endl;
    output_file_B<<"x;y;z"<<std::endl;
    output_file_C<<"x;y;z"<<std::endl;
    
    // for(int i = 0; i<N_STEPS-1; i++){
    //     output_file_A << x_1[0][i] << ";" << x_1[1][i] << ";" << x_1[2][i]<< std::endl;
    //     output_file_B << x_2[0][i] << ";" << x_2[1][i] << ";" << x_2[2][i]<< std::endl;
    //     output_file_C << x_3[0][i] << ";" << x_3[1][i] << ";" << x_3[2][i]<< std::endl;
    // }    
    for(int i = 0; i<N_STEPS-1; i++){
        output_file_A << a_1[0][i] << ";" << a_1[1][i] << ";" << a_1[2][i]<< std::endl;
        output_file_B << a_2[0][i] << ";" << a_2[1][i] << ";" << a_2[2][i]<< std::endl;
        output_file_C << a_3[0][i] << ";" << a_3[1][i] << ";" << a_3[2][i]<< std::endl;
    }    
    output_file_A.close();
    output_file_B.close(); 
    output_file_C.close();

    return 0;
}