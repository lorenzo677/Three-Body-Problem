#include <iostream>
#include <cmath>
#include <array>
#include "integrators.h"

#define DIM 3
#define G 6.67408e-11
#define N_BODIES 3
#define N_STEPS 1500

std::array<float, 3> compute_vector_cm(float mass_1, float mass_2, float mass_3, float *vector_1, float *vector_2, float *vector_3){   
    std::array<float, 3> vector_cm;
    for (int i = 0; i < 3; i++) {   
        vector_cm[i] = ( mass_1 * vector_1[i] + mass_2 * vector_2[i] + mass_3 * vector_3[i] ) / (mass_1 + mass_2 + mass_3);
    }
    return vector_cm; 
}

float distance(float pos_1[DIM], float pos_2[DIM]){
    return sqrt(pow(pos_1[0]-pos_2[0], 2) + pow(pos_1[1]-pos_2[1], 2) + pow(pos_1[2]-pos_2[2], 2));
}

// Funzioni per le accelerazioni //

float acceleration(float mass_1, float mass_2, float pos_1, float pos_2, float pos_3){ 
    //compute the acceleration along one axis of the body 3
    return G * (mass_1 / pow(pos_3-pos_1, 2) + mass_2 / pow(pos_3-pos_2, 2));
}



int main(){
    // Definition of constants (be careful of Unit Measurements)
    float mass_1 = 2e28, mass_2 = 3e30, mass_3 = 1e31;                       // masses
    float x0_1[DIM], x0_2[DIM], x0_3[DIM];                                  // initial positions
    float v0_1[DIM], v0_2[DIM], v0_3[DIM];                                  // initial velocity
    float h = 1;

    float time[N_STEPS], x_1[DIM][N_STEPS], x_2[DIM][N_STEPS], x_3[DIM][N_STEPS], v_1[DIM][N_STEPS], v_2[DIM][N_STEPS], v_3[DIM][N_STEPS], a_1[DIM][N_STEPS], a_2[DIM][N_STEPS], a_3[DIM][N_STEPS];
          

    for (int i=0;i<N_STEPS; i++){time[i]= i * h;}

    x_1[0][0] = x0_1[0];
    x_2[0][0] = x0_2[0];
    x_3[0][0] = x0_3[0];
    v_1[0][0] = v0_1[0];
    v_2[0][0] = v0_2[0];
    v_3[0][0] = v0_3[0];
    x_1[1][0] = x0_1[1];
    x_3[1][0] = x0_2[1];
    v_1[1][0] = x0_3[1];
    v_2[1][0] = v0_1[1];
    x_2[1][0] = v0_2[1];
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
for(int j=0; j<DIM; j++){
    for (int i=0; i<N_STEPS; i++) {
    
    a_1[j][i] = acceleration(mass_2, mass_3, x_2[j][i], x_3[j][i], x_1[j][i]);
    a_2[j][i] = acceleration(mass_1, mass_3, x_1[j][i], x_3[j][i], x_2[j][i]);
    a_3[j][i] = acceleration(mass_1, mass_2, x_1[j][i], x_2[j][i], x_3[j][i]);

    v_1[j][i + 1] = v_1[j][i] + a_1[j][i] * h;
	v_2[j][i + 1] = v_2[j][i] + a_2[j][i] * h;
	v_2[j][i + 1] = v_3[j][i] + a_3[j][i] * h;

    x_1[j][i + 1] = x_1[j][i] + v_1[j][i] * h;
    x_2[j][i + 1] = x_2[j][i] + v_2[j][i] * h;
    x_3[j][i + 1] = x_3[j][i] + v_3[j][i] * h;

}
}

    // diamo delle condizioni iniziali e stampiamo l'array x_1
    // TODO


    //qua sotto per vedere se funziona compute_vector_cm
    float a[3] = {3.0, 2.0, 3.4};
    float b[3] = {2.0, 2.0, 40.4};
    float c[3] = {2.0, 1.0, 4.4};
    std::array<float, 3> cm;
    cm = compute_vector_cm(4.0, 5.4, 6.2, a, b, c);

    std::cout<<cm[0]<<'\t'<<cm[1]<<'\t'<<cm[2]<<std::endl;

    return 0;
}