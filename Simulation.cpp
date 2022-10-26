#include <iostream>
#include <cmath>
#include "integrators.h"

#define DIM 3
#define G 6.67408e-11
#define N_BODIES 3

void compute_vector_cm(float mass_1, float mass_2, float mass_3, float vector_1[3], float vector_2[3], float vector_3[3]){   
    float vector_cm[3];
    for (int i = 0; i < 4; i++) {   
        vector_cm[i] = ( mass_1 * vector_1[i] + mass_2 * vector_2[i] + mass_3 * vector_3[i] ) / (mass_1 + mass_2 + mass_3);
    }
}

float distance(float pos_1[DIM], float pos_2[DIM]){
    return sqrt(pow(pos_1[0]-pos_2[0], 2) + pow(pos_1[1]-pos_2[1], 2) + pow(pos_1[2]-pos_2[2], 2));
}

int main(){
    // Definition of constants (be careful of Unit Measurements)
    float mass_1, mass_2, mass_3;                       // masses
    float x0_1[DIM], x0_2[DIM], x0_3[DIM];              // initial positions




    return 0;
}