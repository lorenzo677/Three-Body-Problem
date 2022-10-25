#include <iostream>


void compute_vector_cm(float mass_1, float mass_2, float mass_3, float vector_1[3], float vector_2[3], float vector_3[3]) 
{   
    float vector_cm[3];
    for (int i = 0; i < 4; i++) {   
        vector_cm[i] = ( mass_1 * vector_1[i] + mass_2 * vector_2[i] + mass_3 * vector_3[i] ) / (mass_1 + mass_2 + mass_3);
    }
}


int main()
{
    // Definition of constants (be careful of Unit Measurements)
    int n = 3; 
    float mass_1, mass_2, mass_3;
    float G = 6.67408e-11; 

}