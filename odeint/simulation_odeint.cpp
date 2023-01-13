#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// Define the state type as a Boost array of double
typedef boost::array< double , 18 > state_type;

// Declare the total energy function
double total_energy( const state_type &x){
    
    double G = 10.0; // gravitational constant
    double m1 = 20.0, m2 = 1.0, m3 = 1.0; // masses of the three bodies
    double r12, r13, r23; // distance between bodies
    double k = 1.0e4; // spring constant
    double L0 = sqrt(8); // rest length of the spring
    double kinetic_energy;
    double potential_energy;

    r12 = sqrt( ( x[0] - x[6] ) * ( x[0] - x[6] ) + ( x[1] - x[7] ) * ( x[1] - x[7] ) + ( x[2] - x[8] ) * ( x[2] - x[8] ) );
    r13 = sqrt( ( x[0] - x[12] ) * ( x[0] - x[12] ) + ( x[1] - x[13] ) * ( x[1] - x[13] ) + ( x[2] - x[14] ) * ( x[2] - x[14] ) );
    r23 = sqrt( ( x[6] - x[12] ) * ( x[6] - x[12] ) + ( x[7] - x[13] ) * ( x[7] - x[13] ) + ( x[8] - x[14] ) * ( x[8] - x[14] ) );

    // kinetic energy
    kinetic_energy = 0.5 * m1 * (x[3] * x[3] + x[4] * x[4] + x[5] * x[5]) +
                     0.5 * m2 * (x[9] * x[9] + x[10] * x[10] + x[11] * x[11]) +
                     0.5 * m3 * (x[15] * x[15] + x[16] * x[16] + x[17] * x[17]);
    
    // gravitational potential energy
    potential_energy  = -G * m1 * m2 / r12 - G * m1 * m3 / r13 - G * m2 * m3 / r23;

    // elastic potential energy
    double dx = x[6] - x[12];
    double dy = x[7] - x[13];
    double dz = x[8] - x[14];
    potential_energy += 0.5 * k * (sqrt( dx * dx + dy * dy + dz * dz) - L0) * (sqrt( dx * dx + dy * dy + dz * dz) - L0);

    return kinetic_energy + potential_energy;
}

// Declare the force function
void three_body_force( const state_type &x , state_type &dxdt , double t )
{
    double G = 10.0; // gravitational constant
    double m1 = 20.0, m2 = 1.0, m3 = 1.0; // masses of the three bodies
    double r12, r13, r23; // distance between bodies
    double k = 1.0e4; // spring constant
    double L0 = sqrt(8); // rest length of the spring

    r12 = sqrt( ( x[0] - x[6] ) * ( x[0] - x[6] ) + ( x[1] - x[7] ) * ( x[1] - x[7] ) + ( x[2] - x[8] ) * ( x[2] - x[8] ) );
    r13 = sqrt( ( x[0] - x[12] ) * ( x[0] - x[12] ) + ( x[1] - x[13] ) * ( x[1] - x[13] ) + ( x[2] - x[14] ) * ( x[2] - x[14] ) );
    r23 = sqrt( ( x[6] - x[12] ) * ( x[6] - x[12] ) + ( x[7] - x[13] ) * ( x[7] - x[13] ) + ( x[8] - x[14] ) * ( x[8] - x[14] ) );

    // First body's position and velocity
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];

    // First body's acceleration
    dxdt[3] = ( G * m2 * ( x[6] - x[0] ) / ( r12 * r12 * r12 ) ) + ( G * m3 * ( x[12] - x[0] ) / ( r13 * r13 * r13 ) );
    dxdt[4] = ( G * m2 * ( x[7] - x[1] ) / ( r12 * r12 * r12 ) ) + ( G * m3 * ( x[13] - x[1] ) / ( r13 * r13 * r13 ) );
    dxdt[5] = ( G * m2 * ( x[8] - x[2] ) / ( r12 * r12 * r12 ) ) + ( G * m3 * ( x[14] - x[2] ) / ( r13 * r13 * r13 ) );

    // Second body's position and velocity
    dxdt[6] = x[9];
    dxdt[7] = x[10];
    dxdt[8] = x[11];

    // Second body's acceleration
    dxdt[9] = ( G * m1 * ( x[0] - x[6] ) / ( r12 * r12 * r12 ) ) + ( G * m3 * ( x[12] - x[6] ) / ( r23 * r23 * r23 ) );
    dxdt[10] = ( G * m1 * ( x[1] - x[7] ) / ( r12 * r12 * r12 ) ) + ( G * m3 * ( x[13] - x[7] ) / ( r23 * r23 * r23 ) );
    dxdt[11] = ( G * m1 * ( x[2] - x[8] ) / ( r12 * r12 * r12 ) ) + ( G * m3 * ( x[14] - x[8] ) / ( r23 * r23 * r23 ) );

    // Third body's position and velocity
    dxdt[12] = x[15];
    dxdt[13] = x[16];
    dxdt[14] = x[17];

    // Third body's acceleration
    dxdt[15] = ( G * m1 * ( x[0] - x[12] ) / ( r13 * r13 * r13 ) ) + ( G * m2 * ( x[6] - x[12] ) / ( r23 * r23 * r23 ) );
    dxdt[16] = ( G * m1 * ( x[1] - x[13] ) / ( r13 * r13 * r13 ) ) + ( G * m2 * ( x[7] - x[13] ) / ( r23 * r23 * r23 ) );
    dxdt[17] = ( G * m1 * ( x[2] - x[14] ) / ( r13 * r13 * r13 ) ) + ( G * m2 * ( x[8] - x[14] ) / ( r23 * r23 * r23 ) );

    // Spring force between second and third bodies
    double dx = x[6] - x[12];
    double dy = x[7] - x[13];
    double dz = x[8] - x[14];
    double spring_force = - k * ( sqrt( dx * dx + dy * dy + dz * dz ) - L0 );

    dxdt[9] += dx * spring_force / r23;
    dxdt[10] += dy * spring_force / r23;
    dxdt[11] += dz * spring_force / r23;
    dxdt[15] += -dx * spring_force / r23;
    dxdt[16] += -dy * spring_force / r23;
    dxdt[17] += -dz * spring_force / r23;

}

// Declare the function that compute the distance of the first body from the line between the second and the third
double distance_from_line(const state_type &x)
{
    double x1 = x[0], y1 = x[1], z1 = x[2]; // position of first body
    double x2 = x[6], y2 = x[7], z2 = x[8]; // position of second body
    double x3 = x[12], y3 = x[13], z3 = x[14]; // position of third body

    double a = (y2 - y3) * z1 - (z2 - z3) * y1 + (z2 * y3 - y2 * z3);
    double b = (z2 - z3) * x1 - (x2 - x3) * z1 + (x2 * z3 - z2 * x3);
    double c = (x2 - x3) * y1 - (y2 - y3) * x1 + (y2 * x3 - x2 * y3);
    double d = sqrt(a * a + b * b + c * c) / sqrt((y2 - y3) * (y2 - y3) + (z2 - z3) * (z2 - z3) + (x2 - x3) * (x2 - x3));
    return d;
}


int main()
{
    // Set up the initial conditions
    //                  x     y     z     vx    vy    vz     x      y     z     vx    vy    vz   x     y    z      
    // state_type x = {{ 10.0 , 0.0 , 2.0 , 0.0 , 0.0 , 0.0 , -10.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0 }};
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , -20.0 , 10.0 , 10.0 , 1.0 , 0.0 , 0.0, -20.0 , 12.0 , 12.0 , -1.0 , 0.0, 0.0 }};
    
    state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 0.0 , -20.0 , 10.0 , 10.0 , 1.0 , 0.0 , 0.0, -20.0 , 12.0 , 12.0 , -1.0 , 0.0, 0.0 }};

    // Declare the energy variable 
    double energy = 0;

    // Declare the distance variable 
    double distance = 0;

    // Set the integration time step
    double dt = 0.002;

    // Set the end time of the simulation
    double tend = 140;

    // Set the output file
    ofstream output("output.csv");

    // Write the headers to the file
    output << "t,x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,x3,y3,z3,vx3,vy3,vz3,energy,distance" << endl;

    controlled_runge_kutta< runge_kutta_dopri5< state_type > > stepper;

    // Integrate the equations of motion
    for (double t = 0; t <= tend; t += dt)
    {
        // Compute the total energy of the system
        energy = total_energy(x);
        
        // Compute the distance of the first body from the line between the second and the third
        distance = distance_from_line(x);

        // Write the current state to the file
        output << t << "," << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << "," << x[5] << "," << x[6] << "," << x[7] << "," << x[8] << "," << x[9] << "," << x[10] << "," << x[11] << "," << x[12] << "," << x[13] << "," << x[14] << "," << x[15] << "," << x[16] << "," << x[17] << "," << energy << "," << distance <<endl;
        
        integrate_adaptive( stepper, three_body_force, x , 0.0 , dt , dt );
    }

    // Close the output file
    output.close();

    return 0;
}