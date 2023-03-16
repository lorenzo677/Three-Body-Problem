#include <iostream>
#include <fstream>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// Define the state type as a Boost array of double
typedef std::array< double , 18 > state_type;

double G = 6.67e-11; // gravitational constant
double m1 = 1.8981e27, m2 = 2.40e22, m3 = 2.40e22; // masses of the giove and europa
// double m1 = 1.98e30, m2 = 1.62e23, m3 = 1.62e23; // masses of the sun and mercury
// double m1 = 6e24, m2 = 3.671e22, m3 = 3.671e22; // masses of the earth and moon

double RR = 1.496e11; // Normalizing distance in km (= 1 AU)
double MM = 6e24; // Normalizing mass
double TT = 365 * 24 * 60 * 60.0; // nomalizing time (1 year)

double FF = (G * MM * MM) / (RR * RR); // Unit force
double EE = FF * RR; // unit energy   

double GG = (MM*G*TT*TT)/(RR*RR*RR);

double mn1 = m1/MM; // Normalized mass1
double mn2 = m2/MM; // Normalized mass2 
double mn3 = m3/MM; // Normalized mass3

double r12, r13, r23; // distance between bodies
double k = 10e4*RR; // spring constant
double L0 = 1e6/RR; // rest length of the spring

// Declare the energy variable 
double energy = 0;

// Set the output file
ofstream output("output_europa_beta_0.04_final.csv");

// Declare the total energy function
double total_energy( const state_type &x){
    
    double kinetic_energy;
    double potential_energy;

    r12 = sqrt( ( x[0] - x[6] ) * ( x[0] - x[6] ) + ( x[1] - x[7] ) * ( x[1] - x[7] ) + ( x[2] - x[8] ) * ( x[2] - x[8] ) );
    r13 = sqrt( ( x[0] - x[12] ) * ( x[0] - x[12] ) + ( x[1] - x[13] ) * ( x[1] - x[13] ) + ( x[2] - x[14] ) * ( x[2] - x[14] ) );
    r23 = sqrt( ( x[6] - x[12] ) * ( x[6] - x[12] ) + ( x[7] - x[13] ) * ( x[7] - x[13] ) + ( x[8] - x[14] ) * ( x[8] - x[14] ) );

    // kinetic energy
    kinetic_energy = 0.5 * mn1 * (x[3] * x[3] + x[4] * x[4] + x[5] * x[5]) +
                     0.5 * mn2 * (x[9] * x[9] + x[10] * x[10] + x[11] * x[11]) +
                     0.5 * mn3 * (x[15] * x[15] + x[16] * x[16] + x[17] * x[17]);
    
    // gravitational potential energy
    potential_energy  = -GG * mn1 * mn2 / r12 - GG * mn1 * mn3 / r13 - GG * mn2 * mn3 / r23;

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
    double beta = 0.004;

    r12 = sqrt( ( x[0] - x[6] ) * ( x[0] - x[6] ) + ( x[1] - x[7] ) * ( x[1] - x[7] ) + ( x[2] - x[8] ) * ( x[2] - x[8] ) );
    r13 = sqrt( ( x[0] - x[12] ) * ( x[0] - x[12] ) + ( x[1] - x[13] ) * ( x[1] - x[13] ) + ( x[2] - x[14] ) * ( x[2] - x[14] ) );
    r23 = sqrt( ( x[6] - x[12] ) * ( x[6] - x[12] ) + ( x[7] - x[13] ) * ( x[7] - x[13] ) + ( x[8] - x[14] ) * ( x[8] - x[14] ) );

    // First body's position and velocity
    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];

    // First body's acceleration
    dxdt[3] = ( GG * mn2 * ( x[6] - x[0] ) / ( r12 * r12 * r12 ) ) + ( GG * mn3 * ( x[12] - x[0] ) / ( r13 * r13 * r13 ) );
    dxdt[4] = ( GG * mn2 * ( x[7] - x[1] ) / ( r12 * r12 * r12 ) ) + ( GG * mn3 * ( x[13] - x[1] ) / ( r13 * r13 * r13 ) );
    dxdt[5] = ( GG * mn2 * ( x[8] - x[2] ) / ( r12 * r12 * r12 ) ) + ( GG * mn3 * ( x[14] - x[2] ) / ( r13 * r13 * r13 ) );

    // Second body's position and velocity
    dxdt[6] = x[9];
    dxdt[7] = x[10];
    dxdt[8] = x[11];

    // Second body's acceleration
    dxdt[9] = ( GG* mn1 * ( x[0] - x[6] ) / ( r12 * r12 * r12 ) ) + ( GG * mn3 * ( x[12] - x[6] ) / ( r23 * r23 * r23 ) );
    dxdt[10] = ( GG * mn1 * ( x[1] - x[7] ) / ( r12 * r12 * r12 ) ) + ( GG * mn3 * ( x[13] - x[7] ) / ( r23 * r23 * r23 ) );
    dxdt[11] = ( GG * mn1 * ( x[2] - x[8] ) / ( r12 * r12 * r12 ) ) + ( GG * mn3 * ( x[14] - x[8] ) / ( r23 * r23 * r23 ) );

    // Third body's position and velocity
    dxdt[12] = x[15];
    dxdt[13] = x[16];
    dxdt[14] = x[17];

    // Third body's acceleration
    dxdt[15] = ( GG * mn1 * ( x[0] - x[12] ) / ( r13 * r13 * r13 ) ) + ( GG * mn2 * ( x[6] - x[12] ) / ( r23 * r23 * r23 ) );
    dxdt[16] = ( GG * mn1 * ( x[1] - x[13] ) / ( r13 * r13 * r13 ) ) + ( GG * mn2 * ( x[7] - x[13] ) / ( r23 * r23 * r23 ) );
    dxdt[17] = ( GG * mn1 * ( x[2] - x[14] ) / ( r13 * r13 * r13 ) ) + ( GG * mn2 * ( x[8] - x[14] ) / ( r23 * r23 * r23 ) );

    // Spring force between second and third bodies
    double dx = x[6] - x[12];
    double dy = x[7] - x[13];
    double dz = x[8] - x[14];
    
    double dvx = x[9] - x[15];
    double dvy = x[10] - x[16];
    double dvz = x[11] - x[17];
    
    double spring_force = - k * ( sqrt( dx * dx + dy * dy + dz * dz ) - L0 ) - beta * sqrt(dvx * dvx + dvy * dvy + dvz * dvz);

    dxdt[9] += dx * spring_force / r23;  // - beta * x[9];
    dxdt[10] += dy * spring_force / r23; // - beta * x[10];
    dxdt[11] += dz * spring_force / r23; // - beta * x[11];
    dxdt[15] += -dx * spring_force / r23; // - beta * x[15];
    dxdt[16] += -dy * spring_force / r23; // - beta * x[16];
    dxdt[17] += -dz * spring_force / r23; // - beta * x[17];

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
int i=0;
template<class Stepper>
struct observer {
    Stepper &m_stepper;
    observer(Stepper &stepper) : m_stepper(stepper) {}
    void operator()(const state_type &x, double t) {
        energy = total_energy(x);
        
        // Compute the distance of the first body from the line between the second and the third
        double distance = distance_from_line(x);

        // Write the current state to the file
        if (i%5000000==0|| i==0){
        output  << t << "," << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << "," << x[5] << "," << x[6] << "," << x[7] << "," << x[8] << "," << x[9] << "," << x[10] << "," << x[11] << "," << x[12] << "," << x[13] << "," << x[14] << "," << x[15] << "," << x[16] << "," << x[17] << "," << energy << "," << distance <<endl;
        }
        i++;
    }
};


int main()
{
    // Set up the initial conditions
    //                  x     y     z     vx    vy    vz     x      y     z     vx    vy    vz   x     y    z      
    // state_type x = {{ 10.0 , 0.0 , 2.0 , 0.0 , 0.0 , 0.0 , -10.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0 , 0.0 , 0.0 , 0.0 , 0.0, 0.0 }};
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 1.0 , 0.0 , -20.0 , 10.0 , 10.0 , 1.0 , 0.0 , 0.0, -20.0 , 12.0 , 12.0 , -1.0 , 0.0, 0.0 }};
    
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 0.0 , -20.0 , 10.0 , 10.0 , 1.0 , 0.1 , 0.5, -20.0 , 12.0 , 12.0 , 1.0 , 0.1, 0.1 }};
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 0.0 , -20.0 , 10.0 , 10.0 , 1.0 , 0.0 , 0.0, -24.0 , 12.0 , 12.0 , -1.0 , 0.0, 0.0 }};
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , -40.0e9/RR , 10.0e9/RR , 40.0e9/RR , -47.0e3*TT/RR , 0.0 , 5.0e3*TT/RR, -40.0e9/RR , 10.5e9/RR , 40.0e9/RR , 47.0e3*TT/RR , 0.0, 0.0 }}; // vere distanze sole mercurio
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , -6.982e10/RR , 0.0/RR , 0.0 , 0.0 , (-38.86e3-3.0256)*TT/RR , 0.0, -6.981e10/RR , 0.0 , 0.0 , 0.0 , (-38.86e3+3.0256)*TT/RR, 0.0 }}; // vere distanze sole mercurio
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 0.0001 , 0.0 , -4.055e8/RR , 0.0 , 0.0 , 0.0 , (-1076)*TT/RR , 0.0, -4.065e8/RR , 0.0 , 0.0 , 0.0 , (-1077)*TT/RR, 0.0 }}; // vere distanze terra luna
    state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , -6.647e8/RR , 0.0 , 0.0 , 0.0 , (-13871)*TT/RR , 0.0, -6.657e8/RR , 0.0 , 0.0 , 0.0 , (-13872)*TT/RR, 0.0 }}; // vere distanze giove europa
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , -20.0 , 30.0 , 10.0 , -1.0 , -1.0 , 1.0, -20.0 , 29.0 , 11.0 , -1.0 , -1.0, 1.0 }}; 
    
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 2.0 , 0.0 , -20.0 , 10.0 , 10.0 , 1.0 , 0.0 , 1.0, -20.0 , 10.2 , 10.2 , -1.0 , 0.0, 1.0 }}; // belle orbite ma nn torna
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.5 , 0.5 , 0.0 , -20.0 , 10.0 , 10.0 , 1.0 , 0.0 , 2.0, -20.0 , 10.2 , 10.2 , -1.0 , 0.0, 2.0 }}; 
    // state_type x = {{ 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , -31.0 , 0.0 , 0.0 , -1.0 , -1.0 , 0.0, -31.3 , 0.3 , 0.0 , -1.0 , -1.0, 0.0 }};  //belle ma nn torna
    // // Declare the energy variable 
    // double energy = 0;

    // // Declare the distance variable 
    // double distance = 0;

    // Set the integration time step
    double dt = 0.06;

    // Set the end time of the simulation
    double tend = 50;

    // // Set the output file
    // ofstream output("output.csv");

    // Write the headers to the file
    output << "t,x1,y1,z1,vx1,vy1,vz1,x2,y2,z2,vx2,vy2,vz2,x3,y3,z3,vx3,vy3,vz3,energy,distance" << endl;

    // scontrolled_runge_kutta< runge_kutta_dopri5< state_type > > stepper;
    typedef runge_kutta_dopri5< state_type > stepper_type;
    auto controlled_stepper = make_controlled( 1.0e-6 , 1.0e-6 , stepper_type());
    observer<decltype(controlled_stepper)> obs(controlled_stepper);

    // Integrate the equations of motion
    // for (double t = 0; t <= tend; t += dt)
    // {
    //     // Compute the total energy of the system
    //     energy = total_energy(x);
        
    //     // Compute the distance of the first body from the line between the second and the third
    //     double distance = distance_from_line(x);

    //     // Write the current state to the file
    //     output << t << ","  << x[0] << "," << x[1] << "," << x[2] << "," << x[3] << "," << x[4] << "," << x[5] << "," << x[6] << "," << x[7] << "," << x[8] << "," << x[9] << "," << x[10] << "," << x[11] << "," << x[12] << "," << x[13] << "," << x[14] << "," << x[15] << "," << x[16] << "," << x[17] << "," << energy << "," << distance <<endl;
        
    //     integrate_adaptive( stepper, three_body_force, x , 0.0 , dt , dt );
        
    // }
   // observer< runge_kutta_dopri5< state_type>> obs(controlled_stepper);
    integrate_adaptive( controlled_stepper, three_body_force, x , 0.0 , tend , dt, obs );
    // Close the output file
    output.close();

    return 0;
}