/*
 Copyright 2010-2012 Karsten Ahnert
 Copyright 2011-2013 Mario Mulansky
 Copyright 2013 Pascal Germroth

 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <iomanip>   // std::setprecision, std::setw
#include <iostream>
#include <vector>
#include <fstream>
#include <boost/numeric/odeint.hpp>

static constexpr double G = 10;

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
        void computeKineticEnergy(){
            energy = 0.5 * m * (pow(v[0],2)+pow(v[1],2)+pow(v[2],2));
        }
};

Planet A(100, -10, 10, -11, -3, 0, 0);   // corpi allineati sull'asse delle x
Planet B(100, 0, 0, 0, 3, 0, 0);
Planet C(0, 10, 14, 12, 0, 0, 0);

 //[ rhs_function
 /* The type of container used to hold the state vector */
typedef std::vector< double > state_type;

const double gam = 0.15;

/* The rhs of x' = f(x) */
void harmonic_oscillator(const state_type& x, state_type& dxdt, const double /* t */)
{

   double mass_A = A.m;
   double mass_B = B.m;
   double mass_C = C.m;
   double posx_A = A.x[0]; 
   double posx_B = B.x[0]; 
   double posx_C = C.x[0]; 
   double posy_A = A.x[1]; 
   double posy_B = B.x[1]; 
   double posy_C = C.x[1]; 
   double posz_A = A.x[2]; 
   double posz_B = B.x[2]; 
   double posz_C = C.x[2];

   dxdt[0] = x[0];
   dxdt[1] = (-1 * G * (mass_C * (posx_C-x[1]) / pow(sqrt(pow(x[1]-posx_C,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posx_B-x[1]) / pow(sqrt(pow(x[1]-posx_B,2)+pow(posy_A-posy_B,2)+pow(posz_A-posz_B,2)), 3)));
   dxdt[2] = x[2];
   dxdt[3] = (-1 * G * (mass_C * (posy_C-x[3]) / pow(sqrt(pow(posx_C-posx_A,2)+pow(x[3]-posy_C,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posy_B-x[3]) / pow(sqrt(pow(posx_A-posx_B,2)+pow(x[3]-posy_B,2)+pow(posz_A-posz_B,2)), 3)));
   dxdt[4] = x[4];
   dxdt[5] = (-1 * G * (mass_C * (posz_C-x[5]) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(x[5]-posz_C,2)), 3) + mass_B * (posz_B-x[5]) / pow(sqrt(pow(posx_A-posx_B,2)+pow(posy_A-posy_B,2)+pow(x[5]-posz_B,2)), 3)));

   dxdt[6] = x[6];
   dxdt[7] = (-1 * G * (mass_C * (posx_C-x[7]) / pow(sqrt(pow(x[7]-posx_C,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_A * (posx_A-x[7]) / pow(sqrt(pow(x[7]-posx_A,2)+pow(posy_A-posy_B,2)+pow(posz_A-posz_B,2)), 3)));
   dxdt[8] = x[8];
   dxdt[9] = (-1 * G * (mass_C * (posy_C-x[9]) / pow(sqrt(pow(posx_C-posx_A,2)+pow(x[9]-posy_C,2)+pow(posz_C-posz_A,2)), 3) + mass_A * (posy_A-x[9]) / pow(sqrt(pow(posx_A-posx_B,2)+pow(x[9]-posy_B,2)+pow(posz_A-posz_B,2)), 3)));
   dxdt[10] = x[10];
   dxdt[11] = (-1 * G * (mass_C * (posz_C-x[11]) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(x[11]-posz_C,2)), 3) + mass_A * (posz_A-x[11]) / pow(sqrt(pow(posx_A-posx_B,2)+pow(posy_A-posy_B,2)+pow(x[11]-posz_A,2)), 3)));

   dxdt[12] = x[12];
   dxdt[13] = (-1 * G * (mass_A * (posx_A-x[13]) / pow(sqrt(pow(x[13]-posx_A,2)+pow(posy_C-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posx_B-x[13]) / pow(sqrt(pow(x[13]-posx_B,2)+pow(posy_C-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
   dxdt[14] = x[14];
   dxdt[15] = (-1 * G * (mass_A * (posy_A-x[15]) / pow(sqrt(pow(posx_C-posx_A,2)+pow(x[15]-posy_A,2)+pow(posz_C-posz_A,2)), 3) + mass_B * (posy_B-x[15]) / pow(sqrt(pow(posx_C-posx_B,2)+pow(x[15]-posy_B,2)+pow(posz_C-posz_B,2)), 3)));
   dxdt[16] = x[16];
   dxdt[17] = (-1 * G * (mass_A * (posz_A-x[17]) / pow(sqrt(pow(posx_C-posx_A,2)+pow(posy_C-posy_A,2)+pow(x[17]-posz_A,2)), 3) + mass_B * (posz_B-x[17]) / pow(sqrt(pow(posx_C-posx_B,2)+pow(posy_C-posy_B,2)+pow(x[17]-posz_B,2)), 3)));

}
//]



//[ rhs_class
/* The rhs of x' = f(x) defined as a class */
class harm_osc {

   double m_gam;

public:
   harm_osc(double gam) : m_gam(gam) { }

   void operator() (const state_type& x, state_type& dxdt, const double /* t */)
   {
      dxdt[0] = x[1];
      dxdt[1] = -x[0] - m_gam * x[1];
   }
};
//]


//[ integrate_observer
struct push_back_state_and_time
{
   std::vector< state_type >& m_states;
   std::vector< double >& m_times;

   push_back_state_and_time(std::vector< state_type >& states, std::vector< double >& times)
      : m_states(states), m_times(times) { }

   void operator()(const state_type& x, double t)
   {
      m_states.push_back(x);
      m_times.push_back(t);
   }
};
//]
std::ofstream output_file_A("positions_A_odeint.csv");
std::ofstream output_file_B("positions_B_odeint.csv");
std::ofstream output_file_C("positions_C_odeint.csv");
struct write_state
{
   void operator()(const state_type& x) const
   {
      std::cout << x[0] << "\t" << x[1] << "\n";
   }
};
 // Observer, prints time and state when called (during integration)
 void my_observer( const state_type &x, const double t ){
    A.v[0] = x[0];  // positionA x
    A.x[0] = x[1];  
    A.v[1] = x[2]; // positionA y
    A.x[1] = x[3]; 
    A.v[2] = x[4]; // positionA z
    A.x[2] = x[5];
        
    B.v[0] = x[6];  // positionB x
    B.x[0] = x[7];  
    B.v[1] = x[8];  // positionB y
    B.x[1] = x[9]; 
    B.v[2] = x[10];  // positionB z
    B.x[2] = x[11];
    
    C.v[0] = x[12];  // positionC x
    C.x[0] = x[13];  
    C.v[1] = x[14];  // positionC y
    C.x[1] = x[15]; 
    C.v[2] = x[16];  // positionC z
    C.x[2] = x[17];
    output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
    output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
    output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;

    //std::cout  << t << "   " << x[0] <<"\t"<<x[1]<<"\t"<<x[2]<< std::endl;   
    }

int main( )
{
   using namespace std;
   using namespace boost::numeric::odeint;

   

   //[ state_initialization
   state_type x(18);

   runge_kutta4< state_type > stepper;
   x[0] = A.v[0];  // positionA x
   x[1] = A.x[0];  
   x[2] = A.v[1]; // positionA y
   x[3] = A.x[1]; 
   x[4] = A.v[2]; // positionA z
   x[5] = A.x[2];
    
   x[6] = B.v[0];  // positionB x
   x[7] = B.x[0];  
   x[8] = B.v[1];  // positionB y
   x[9] = B.x[1]; 
   x[10] = B.v[2];  // positionB z
   x[11] = B.x[2];
   
   x[12] = C.v[0];  // positionC x
   x[13] = C.x[0];  
   x[14] = C.v[1];  // positionC y
   x[15] = C.x[1]; 
   x[16] = C.v[2];  // positionC z
   x[17] = C.x[2];
    // std::ofstream output_file_A("positions_A_odeint.csv");
    // std::ofstream output_file_B("positions_B_odeint.csv");
    // std::ofstream output_file_C("positions_C_odeint.csv");
    output_file_A<<"x;y;z"<<std::endl;
    output_file_B<<"x;y;z"<<std::endl;
    output_file_C<<"x;y;z"<<std::endl;
    output_file_A << A.x[0] << ";" << A.x[1] << ";" << A.x[2]<< std::endl;
    output_file_B << B.x[0] << ";" << B.x[1] << ";" << B.x[2]<< std::endl;
    output_file_C << C.x[0] << ";" << C.x[1] << ";" << C.x[2]<< std::endl;
   integrate_const(stepper, harmonic_oscillator, x, 0.0, 500.0, 0.008, my_observer);
      output_file_A.close();
    output_file_B.close(); 
    output_file_C.close();
}