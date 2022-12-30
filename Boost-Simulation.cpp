#include <boost/numeric/odeint.hpp>
#include <array>

using namespace std;
using namespace boost::numeric::odeint;

// Definire le costanti del sistema
const double G = 6.67408e-11; // Costante gravitazionale
const double M1 = 1.989e30; // Massa del corpo 1
const double M2 = 5.972e24; // Massa del corpo 2
const double M3 = 7.347e22; // Massa del corpo 3

// Definire le variabili del sistema
array<double, 3> x1, x2, x3; // Posizioni dei tre corpi
array<double, 3> v1, v2, v3; // Velocità dei tre corpi

// Definire le funzioni di derivata per il sistema di equazioni differenziali
void derivs(const array<double, 6> &y, array<double, 6> &dydt, double t)
{
// Calcolare le derivate delle posizioni dei tre corpi
dydt[0] = y[3];
dydt[1] = y[4];
dydt[2] = y[5];

// Calcolare le derivate delle velocità dei tre corpi
dydt[3] = (G * M2 * (y[0] - y[3])) / pow(pow(y[0] - y[3], 2) + pow(y[1] - y[4], 2) + pow(y[2] - y[5], 2), 1.5);
dydt[4] = (G * M2 * (y[1] - y[4])) / pow(pow(y[0] - y[3], 2) + pow(y[1] - y[4], 2) + pow(y[2] - y[5], 2), 1.5);
dydt[5] = (G * M2 * (y[2] - y[5])) / pow(pow(y[0] - y[3], 2) + pow(y[1] - y[4], 2) + pow(y[2] - y[5], 2), 1.5);


}

int main()
{
// Definire l'integratore
runge_kutta4<array<double, 6>> integrator;

// Definire il vettore di stato iniziale
array<double, 6> y = {x1[0], x1[1], x1[2], v1[0], v1[1], v1[2]};

// Eseguire l'integrazione per un intervallo di tempo specificato
integrator.do_step(derivs, y, 0.0, 1.0);

// Stampa dei risultati
cout << "Posizione finale del corpo 1: (" << y[0] << ", " << y

}