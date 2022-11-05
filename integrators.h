#include <iostream>
#include <functional>

float RK4(float x0, float y0, float h, float (*func)(float, float, float, float, float, float, float), float m1, float m2, float p1, float p2){
    // Finds value of y for a given x using step size h
    // and initial value y0 at x0.
    float k1, k2, k3, k4;

    k1 = h * func(x0, y0, m1, m2, p1, p2, x0);
    k2 = h * func(x0 + 0.5 * h, y0 + 0.5 * k1, m1, m2, p1, p2, x0);
    k3 = h * func(x0 + 0.5 * h, y0 + 0.5 * k2, m1, m2, p1, p2, x0);
    k4 = h * func(x0 + h, y0 + k3, m1, m2, p1, p2, x0);
  
    // Update next value of y
    y0 = y0 + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);;
  
    // Update next value of x
    x0 = x0 + h;

    return y0;
}
