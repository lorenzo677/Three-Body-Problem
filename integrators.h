#include <iostream>
#include <functional>

float RK4(float x0, float y0, float x, float h, float (*func)(float, float)){
    // Finds value of y for a given x using step size h
    // and initial value y0 at x0.
    float k1, k2, k3, k4;

    int n = (int)((x - x0) / h);
    float y = y0;

    for (int i=0; i<n; i++){
        k1 = h * func(x0, y);
        k2 = h * func(x0 + 0.5 * h, y + 0.5 * k1);
        k3 = h * func(x0 + 0.5 * h, y + 0.5 * k2);
        k4 = h * func(x0 + h, y + k3);
  
        // Update next value of y
        y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4);;
  
        // Update next value of x
        x0 = x0 + h;
    }
    return y;
}