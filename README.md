# Three-Body-Problem

These scripts are useful to integrate the three-body-problem. in addition to the classical configuration, it is possible to introduce a spring between two of the three bodies. This allow, for example, to simulate the Moon-Earth locked system, showing that after a certain time, the moon shows to the Earth always the same face.

## File explication
The file `Simulation.cpp` is the script to integrate the system with Leaprfrog, Symplectic Euler and Runge Kutta 4.
Instead in the cartel `odeint`, in the file `simulation_odeint.cpp` is written the code to integrate the system with the `boost` library


![](https://github.com/lorenzo677/Three-Body-Problem/blob/main/animation.gif)
