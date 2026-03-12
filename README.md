# shallow_water_dycore

This program solves the (inviscid) linearized shallow water equations on a sphere.

Spatial derivatives are calculated using second-order centered finite differences. Tendencies (time derivative quantities) are calculated using SSPRK3 (Strong Stability Preserving Runge-Kutta 3). The staggered grid layout follows the Arakawa C-grid.