# shallow_water_dycore

This program solves the (inviscid) linearized shallow water equations on a sphere.

Spatial derivatives are calculated using centered finite differences. Tendencies (time derivative quantities) are calculated using a forward Euler scheme. The grid layout follows the Arakawa C-grid. 