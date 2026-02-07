# Governing differential equations

## Eulerian fluid phase

The fluid flow is governed by the nondimensional incompressible Navier–Stokes equations
$$
\frac{\partial u_i}{\partial t} + \frac{\partial u_i u_j}{\partial x_j} = -\frac{1}{\rho_f}\frac{\partial p}{\partial x_i} + \nu \frac{\partial u_i}{\partial x_j \partial x_j} + a_i + \tilde{a_i},
$$
subject to the divergence free condition
$$
\frac{\partial u_i}{\partial x_i} = 0,
$$
with $u_i \equiv (u_x, u_y, u_z)$ being the velocity, $\rho_f$ being the density of the fluid, $p$ being the pressure, $\nu$ the fluid kinematic viscosity, $a_i$ being the source acceleration term used to model the interaction between the phases, and $\tilde{a_i}$ being any additional source terms introduced by the user. Note that the particle volume fraction is not taken into account under the assumption that it is very small (in the case of a sparse/dilute distribution of particles much smaller than cell sizes).

## Lagrangian particle phase

The point particles are governed by the Maxey–Riley equation, 
although only the steady-state drag and gravity are implemented at the moment. Therefore, the point particles are governed by
$$
\frac{dU_i}{dt} = \frac{3}{4} \frac{C_d}{d_p} \frac{\rho_p}{\rho_f} \left(\tilde{u_i} - U_i\right) \left|\tilde{u_i} - U_i\right| + g_i\left(1 - \frac{\rho_f}{\rho_p}\right),
$$
with $U_i$ being the particle velocity, $C_d$ being the steady state drag coefficient, $\tilde{u_i}$ being the undisturbed fluid velocity at the location of the particle, $\rho_p$ being the density of the particle, $d_p$ being the particle diameter, and $g_i$ being the acceleration due to gravity.

## Particle-Fluid coupling

The source acceleration on the fluid $a_i$ is modelled as a sum of contributions of point forces from particles as
$$
a_i = \Sigma \left( \frac{\rho_f}{\rho_p} \frac{dU_i}{dt} \delta_D(X_i) \right)
$$
with $X_i$ being the location of the particle, and $\delta_D$ being the Dirac delta function. 