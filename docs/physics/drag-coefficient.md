# Drag coefficient models

The steady state drag coefficient in [governing differential equations](governing-equations.md) is modelled as the drag coefficient of a sphere of a given diameter subject to free stream flow. Several drag coefficient models can be found in literature, but three of these are implemented in the code currently. All of these drag coefficients are given as functions of the particle Reynolds number $Re_p \equiv U_\infty d_p/ \nu$ with $U_\infty$ being the free stream velocity, $d_p$ being the particle diameter, and $\nu$ being the kinematic viscosity of air.

## Stokes

The Stokes drag coefficient model derived from the linear drag on a sphere in Stokes flow given by
$$
C_d(Re_p) = \frac{24}{Re_p},
$$
is only valid for small Reynolds numbers in the Stokes regime, typically for $Re_p \ll 1$.

## Schiller-Naumann

The model proposed by Schiller and Naumann [1] given by
$$
C_D = 
\begin{cases} 
\frac{24}{Re}(1 + 0.15 Re^{0.687}) & \text{for} \quad Re_p \leq 1000 \\
0.44 & \text{for} \quad Re_p > 1000 
\end{cases}
$$
works well for particles in gases, especially at moderate particle Reynolds numbers.

## Morsi-Alexander

Another empirical formulation by Morsi and Alexander [2] applicable for a large range of particle Reynolds numbers is given by
$$
C_d(Re_p) = c_1 + \frac{c_2}{Re_p} + \frac{c_2}{{Re_p}^2},
$$
with
$$
\begin{align*}
c_1 &= 0,      \quad & c_2 &= 24,       \quad & c_3 &= 0       \quad & \text{for} \quad & Re_p < 0.1 \\
c_1 &= 3.690,  \quad & c_2 &= 22.73,    \quad & c_3 &= 0.0903  \quad & \text{for} \quad & 0.1 \le Re_p < 1 \\
c_1 &= 1.222,  \quad & c_2 &= 29.1667,  \quad & c_3 &= -3.8889 \quad & \text{for} \quad & 1 \le Re_p < 10 \\
c_1 &= 0.6167, \quad & c_2 &= 46.50,    \quad & c_3 &= -116.67 \quad & \text{for} \quad & 10 \le Re_p < 100 \\
c_1 &= 0.3644, \quad & c_2 &= 98.33,    \quad & c_3 &= -2778   \quad & \text{for} \quad & 100 \le Re_p < 1000 \\
c_1 &= 0.357,  \quad & c_2 &= 148.62,   \quad & c_3 &= -47500  \quad & \text{for} \quad & 1000 \le Re_p < 5000 \\
c_1 &= 0.46,   \quad & c_2 &= -490.546, \quad & c_3 &= 578700  \quad & \text{for} \quad & 5000 \le Re_p < 10000 \\
c_1 &= 0.5191, \quad & c_2 &= -1662.5,  \quad & c_3 &= 5416700 \quad & \text{for} \quad & Re_p \ge 10000
\end{align*}
$$
