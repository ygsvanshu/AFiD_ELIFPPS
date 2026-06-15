# `solver.in`

This file contains the main input parameters required to run the fluid solver. This input file is structured to contain a three digit integer followed by relevant values in single quotes. The three-digit number identifies the input variable(s) being read in, and anything in quotes is read in as the input variable(s).

## Initialization options (100-199)

- `NREAD`: Option to restart simulation from a "continua" checkpoint or to start afresh
	- `NREAD = y`: Read in the continua checkpoint and continue calculation
	- `NREAD = n`: Restart calculation from beginning

## Grid and domain specification (200-299)

- Number of grid cells of the Cartesian grid
	> **Note**  
	> We supply the number of grid cells The number of grid points is always one greater than the number of grid cells
	- `NXM`: Number of grid cells along x-axis
	- `NYM`: Number of grid cells along y-axis
	- `NZM`: Number of grid cells along z-axis
	> **Note**  
	> 2D simulations can be run by setting `NZM = 1` However, setting `NYM = 1` will result in the code crashing

- Dimensions of the domain
	- `XLEN`: Dimension of the domain along x-axis
	- `YLEN`: Dimension of the domain along y-axis
	- `ZLEN`: Dimension of the domain along z-axis

- Stretched grid parameters
	- `STRAXS`: Grid stretching axis
		- `STRAXS = x`: Stretched grid in x-direction
		- `STRAXS = y`: Stretched grid in y-direction
		- `STRAXS = z`: Stretched grid in z-direction
		- `STRAXS = n`: Uniform grid in all directions
	- `STRTYP`: Grid point distribution scheme for stretched grid
		- `STRTYP = 1`: Uniform grid with no grid stretching
		- `STRTYP = 2`: Hyperbolic tangent type clustering
		- `STRTYP = 3`: Clipped Chebyshev type clustering
		- `STRTYP = 4`: Clustering using the error function
		- Any other value of `STRTYP`: Read in custom grid from `sgridin`
			> **Note**  
			> If `STRTYP` is not set to a value between 1 and 4, and no `sgridin` file is found in `Inputs/`, then the code crashes
	- `STRVAL`: Grid stretching parameter
		> **Note**  
		> This value should be a positive floating point number Lower values typically correspond to more stretched grid

## Runtime limits (300-399)

- `NTST`: Maximum number of time steps before automatic termination of simulation run

- `TMAX`: Maximum dimensionless physical time before automatic termination of simulation run

- `WALLTIMEMAX`: Maximum wall clock time in seconds before automatic termination of simulation run

## Numerical parameters (400-499)

- `NSST`: Parameter to select the time stepping scheme
	- `NSST = 1`: 2nd order Adams–Bashforth scheme
	- `NSST = 3`: 3rd order Runge–Kutta scheme

- `IDTV`: Parameter to enable or disable dynamic (adaptive) time step
	- `IDTV = y`: Dynamic (adaptive) time step enabled
	- `IDTV = n`: Time step value fixed by the parameter `DT`

- `DT`: Parameter to pass in first time step
	- `IDTV = y`: `DT` is only used for the first time step
	- `IDTV = n`: `DT` is used for all time steps
	> **Note**  
	> At the beginning of the simulation, no information about the previous time step is available, which reduces the stability margin of the AB2 or RK3 time stepping scheme Therefore, a smaller time step of 10E-8 is advised for the startup

- Maximum and minimum limits for dynamic time stepping
	- `DTMIN`: Minimum value permitted for dynamic time step
	- `DTMAX`: Maximum value permitted for dynamic time step
	> **Note**  
	> When the dynamic time stepping is turned on by setting parameter `IDTV = y`, the simulation will try to maximize each time step within the stability margin dictated by constraints arising from the governing differential equations Recommended default values for these parameters are `DTMIN = 10E-8` and `DTMAX = 5E-03`

- Stability parameters for dynamic time stepping
	- `CFLMAX`: Maximum permitted CFL value The dynamic time step is adjusted such that the CFL number always remains less than the supplied value
		> **Note**  
		> For the second order Adams–Bashforth scheme (`NSST = 1`), the theoretical maximum CFL for stable time stepping is 05 while for the third order Runge–Kutta scheme (`NSST = 3`), the theoretical maximum CFL for 

- `RESID`: Maximum permitted divergence If the divergence goes higher than this value, the simulation automatically terminates with an error

- `VLIM`: Maximum permitted value of velocity components
	> **Note**  
	> Since the equations solved by the code are non-dimensional, the reason for a large velocity is often a blow-up caused by a departure from the stability margin of the governing equations

- `CVEL`: Parameter that controls the non-reflecting outflow advection speed.
	> **Note**  
	> A value of 0 would correspond to a static condition fixed at initial state, and a value of infinity would correspond to a zero gradient Neumann condition. Typically lower values (less than 1.0) are more numerically stable.

- Perturbed initial condition
    - `EPSNUM`: Parameter that determines the number of modes of sinusoidal perturbation functions added to the initial condition
    - `EPS`: Parameter that determines the magnitude of perturbation in the scalar initial condition (see [initial conditions]() for further information on the initial conditions)

## Physical parameters (500-599)

- `REY`: Characteristic Reynolds number of the flow (see [Nondimensionalization](../../physics/nondimensionalization.md) for more information)

## Output parameters (600-699)

- `NOUT`: Number of intervals between global statistics calculation and printout to stdout.

- Spatially averaged 1D profiles 
	- `1DON`: Saving spatially averaged profiles
		- `1DON = y`: Saving enabled
		- `1DON = n`: Saving disabled
	- `1DST`: Save start non-dimensional physical time 
	- `1DFQ`: Non-dimensional physical time Interval between saves

- Spatially averaged 2D slices 
	- `2DON`: Saving spatially averaged slices
		- `2DON = y`: Saving enabled
		- `2DON = n`: Saving disabled
	- `2DST`: Save start non-dimensional physical time 
	- `2DFQ`: Non-dimensional physical time Interval between saves

- 3D Flow-field snapshots
	- `3DON`: Saving spatially averaged profiles
		- `3DON = y`: Saving enabled
		- `3DON = n`: Saving disabled
	- `3DST`: Save start non-dimensional physical time 
	- `3DFQ`: Non-dimensional physical time Interval between saves