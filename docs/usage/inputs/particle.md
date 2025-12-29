# `particle.in`

This file contains the main input parameters required to run the particle solver. This input file is structured to contain a three digit integer followed by relevant values in single quotes. The three-digit number identifies the input variable(s) being read in, and anything in quotes is read in as the input variable(s).

In the absence of this file, any functionality pertaining to the particle solver will be turned off and only uncoupled the fluid flow will be solved. This is useful in case one wishes to solve just for the fluid phase.

## Drag model (100-199)

- `DMOD`: Steady state drag coefficient model (see [Drag coefficient models](../../physics/drag-coefficient.md))
    - `DMOD = 0`: Stokes
    - `DMOD = 1`: Schiller-Naumann
    - `DMOD = 2`: Morsi-Alexander

- `SCOR`: Slip correction 
    - `SCOR = y`: Slip correction enabled
    - `DMOD = n`: Slip correction disabled
    
## Coupling (200-299)

- Coupling force multiplier between Eulerian and Lagrangian phases (i.e., between fluids and particles)
    - `E->L`: Coupling force multiplier from Eulerian phase to Lagrangian phase (i.e., forces exerted by the fluid on the particle)
    - `L->2`: Coupling force multiplier from Lagrangian phase to Eulerian phase (i.e., forces exerted by the particle on the fluid)

## Gravity (300-399)

- Components of acceleration due to gravity
    - `GRVX`: Component of gravitational acceleration vector along x-axis
    - `GRVY`: Component of gravitational acceleration vector along y-axis
    - `GRVZ`: Component of gravitational acceleration vector along z-axis

## Numerical stability (400-499)

- `STOL`: Nondimensinal physical time tolerance/granularity for detecting a new spawn event

- `CLIM`: Maximum prescribed limit for the CFL analogue of particles

- `TLIM`: Maximum prescribed limit for the ratio of adaptive timestep to the particle response time

## Outputs (500-599)

- Saving particle details during events where they exit the domain
    - `EXON`: Saving particle exit events
		- `EXON = y`: Saving enabled
		- `EXON = n`: Saving disabled
	- `EXST`: Save start non-dimensional physical time 
	- `EXFQ`: Non-dimensional physical time Interval between saves

- Saving snapshots of all particles in the domain
    - `LPON`: Saving particle exit events
		- `LPON = y`: Saving enabled
		- `LPON = n`: Saving disabled
	- `LPST`: Save start non-dimensional physical time 
	- `LPFQ`: Non-dimensional physical time Interval between saves