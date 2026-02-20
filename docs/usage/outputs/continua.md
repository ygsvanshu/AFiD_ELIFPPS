# `continua.h5`

The binary HDF5 file `continua.h5` is a checkpoint that contains all the information that is necessary to restart a simulation. When stored as a [checkpoint](../inputs/solver#output-parameters-600-699), it can be found directly in the run directory. When saved as [3D flow snapshots](../inputs/solver#output-parameters-600-699), a suffix for the nondimensional flow time is added (e.g. `continua_00000100.h5`) and is stored in the `Results` directory. The keys in the output file are explained below.

## General keys

- `nxm`: Number of grid cells in x-direction, datatype: `int32`, shape: `[1]`
- `nym`: Number of grid cells in y-direction, datatype: `int32`, shape: `[1]`
- `nzm`: Number of grid cells in z-direction, datatype: `int32`, shape: `[1]`

- `straxs`: Stretched axis (1 for x-axis, 2 for y-axis, 3 for z-axis), datatype: `int32`, shape: `[1]`

- `perx`: Periodic in x-axis (0 is False and 1 is True), datatype: `int32`, shape: `[1]` 
- `pery`: Periodic in y-axis (0 is False and 1 is True), datatype: `int32`, shape: `[1]` 
- `perz`: Periodic in z-axis (0 is False and 1 is True), datatype: `int32`, shape: `[1]` 

- `xlen`: Domain extent along x-axis (nondimensional), datatype: `float64`, shape: `[1]`
- `ylen`: Domain extent along y-axis (nondimensional), datatype: `float64`, shape: `[1]`
- `zlen`: Domain extent along z-axis (nondimensional), datatype: `float64`, shape: `[1]`

- `xc`: Grid node (face center) coordinates for x-axis (nondimensional), datatype: `float64`, shape: `[nx]`
- `xm`: Grid cell (cell center) coordinates for x-axis (nondimensional), datatype: `float64`, shape: `[nxm]` 
- `yc`: Grid node (face center) coordinates for y-axis (nondimensional), datatype: `float64`, shape: `[ny]`
- `ym`: Grid cell (cell center) coordinates for y-axis (nondimensional), datatype: `float64`, shape: `[nym]` 
- `zc`: Grid node (face center) coordinates for z-axis (nondimensional), datatype: `float64`, shape: `[nz]`
- `zm`: Grid cell (cell center) coordinates for z-axis (nondimensional), datatype: `float64`, shape: `[nzm]` 

- `time`: Nondimensional flow time at which the `continua.h5` file was stored, datatype: `float64`, shape: `[1]`  
- `dt`: Nondimensional time step after which the `continua.h5` file was stored, datatype: `float64`, shape: `[1]` 
- `rey`: Reynolds number, datatype: `float64`, shape: `[1]`

- `vx`: 3D array of nondimensional velocity component along x-axis, datatype: `float64`, shape: `[nzm,nym,nx]`
- `vy`: 3D array of nondimensional velocity component along y-axis, datatype: `float64`, shape: `[nzm,ny,nxm]`
- `vz`: 3D array of nondimensional velocity component along z-axis, datatype: `float64`, shape: `[nz,nym,nxm]`
- `pr`: 3D array of nondimensional pressure (pressure coefficient), datatype: `float64`, shape: `[nzm,nym,nxm]` 

## Keys pertaining to particle solver

- `lpp_dmod`: [Drag coefficient model](../../physics/drag-coefficient.md) for particle solver, datatype: `int32`, shape: `[1]` [(0 for Stokes, 1 for Schiller-Naumann, 2 for Morsi-Alexander)](../inputs/particle.md#drag-model-100-199), datatype: `int32`, shape: `[1]`
- `lpp_scor`: [Slip correction](../../numerics/slip-correction.md) [(0 is False and 1 is True)](../inputs/particle.md#coupling-200-299), datatype: `int32`, shape: `[1]`
- `lpp_grav`: [Acceleration due to gravity](../inputs/particle.md#gravity-300-399), datatype: `float64`, shape: `[3]`
- `e2l_mult`: [Eulerian -> Lagrangian coupling multiplier](../inputs/particle.md#coupling-200-299), datatype: `float64`, shape: `[1]`
- `l2e_mult`: [Lagrangian -> Eulerian coupling multiplier](../inputs/particle.md#coupling-200-299), datatype: `float64`, shape: `[1]`
- `lpp_spwn`: Total number of particles spawned until the `continua.h5` file was stored, datatype: `int32`, shape: `[1]`
- `lpp_exit`: Total number of particles that exited the domain until the `continua.h5` file was stored, datatype: `int32`, shape: `[1]` 

- `lpp_num`: Number of particles in the domain at the instant the `continua.h5` file was stored, datatype: `int32`, shape: `[1]`

- `src_idx`: Source indices of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `int32`, shape: `[lpp_num]`
- `lpp_lft`: Life time since spawn of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[lpp_num]`
- `lpp_dia`: Diameters of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[lpp_num]`
- `lpp_den`: Densities of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[lpp_num]`
- `lpp_rey`: Particle reynolds numbers of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[lpp_num]`
- `lpp_pos`: Positions of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[3,lpp_num]`
- `lpp_vel`: Velocities of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[3,lpp_num]`
- `acc_old`: Accelerations at previous substep of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[3,lpp_num]`
- `acc_now`: Accelerations at current substep of all particles in the domain at the instant the `continua.h5` file was stored, datatype: `float64`, shape: `[3,lpp_num]`