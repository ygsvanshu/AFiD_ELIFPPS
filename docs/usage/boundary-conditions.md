# Defining boundary conditions

## Predefined boundary condition configurations

The code provides a set of boundary conditions

## User-defined boundary condition configurations

The code also provides flexibility to the slightly more savvy user (albeit in a slightly involved process) to create custom boundary conditions if so desired. For this, a basic understanding of Fortran90 and its syntax is required.

The boundary conditions are implemented in source files named `BoundaryConditions.F90` in individual directories at [`/src/BoundaryConditions`](../../src/BoundaryConditions/), each for a different configuration. 

If needed, one can implement their own configuration of boundary conditions applied on the six faces of the box domain as outlined below

1. Make a clone/copy of the code
2. Create a directory `MyCustomBC` in [`/src/BoundaryConditions`](../../src/BoundaryConditions/)
3. Copy the file [`BoundaryConditions.F90`](../../src/BoundaryConditions/LidDrivenCavity/BoundaryConditions.F90) from [`/src/BoundaryConditions/LidDrivenCavity`](../../src/BoundaryConditions/LidDrivenCavity) to `/src/BoundaryConditions/MyCustomBC`
4. The file `/src/BoundaryConditions/MyCustomBC/BoundaryConditions.F90` can now be edited as required. A brief note is provided in the comments at the top of the file on how to use the quantities available in the code to modify the file.
5. Modifications for sinusoidal "pump" with open boundaries
6. Then, recompile the code with
    ```bash
    make BOUDIR=MyCustomBC
    ```
7. Run the case (see [getting started](../getting-started.md))
8. Plot data of `Avg_Vy` against data of `Time` from `Results/global.out`. It should be possible to observe a preferential flow in the y-direction.