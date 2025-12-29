# Getting Started

AFiD_ELIFPPS is primarily designed to run on high performance computing clusters, but the examples provided can typically be run on modern hardware (laptops/desktops) to get started.

## Prerequisites

The following need to be installed before to successfully build and run AFiD_ELIFPPS

- A Fortran 90 compiler (preferably GNU gfortran)
- LAPACK (or an equivalent package such as Intel MKL)
- MPI implementation (preferably OpenMPI)
- FFTW
- HDF5

## Building the solver

To compile the source code into an executable, a default [Makefile](../Makefile) is provided, but users might want to modify this to their use case, especially since recompiling is necessary when the boundary conditions or initial conditions are changed.

The makefile provided relies on the parallel HDF5-MPI wrapper of the Fortran 90 compiler `h5pfc`. 

To compile, simply navigate to the directory where the source code has been cloned and type
```bash
make
```
This should create the executable `afid`. 


> **Note**   
> For running some example cases, the invocation to compile is slightly different with additional environment/shell variables passed as arguments to the `make` command. For more details, refer to the individual `README` files in the respective example folders.


## Running examples

It is strongly suggested to create the environment variables
- `AFID_PATH` with the explicit path to the location where the executable `afid` is stored. This can be achieved by 
    ```bash
    export AFID_PATH=<path to afid executable>
    ```
- `AFID_DATA` with the explicit path to the location to the directory [`data`](../data) where the file [`slip_correction.dat`](../data/slip_correction.dat) is stored. This can be achieved by 
    ```bash
    export AFID_DATA=<path to data directory containing the file slip_correction.dat>
    ```
    Although this is only strictly necessary when using slip correction for the particle model, it is still recommended to always set this variable as a good practice. In absence of this variable, the solver will fail when using slip correction with the particle model. 

Clone the example [`ParticleDeceleration`](../examples/ParticleDeceleration) and navigate to the directory. Ensure that it contains the directory `Inputs` containing input files.

Running the simulation can be achieved by invoking
```bash
mpirun -np <number of processes> $(AFID_PATH)/afid <number MPI pencil rows> <number of MPI pencil columns>
```
Note that `mpiexec` can also be used in place of `mpirun`. When using a scheduler like SLURM, `srun` may also be used in the batch file if it has been set up appropriately.

After a successful run, the stdout stream to the terminal should read 
```
time greater than tmax
continuation updated
```
A directory `Results` containing the outputs/results of the simulation should be generated.

For more details on how to use the solver for your own application, refer to the [usage](usage/index.md) manual.