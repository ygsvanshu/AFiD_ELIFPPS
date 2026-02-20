# `particle_exit.h5`

The binary HDF5 file `particle_exit.h5` consists of information of all the particles that exited the domain. To keep a consistent track of the case information for each run in the case of continued simulations, the file is designed as follows.

- `runs`: The number of runs, datatype: `int32`, shape: `[1]`
- `info`: The case information corresponding to all runs
    - run index: Run indices from 1 to `runs`
        - `stt_indx`: Start index of selected run index for arrays in `data`, datatype: `int32`, shape: `[1]`
        - `lpp_dmod`: [Drag coefficient model](../../physics/drag-coefficient.md) for the selected run index, datatype: `int32`, shape: `[1]` [(0 for Stokes, 1 for Schiller-Naumann, 2 for Morsi-Alexander)](../inputs/particle.md#drag-model-100-199), datatype: `int32`, shape: `[1]`
        - `lpp_scor`: [Slip correction](../../numerics/slip-correction.md) for the selected run index [(0 is False and 1 is True)](../inputs/particle.md#coupling-200-299), datatype: `int32`, shape: `[1]`
        - `lpp_grav`: [Acceleration due to gravity](../inputs/particle.md#gravity-300-399) for the selected run index, datatype: `float64`, shape: `[3]`
        - `e2l_mult`: [Eulerian -> Lagrangian coupling multiplier](../inputs/particle.md#coupling-200-299) for the selected run index, datatype: `float64`, shape: `[1]`
        - `l2e_mult`: [Lagrangian -> Eulerian coupling multiplier](../inputs/particle.md#coupling-200-299) for the selected run index, datatype: `float64`, shape: `[1]`
        - `lpp_srcs`: List of particle sources that were used for the selected run index
            - `src_num`: Number of sources present during the selected run index, datatype: `int32`, shape: `[1]`
            - `src_idx`: Indices of sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
            - `src_sta`: Spawn start times of sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
            - `src_end`: Spawn end times of sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
            - `src_frq`: Spawn frequencies of sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
            - `src_dia`: Diameters of particles spawned by sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
            - `src_den`: Densities of particles of sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
            - `src_pos`: Positions (nondimensional) of sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
            - `src_vel`: Velocities (nondimensional) of sources present during the selected run index, datatype: `float64`, shape: `[src_num]`
- `data`: The data of exited particles
    - `pex_num`: Total number of exited particles, datatype: `int32`, shape: `[1]`
    - `src_idx`: Source indices of exited particles, datatype: `int32`, shape: `[pex_num]`
    - `pex_pln`: Exit plane (`x`, `y`, or `z`) of exited particles, datatype: `str` of length 1, shape: `[pex_num]`
    - `pex_lft`: Life time (nondimensional) since spawn for exited particles, datatype: `float64`, shape: `[pex_num]`
    - `pex_eft`: Flow time (nondimensional) at the exit event for exited particles, datatype: `float64`, shape: `[pex_num]`
    - `pex_pos`: Position (nondimensional) at the exit event for exited particles, datatype: `float64`, shape: `[3,pex_num]`
    - `pex_vel`: Velocity (nondimensional) at the exit event for exited particles, datatype: `float64`, shape: `[3,pex_num]`
