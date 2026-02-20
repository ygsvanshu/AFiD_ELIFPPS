# `particle_history.h5`

The binary HDF5 file `particle_history.h5` consists of a series of snapshots containing information about all the particles within the domain at the time instant the snapshot was recorded. To keep a consistent track of the case information for each run in the case of continued simulations, the file is designed as follows.

- `runs`: The number of runs, datatype: `int32`, shape: `[1]`
- `steps`: The total number of snapshots recorded in the file
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
- `data`: Data of particle snapshots
    - step: Selected snapshot ranging from 1 to `steps`
        - `time`: Nondimensional flow time when the snapshot was recorded
        - `dt`: Nondimensional time step when the snapshot was recorded
        - `lpp_num`: Number of particles in the domain at the instant the snapshot was recorded, datatype: `int32`, shape: `[1]`
        - `src_idx`: Source indices of all particles in the domain at the instant the snapshot was recorded, datatype: `int32`, shape: `[lpp_num]`
        - `lpp_lft`: Life time since spawn of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[lpp_num]`
        - `lpp_dia`: Diameters of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[lpp_num]`
        - `lpp_den`: Densities of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[lpp_num]`
        - `lpp_rey`: Particle reynolds numbers of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[lpp_num]`
        - `lpp_pos`: Positions of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[3,lpp_num]`
        - `lpp_vel`: Velocities of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[3,lpp_num]`
        - `acc_old`: Accelerations at previous substep of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[3,lpp_num]`
        - `acc_now`: Accelerations at current substep of all particles in the domain at the instant the snapshot was recorded, datatype: `float64`, shape: `[3,lpp_num]`