# `movie_slices_*.in`

The files `movie_slices_x.in`, `movie_slices_y.in`, and `movie_slices_z.in` supply inputs about what slice information to store for post-processing (e.g., making movies of visualizations). The file format is essentially an ordered table as follows. 

- The header (if present) should start with a '`#`'

- Each row (line) starts with the ***nondimensional*** location of the slice and the subsequent columns state whether certain quantities need be stored on disk. Valid inputs for each of the quantities are `T` for true / `F` for false. The fields are further explained below
    - `Location`: The row (line) starts with the ***nondimensional*** location of a single slice plane
    - `SaveVx`: Whether to save the x-component of ***nondimensional*** fluid velocity at the slice
    - `SaveVy`: Whether to save the y-component of ***nondimensional*** fluid velocity at the slice
    - `SaveVz`: Whether to save the z-component of ***nondimensional*** fluid velocity at the slice
    - `SavePr`: Whether to save the x-component of ***nondimensional*** pressure (pressure coefficient) at the slice
    - `SaveOx`: Whether to save the x-component of ***nondimensional*** fluid vorticity at the slice
    - `SaveOy`: Whether to save the y-component of ***nondimensional*** fluid vorticity at the slice
    - `SaveOz`: Whether to save the z-component of ***nondimensional*** fluid vorticity at the slice

For an example, see the corresponding [particle_sources.in](../../../examples/LidDrivenCavity/Inputs/movie_slices_z.in) file in the [lid driven cavity example case](../../../examples/LidDrivenCavity)