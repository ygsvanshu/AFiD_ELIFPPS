# Inputs

AFiD_ELIPPS inputs consists of required and optional ascii plain-text files all placed in the directory `Inputs/`. The only mandatory input file is `solver.in`. If the solver doesn't find it, the code execution immediately stops with `MPI_Abort`. Conditionally optional input files include `sgrid.in`, `particle.in`, `particle_sources.in`, `movie_slices_x.in`, `movie_slices_y.in`, and `movie_slices_z.in`. These files will be explained in detail in the following links below.

- [`solver.in`](solver.md)
- [`sgrid.in`](sgrid.md)
- [`movie_slices_x.in`](movie_slices.md)
- [`movie_slices_y.in`](movie_slices.md)
- [`movie_slices_z.in`](movie_slices.md)
- [`particle.in`](particle.md)
- [`particle_sources.in`](particle_sources.md)