# Post processing using python

The files written out by the solver can easily be post-processed using python. Here we make use of well-known libraries like [`numpy`](https://numpy.org/), [`scipy`](https://scipy.org/), [`h5py`](https://www.h5py.org/), and [`matplotlib`](https://matplotlib.org/) to read in the data, process it, and to create plots.

> **Note**  
> For more interactive plotting, especially of 2D slices saved in [`Results/movie_slices_*.h5`](./outputs/movie_slices.md) using a 3D perspective view, the library [`plotly`]() is recommended. For visualizations of 3D flow-fields saved in [`Results/continua_*.h5`](./outputs/continua.md), [Paraview](https://www.paraview.org/) is recommended.

## Processing ascii files with `.out` extension

- [global.out](global.md)
- [particle_stats.out](particle_stats.md)

## Processing HDF5 binary files with `.h5` extension

- [continua.h5](continua.md)
- [movie_slices_*.h5](movie_slices.md)
- [particle_exit.h5](particle_exit.md)
- [particle_history.h5](particle_history.md)