# `sgrid.in`

This file consists of a list of coordinates for a custom user generated grid when the input switch `STRTYP` is set to any value other than `1`, `2`, `3`, or `4`.

> **Note**  
> The coordinates listed must be in increasing order, the first coordinate listed in `sgrid.in` must be zero, and the last coordinate listed in `sgrid.in` must be identical to the length of the corresponding axis in `solver.in`, and the number of coordinates must be one larger than the corresponding grid cell count in `solver.in` (for example, if `STRAXS` is set to `y`, then the last coordinate in `sgrid.in` must match the value set for `YLEN`, and the total number of coordinates listed must be equal to one greater than the value set for `NYM`)

The format of the file is as follows

- The header (if present) should start with a '`#`'

- Each line consists of a single ***nondimensional*** location/coordinate of grid point of the stretched grid