# Limitations and Roadmap

This document summarizes known limitations of the current
AFiD_ELIFPPS solver and outlines planned developments.

This is a **non-binding roadmap** and subject to change.

---

## Current Limitations

### Physics

- No implementation of Robin boundary conditions
- No temperature field to model heat transfer or buoyancy effects
- Particle–particle interactions (collisions) not yet implemented
- Particle-boundary interactions are not yet implemented, and the default behaviour is to assume an exit event
- Limited range of particle Reynolds numbers, and grid aspect ratios, and ratio of particle diameter to grid cell dimension (see header of the [slip correction data](../data/slip_correction.dat) for parameter ranges)
- Only steady-state drag forces are modelled. This approximation works well for heavy particles (with density much larger than the fluid)
- No model for lift forces


### Numerics

- No immersed boundary implementation for solid geometries

### Parallelization

- Currently only supports manual distribution of MPI processes in pencil rows and columns passed as parameters to the mpirun/mpiexec command. Auto tune mode may fail.
- No load balancing for particles between different processes

---

## Planned Near-Term Features

### Physics (Planned)

- Particle-particle interactions 
- Particle-boundary interactions
- Modelling all terms of Maxey-Riley equation (drag + lift)

### Numerics (Planned)

- Immersed boundary implementation

### Parallelization

- Custom subroutine to safely perform auto-tuning for best performance (including primitive particle load estimation)

---

## How to Contribute

Contributions addressing items in the roadmap are welcome.
Please open an issue to discuss implementation details before
submitting major changes.