# `particle_sources.in`

This file essentially contains a table with details of particle sources in the following format

- The header (if present) should start with a '`#`'

- Each row (line) consists of ***nondimensional*** information about a single source

- The following ordered columns in each row give details about the source 
    - `StartTime`: Particle injection start time
    - `EndTime`: Particle injection stop time
    - `Frequency`: Particle injection frequency 
    - `Diameter`: Diameter of particles injected by the source
    - `Density`: Relative density of particles injected by the source
    - `PosX`: Source location x-coordinate
    - `PosY`: Source location y-coordinate
    - `PosZ`: Source location z-coordinate
    - `VelX`: x-component of injected particle velocity 
    - `VelY`: y-component of injected particle velocity 
    - `VelZ`: z-component of injected particle velocity 

- For an example, see the corresponding [particle_sources.in](../../../examples/ParticleDeceleration/Inputs/particle_sources.in) file in the [decelerating particle example case](../../../examples/ParticleDeceleration)