# `particle_stats.h5`

The ascii file consists of global statistics of particles as a time series. The columns of the data are explained as follows.

1.  `Time`: Nondimensional flow time

2.  `Injected`: Number of particles spawned/injected into the domain
3.  `Active`: Number of particles currently in the domain
4.  `Exited`: Number of particles that have exited out of the domain

    > **Note**  
	> The value of `Injected` must be equal to the sum of values of `Active` and `Exited`. Usually this is the case. In case the sum doesn't hold true, please make a bug report.

5.  `Re_Max_Vx`: Particle Reynolds number corresponding to absolute maximum of x-component of particle velocity
6.  `Re_Max_Vy`: Particle Reynolds number corresponding to absolute maximum of y-component of particle velocity
7.  `Re_Max_Vz`: Particle Reynolds number corresponding to absolute maximum of z-component of particle velocity

8.  `Re_Avg_Vx`: Particle Reynolds number corresponding to the average x-component of particle velocity of all particles in the domain
9.  `Re_Avg_Vy`: Particle Reynolds number corresponding to the average y-component of particle velocity of all particles in the domain
10. `Re_Avg_Vz`: Particle Reynolds number corresponding to the average z-component of particle velocity of all particles in the domain

11. `Re_RMS_Vx`: Particle Reynolds number corresponding to the root mean squared (RMS) x-component of particle velocity of all particles in the domain
12. `Re_RMS_Vy`: Particle Reynolds number corresponding to the root mean squared (RMS) y-component of particle velocity of all particles in the domain
13. `Re_RMS_Vz`: Particle Reynolds number corresponding to the root mean squared (RMS) z-component of particle velocity of all particles in the domain

14. `BF_Max_Vx`: Absolute maximum of x-component of body force applied in the grid cells of the domain 
15. `BF_Max_Vy`: Absolute maximum of y-component of body force applied in the grid cells of the domain 
16. `BF_Max_Vz`: Absolute maximum of z-component of body force applied in the grid cells of the domain 

17. `BF_Avg_Vx`: Average value of x-component of body force applied in the grid cells of the domain
18. `BF_Avg_Vy`: Average value of y-component of body force applied in the grid cells of the domain
19. `BF_Avg_Vz`: Average value of z-component of body force applied in the grid cells of the domain

20. `BF_RMS_Vx`: Root mean squared (RMS) value of x-component of body force applied in the grid cells of the domain
21. `BF_RMS_Vy`: Root mean squared (RMS) value of y-component of body force applied in the grid cells of the domain
22. `BF_RMS_Vz`: Root mean squared (RMS) value of z-component of body force applied in the grid cells of the domain