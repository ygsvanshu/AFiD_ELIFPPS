# Post-processing [`particle_stats.out`](../outputs/particle_stats.md)

First, the required libraries (in this case [`numpy`](https://numpy.org/) and [`matplotlib`](https://matplotlib.org/)) are imported.

```python

import numpy as np
import matplotlib.pyplot as plt

```

Since [`particle_stats.out`](../outputs/particle_stats.md) is an ascii output file written in a space delimited tabular fashion, the data can directly be read in using [`numpy.genfromtxt`](https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html) as a 2D array as follows.

```python

# Set datatypes
datatypes = ["f8"]  # (1)!
datatypes = datatypes + ["i4" for i in range(3)] # (2)!
datatypes = datatypes + ["f8" for i in range(18)] # (3)!

# Read data from the .out file, skip the string header in the first line
data = np.genfromtxt("Results/particle_stats.out", skip_header=1, dtype=datatypes)

```

1. For `float64` datatype [`Time`](../outputs/particle_stats.md) at column 0
2. For `int32` dataypes [`Injected`](../outputs/particle_stats.md), [`Active`](../outputs/particle_stats.md), [`Exited`](../outputs/particle_stats.md) at columns 1, 2, 3, respectively.
3. For `float64` [datatypes of maximum, average, and RMS particle Reynolds numbers and body forces](../outputs/particle_stats.md) at columns 4-22

Check that value of [`Injected`](../outputs/particle_stats.md) must be equal to the sum of values of [`Active`](../outputs/particle_stats.md) and [`Exited`](../outputs/particle_stats.md). If not, print out a warning.

```python

result = data[:,1] - data[:,2] - data[:,3] # (1)!
if (np.max(np.abs(result)) != 0): 
    print("Warning! Inconsistency in particle accounting!")

```

2. `data[:,1]` is subarray for [`Injected`](../outputs/particle_stats.md), `data[:,2]` is subarray for [`Active`](../outputs/particle_stats.md), and `data[:,3]` is subarray for [`Exited`](../outputs/particle_stats.md) at columns 1, 2, 3, respectively.

Now one can use the data to plot the maximum, average, and RMS values of particle velocities and body forces using [`matplotlib`](https://matplotlib.org/). First the figure and subplots are set up with certain figure size and layout.

```python

# Make plots using matplotlib

# Set up figure, subplots, size and layout
fig,axs = plt.subplots(2,3, figsize=(8,5), constrained_layout=True)

```

Now one can use the data to plot the maximum, average, and RMS values of particle velocities.

```python

# Select subplot at row=0, column=0
axc = axs[0,0]
# Set axis title
axs.set_title("Maximum velocity")
# Plot the maximum velocity data
axc.plot(data[:,0], data[:,4], label="Vx")  # (1)!
axc.plot(data[:,0], data[:,5], label="Vy")  # (2)!
axc.plot(data[:,0], data[:,6], label="Vz")  # (3)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional velocity")
# Make legend
axc.legend()

# Select subplot at row=0, column=1
axc = axs[0,1]
# Set axis title
axs.set_title("Average velocity")
# Plot the average velocity data
axc.plot(data[:,0], data[:,7], label="Vx")  # (4)!
axc.plot(data[:,0], data[:,8], label="Vy")  # (5)!
axc.plot(data[:,0], data[:,9], label="Vz")  # (6)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional velocity")
# Make legend
axc.legend()

# Select subplot at row=0, column=2
axc = axs[0,2]
# Set axis title
axs.set_title("Root mean squared (RMS) velocity")
# Plot the root mean squared (RMS) velocity data
axc.plot(data[:,0], data[:,10], label="Vx") # (7)!
axc.plot(data[:,0], data[:,11], label="Vy") # (8)!
axc.plot(data[:,0], data[:,12], label="Vz") # (9)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional velocity")
# Make legend
axc.legend()

```

1.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,4]` is subarray [`Re_Max_Vx`](../outputs/particle_stats.md) at column 4
2.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,5]` is subarray [`Re_Max_Vy`](../outputs/particle_stats.md) at column 5
3.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,6]` is subarray [`Re_Max_Vz`](../outputs/particle_stats.md) at column 6
4.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,7]` is subarray [`Re_Avg_Vx`](../outputs/particle_stats.md) at column 7
5.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,8]` is subarray [`Re_Avg_Vy`](../outputs/particle_stats.md) at column 8
6.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,9]` is subarray [`Re_Avg_Vz`](../outputs/particle_stats.md) at column 9
7.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,10]` is subarray [`Re_RMS_Vx`](../outputs/particle_stats.md) at column 10
8.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,11]` is subarray [`Re_RMS_Vy`](../outputs/particle_stats.md) at column 11
9.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,12]` is subarray [`Re_RMS_Vz`](../outputs/particle_stats.md) at column 12

Now one can use the data to plot the maximum, average, and RMS values of body forces applied on the fluid.

```python

# Select subplot at row=1, column=0
axc = axs[1,0]
# Set axis title
axs.set_title("Maximum body force")
# Plot the maximum body force data
axc.plot(data[:,0], data[:,13], label="Bfx")  # (1)!
axc.plot(data[:,0], data[:,14], label="Bfy")  # (2)!
axc.plot(data[:,0], data[:,15], label="Bfz")  # (3)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional body force")
# Make legend
axc.legend()

# Select subplot at row=1, column=1
axc = axs[1,1]
# Set axis title
axs.set_title("Average body force")
# Plot the average body force data
axc.plot(data[:,0], data[:,16], label="Bfx")  # (4)!
axc.plot(data[:,0], data[:,17], label="Bfy")  # (5)!
axc.plot(data[:,0], data[:,18], label="Bfz")  # (6)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional body force")
# Make legend
axc.legend()

# Select subplot at row=1, column=2
axc = axs[1,2]
# Set axis title
axs.set_title("Root mean squared (RMS) body force")
# Plot the root mean squared (RMS) body force data
axc.plot(data[:,0], data[:,19], label="Bfx") # (7)!
axc.plot(data[:,0], data[:,20], label="Bfy") # (8)!
axc.plot(data[:,0], data[:,21], label="Bfz") # (9)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional body force")
# Make legend
axc.legend()

```

1.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,13]` is subarray [`BF_Max_Vx`](../outputs/particle_stats.md) at column 13
2.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,14]` is subarray [`BF_Max_Vy`](../outputs/particle_stats.md) at column 14
3.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,15]` is subarray [`BF_Max_Vz`](../outputs/particle_stats.md) at column 15
4.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,16]` is subarray [`BF_Avg_Vx`](../outputs/particle_stats.md) at column 16
5.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,17]` is subarray [`BF_Avg_Vy`](../outputs/particle_stats.md) at column 17
6.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,18]` is subarray [`BF_Avg_Vz`](../outputs/particle_stats.md) at column 18
7.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,19]` is subarray [`BF_RMS_Vx`](../outputs/particle_stats.md) at column 19
8.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,20]` is subarray [`BF_RMS_Vy`](../outputs/particle_stats.md) at column 20
9.  `data[:,0]` is subarray [`Time`](../outputs/particle_stats.md) at column 0 and `data[:,21]` is subarray [`BF_RMS_Vz`](../outputs/particle_stats.md) at column 21


Save the figure to file.

```python

# Save figure
plt.savefig("particle_stats.png", dpi=300, bbox_inches="tight")

```

Finally close the figure as a good practice.

```python

# Close figure
fig.close()

```