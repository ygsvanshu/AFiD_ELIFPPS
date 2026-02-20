# Post-processing [`global.out`](../outputs/global.md)

First, the required libraries (in this case [`numpy`](https://numpy.org/) and [`matplotlib`](https://matplotlib.org/)) are imported.

```python

import numpy as np
import matplotlib.pyplot as plt

```

Since [`global.out`](../outputs/global.md) is an ascii output file written in a space delimited tabular fashion, the data can directly be read in using [`numpy.genfromtxt`](https://numpy.org/doc/stable/reference/generated/numpy.genfromtxt.html) as a 2D array as follows.

```python

# Read data from the .out file, skip the string header in the first line
data = np.genfromtxt("Results/global.out", skip_header=1, dtype=np.float64)

```

Now one can use the data to plot the maximum, average, and RMS values of flow velocities using [`matplotlib`](https://matplotlib.org/). First the figure and subplots are set up with certain figure size and layout.

```python

# Make plots using matplotlib

# Set up figure, subplots, size and layout
fig,axs = plt.subplots(1,3, figsize=(8,3), constrained_layout=True)

```

Now one can use the data to plot the maximum, average, and RMS values of flow velocities.

```python

# Select first subplot
axc = axs[0]
# Set axis title
axs.set_title("Maximum velocity")
# Plot the maximum velocity data
axc.plot(data[:,0], data[:,1], label="Vx")  # (1)!
axc.plot(data[:,0], data[:,2], label="Vy")  # (2)!
axc.plot(data[:,0], data[:,3], label="Vz")  # (3)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional velocity")
# Make legend
axc.legend()

# Select second subplot
axc = axs[0]
# Set axis title
axs.set_title("Average velocity")
# Plot the average velocity data
axc.plot(data[:,0], data[:,4], label="Vx")  # (4)!
axc.plot(data[:,0], data[:,5], label="Vy")  # (5)!
axc.plot(data[:,0], data[:,6], label="Vz")  # (6)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional velocity")
# Make legend
axc.legend()

# Select first subplot
axc = axs[0]
# Set axis title
axs.set_title("Root mean squared (RMS) velocity")
# Plot the root mean squared (RMS) velocity data
axc.plot(data[:,0], data[:,7], label="Vx") # (7)!
axc.plot(data[:,0], data[:,8], label="Vy") # (8)!
axc.plot(data[:,0], data[:,9], label="Vz") # (9)!
axc.plot(data[:,0], data[:,10], label="Magnitude") # (10)!
# Write axis labels
axc.set_xlabel("Nondimensional flow time")
axc.set_ylabel("Nondimensional velocity")
# Make legend
axc.legend()

```

1.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,1]` is subarray [`Max_Vx`](../outputs/global.md) at column 1
2.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,2]` is subarray [`Max_Vy`](../outputs/global.md) at column 2
3.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,3]` is subarray [`Max_Vz`](../outputs/global.md) at column 3
4.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,1]` is subarray [`Avg_Vx`](../outputs/global.md) at column 4
5.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,2]` is subarray [`Avg_Vy`](../outputs/global.md) at column 5
6.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,3]` is subarray [`Avg_Vz`](../outputs/global.md) at column 6
7.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,1]` is subarray [`RMS_Vx`](../outputs/global.md) at column 7
8.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,2]` is subarray [`RMS_Vy`](../outputs/global.md) at column 8
9.  `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,3]` is subarray [`RMS_Vz`](../outputs/global.md) at column 9
10. `data[:,0]` is subarray [`Time`](../outputs/global.md) at column 0 and `data[:,4]` is subarray [`RMS_VxVyVz`](../outputs/global.md) at column 10

Save the figure to file.

```python

# Save figure
plt.savefig("global.png", dpi=300, bbox_inches="tight")

```

Finally close the figure as a good practice.

```python

# Close figure
fig.close()

```