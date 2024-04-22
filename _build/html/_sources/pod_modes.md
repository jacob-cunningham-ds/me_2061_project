# POD Modes

Compute the proper orthogonal decomposition (POD) modes for the data matrix `T`. Plot the first 4 POD modes. Plot the first 50 singular values in a semilogy figure. How many modes are required to capture 99% of the energy? Let that number be denoted by $r$ and use $r$ POD modes in the rest of this project.

MATLAB code was used to generate the first four POD modes shown in {numref}`pod_mode_1` through {numref}`pod_mode_4`, plots the first 50 singular values in {numref}`first_50_svs`, finds the number of modes to capture 99% of the total energy, and plots the cumulative energy as a percentage in {numref}`cumsum`.

## First 4 POD Modes

{numref}`pod_mode_1` through {numref}`pod_mode_4` show the first four POD modes.

```{figure} img/pod_mode_1.png
---
name: pod_mode_1
---
POD mode 1
```

```{figure} img/pod_mode_2.png
---
name: pod_mode_2
---
POD mode 2
```

```{figure} img/pod_mode_3.png
---
name: pod_mode_3
---
POD mode 3
```

```{figure} img/pod_mode_4.png
---
name: pod_mode_4
---
POD mode 4
```

## First 50 Singular Values

{numref}`first_50_svs` shows the first 50 singular values.

```{figure} img/first_50_svs.png
---
name: first_50_svs
---
First 50 singular values
```

## Modes Required to Capture 99% of the Energy

The number of modes to required to capture 99% of the energy is **6 modes**.  

The cumulative energy as a percentage is shown in {numref}`cumsum`.

```{figure} img/cumsum.png
---
name: cumsum
---
Cumulative energy as a percentage
```