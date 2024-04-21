# Data

Generate a data matrix by solving the full order model (FOM) for $N_{s} = 11$ samples for the following choices of heat sources: $\theta = $`linspace(180, 0, Ns)`. Store the resulting data in a matrix named `T`. The size of this matrix must be $N \times (N_{t}N_{s})$, where $n$ is the number of FEM nodes. Plot the last snapshot of the temperature for each of the cases.

MATLAB code was used to generate the matrix `T` and plots the last snaphshot of the temperature for each of the cases as shown in {numref}`fom_theta_0` through {numref}`fom_theta_180`.

## T Matrix

The size of the `T` matrix is $4338 \times 2211$. `N` is $4338$, `ntime` is $201$, `Ns` is $11$, and `ntime` $\times$ `Ns` is $2211$.  Therefore `T` has the correct dimensions.

## Plots

```{figure} img/fom_theta_0.png
---
name: fom_theta_0
---
Temperature profile for last time snapshot at $\theta = 0^{o}$
```

```{figure} img/fom_theta_18.png
---
name: fom_theta_18
---
Temperature profile for last time snapshot at $\theta = 18^{o}$
```

```{figure} img/fom_theta_36.png
---
name: fom_theta_36
---
Temperature profile for last time snapshot at $\theta = 36^{o}$
```

```{figure} img/fom_theta_54.png
---
name: fom_theta_54
---
Temperature profile for last time snapshot at $\theta = 54^{o}$
```

```{figure} img/fom_theta_72.png
---
name: fom_theta_72
---
Temperature profile for last time snapshot at $\theta = 72^{o}$
```

```{figure} img/fom_theta_90.png
---
name: fom_theta_90
---
Temperature profile for last time snapshot at $\theta = 90^{o}$
```

```{figure} img/fom_theta_108.png
---
name: fom_theta_108
---
Temperature profile for last time snapshot at $\theta = 108^{o}$
```

```{figure} img/fom_theta_126.png
---
name: fom_theta_126
---
Temperature profile for last time snapshot at $\theta = 126^{o}$
```

```{figure} img/fom_theta_144.png
---
name: fom_theta_144
---
Temperature profile for last time snapshot at $\theta = 144^{o}$
```

```{figure} img/fom_theta_162.png
---
name: fom_theta_162
---
Temperature profile for last time snapshot at $\theta = 162^{o}$
```

```{figure} img/fom_theta_180.png
---
name: fom_theta_180
---
Temperature profile for last time snapshot at $\theta = 180^{o}$
```