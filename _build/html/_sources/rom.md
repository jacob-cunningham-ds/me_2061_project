# ROM

Build the reduced order model (ROM) low-rank matrix and vector `Kr` and `fr`, where `Kr` $\in \mathbb{R}^{r \times r}$ is the reduced stiffness matrix and `fr` $\in \mathbb{R}^{r \times 1}$ for any given parameter $\theta$. 

To solve the ROM you must first complete the function `ROM_rhs` and then use `SolveROM`. Check the validity of the ROM by solving the ROM for the last sample $\theta = 0$. 

Compare the FOM and ROM temperature contours at the last snapshot. Also, compare the FOM and ROM temperature profile on the top half of the cylinder at times $t = 0.01, 0.03, 0.05$. Use the function `ExtractCylinderT` to extract the temperature on the top half of the cylinder surface ($0 \leq \theta \leq 180$). An example of using this function is shown in `main.m`.

MATLAB code was used to build the ROM, check the validity of the ROM by solving the ROM for the last sample $\theta = 0$ as shown in {numref}`fom_rom_theta_0`, compare the FOM and ROM temperature contours at the last snapshot as shown in {numref}`fom_rom_theta_0` through {numref}`fom_rom_theta_180`, and compare the FOM and ROM temperature profile on the top half of the cylinder at times $t = 0.01, 0.03, 0.05$ in {numref}`cyl_theta_0` through {numref}`cyl_theta_180`.

## Comparing FOM and ROM Last Snapshots

```{note}
The ROM was constructed from a non-dimensionalized form of the FOM.
```

The ROM captures the relative pattern of the temperature distribution but not the scale.  This could be due to the fact that the ROM was constructed from a non-dimensionalized form of the FOM.

The ROM took roughly 0.0082 seconds to compute the result at $\theta = 0$ whereas the FOM took roughly 1.91 seconds to compute the result at $\theta = 0$.  The ROM is over 200 times faster than the FOM at $\theta = 0$.  This trend is expected to be true for the other samples.

```{figure} img/fom_rom_theta_0.png
---
name: fom_rom_theta_0
---
Temperature profile for last time snapshot at $\theta = 0^{o}$
```

```{figure} img/fom_rom_theta_18.png
---
name: fom_rom_theta_18
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 18^{o}$
```

```{figure} img/fom_rom_theta_36.png
---
name: fom_rom_theta_36
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 36^{o}$
```

```{figure} img/fom_rom_theta_54.png
---
name: fom_rom_theta_54
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 54^{o}$
```
```{figure} img/fom_rom_theta_72.png
---
name: fom_rom_theta_72
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 72^{o}$
```

```{figure} img/fom_rom_theta_90.png
---
name: fom_rom_theta_90
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 90^{o}$
```

```{figure} img/fom_rom_theta_108.png
---
name: fom_rom_theta_108
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 108^{o}$
```

```{figure} img/fom_rom_theta_126.png
---
name: fom_rom_theta_126
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 126^{o}$
```

```{figure} img/fom_rom_theta_144.png
---
name: fom_rom_theta_144
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 144^{o}$
```

```{figure} img/fom_rom_theta_162.png
---
name: fom_rom_theta_162
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 162^{o}$
```

```{figure} img/fom_rom_theta_180.png
---
name: fom_rom_theta_180
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 180^{o}$
```

## Comparing FOM and ROM Cylinders

The ROM captures the relative pattern of the temperature distribution but not the scale.  This could be due to the fact that the ROM was constructed from a non-dimensionalized form of the FOM.

```{figure} img/cyl_theta_0.png
---
name: cyl_theta_0
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 0^{o}$
```

```{figure} img/cyl_theta_18.png
---
name: cyl_theta_18
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 18^{o}$
```

```{figure} img/cyl_theta_36.png
---
name: cyl_theta_36
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 36^{o}$
```

```{figure} img/cyl_theta_54.png
---
name: cyl_theta_54
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 54^{o}$
```

```{figure} img/cyl_theta_72.png
---
name: cyl_theta_72
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 72^{o}$
```

```{figure} img/cyl_theta_90.png
---
name: cyl_theta_90
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 90^{o}$
```

```{figure} img/cyl_theta_108.png
---
name: cyl_theta_108
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 108^{o}$
```

```{figure} img/cyl_theta_126.png
---
name: cyl_theta_126
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 126^{o}$
```

```{figure} img/cyl_theta_144.png
---
name: cyl_theta_144
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 144^{o}$
```

```{figure} img/cyl_theta_162.png
---
name: cyl_theta_162
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 162^{o}$
```

```{figure} img/cyl_theta_180.png
---
name: cyl_theta_180
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 180^{o}$
```