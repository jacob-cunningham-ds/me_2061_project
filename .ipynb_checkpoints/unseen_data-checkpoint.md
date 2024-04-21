# Performance of ROM for Unseen Data

Use the ROM to solve for a new heat source parameter at $\theta = 170$. Note that a heat source at this angle is not included in the training samples. You need to update `fr` according to the new heat source distribution. To this end, you can call the function `[model, FEM_M, FEM_K, FEM_F] = GetFEMMatmodel(xc, yc, model)` to update `FEM_F` by passing the new `xc` and `yc`.

Compare the FOM and ROM temperature contours for the last snapshot. Also, compare the FOM and ROM temperature profiles on the top half of the cylinder at times $t = 0.01, 0.03, 0.05$.

MATLAB code was used to solve for a new heat source at $\theta = 170$, compare the FOM and ROM tempearature contours for the last snapshot as shown in Figure {numref}`fom_rom_theta_170`, and compare the FOM and ROM temperature profiles on the top half of the cylinder at times $t = 0.01, 0.03, 0.05$ as shown in Figure {numref}`cyl_theta_170`.

## Comparing FOM and ROM Last Snapshots

```{figure} img/fom_rom_theta_170.png
---
name: fom_rom_theta_170
---
Comparison of the FOM solution and ROM solution for the last time snapshot at $\theta = 170^{o}$
```

## Comparing FOM and ROM Cylinders

```{figure} img/cyl_theta_170.png
---
name: cyl_theta_170
---
Comparison of the FOM solution and ROM solution for the cylinder temperature profile at $\theta = 170^{o}$
```