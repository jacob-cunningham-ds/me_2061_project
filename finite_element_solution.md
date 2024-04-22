# Finite Element Solution

For the purpose of this project, the equations were solved using finite element method (FEM) using MATLAB PED solver toolbox. The FEM mesh is genererated in MATLAB and it is stored in the variable `msh`. 

You can control the size of the mesh by changing the variable `Hmax`. You need to report all of your results with `Hmax=0.05`. However, for debugging purposes you may use a larger value for `Hmax`. 

The function `SolveFOM` solves the finite element problem and returns the solution snapshots in the interval of $0 \leq t \leq T_{f}$ at $N_{t} = 201$ uniformally distributed time snapshots, where $T_{f} = 0.05$.

## Geometry and Mesh Generation

MATLAB code was used to generate the geometry and mesh for this project shown in {numref}`geometry_build` and {numref}`mesh`.

```{figure} img/geometry_build.png
---
name: geometry_build
---
Resultant geometry
```

```{figure} img/mesh.png
---
name: mesh
---
Resultant mesh
```