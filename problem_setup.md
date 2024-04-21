# Problem Setup

The objective of this project is to build a parametric reduced order model (ROM) for a transient heat conduction problem in a plate with a cylindrical hole.  The geometry of the problem is shown in {numref}`geometry`.

```{figure} img/geometry.png
---
name: geometry
---
Schematic of the transient heat conduction problem
```
The cylinder is centered at $(x_{0}, y_{0}) = (0.5, -0.25)$ and its radius is $r_{0} = 0.25$. The governing equation is given by Equation {eq}`gov_pde`:

```{math}
:label: gov_pde
\frac{\partial T}{\partial t} = \kappa \bigg( \frac{\partial^{2} T}{\partial x^{2}} + \frac{\partial^{2} T}{\partial y^{2}} \bigg) + f(x, y)
```

where $(x, y)$ denote the spatial coordinates, $T$ is temperature, $t$ is time, $f(x,y)$ is the heat source, and $\kappa = 1$ is the thermal conductivity. The above equation is assumed to be non-dimensionalized. The cylinder is insulated and as a result the boundary condition is $\partial T / \partial n = 0$ at the surface of the cylinder, where $n$ is the direction normal to the cylinder surface. The boundary condition at all other edges are $T = 0$. The initial condition is $T = 0$. The heat source is localized and its maximum location is parameterized by Equation {eq}`param_eqn`:

```{math}
:label: param_eqn
f(x, y) = 10^{4}\exp{ \bigg( - \frac{\big( x - x_{c} \big)^{2} + \big( y - y_{c} \big)^{2}}{0.05} \bigg)}
```

where $(x_{c}, y_{c})$ are parameters that specify the center of the heat source. In this project $(x_{c}, y_{c})$ can be anywhere on a circle with the radius $R=0.5$ with the same center as the hole. Therefore, the angle $\theta$ (shown in {numref}`geometry`) parametrizes the location of the heat source as in Equation {eq}`x_param_eqn` and Equation {eq}`y_param_eqn`:

```{math}
:label: x_param_eqn
x_{c} = R\cos{\theta} + x_{0}
```
```{math}
:label: y_param_eqn
y_{c} = R\sin{\theta} + y_{0}
```

The possible location of the heat source is shown by the red line in {numref}`geometry`.