# Data

Generate a data matrix by solving the full order model (FOM) for $N_{s} = 11$ samples for the following choices of heat sources: $\theta = $`linspace(180, 0, Ns)`. Store the resulting data in a matrix named `T`. The size of this matrix must be $N \times (N_{t}N_{s})$, where $n$ is the number of FEM nodes. Plot the last snapshot of the temperature for each of the cases.

## Code

The following MATLAB code was used to generate the matrix `T` and plots the last snaphshot of the temperature for each of the cases as shown in {numref}`fom_theta_0` through {numref}`fom_theta_180`.

```{code} matlab
:number-lines: true

%% ------ Time discretization parameters ---------- %%
ntime = 201;
startTime = 0;
endTime = .05;
t = linspace(startTime,endTime,ntime);

%% ------ Initial Condition ---------- %%
T0 = IC(P(1,:), P(2,:));

%% ------ Part 1: Data ---------- %%

% Intialize arrays
T = [];
F = [];

% Define samples
Ns = 11;
theta = linspace(180,0,Ns);

% Determine the number of nodes and initalize Temperature matrix
N = size(P, 2);
T = zeros(N, ntime * Ns);

for n=1:Ns
    % Calculate the coordinates of the heat source based on theta
    xc(n) = 0.5 * cosd(theta(n)) + 0.5;  % x-coordinate of heat source
    yc(n) = 0.5 * sind(theta(n)) - 0.25;  % y-coordinate of heat source
    
    % Update the PDE model with new coefficients reflecting the current heat source location
    [model, FEM_M, FEM_K, FEM_F] = GetFEMMatmodel(xc(n), yc(n), model);
    
    % Solve the FOM
    [t, Tn] = SolveFOM(T0, t, xc(n), yc(n));
    
    % Store the results in T for each time step and each theta
    T(:, (n-1)*ntime + 1 : n*ntime) = Tn;
    
    % Plot the last temperature snapshot for the current value of theta
    figure;
    pdeplot(model, 'XYData', Tn(:, end), 'ZData', Tn(:, end), 'Contour', 'on', 'ColorMap', 'jet');
    title(sprintf('Temperature Distribution at $\\theta = %g^\\circ$', theta(n)));
    xlabel('$x$');
    ylabel('$y$');
    view([0 0 1])
    colorbar;
    axis equal;
    drawnow;
end
```

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