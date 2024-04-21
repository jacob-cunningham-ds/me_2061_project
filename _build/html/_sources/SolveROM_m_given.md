# SolveROM.m 

The following MATLAB code was provided for the function `SolveROM.m`.

```{code} matlab
:number-lines: true
function [t,y]=SolveROM(y0,t)
odesolve_tol = 1e-6;
options = odeset('AbsTol',odesolve_tol ,'RelTol',odesolve_tol ,'Stats','on');
disp('Solving ROM');
tic
[t , y] = ode45(@ROM_rhs,t,y0,options );
toc
```