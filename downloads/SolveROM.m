% Solves the ROM
function [t,y]=SolveROM(y0,t)
odesolve_tol = 1e-6; % Set the tolerance level

% Configure the options for the ODE solver
options = odeset('AbsTol',odesolve_tol ,'RelTol',odesolve_tol ,'Stats','on');
disp('Solving ROM');
tic
[t , y] = ode45(@ROM_rhs,t,y0,options ); % Solve the ODEs
toc