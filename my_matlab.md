# Source Code

The following MATLAB code was used for this project.  The MATLAB files, and associated files to build this Jupyter Book, can be retreived from this [GitHub repository](https://github.com/jacob-cunningham-ds/me_2061_project).

## main.m

```{code} matlab
:number-lines: true
% Clear the workspace
clear,clc,close all;

% Set the interpreter
set(0,'defaulttextinterpreter','latex')

% Declare global variables
global model FEM_M FEM_K FEM_F U Kr fr EID;

% Create a pde model container
model = createpde(1);

%% ------ Geometry ---------- %%

% Define the plate dimensions
R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';

% Define the inner circle dimensions
C1 = [1,0.5,-0.25,0.25]';

% Pad with zeros to match the length of R1
C1 = [C1;zeros(length(R1) - length(C1),1)];

% Concatenate the geometry
gm = [R1,C1];

% Subtract the cylinder set from the rectangle set
sf = 'R1-C1';

% Capture the names
ns = char('R1','C1');
ns = ns';

% Construct the geometry
g = decsg(gm,sf,ns);

% Import the geometry into the PDE model
geometryFromEdges(model,g);

% Plot the geometry
pdegplot(model,'EdgeLabels','on')
axis equal
xlim([-1.1,1.1])

% Generate the mesh
msh = generateMesh(model,'GeometricOrder','quadratic','Hmax',.05);

% Open a new figure
figure

% Plot the mesh
pdeplot(msh); hold on
pdegplot(model,'EdgeLabels','on')
set(gca,'Fontsize',15);
xlabel('$x$')
ylabel('$y$')
drawnow

% Convert the mesh to PET format
[P,E,~] = meshToPet( msh );

% Identify the nodes on the edges
EID = findNodes(msh,'region','Edge',[1:4]);

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
    drawnow;
end

%% ------ Part 2: POD Modes ---------- %%
C = T'*FEM_M*T;
C = (C+C')/2;

% Compute SVD
[U, S, V] = svd(T, 'econ');

% Loop through the first four modes and plot each
for i = 1:4

    % Creat a new figure
    figure;

    % Reshape the ith column of U back into the spatial domain dimensions
    mode_shape = reshape(U(:, i), N, []);
    
    % Use pdeplot to plot
    pdeplot(model, 'XYData', mode_shape(:, end), 'Contour', 'on', 'ColorMap', 'jet');
    title(sprintf('POD Mode %d', i));
    xlabel('$x$');
    ylabel('$y$');
    colorbar;
    drawnow;
end

% Extract the diagonal elements of S
singular_values = diag(S);

% Create a semilogy plot of the first 50 singular values
figure;
semilogy(1:50, singular_values(1:50), 'o-');
title('First 50 Singular Values');
xlabel('Index');
ylabel('Singular Value');
grid on;

% Square the singular values to get the energy of each mode
energy_per_mode = singular_values.^2;

% Compute the cumulative sum of energies
cumulative_energy = cumsum(energy_per_mode);

% Normalize by the total sum to get the fraction or percentage of total energy
cumulative_energy_percentage = cumulative_energy / sum(energy_per_mode);

% Find how many modes are needed to capture at least 99% of the energy
modes_needed_99_percent = find(cumulative_energy_percentage >= 0.99, 1);

% Display the number of modes required
fprintf('Number of modes required to capture 99%% of the energy: %d\n', modes_needed_99_percent);

% Plot the cumulative energy as a percentage
figure;
plot(cumulative_energy_percentage, 'o-');
title('Cumulative Energy Captured by Singular Values');
xlabel('Number of Modes');
ylabel('Cumulative Energy Percentage');
grid on;
axis tight;

%% ------ Part 3: ROM ---------- %%

% Reduced basis containing the dominant modes
Ur = U(:, 1:modes_needed_99_percent);

% Initialize storage for ROM solutions
Y_full = zeros(N, ntime * Ns);

% Loop through all theta samples
for n = 1:Ns
    % Calculate the coordinates of the heat source based on theta
    xc = 0.5 * cosd(theta(n)) + 0.5;  % x-coordinate of heat source
    yc = 0.5 * sind(theta(n)) - 0.25;  % y-coordinate of heat source

    % Update the PDE model coefficients reflecting the current heat source location
    [model, FEM_M, FEM_K, FEM_F] = GetFEMMatmodel(xc, yc, model);

    % Reduce FEM_K and FEM_F based on the current reduced basis Ur
    Kr = Ur' * FEM_K * Ur;
    fr = Ur' * FEM_F;

    % Project initial conditions into the reduced space
    y0 = Ur' * T0';

    % Solve the ROM
    [t, y] = SolveROM(y0, t);

    % Map reduced state back to the full state for plotting
    Y_full(:, (n-1)*ntime + 1:n*ntime) = Ur * y';
end

%% ------ Compare ROM and FOM ---------- %%

% Number of theta samples
Ns = length(theta);

% Set up a counter for subplot indexing
subplotCounter = 1;

% Loop through all sampled values of theta
for n = 1:Ns
    
    % Create a new figure
    figure;

    % Calculate the index for the solution corresponding to the current theta
    currentIndex = (n-1)*ntime + 1 : n*ntime;
    
    % Reshape the last snapshot from FOM for spatial dimensions
    FOM_current_snapshot = reshape(T(:, currentIndex(end)), N, []);
    
    % Reshape the last snapshot from ROM for spatial dimensions
    ROM_current_snapshot = reshape(Y_full(:, currentIndex(end)), N, []);
    
    % Subplot for the FOM solution
    subplot(2, 1, 1);
    pdeplot(model, 'XYData', FOM_current_snapshot(:, end), 'Contour', 'on', 'ColorMap', 'jet');
    title(sprintf('FOM at $\\theta = %d^{\\circ}$', theta(n)));
    xlabel('$x$');
    ylabel('$y$');
    colorbar;
    drawnow;
    
    % Subplot for the ROM solution
    subplot(2, 1, 2);
    pdeplot(model, 'XYData', ROM_current_snapshot(:, end), 'Contour', 'on', 'ColorMap', 'jet');
    title(sprintf('ROM at $\\theta = %d^{\\circ}$', theta(n)));
    xlabel('$x$');
    ylabel('$y$');
    colorbar;
    drawnow;

    % Add a general title to the figure
    sgtitle('Comparison of FOM and ROM Solutions');
end
%% ------ ExtractCylinderT ---------- %%

% Time points for comparison
time_points = [0.01, 0.03, 0.05];

% Loop through all sampled values of theta
for theta_idx = 1:length(theta)
    
    % Create a new figure for each theta
    figure;
    
    % Loop through the specified time points for comparison
    for time_idx = 1:length(time_points)
        % Find the index of the current time point in the t array
        index_in_t = arrayfun(@(tp) find(min(abs(t-tp))==abs(t-tp),1), time_points);

        % Extract FOM temperature on the cylinder at the current time index
        [T_cyl_FOM, theta_cyl_FOM, ~] = ExtractCylinderT(msh, T(:, (theta_idx-1)*ntime + index_in_t));
        valid_idx_FOM = theta_cyl_FOM <= 180;  % Ensure theta is in the range [0, 180]

        % Extract ROM temperature on the cylinder at the current time index
        [T_cyl_ROM, theta_cyl_ROM, ~] = ExtractCylinderT(msh, Y_full(:, (theta_idx-1)*ntime + index_in_t));
        valid_idx_ROM = theta_cyl_ROM <= 180;  % Ensure theta is in the range [0, 180]

        % Subplot for the temperature profile at the current time point
        subplot(length(time_points), 2, 2*(time_idx-1)+1);
        
        % Plot FOM temperature profile
        plot(theta_cyl_FOM(valid_idx_FOM), T_cyl_FOM(valid_idx_FOM), 'b-', 'LineWidth', 2);

        % Annotate the plot
        title(sprintf('FOM Temperature at $t = %.2fs$', time_points(time_idx)));
        xlabel('$\theta$ (degrees)');
        ylabel('Temperature');
        xlim([0, 180]);
        grid on;

        % Subplot for the ROM temperature profile at the current time point
        subplot(length(time_points), 2, 2*(time_idx-1)+2);

        % Plot ROM temperature profile
        plot(theta_cyl_ROM(valid_idx_ROM), T_cyl_ROM(valid_idx_ROM), 'r--', 'LineWidth', 2);
        
        % Annotate the plot
        title(sprintf('ROM Temperature at $t = %.2fs$', time_points(time_idx)));
        xlabel('$\theta$ (degrees)');
        ylabel('Temperature');
        xlim([0, 180]);
        grid on;

    end

    % Add a general title to the figure
    sgtitle(sprintf('Temperature on the Top Half of the Cylinder Surface at $\\theta = %d^{\\circ}$', theta(theta_idx)));
end

%% ------ Part 4: Performance of ROM for unseen data ---------- %%

% New heat source angle
theta_new = 170;

% Calculate new heat source coordinates
xc_new = 0.5 * cosd(theta_new) + 0.5;
yc_new = 0.5 * sind(theta_new) - 0.25;

% Update FEM_F for the new heat source distribution
[model, ~, ~, FEM_F_new] = GetFEMMatmodel(xc_new, yc_new, model);
FEM_F = FEM_F_new;

% Solve the FOM with the new heat source force vector
[~, ~] = SolveFOM(T0, t, xc_new, yc_new);

% Reduce the force vector for the new heat source
fr_new = Ur' * FEM_F_new;
fr = fr_new;

% Solve ROM with the new reduced force vector
[~, y_ROM_new] = SolveROM(y0, t);

% Map reduced state back to the full state for the new parameter
Y_full_new = Ur * y_ROM_new';

%% ------ Plot FOM and ROM temperature contours for the last snapshot ------ %%

% Solve the FOM for the new parameter
[~, Tn_new] = SolveFOM(T0, t, xc_new, yc_new);

% Extract the last snapshot for both FOM and ROM
FOM_last_new = reshape(Tn_new(:, end), N, []);
ROM_last_new = reshape(Y_full_new(:, end), N, []);

% Create figure for comparison
figure;

% Subplot for the FOM solution
subplot(2, 1, 1);
pdeplot(model, 'XYData', FOM_last_new, 'Contour', 'on', 'ColorMap', 'jet');
title(sprintf('FOM Temperature Distribution at $\\theta = %d^{\\circ}$', theta_new));
xlabel('$x$');
ylabel('$y$');
colorbar;
drawnow;

% Subplot for the ROM solution
subplot(2, 1, 2);
pdeplot(model, 'XYData', ROM_last_new, 'Contour', 'on', 'ColorMap', 'jet');
title(sprintf('ROM Temperature Distribution at $\\theta = %d^{\\circ}$', theta_new));
xlabel('$x$');
ylabel('$y$');
colorbar;
drawnow;

% Add a figure title
sgtitle('Comparison of FOM and ROM Temperature Distributions for Unseen Heat Source');

%% ------ Plot FOM and ROM cylinder temperature profiles ------ %%
% Create a figure to plot the temperature profiles
figure;

% Loop over the selected time points
for i = 1:length(time_points)
    % Find the index of the current time point in the t array
    index_in_t = arrayfun(@(tp) find(min(abs(t-tp))==abs(t-tp),1), time_points);

    % Extract FOM temperature on the cylinder at the current time index
    [T_cyl_FOM, theta_cyl_FOM, ~] = ExtractCylinderT(msh, Tn_new(:, index_in_t));

    % Extract ROM temperature on the cylinder at the current time index
    [T_cyl_ROM, theta_cyl_ROM, ~] = ExtractCylinderT(msh, Y_full_new(:, index_in_t));

    % Ensure theta is in the range [0, 180] for the top half of the cylinder
    valid_idx_FOM = theta_cyl_FOM <= 180;
    valid_idx_ROM = theta_cyl_ROM <= 180;
    
    % Subplot for the FOM temperature profile at the current time point
    subplot(length(time_points), 2, 2*i-1); % Subplots on the left for FOM
    
    % Plot FOM temperature profile
    plot(theta_cyl_FOM(valid_idx_FOM), T_cyl_FOM(valid_idx_FOM), 'b-', 'LineWidth', 2);
    
    % Annotate the plot
    title(sprintf('FOM Temperature at $t = %.2fs$', time_points(i)), 'Interpreter', 'latex');
    xlabel('$\theta$ (degrees)', 'Interpreter', 'latex');
    ylabel('Temperature', 'Interpreter', 'latex');
    xlim([0, 180]);
    grid on;
    
    % Subplot for the ROM temperature profile at the current time point
    subplot(length(time_points), 2, 2*i);
    
    % Plot ROM temperature profile
    plot(theta_cyl_ROM(valid_idx_ROM), T_cyl_ROM(valid_idx_ROM), 'r--', 'LineWidth', 2);
    
    % Annotate the plot
    title(sprintf('ROM Temperature at $t = %.2fs$', time_points(i)), 'Interpreter', 'latex');
    xlabel('$\theta$ (degrees)', 'Interpreter', 'latex');
    ylabel('Temperature', 'Interpreter', 'latex');
    xlim([0, 180]);
    grid on;
end

% Add a general title to the figure
sgtitle(sprintf('Temperature on the Top Half of the Cylinder Surface at $\\theta = %d^{\\circ}$', theta_new));

%% ------ Part 5: Field Reconstruction ---------- %%
% New heat source coordinates
theta_new = 15;
xc_new = 0.5 * cosd(theta_new) + 0.5;
yc_new = 0.5 * sind(theta_new) - 0.25;

% Update FEM_F for the new heat source distribution
[model, ~, ~, FEM_F_new] = GetFEMMatmodel(xc_new, yc_new, model);
FEM_F = FEM_F_new;

% Solve FOM for the new parameter
[t_FOM_new, Tn_new] = SolveFOM(T0, t, xc_new, yc_new);

% Use the function ExtractCylinderT to extract the temperature
[~, ~, Id_cyl] = ExtractCylinderT(msh, Tn_new(:, end));

% Get the temperature at the cylinder's top surface at the final time
T_cyl_FOM = Tn_new(Id_cyl, end);

% POD modes on cylinder surface
Ur_surface = Ur(Id_cyl, :);

% Solve a least squares problem to find the coefficients of the POD modes
coefficients = Ur_surface \ T_cyl_FOM;

% Reconstruct the temperature field
T_reconstructed = Ur * coefficients;

% Reshape T_reconstructed for plotting
T_reconstructed_reshaped = reshape(T_reconstructed, size(Tn_new, 1), []);

% Calculate the error contours
error_field = Tn_new(:, end) - T_reconstructed_reshaped(:, end);

% Create figure
figure;

% Subplot for the true temperature field
subplot(3, 1, 1)

% Plot the true temperature field
pdeplot(model, 'XYData', Tn_new(:, end), 'Contour', 'on', 'ColorMap', 'jet');
title('True Temperature Distribution');
xlabel('x');
ylabel('y');
colorbar;
drawnow;

% Subplot for the reconstructed temperature field
subplot(3, 1, 2)

% Plot the reconstructed temperature field
pdeplot(model, 'XYData', T_reconstructed_reshaped(:, end), 'Contour', 'on', 'ColorMap', 'jet');
title('Reconstructed Temperature Distribution');
xlabel('x');
ylabel('y');
colorbar;
drawnow;

% Subplot for the reconstructed temperature field
subplot(3, 1, 3)

% Plot contours of the error
pdeplot(model, 'XYData', error_field, 'Contour', 'on', 'ColorMap', 'jet');
title('Error Contours');
xlabel('x');
ylabel('y');
colorbar;
drawnow;

% Add a general title to the figure
sgtitle(sprintf('Comparison of the Temperature Fields at $\\theta = %d^{\\circ}$', theta_new));
```

## IC.m

```{code} matlab
:number-lines: true
% Specifies the Initial Condition of the PDE
function  f = IC(x,y)
nr = length(x); % Find the number of nodes in the mesh
f = zeros(1, nr); % Set the temperature at these nodes to 0
end
```

## GetFEMMatmodel.m

```{code} matlab
:number-lines: true
% Get the FEM matrices for the PDE model
function [model,FEM_M,FEM_K,FEM_F] = GetFEMMatmodel(xc,yc,model)
global EID

% Heat source function
f = @(location,state)10000*exp(-((location.x - xc).^2+(location.y-yc).^2)/0.05);

% Coefficients for the PDE
specifyCoefficients(model,'m',0,'d',1,'c',1,'a',0,'f',f);

% Get the FEM matrices from the model
model_FEM_matrices = assembleFEMatrices(model);
FEM_M = model_FEM_matrices.M; % Mass matrix
FEM_K = model_FEM_matrices.K; % Stiffness matrix
FEM_F = model_FEM_matrices.F; % Load vector

% Enforce boundary conditions
FEM_M(EID,:)= 0 ;
FEM_M(:,EID)= 0 ;
FEM_M(EID,EID)= eye(length(EID));
```

## SolveFOM.m

```{code} matlab
:number-lines: true
% Solves FOM
function [t,T]=SolveFOM(u0,t,xc,yc)
global model FEM_M FEM_K FEM_F

odesolve_tol = 1e-6; % Set the tolerance level

% Configure the options for the ODE solver
options = odeset('Mass',FEM_M ,'AbsTol',odesolve_tol ,'RelTol',odesolve_tol ,'Stats','on');
disp('Solving FOM');
tic
[t , T] = ode45(@FOM_rhs,t,u0',options ); % Solve the ODEs
toc
T = T'; % Tranpose so each column is a solution at a point in time 
```

## FOM_rhs.m

```{code} matlab
:number-lines: true
% Defines the RHS for the ODEs
function Yout = FOM_rhs(t,Y)
global  FEM_K  FEM_F EID
Yout = -(FEM_K*Y) + FEM_F; % Store the linear combination
Yout(EID') = 0; % Apply boundary conditions
end
```

## SolveROM.m

```{code} matlab
:number-lines: true
% Solves the ROM
function [t,y]=SolveROM(y0,t)
odesolve_tol = 1e-6; % Set the tolerance level

% Configure the options for the ODE solver
options = odeset('AbsTol',odesolve_tol ,'RelTol',odesolve_tol ,'Stats','on');
disp('Solving ROM');
tic
[t , y] = ode45(@ROM_rhs,t,y0,options ); % Solve the ODEs
toc
```

## ROM_rhs.m

```{code} matlab
:number-lines: true
% Defines the RHS for the ODEs
function Yout = ROM_rhs(t,Y)
global Kr fr
Yout = -Kr * Y + fr; % Store the linear combination
end
```

## ExtractCylinderT.m

```{code} matlab
:number-lines: true
% Extract temperature data along the cylinder
function [T_cyl,theta_cyl,EID_cyl] = ExtractCylinderT(msh,T)

% Find node ids of the elements that lie on the specified edges
EID_cyl = findNodes(msh,'region','Edge',[8 7]);

% Extract node positions and elements
[P,E,~] = meshToPet( msh );

% Define the center and radius of the cylinder
x0 = 0.5;
r_cyl = 0.25;

% Get the x-coord of the nodes on the cylinder boundary
x_cyl = P(1,EID_cyl);

% Conver to angular position
theta_cyl = acosd((x_cyl-x0)/r_cyl);

% Sort the angles
[theta_cyl,I] = sort(theta_cyl);

% Get the temperature
T_cyl = T(EID_cyl(I),:);
```