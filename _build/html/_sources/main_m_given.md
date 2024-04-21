# main.m

The following MATLAB code was provided for the script `main.m`.

```{code} matlab
:number-lines: true
clear,clc,close all;
set(0,'defaulttextinterpreter','latex')
global model FEM_M FEM_K FEM_F U Kr fr EID;
model = createpde(1);

%% ------ Geometry ---------- %%
R1 = [3,4,-1,1,1,-1,0.5,0.5,-0.75,-0.75]';
C1 = [1,0.5,-0.25,0.25]';
C1 = [C1;zeros(length(R1) - length(C1),1)];
gm = [R1,C1];
sf = 'R1-C1';
ns = char('R1','C1');
ns = ns';
g = decsg(gm,sf,ns);
geometryFromEdges(model,g);
pdegplot(model,'EdgeLabels','on')
axis equal
xlim([-1.1,1.1])

%% ------ Mesh ---------- %%
msh = generateMesh(model,'GeometricOrder','quadratic','Hmax',.05);
% figure(1)
% pdeplot(msh); hold on
% pdegplot(model,'EdgeLabels','on')
set(gca,'Fontsize',15);
xlabel('$x$')
ylabel('$y$')
drawnow
[P,E,~] = meshToPet( msh );
EID = findNodes(msh,'region','Edge',[1:4]);

%% ------ Time discretization parameters ---------- %%
ntime = 201;
startTime = 0;
endTime = .05;
t = linspace(startTime,endTime,ntime);

%% ------ Initial Condition ---------- %%
T0 = IC(P(1,:), P(2,:));

%% ------ Part 1: Data ---------- %%
T = [];
F = [];
Ns = 11;
theta = linspace(180,0,Ns);
% ---- Sample lines for running FOM and ploting temperature fields
% -------------------------------------------------------------------------
n=1;
xc(n) = 0.5*cosd(theta(n))+0.5;
yc(n) = 0.5*sind(theta(n))-0.25;
[model,FEM_M,FEM_K,FEM_F] = GetFEMMatmodel(xc(n),yc(n),model);
[t , Tn] = SolveFOM(T0,t,xc(n),yc(n));
figure; pdeplot(model ,'XYData',Tn(:,end),'ZData',Tn(:,end),'Contour','on','ColorMap','jet'); view([0 0 1]);drawnow
% -------------Edn of sample lines ----------------------------------------
for n=1:Ns
    % complete this loop
end

%% ------ Part 2: POD Modes ---------- %%
C = T'*FEM_M*T;
C = (C+C')/2;
% Complete the rest of part 3 here ...
% To plot the POD modes use the pdeplot as shown above

%% ------ Part 3: ROM ---------- %%
% Build Kr and fr here
% Use SolveROM to solve the ROM
%[t,y] = SolveROM(y0,t); y0 is the inital condition for ROM

%% ------ Part 4: Performance of ROM for unseen data ---------- %%
theta_new= 170;
xc = 0.5*cosd(theta_new)+0.5;
yc = 0.5*sind(theta_new)-0.25;
% update FEM_M
%  Solve FOM for the new parameter --------------
%  Solve ROM for the new parameter --------------
% ------------- Example of using ExtractCylinderT -----------
figure
for i=1:length(t)
    [T_cyl_FOM, theta_cyl,~] = ExtractCylinderT(msh,Tn(:,i));
    plot(theta_cyl,T_cyl_FOM); 
    ylim([0 180])
    drawnow
end

%% ------ Part 5: Field Reconstruction ---------- %%
theta_new= 15;
xc = 0.5*cosd(theta_new)+0.5;
yc = 0.5*sind(theta_new)-0.25;
% update FEM_M
% Solve FOM for the new parameter --------------
% [t,Tn] = SolveFOM(...) 
[~, ~,Id_cyl] = ExtractCylinderT(msh,Tn(:,end));
% Id_cyl: is the ID of the nodes at top surface of the cylinder
T_cyl_FOM = Tn(Id_cyl,end);
% Solve a least square problem to find the coefficent of the POD modes
% using only the surafce measurements 
```