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
