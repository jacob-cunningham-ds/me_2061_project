% Defines the RHS for the ODEs
function Yout = FOM_rhs(t,Y)
global  FEM_K  FEM_F EID
Yout = -(FEM_K*Y) + FEM_F; % Store the linear combination
Yout(EID') = 0; % Apply boundary conditions
end