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