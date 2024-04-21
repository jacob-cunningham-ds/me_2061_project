% Specifies the Initial Condition of the PDE
function  f = IC(x,y)
nr = length(x); % Find the number of nodes in the mesh
f = zeros(1, nr); % Set the temperature at these nodes to 0
end