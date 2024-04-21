% Defines the RHS for the ODEs
function Yout = ROM_rhs(t,Y)
global Kr fr
Yout = -Kr * Y + fr; % Store the linear combination
end