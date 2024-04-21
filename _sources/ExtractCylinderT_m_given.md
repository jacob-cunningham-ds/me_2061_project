# ExtractCylinderT.m 

The following MATLAB code was provided for the function `ExtractCylinderT.m`.

```{code} matlab
:number-lines: true
function [T_cyl,theta_cyl,EID_cyl] = ExtractCylinderT(msh,T)
EID_cyl = findNodes(msh,'region','Edge',[8 7]);
[P,E,~] = meshToPet( msh );
x0 = 0.5;
r_cyl = 0.25;
x_cyl = P(1,EID_cyl);
theta_cyl = acosd((x_cyl-x0)/r_cyl);
[theta_cyl,I] = sort(theta_cyl);
T_cyl = T(EID_cyl(I),:);
```