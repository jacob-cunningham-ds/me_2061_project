# FOM_rhs.m 

The following MATLAB code was provided for the function `f.m`.

```{code} matlab
:number-lines: true
function Yout = FOM_rhs(t,Y)
global  FEM_K  FEM_F EID
Yout = -(FEM_K*Y) + FEM_F;
Yout(EID') = 0;
end
```