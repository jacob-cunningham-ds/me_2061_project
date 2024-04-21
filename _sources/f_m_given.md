# f.m 

The following MATLAB code was provided for the function `f.m`.

```{code} matlab
:number-lines: true
function f(location,state)
global xc yc
f = 10000*exp(-((location.x - xc).^2+(location.y-yc).^2)/0.05);
```