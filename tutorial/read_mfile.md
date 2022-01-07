- read_mfile

  <span style="color: orange">**function [face,vertex,extra] = read_mfile(filename)**</span>

```matlab
% input: file name like A.mat
% output are face, vertex and extra value after calculation

% example:
[F,V,E]=read_mfile('A.mat');
% you can use F(face), V(vertex), E(Extra value) separately.

% Then, use thses parameters in Figure class
param.face = F;
param.vertex = V;
param.EValue = E;
```

