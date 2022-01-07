%%		
% Description  -- function fn = laplacian_smooth(uv,f, W)
%		make a grid laplacian smooth
%
% Parameter(s): 
%		uv[Nx2 double array]    --  uv coordinate in plane
%		f[Nxd double array]     --  feature vector in the plane 
%		W[Nx1 double array]     --  
% return: 
%		fn[Nxd double]          --  the result of grid smooth
%		s[1x1 double]           --  the optimal smooth scale
%		
%%
function fn = laplacian_smooth(uv,f, W)
if ~exist('W','var')
    W = ones(size(f,1),1);
end

Mres = 200; Nres = 200;
x = {linspace(-1,1,Mres),linspace(-1,1,Nres)};
[xx,yy] = ndgrid(x{1},x{2});
grid_vec = [xx(:) yy(:)];
id0 = xx.^2+yy.^2 >1; % take r>1 points



Fw = scatteredInterpolant(uv(:,1), uv(:,2), W);
W = reshape(Fw(grid_vec),[Mres,Nres]);
W(id0)=0;
W(W<0)=0;
    
for i = 1:size(f,2)
    % Create   data on grid, ecc,R2, ang
    feature = f(:,i);
    Fr = scatteredInterpolant(uv(:,1), uv(:,2), feature);
    f_on_grid = reshape(Fr(grid_vec),[Mres,Nres]);
%     f_on_grid(id0)=0;
    
    %  Do smoothing
    smooth_f_on_grid = smoothn(f_on_grid,W);
    fn(:,i) = griddata(xx, yy, smooth_f_on_grid, uv(:,1), uv(:,2));
end

end