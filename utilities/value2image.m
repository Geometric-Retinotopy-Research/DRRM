%% Description  --  value2image(F,V, vis, rect, fn)
%		convert values to image, and save image to file path.		
% 
% Parameter(s): 
%     F[double array]   -- connectivity of mesh
%     V[double array]   -- vertex of mesh
%	  vis               -- v value
%     rect              -- location and size infomation
%     fn                -- file name to save
% 
function value2image(F,V, vis, rect, fn)
global N;
minvis = min(vis);

vis = vis -minvis;
maxvis = max(vis);
vis = vis/maxvis;

[xx,yy]=meshgrid(linspace(rect(1), rect(2),N), linspace(rect(3), rect(4),N));
Fi = scatteredInterpolant(V(:,1), V(:,2), vis);
zz = Fi(xx(:),yy(:)); 
cx = (rect(1)+rect(2))/2;
cy = (rect(3)+rect(4))/2;
r = (rect(2)-rect(1))/2;
zz((xx-cx).^2 + (yy-cy).^2 >r^2)=0;
img=reshape(zz, N,N);
imwrite(img,fn);
end

