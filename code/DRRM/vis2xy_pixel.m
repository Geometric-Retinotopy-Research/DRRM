%% 
% Description  --  function index = vertex_search(uv,vertex)
%		recovery from hat_vis to xy in pixel		
% 
% Parameter(s): 
%     hat_vis[double array]     -- vis data
% Return:
%     x[double]     -- x value after recovery
%     y[double]     -- y value after recovery
%% 
function [x,y] = vis2xy_pixel(hat_vis)

% from hat_vis to new x, y in visual image pixel
ecc = hat_vis(:,1)/8*100;
anginrad = hat_vis(:,2)/180*pi;  

x = ecc.*cos(anginrad) + 100.5;
y = 100.5 - ecc.*sin(anginrad);

end

