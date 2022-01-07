%% 
% Description  --  function index = vertex_search(uv,vertex)
%		convert values to image, and save image to file path.		
% 
% Parameter(s): 
%     uv[double array]     -- connectivity of mesh
%     vertex[double array]     -- vertex of mesh
% Return:
%     index[int]     -- the index of the search vertex.
%% 
function index = vertex_search(uv,vertex)

% This function searches vertex indexes are nearest to the input XY / XYZ.
%
% Inputs
% uv: (X,Y)/(X,Y,Z) coordinates of ponits in the form of k x 2  
%
% Output :
%	index : k x 1 vertex indexes of the points  

k = size(uv,1); n = size(vertex,1); 
v = vertex'; 
[~,index] = min(sqrt((repmat(v(1,:),k,1)-repmat(uv(:,1),1,n)).^2 +...
    (repmat(v(2,:),k,1)-repmat(uv(:,2),1,n)).^2),[],2);
    

end