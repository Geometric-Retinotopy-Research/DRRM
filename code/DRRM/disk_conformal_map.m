%% Description  -- disk_conformal_map(face,vertex)
%     nake a conformal map from a disk plot
%
%% parameter(s): 
%      face[double array]  -- connectivity of mesh
%      vertex[double array]  -- vertex of mesh	
%% return: 
%      uv[double array]  -- coordinate of mesh
% 
%% 
function uv = disk_conformal_map(face,vertex)
uv = disk_harmonic_map(face,vertex);
bd = compute_bd(face);
for i=1:1
    mu = compute_bc(face,uv,vertex);
    [uv,muc] = linear_beltrami_solver(face,uv,mu,bd,uv(bd,:));
    fprintf('max(abs(muc))=%f\n',max(abs(muc)));
end

 