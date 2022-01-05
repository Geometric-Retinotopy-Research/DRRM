%%
% Description  -- function hcpsub = load_hcp_subject_32k(fn)
%		load the subject content from HCP 32K data file to result.
%
%%
function hcpsub = load_hcp_subject_32k(fn)
res = load(fn);
hcpsub = res.mesh_sub;

% 
% ang = cos(hcpsub.lh.pRF(:,1)/180*pi);
% ecc = hcpsub.lh.pRF(:,1);
% hcpsub.lh.vis = [ecc.*cos(ang) ecc.*sin(ang) hcpsub.lh.pRF(:,3:end)];
% 
% 
% ang = cos(hcpsub.rh.pRF(:,1)/180*pi);
% ecc = hcpsub.rh.pRF(:,1);
% hcpsub.rh.vis = [ecc.*cos(ang) ecc.*sin(ang) hcpsub.rh.pRF(:,3:end)];