%%
% Description  -- function hcpsub = load_hcp_subject_32k(fn)
%		load the subject content from HCP 32K data file to result.
%
%%
function hcpsub = load_hcp_subject_32k(fn)
res = load(fn);
hcpsub = res.mesh_sub;
 