
function [Bay_ecc32k,Bay_ang32k, Bay_sigma32k, Bay_varea32k] = load_benson_inference_cvt_to_32k(subj_nm, lr)
addpath('freesurferReader')
% subj_nm = '100610';
% lr = 'lh'; % hemisphere lh for left and rh for right

% Load Benson's Bayesian Inference
Bay_ecc = squeeze(fs_load_mgh(sprintf('../data/Bayesian_inferred_maps/%s.%s.inf-MSMAll_eccen.native32k.mgz', subj_nm, lr))); 
Bay_ang = squeeze(fs_load_mgh(sprintf('../data/Bayesian_inferred_maps/%s.%s.inf-MSMAll_angle.native32k.mgz',subj_nm, lr))); 
Bay_sigma = squeeze(fs_load_mgh(sprintf('../data/Bayesian_inferred_maps/%s.%s.inf-MSMAll_sigma.native32k.mgz',subj_nm, lr))); 
Bay_varea = squeeze(fs_load_mgh(sprintf('../data/Bayesian_inferred_maps/%s.%s.inf-MSMAll_varea.native32k.mgz',subj_nm, lr))); 
 
if lr(1)=='l'
    Bay_ang  = 90-Bay_ang;
else
    Bay_ang = Bay_ang+90;
    
end

% Load HCP 32k mesh and prf results
fn = sprintf('../data/HCP32Kmat/%s_32K_pRF_fMRI', subj_nm);
hcpsub_32k = load_hcp_subject_32k(fn);

% Load Native Freesurfer Results
Freesurfer_Result = freesurfer_read_subject(sprintf('../data/HCP_FS_Structure/%s/T1w/%s',subj_nm,subj_nm));
BensonData_on_fslr = fsaverage2fslr(Freesurfer_Result.(lr).spherereg.vertices,lr);
 

%% interpretate the inference to 32k mesh by  
uvw_s = BensonData_on_fslr;
uvw_q = hcpsub_32k.(lr).sphere.vertices;
% 
% % convert to xy coordiante to do interpretate and convert back
% xy = unit_cvt([Bay_ecc Bay_ang],'p2c');
% xy32k  =  sphere_interp(uvw_s, xy, uvw_q); 
% ep32k = unit_cvt(xy32k,'c2p'); 
% Bay_ecc32k = ep32k(:,1);
% Bay_ang32k = ep32k(:,2); 

% interepretate to 32k 
Bay_ecc32k = sphere_interp(uvw_s, Bay_ecc, uvw_q);
Bay_ang32k = sphere_interp(uvw_s, Bay_ang, uvw_q);

Bay_sigma32k  =  sphere_interp(uvw_s, Bay_sigma, uvw_q);
Bay_varea32k  = sphere_interp(uvw_s, Bay_varea, uvw_q);



end
