%% Description  -- function subjectout = load_studyforrest_subject_mat(sid)
%		load study studyforrest subject file with id 
%
% Parameter(s): 
%	     sid [int]  --  id with file name, like 'sub-01.mat'  
%
% Return: 
%      subjectout[struct]  -- suject struture, including information of
%      sunject file
%%
function subjectout = load_studyforrest_subject_mat(sid)
fn = sprintf('../data/Studyforrest/sub-%02d.mat',sid);
submat = load(fn); 
subjectout = submat.subject;
mask = get_59kmesh_mask(submat.subject);
N = 59292;
fns = fieldnames(submat.results);
for i=1:length(fns)
    fn = fns{i};
    dn = submat.results.(fn);
    if iscell(dn) || isempty(dn) ||isstruct(dn)
        continue
    end
    if length(size(dn))>2
        dn = squeeze(dn)';
    end
    datacols = size(dn,2);
    feature =  zeros(N*2,datacols);
    feature(find(mask),:) = dn;
    
    subjectout.lh.(fn) = feature(1:N,:);
    subjectout.rh.(fn) = feature(N+1:N*2,:);
    
end

for lr = {'lh','rh'}
    lr = lr{1};
    ret_scan_names = {'con','exp','clw', 'ccw'};
    fmricat = [];
    for i=1:length(ret_scan_names)
        fn = ['fmri_' ret_scan_names{i}];
        fmricat = [fmricat subjectout.(lr).(fn)];
        subjectout.(lr).(fn) = [];
    end
    subjectout.(lr).fMRI  = fmricat;
    
    
    subjectout.(lr).pRF = [subjectout.(lr).ang subjectout.(lr).ecc, ...
                            subjectout.(lr).gain subjectout.(lr).meanvol, ...
                            subjectout.(lr).R2, subjectout.(lr).rfsize];
                        
                        
    mesh_names = {'sphere','white','pial', 'inflated'};                     
    for i=1:length(mesh_names)
        fn = mesh_names{i};         
        meshstr.faces = subjectout.(lr).faces;
        meshstr.vertices = subjectout.(lr).(fn);
        subjectout.(lr).(fn) = meshstr;
    end
    
end
 
 
end