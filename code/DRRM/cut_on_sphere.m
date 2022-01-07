%% Description  -- function sf_subj_cut = cut_on_sphere(subj)
% 
%		Cut the subject by the boundary defined on the fslr sphere
% Parameter(s): 
%		subj[struct]  --  subject needs to cut.
%
% Return: 
%		sf_subj_cut[struct]  -- subject with cut handle, including left and
%		right himesphere
%
function sf_subj_cut = cut_on_sphere(subj)

load('HCPTemplate.mat','template');
lrs = ['lh'; 'rh'];
for i=1:2
    lr = lrs(i,:);
    
    % Process 
    hemi = subj.(lr); 
 
    Vs = fsaverage2fslr(double(hemi.sphere.vertices),lr); % fslr sphere
   
    % stereographic project
    xs = Vs(:,1)./(100-Vs(:,3));
    ys = Vs(:,2)./(100-Vs(:,3));
    
    hemi_template = template.(lr);
    in=inpolygon(xs, ys, hemi_template.bd_suv(:,1),hemi_template.bd_suv(:,2));
    
    [Fcut,~,father] = gf_remove_mesh_vertices(hemi.pial.faces,hemi.pial.vertices, ~in); % get a list of vertex id to delete
    
    

    sf_subj_cut.(lr).father = father;
    sf_subj_cut.(lr).pial = hemi.pial.vertices(father,:);
    sf_subj_cut.(lr).pRF = hemi.pRF(father,:);
    sxsy2u=scatteredInterpolant(hemi_template.suv,hemi_template.uv(:,1));
    sxsy2v=scatteredInterpolant(hemi_template.suv,hemi_template.uv(:,2));
    sf_subj_cut.(lr).uv =  [sxsy2u(xs(father),ys(father)) sxsy2v(xs(father),ys(father))];
    
    Fcut =  delaunay(sf_subj_cut.(lr).uv);
    
    sf_subj_cut.(lr).F = Fcut;
    
    if iscell(hemi.fMRI) 
        for j=1:length(hemi.fMRI)
            sf_subj_cut.(lr).fMRI{j} = hemi.fMRI{j}(father,:);
        end
    else
         sf_subj_cut.(lr).fMRI = hemi.fMRI(father,:);
    end
    
    bd = compute_bd(Fcut);
    bdnorm = sqrt(sf_subj_cut.(lr).uv(bd,1).^2 + sf_subj_cut.(lr).uv(bd,2).^2);
    sf_subj_cut.(lr).uv(bd,:) = sf_subj_cut.(lr).uv(bd,:)./bdnorm;
    
    try
        sf_subj_cut.(lr).atlas_hcp = hemi.atlas_hcp(father,:);
        sf_subj_cut.(lr).atlas_wang	 = hemi.atlas_wang(father,:);
    catch
        
    end
    
end

end