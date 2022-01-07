% get_59kmesh_mask: find a mask to keep a region of interest, the roi is
% defined in HCPTemplate.mat
% remove vertices with given ids to remove

%% Syntax
% mask = get_59kmesh_mask(subject)

%% Description
%remove the vertex and face by id2del, also get each vertex fatehr
%% Contribution
%  Author : Yanshuai Tu
%  Created: 2019/2/6
%  Revised: 2022/1/5
%
%  Copyright @ Geometry Systems Laboratory
%  School of Computing, Informatics, and Decision Systems Engineering, ASU
%  http://gsl.lab.asu.edu/

function mask = get_59kmesh_mask(subject)

    hcp_subj_cut = find_roi_on_sphere(subject);

    offset = size(subject.lh.pial,1);
    mask = zeros(2*offset,1);
    mask([hcp_subj_cut.lh.father;hcp_subj_cut.rh.father+offset])=1;
end


function hcp_subj_cut = find_roi_on_sphere(sf_subj)
    % Process left/right 
    for lr={'lh','rh'}
        lr =lr{1};
        hemi = sf_subj.(lr);

        Vlr =fsaverage2fslr(hemi.sphere,lr);

        Vs = double(Vlr/100);
        % stereographic project
        xs = Vs(:,1)./(1-Vs(:,3));
        ys = Vs(:,2)./(1-Vs(:,3));

        tdata = load('HCPTemplate.mat');
        hemi_t = tdata.template.(lr);
        in=inpolygon(xs, ys,  hemi_t.bd_suv(:,1), hemi_t.bd_suv(:,2));
        [Fcut,~,father] = gf_remove_mesh_vertices(hemi.faces,hemi.pial, ~in); % get a list of vertex id to delete
        hcp_subj_cut.(lr).father =  father;
    end



end