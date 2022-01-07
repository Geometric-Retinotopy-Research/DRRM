classdef GenDataManager
    methods(Static)
%%
% Description  -- result = getTemplateTxt()  
%		load the txt content from .txt template file to result.
%
% return: 
%      result[double array]  -- file content
%%
       function result = getTemplateTxt()
            if ~exist([GenConsts.kTXTURL, 'template_landmark.txt'],'file')
                xy = ginput(21);
                save([GenConsts.kTXTURL, 'template_landmark.txt'], '-ascii'); 
            else
                xy = load([GenConsts.kTXTURL, 'template_landmark.txt']);
            end
            result = xy;
	   end

%%
% Description  -- function result = loadSubject(filePath)
%		load the subject content from HCP 32K data file to result.
%
% return: 
%      result[double array]  -- file content
%%
	   function result = loadSubject(filePath)
		   subject =load_hcp_subject_32k(filePath);
		   % set the fovea index manually 
		   subject.lh.fovid = 23686;
		   subject.rh.fovid = 23778;
		   result = subject;
	   end
%% Description  -- function [face,vertex,extra] = read_mfile(filename)
%     load benson's template and handle data in order to use in the current
%     project.
%
%% parameter(s): 
%		lr [String]  -- the type of the current template need
%		hemi_cut[strcut] -- hemisphere cut 	
%
%% return: 
%      template[struct]  -- template including benson's template data.
%      b_template[struct]  -- original benson's template struct
% 
%%   
	   function [template, b_template] = loadBensonAtlasTempelate(lr, hemi_cut)
		   benson_template = load_benson_atlas(lr); 
		   b_template = benson_template;
		   Vlr = fsaverage2fslr(benson_template.V, lr);
		   uv_stereo = [Vlr(:,1)./(100-Vlr(:,3)) Vlr(:,2)./(100-Vlr(:,3)) ];
		   in = inpolygon(uv_stereo(:,1), uv_stereo(:,2), hemi_cut.bd_suv(:,1), hemi_cut.bd_suv(:,2));

		   % fit a map from uv_stereo to
		   uv_benson = zeros(length(find(in)),2);
		   for j=1:2
			   Fi{j} = scatteredInterpolant(double(hemi_cut.suv), double(hemi_cut.uv(:,j)));
			   uv_benson(:,j) = Fi{j}(uv_stereo(in,:));
		   end

			[Ftemplate, ~,father_benson]= gf_remove_mesh_vertices(benson_template.F, benson_template.V, ~in);

			template.F = Ftemplate;
			template.V = Vlr(father_benson,:);
			template.uv = uv_benson;
			template.suv = uv_stereo(father_benson,:);
			template.bd_suv = template.suv(compute_bd(Ftemplate),:);
			template.bd_suv_description = 'bd_suv is the boundary in space of stereographic uv form fs_LR ';
			template.suv2uv = @(suv) [Fi{1}(suv) Fi{2}(suv)];
			template.ecc = benson_template.ecc(father_benson);
			template.ang = benson_template.ang(father_benson);
			template.sigma = benson_template.sigma(father_benson);
			template.varea = benson_template.varea(father_benson);
	   end
	   
    end
end