%%  Reproduce Figure 5. 

% The Retinotopic Template (left hemisphere). 
% (a) decoded polar angles of the group-average HCP retinotopic map in the disk domain; 
% (b) decoded eccentricities of the group-average HCP retinotopic map in the disk domain, 
% (c) BW’s retinotopic model in the 2D domain, 
% (d) the final template, 
% (e) the template on the fsLR sphere.

clear;clc;close all;
% generate the atlas 
%% Load / Cut  the average data surface  
sid = '999999';
filePath = sprintf('%s/%s_32k_pRF_fMRI.mat',GenConsts.kDataHCP32URL, sid);
subject = GenDataManager.loadSubject(filePath);
plot_subject(subject);
%% Gen Fig 5ab, and save template combined with left and right hemisphere
[template.lh, benson_template_lh] = genFigure(subject, 'lh', 1);
[template.rh, benson_template_rh] = genFigure(subject, 'rh', 2);
template.description ='The template is defined on the fs_LR sphere.  fslr==hcp sphere   fsaverage == freesurfer';
save([GenConsts.kTemplatesURL,'HCPTemplate'],'template');
% now generate Figure 5c
load([GenConsts.kTemplatesURL,'HCPTemplate'],'template');

%% create uv coordinate from template generated from above
themi = template.lh;  
bd= compute_bd(themi.F);
uv = themi.uv;  
rbd = sqrt(uv(bd,1).^2 + uv(bd,2).^2 );
uv(bd,1) = uv(bd,1) ./ rbd;     
uv(bd,2) = uv(bd,2) ./ rbd;

%% generate figures
% plot ang figure
param5c.face = themi.F;
param5c.vertex = uv;
param5c.EValue = themi;
param5c.renderColor.value = themi.ang;
param5c.renderColor.valueType = 'lh';
param5c.renderColor.ecc = 0.7;
GenFigure5c = GenFigureAtlas(param5c);
GenFigure5c.draw().save('Fig5_4');

% plot 3D figure with overlappig
%%
figure;
color = benson_template_rh.V*0 +0.8;
plot_surf(benson_template_rh.F, benson_template_rh.V*0.99, color); hold on;

param5c1_a.face = themi.F;
param5c1_a.vertex = themi.V;
param5c1_a.EValue = themi;

GenFigure5c1 = GenFigureAtlas(param5c1_a);
GenFigure5c1.subPlot('lh', themi.ang).setConfigurationExtra().setLightConfigurationExtra().save('Fig5_5');


%% Description  -- function [template, benson_template] = genFigure(subject, lr, i)
%		generate figures with benson's template. This function will make
%		some figures: including subject,including lr(left and right hemisphere) and i(in order to compute):
%			1. plot angle figure
% 			2. plot ecc figure
%           3. plot figure with use uv coordinate ues benson's template
%			4. plot figure with register Benson's model to the average subject 
% parameter: 
%     subject [struct]  -- data of processing of plot
%     lr[String]        -- type of plot, left or right hemisphere
%	  i                 -- extra value, just use to compute.
%
% return: 
%      template[struct]  -- new template with data handle from benson's
%      template.
%
%	   benson_template   -- original benson's template throught load
%	   functuion.
%
function [template, benson_template] = genFigure(subject, lr, i)
	rotangs = [2.6-pi 3.8]; % for clear viewer, does not affect results
    hemi = subject.(lr);
    Fpial = double(hemi.pial_MSMAll.faces);
    Vpial = double(hemi.pial_MSMAll.vertices);
    % perform a fast march to get region to cut
    D = perform_fast_marching_mesh(Vpial,Fpial, hemi.fovid);
    ind2del  = find(D > 100); % find the index to delete
    % get a list of vertex id to delete
    [Fcut,Vcut,father] = remove_mesh_vertices(Fpial,Vpial, ind2del);
    % map to disk,
    uv = disk_conformal_map(Fcut, Vcut); % harmonic res of cutted
    rotangi = rotangs(i);
    uv = uv* [cos(rotangi) -sin(rotangi);  sin(rotangi) cos(rotangi)]';%  rotate
    hemi_cut.uv = uv;
    hemi_cut.F = Fcut;
    
%% show the polar ang
	paramPA.face = Fcut;
	paramPA.vertex = hemi_cut.uv;
	paramPA.renderColor.value = hemi.pRF(father,1);
	paramPA.renderColor.valueType = lr;
	paramPA.renderColor.maxValue = -1;
	genFigureAtlasPA = GenFigureAtlas(paramPA);
	genFigureAtlasPA.draw();
	if (i == 1)
		genFigureAtlasPA.save('Fig5_1');
	end
	
	
%% show ecc
	paramECC.face = hemi_cut.F;
	paramECC.vertex = hemi_cut.uv;
	paramECC.renderColor.value = hemi.pRF(father,2);
	paramECC.renderColor.valueType = 'ecc';	
	paramECC.renderColor.maxValue = 8;
	genFigureAtlasECC = GenFigureAtlas(paramECC);
	genFigureAtlasECC.draw();
	if (i == 1)
		genFigureAtlasECC.save('Fig5_2');
	end

	
	
%% save the boundary of stereographic projection
    bd = compute_bd(Fcut); % compute the boundary of cutted
    Vtmp = hemi.sphere.vertices(father,:); % This sphere (32k) is on fslr space
     
    hemi_cut.suv =  [Vtmp(:,1)./(100-Vtmp(:,3)) Vtmp(:,2)./(100-Vtmp(:,3))]; % sub is the projected uv by fslr sphere
    hemi_cut.bd_suv = hemi_cut.suv(bd,:);
    hemi_cut.ecc = hemi.pRF(father,2);
    hemi_cut.ang = hemi.pRF(father,1);    
       
% Load Benson's template and plot left or right hemisphere with uv
% coordinate
	[template,b_template] = GenDataManager.loadBensonAtlasTempelate(lr, hemi_cut); 
	benson_template = b_template;
	bd= compute_bd(template.F);
	uv = template.uv;  
	rbd = sqrt(uv(bd,1).^2 + uv(bd,2).^2 );
    uv(bd,1) = uv(bd,1) ./ rbd;     
	uv(bd,2) = uv(bd,2) ./ rbd;
	
	paramTMP.face = template.F;
	paramTMP.vertex = uv;
	paramTMP.EValue = template;
	paramTMP.renderColor.value = template.ang;
	paramTMP.renderColor.valueType = lr;
	paramTMP.renderColor.ecc = 0.7;
	genFigureAtlasTMP = GenFigureAtlas(paramTMP);
	genFigureAtlasTMP.draw();
	if (i == 1)
		genFigureAtlasTMP.save('Fig5_3');
	end
    
%% now register Benson's model to the average subject 
    themi = template;
    static = hemi_cut; 
    static.feature = [static.ecc, static.ang/180*pi];
    moving = themi;
    moving.feature = [moving.ecc, moving.ang/180*pi];
    moving.F = delaunay(moving.uv);
    moving.R2 = ones(size(moving.uv,1),1);
    
    [map,map_mu]  = disk_registration(static, moving, [], []);
    % update to new uv     
    template.uv_ben = template.uv;
 
    bd = compute_bd(template.F);
    bdnorm =  sqrt(map(bd,1).^2 + map(bd,2).^2);
    
    template.uv = map;
    
	paramBen.face = template.F;
	paramBen.vertex = map;
	paramBen.renderColor.value = template.ang;
	paramBen.renderColor.valueType = lr;
	genFigureAtlasBen = GenFigureAtlas(paramBen);
    genFigureAtlasBen.draw();
end


 