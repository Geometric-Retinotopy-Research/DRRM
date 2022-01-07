% Gen Table 6
clear;clc;close;
% no longer used


dbstop if error
%% Load Template and Subject

global template_hemi hemi_cut hcpsub sid  N debug

debug =0;
% Load subject
dirinfo = dir([GenConsts.kDataHCP32URL, '*_32k_pRF_fMRI.mat']); 

global DataBus
DataBus = [];
mkdir([GenConsts.kTablesURL, 'tmp_table6'])
for i_sid = 1:20
    sid = dirinfo(i_sid).name;

	matfn = sprintf('%stmp_table6/running_sub_%s.mat', GenConsts.kTablesURL, sid);
    if exist(matfn)
        loadstruct = load(matfn);
        DataBus = loadstruct.DataBus;
        continue
    end

    for lr = {'lh', 'rh'}
        lr = lr{1};
        % load subject
        fn = [GenConsts.kDataHCP32URL, sid];
        N = 100;
        rect = [-1 1 -1 1];

        hcpsub = load_hcp_subject_32k(fn);
        %%
        hcp_subj_cut = cut_on_sphere_32k(hcpsub);
        hemi_cut = hcp_subj_cut.(lr);

        % load template
		load([GenConsts.kTemplatesURL,'HCPTemplate'],'template');
        template_hemi = template.(lr);

        %% prepare data
        moving.F = hemi_cut.F;
        moving.uv = hemi_cut.uv;
        moving.label = hemi_cut.atlas_hcp;
        % moving.feature = unit_cvt(hemi_cut.pRF(:,[2 1]),'p2c');
        moving.feature = double([hemi_cut.pRF_half1(:,2) hemi_cut.pRF_half1(:,1)/180*pi]);
        moving.R2 = double(hemi_cut.pRF_half1(:,5));

        static.F = template_hemi.F;
        static.uv = template_hemi.uv;
        static.label = template_hemi.varea;
        % static.feature = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
        static.feature = double([template_hemi.ecc template_hemi.ang/180*pi]);
        %
        bd = compute_bd(moving.F);
        bd_src = moving.uv(bd,:);
        bd_target = moving.uv(bd,:);

      

        %%  Raw data  
        evaulate_half(moving.uv, moving.feature) % The databus will be udpated within evaulate function

       

         %% Benson's data
        tic
        t1 = toc;
        [ecc,ang, sigma, varea] = load_benson_inference_cvt_to_32k(sid(1:6), lr);
         % Benson result:  we need different evaulation function 
        evaulate_inference(ecc, ang, sigma, varea); 
         
        DataBus(end, 13) = DataBus(end, 2)/size(moving.F,1);
         
        %% Proposed without  landmark
        tic
        map  = disk_registration(static, moving, [], []);
        t1 = toc;
        evaulate_half(map, moving.feature)
        % %% Proposed with  landmark
        % map  = disk_registration(static, moving, landmark, ptdst);
        % evaulate(map, moving.feature)


    end
    save(matfn,'DataBus')
    close all
end


save('Table6','DataBus')
%%

load('Table6','DataBus')
clc
for i=[1 3]
    meandata = mean(DataBus(i:3:end,[  4  11 ]));  
    fprintf(' %1.5f\t%1.3f\n', meandata)
end


%% evaulate metric and add to DataBus
function evaulate_half(map, moving_feature)

global DataBus template_hemi hemi_cut


template_vis = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
visxfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,1)); % rampa
visyfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,2));
wrapfeature =unit_cvt([visxfun(map), visyfun(map)],'c2p');

interested_id = find(hemi_cut.atlas_wang>=1 & hemi_cut.atlas_wang<=17  & ...
    wrapfeature(:,1)<8 & hemi_cut.atlas_wang~=11 &  hemi_cut.atlas_wang~=10 );
fMRI =[];
for i=1:6
    fMRI(:,(i-1)*150+1:i*150) =  hemi_cut.fMRI{i}(interested_id,151:300);
    
end

pRF0 = hemi_cut.pRF_half1(interested_id,:);
pRF =  pRF0;

[pRF0(:,1),pRF0(:,2)] = vis2xy_pixel([pRF0(:,2), pRF0(:,1)]);
[pRF(:,1),pRF(:,2)] = vis2xy_pixel([wrapfeature(interested_id,1)  wrapfeature(interested_id,2)]);

[NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new, pc_raw, pc_new]= compute_fMRI_metrics_half(fMRI, pRF,pRF0);

roiFflag = (ismember(hemi_cut.F(:,1), interested_id)| ismember(hemi_cut.F(:,2), interested_id) | ismember(hemi_cut.F(:,3), interested_id));
roiF = hemi_cut.F(roiFflag,:);
flipnum = length(find(abs(compute_bc(roiF, hemi_cut.uv, map))>1.05));


mov_vis = unit_cvt([moving_feature(:,1)  moving_feature(:,2)/pi*180],'p2c');
visualdiff = nanmean(vecnorm([visxfun(map(interested_id,:)), visyfun(map(interested_id,:))]'-mov_vis(interested_id,:)'));

mu = abs(compute_bc(hemi_cut.F, hemi_cut.uv, map));
mumean = mean(abs(mu));
DataBus = [DataBus; visualdiff, flipnum, NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new , mumean, pc_raw, pc_new 0 0]; % last 0 is allocate for time


end

 

%%

%%
function evaulate_inference( BS_ecc32k, BS_ang32k, Bs_sigma32k, BS_varea)

global DataBus hemi_cut  template_hemi

map = hemi_cut.uv;
template_vis = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
visxfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,1)); % rampa
visyfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,2));
wrapfeature =unit_cvt([visxfun(map), visyfun(map)],'c2p');
 
% Benson's 32k is for hemisphre, but we only need a portion, which can be
% get by id_in_cut_disk. Lucikly, we have the cut father index for disk. 
id_in_cut_disk = hemi_cut.father;
BS_ecc_cut = BS_ecc32k(id_in_cut_disk);
BS_ang_cut = BS_ang32k(id_in_cut_disk);
BS_sigma_cut = Bs_sigma32k(id_in_cut_disk);

% Even in the disk, we only interested for specific visual area (selected by
% atlas wang >1 <17 not 11, 12) with ecc <8
interested_id = find( hemi_cut.atlas_wang>=1 & hemi_cut.atlas_wang<=17  & ...
                wrapfeature(:,1)<8 & BS_ecc_cut~=0 & hemi_cut.atlas_wang~=11 &...
                hemi_cut.atlas_wang~=10 );

% now compute the visual coordinate distance of template and benson, which
% is a metric about the registration
flipnum = NaN;
BS_vis = unit_cvt([BS_ecc_cut(interested_id)  BS_ang_cut(interested_id)],'p2c');
pRF_vis = unit_cvt(hemi_cut.pRF_half1(interested_id, [2 1]),'p2c');
visualdiff = nanmean(vecnorm(pRF_vis'-BS_vis'));
mumean = NaN;


fMRI =[];
for i=1:6
    fMRI(:,(i-1)*150+1:i*150) =  hemi_cut.fMRI{i}(interested_id,151:300);    
end
pRF0 = hemi_cut.pRF_half1(interested_id,:);
pRF =  pRF0;
[pRF0(:,1),pRF0(:,2)] = vis2xy_pixel([pRF0(:,2), pRF0(:,1)]);
[pRF(:,1), pRF(:,2) ] = vis2xy_pixel([BS_ecc_cut(interested_id), BS_ang_cut(interested_id)]);
pRF(:,6) = BS_sigma_cut(interested_id); 
% accurate
[NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new, pc_raw, pc_new]= compute_fMRI_metrics_half(fMRI, pRF,pRF0);


DataBus = [DataBus; visualdiff, flipnum, NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new , mumean, pc_raw, pc_new NaN NaN]; % last 0 is allocate for time

end 
 