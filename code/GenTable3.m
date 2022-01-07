%% Reproduce Table 3
% Table 3. Average performance of DRRM-registered retinotopic maps of all (N=180) 
% observers in the HCP retinotopy dataset relative to “Raw” retinotopic maps (7T).
% The DataBus stores every evaulated metric; we display at the end


clear;clc;close all;
% Gen Table 2
dbstop if error
%% Load Template and Subject
global DataBus
DataBus = [];
global template_hemi hemi_cut hcpsub sid
dirinfo = dir([GenConsts.kDataMeshURL,'*lh.m']);
mkdir([GenConsts.kTablesURL, 'tmp_table3'])
for iter_id=1:length(dirinfo) % loop through subjects
    % Load subject
    sid = dirinfo(iter_id).name(1:end-4);
    matfn = sprintf('%stmp_table3/running_sub_%s.mat', GenConsts.kTablesURL, sid);
    iter_id
    if exist(matfn,'file')
        DataBus = load(matfn).DataBus;        
        continue
    end    iter_id       
        
    
    for lr = {'lh', 'rh'}
        lr = lr{1};
        % load subject
        fn = [GenConsts.kDataHCP32URL, sid, '_32k_pRF_fMRI.mat'];
        hcpsub = load_hcp_subject_32k(fn);
        figure;plot_surf(hcpsub.lh.sphere.faces, hcpsub.lh.sphere.vertices, prf_value_2_color('ecc',hcpsub.lh.pRF(:,2),8)); view(24,10)
        %%
        hcp_subj_cut = cut_on_sphere_32k(hcpsub);
        hemi_cut = hcp_subj_cut.(lr);
        
        % load template
		load([GenConsts.kTemplatesURL,'HCPTemplate'],'template');
        template_hemi = template.(lr);
        
        %% prepare data
        moving.F = hemi_cut.F;
        moving.uv = hemi_cut.uv;
       
        
        % moving.feature = unit_cvt(hemi_cut.pRF(:,[2 1]),'p2c');
        moving.feature = [hemi_cut.pRF(:,2) hemi_cut.pRF(:,1)/180*pi];
        moving.R2 = hemi_cut.pRF(:,5);        
       
        
        static.F = template_hemi.F;
        static.uv = template_hemi.uv;
        % static.feature = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
        static.feature = [template_hemi.ecc template_hemi.ang/180*pi];
         
        bd = compute_bd(moving.F);
        bd_src = moving.uv(bd,:);
        bd_target = moving.uv(bd,:);
        
        %% evaulate without reg.
        evaulate( moving.uv, moving.feature)        
        
        %% Proposed 
       
        [map,map_mu]  = disk_registration(static, moving, [], []);        
       
        evaulate(map, moving.feature)
    end
  
    save(matfn,'DataBus')
end
save('Table3','DataBus')


%%

clc
for i=1:10
    fprintf('%1.3f\t%d\t%1.5f\t%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',...
        DataBus((i-1)*2+1,1),     DataBus((i-1)*2+2,2), ...
        DataBus((i-1)*2+1,4),     DataBus((i-1)*2+2,4), ...
        DataBus((i-1)*2+1,6),     DataBus((i-1)*2+2,6), ...
        DataBus((i-1)*2+1,8),     DataBus((i-1)*2+2,8))
    
end
ids = 1:2:length(dirinfo);
fprintf('%1.3f\t%d\t%1.5f\t%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',...
    mean(DataBus((ids-1)*2+1,1)),  mean(DataBus((ids-1)*2+2,2)), ...
    mean(DataBus((ids-1)*2+1,4)),      mean(DataBus((ids-1)*2+2,4)), ...
    mean(DataBus((ids-1)*2+1,6)),     mean( DataBus((ids-1)*2+2,6)), ...
    mean( DataBus((ids-1)*2+1,8)),      mean(DataBus((ids-1)*2+2,8)))

ids = 2:2:length(dirinfo);
fprintf('%1.3f\t%d\t%1.5f\t%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',...
    mean(DataBus((ids-1)*2+1,1)),  mean(DataBus((ids-1)*2+2,2)), ...
    mean(DataBus((ids-1)*2+1,4)),      mean(DataBus((ids-1)*2+2,4)), ...
    mean(DataBus((ids-1)*2+1,6)),     mean( DataBus((ids-1)*2+2,6)), ...
    mean( DataBus((ids-1)*2+1,8)),      mean(DataBus((ids-1)*2+2,8)))

% stat
ids = 1:2*length(dirinfo);
length(find((DataBus((ids-1)*2+2,8) -  DataBus((ids-1)*2+1,8))<=0))
length(find((DataBus((ids-1)*2+2,8) -  DataBus((ids-1)*2+1,8))>0))
%% evaulate metric and add to DataBus   
function evaulate(map, moving_feature)

    global DataBus template_hemi hemi_cut 


    template_vis = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
    visxfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,1)); % rampa
    visyfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,2));
    wrapfeature =unit_cvt([visxfun(map), visyfun(map)],'c2p');




    interested_id = find(hemi_cut.atlas_wang>=1 & hemi_cut.atlas_wang<=17  & ...
                            wrapfeature(:,1)<8 & hemi_cut.atlas_wang~=11 &  hemi_cut.atlas_wang~=10 );
    % interested_id = interested_id(1:20:end); % for speed


    fMRI =[]; 
    for i=1:6
        fMRI(:,(i-1)*300+1:i*300) =  hemi_cut.fMRI{i}(interested_id,:);

    end



    pRF0 = hemi_cut.pRF(interested_id,:);
    pRF =  pRF0;

    [pRF0(:,1),pRF0(:,2)] = vis2xy_pixel([pRF0(:,2), pRF0(:,1)]);
    [pRF(:,1),pRF(:,2)] = vis2xy_pixel([wrapfeature(interested_id,1)  wrapfeature(interested_id,2)]);

    [NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new]= compute_fMRI_metrics(fMRI, pRF,pRF0);

    flipnum = length(find(abs(compute_bc(hemi_cut.F, hemi_cut.uv, map))>1.5));


    mov_vis = unit_cvt([moving_feature(:,1)  moving_feature(:,2)/pi*180],'p2c'); 
    visualdiff = nanmean(vecnorm([visxfun(map(interested_id,:)), visyfun(map(interested_id,:))]'-mov_vis(interested_id,:)'));

    DataBus = [DataBus; visualdiff, flipnum, NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new  ];
end

 
