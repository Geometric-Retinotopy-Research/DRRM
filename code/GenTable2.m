clear;clc;close all;
% Gen Table 2
dbstop if error
%% Load Template and Subject

global template_hemi hemi_cut hcpsub sid  N debug

debug =0;
% Load subject

dirinfo = dir([GenConsts.kDataHCP32URL,'*_32k_pRF_fMRI.mat']);

global DataBus
DataBus = [];
mkdir([GenConsts.kTablesURL, 'tmp_table2'])
for i_sid = 1:20
    sid = dirinfo(i_sid).name;
    
    matfn = sprintf('%stmp_table2/running_sub_%s.mat', GenConsts.kTablesURL, sid);
    if exist(matfn)
        continue
    end
        
    for lr = {'lh', 'rh'}
        lr = lr{1};
        % load subject
        fn = [GenConsts.kDataHCP32URL, sid];
        N = 100;
        rect = [-1 1 -1 1];
        
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
        moving.label = hemi_cut.atlas_hcp;
        % moving.feature = unit_cvt(hemi_cut.pRF(:,[2 1]),'p2c');
        moving.feature = [hemi_cut.pRF(:,2) hemi_cut.pRF(:,1)/180*pi];
        moving.R2 = hemi_cut.pRF(:,5);
        
        static.F = template_hemi.F;
        static.uv = template_hemi.uv;
        static.label = template_hemi.varea;
        % static.feature = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
        static.feature = [template_hemi.ecc template_hemi.ang/180*pi];
        
        
        %
        bd = compute_bd(moving.F);
        bd_src = moving.uv(bd,:);
        bd_target = moving.uv(bd,:);
        
        %
        [landmark, target] = detect_landmark_by_R2(static, moving);
        ptdst = target;
        ptsrc = moving.uv(landmark,:);
        % ptsrc = [];
        % ptdst = [];
        % generate image for D-Demos
        x_t_fn = [GenConsts.kImageURL,'x_t.png'];
        y_t_fn = [GenConsts.kImageURL,'y_t.png'];
        
        x_s_fn = [GenConsts.kImageURL, 'ang_s.png'];
        y_s_fn = [GenConsts.kImageURL, 'ecc_s.png'];
        
        value2image(static.F, static.uv, static.feature(:,1), rect, x_t_fn);
        value2image(static.F, static.uv, static.feature(:,2), rect, y_t_fn);
        
        value2image(moving.F, moving.uv, moving.feature(:,1), rect, x_s_fn);
        value2image(moving.F, moving.uv, moving.feature(:,2), rect, y_s_fn);
        
        y_t = im2double(imread(y_t_fn));
        x_t = im2double(imread(x_t_fn));
        y_s = im2double(imread(y_s_fn));
        x_s = im2double(imread(x_s_fn));
        
        %% Register by methods
        
        % Raw data
        evaulate(moving.uv, moving.feature) % The databus will be udpated within evaulate function
        
        
        %% TPS: It use landmark
        tic;
        disp_tps = tpswarp2d(moving.uv, [ptsrc; bd_src], [ptdst; bd_target]);
        t1 = toc;
        evaulate(moving.uv+ disp_tps, moving.feature)
        DataBus(end, 12) = t1;
        DataBus(end, 13) = DataBus(end, 2)/size(moving.F,1);
        
        %% Bayessian
        tic
        
        disp_bayessian = bayessian_reg(moving.F, moving.uv, moving.feature, ...
            static.uv, static.feature, [ptsrc; bd_src], [ptdst ; bd_target]);
        t1 = toc;
        evaulate(moving.uv + disp_bayessian, moving.feature)
        DataBus(end, 12) = t1;
        DataBus(end, 13) = DataBus(end, 2)/size(moving.F,1);
        %% Log Demo
        tic
        res = logdemo(x_s, x_t);
        for iter = 1:50
            res = logdemo(x_s, x_t,  res.vx, res.vy);
            res = logdemo(y_s, y_t, res.vx, res.vy);
        end
        t1 = toc;
        evaulate( moving.uv + res2vxvy(res, moving.uv), moving.feature)
        DataBus(end, 12) = t1;
        DataBus(end, 13) = DataBus(end, 2)/size(moving.F,1);
        %% Proposed without  landmark
        tic
        
        map  = disk_registration(static, moving, [], []);
        t1 = toc;
        
        evaulate(map, moving.feature)
        % %% Proposed with  landmark
        % map  = disk_registration(static, moving, landmark, ptdst);
        % evaulate(map, moving.feature)
        DataBus(end, 12) = t1;
        DataBus(end, 13) = DataBus(end, 2)/size(moving.F,1);
        
    end
    save(matfn,'DataBus')
end


save('Table2','DataBus')

%% 
load('Table2','DataBus')
clc
 


 for i=1:5
     meandata = mean(DataBus(i:5:end,[1 2 13 4 7 11 12]));
     if meandata(2)>0
         meandata(2) = round(meandata(2));
     end
     if meandata(6)<0
         meandata(6) = 0.001;
     end
     meandata(3)  = meandata(3) *100;
fprintf('%1.5f\t%d\t%1.2f\t%1.5f\t%1.3f\t%1.3f\t%1.2f\n', meandata)
end


% [visualdiff, flipnum, NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new , mumean, pc_raw, pc_new, time]

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
 
[NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new, pc_raw, pc_new]= compute_fMRI_metrics(fMRI, pRF,pRF0);
 

flipnum = length(find(abs(compute_bc(hemi_cut.F, hemi_cut.uv, map))>1.001));


mov_vis = unit_cvt([moving_feature(:,1)  moving_feature(:,2)/pi*180],'p2c'); 
visualdiff = nanmean(vecnorm([visxfun(map(interested_id,:)), visyfun(map(interested_id,:))]'-mov_vis(interested_id,:)'));

mu = abs(compute_bc(hemi_cut.F, hemi_cut.uv, map));
mumean = mean(abs(mu));
DataBus = [DataBus; visualdiff, flipnum, NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new , mumean, pc_raw, pc_new 0 0]; % last 0 is allocate for time 
end


      

%%
function vxvy = res2vxvy(res, V)
N =100;
[x,y]=meshgrid(1:N,1:N);
vx = res.vx(:,:, end);
vy = res.vy(:,:, end);
Fx = scatteredInterpolant(x(:),y(:), vx(:));
Fy = scatteredInterpolant(x(:),y(:), vy(:));
vxvy = [Fx(V) Fy(V)];
end


%%
function  [landmark, target] = detect_landmark_by_R2(static, moving)
landmark = find(moving.R2 > 25); 
for i=1:length(landmark)
    li = landmark(i);
    
    dis = norm((moving.feature(li,:) - static.feature)');
    dis(static.label~=moving.label(li)) = inf;
    [~, mid] = min(dis);          
    target(i,:) = static.uv(mid,:);
end



end  