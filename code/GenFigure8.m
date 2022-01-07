% Reproduce figure8: boundary delination
% Figure 8. The visual area inferred by various methods: (a) MSMALL,
% (b) TPS, (c) Benson’s Method, (d) D-Demos, and (e) DRRM. 

% The DataBus stores every evaulated metric; we display at the end


clear;clc;close all;
% Gen Table 2
dbstop if error
%% Load Template and Subject

global template_hemi hemi_cut hcpsub sid  N 

% Load subject

dirinfo = dir([GenConsts.kDataHCP32URL,'*_32k_pRF_fMRI.mat']);


i_sid = 3;
sid = dirinfo(i_sid).name;


% for lr = {'lh', 'rh'}
%     lr = lr{1};
lr = 'lh';
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
global template_hemi
template_hemi = template.(lr);

%% prepare data
moving.F = hemi_cut.F;
moving.uv = hemi_cut.uv;
moving.label = hemi_cut.atlas_hcp;
% moving.feature = unit_cvt(hemi_cut.pRF(:,[2 1]),'p2c');
moving.feature = [hemi_cut.pRF(:,2) hemi_cut.pRF(:,1)];
moving.R2 = hemi_cut.pRF(:,5);
global static
static.F = template_hemi.F;
static.uv = template_hemi.uv;
static.label = template_hemi.varea;
% static.feature = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
static.feature = [template_hemi.ecc template_hemi.ang];


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
x_t_fn = [GenConsts.kImageURL, 'x_t.png'];
y_t_fn = [GenConsts.kImageURL, 'y_t.png'];

x_s_fn = [GenConsts.kImageURL,'ang_s.png'];
y_s_fn = [GenConsts.kImageURL,'ecc_s.png'];

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
render_boundary(moving.F,moving.uv, moving.feature).save('Fig8_1'); % The databus will be udpated within evaulate function

%% TPS: It use landmark
tic;
disp_tps = tpswarp2d(moving.uv, [ptsrc; bd_src], [ptdst; bd_target]);
t1 = toc;
render_boundary(moving.F,moving.uv+ disp_tps, moving.feature).save('Fig8_2');

%% Bayessian
tic

disp_bayessian = bayessian_reg(moving.F, moving.uv, moving.feature, ...
    static.uv, static.feature, [ptsrc; bd_src], [ptdst ; bd_target]);
t1 = toc;
render_boundary(moving.F,moving.uv + disp_bayessian, moving.feature).save('Fig8_3');

%% Log Demo
tic
res = logdemo(x_s, x_t);
for iter = 1:50
    res = logdemo(x_s, x_t,  res.vx, res.vy);
    res = logdemo(y_s, y_t, res.vx, res.vy);
end
t1 = toc;
render_boundary(moving.F, moving.uv + res2vxvy(res, moving.uv), moving.feature).save('Fig8_4');

%% Proposed without  landmark
tic
map  = disk_registration(static, moving, [], []);
t1 = toc;
%
render_boundary(moving.F, map, moving.feature).save('Fig8_5')

%% Description  -- function render_boundary(F, map, feature)
%		render boundary with features
%
% Parameter(s): 
%     F[double array]       -- connectivity of mesh
%	  map[double array]     -- vertex of mesh
%     feature[double array] -- feature array use to compute render color
% 
%% plot boundary of the result
function genFigureCB = render_boundary(F, map, feature)
	close all;
	global static   template_hemi
	param.face = F;
	param.vertex = map;
	param.renderColor.value = feature(:,2);
	param.renderColor.valueType = 'lh';
	genFigureCB = GenFigureComBoundary(param); 
	genFigureCB.draw().subPlot(static, template_hemi);
end


%% Description  -- function [landmark, target] = detect_landmark_by_R2(static, moving)
%		convert res data to vxvy coordinate.
function vxvy = res2vxvy(res, V)
    N =100;
    [x,y]=meshgrid(1:N,1:N);
    vx = res.vx(:,:, end);
    vy = res.vy(:,:, end);
    Fx = scatteredInterpolant(x(:),y(:), vx(:));
    Fy = scatteredInterpolant(x(:),y(:), vy(:));
    vxvy = [Fx(V) Fy(V)];
end

%% Description  -- function [landmark, target] = detect_landmark_by_R2(static, moving)
%		detect the landmark of by R2 data
%
% Parameter(s): 
%		static[struct]         -- the static mesh structure 
%		moving[struct]         -- the moving mesh structure
% Return: 
%		landmark[double array] -- landmark including boundary index
%		target[double array]   -- landmark target 
% 
%%
function [landmark, target] = detect_landmark_by_R2(static, moving)
landmark = find(moving.R2 > 25);
for i=1:length(landmark)
    li = landmark(i);
    
    dis = norm((moving.feature(li,:) - static.feature)');
    dis(static.label~=moving.label(li)) = inf;
    [~, mid] = min(dis);
    target(i,:) = static.uv(mid,:);
end
end



