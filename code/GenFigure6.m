%% Reproduce  Figure 6 
% Figure 6. Before and after registration for the first observer: (a) the eccentricity of first subject, (b) polar angle of the first subject, (c) registered polar angle, and (d) registered eccentricity. In (c)(d), data with eccentricity >8° are removed for clear comparison (since subjects’ max eccentricity is 8°).
close all;clc;clear;
%% load subject and show it 
sid = '100610';
lr = 'lh';
fn = sprintf('%s%s_32k_pRF_fMRI.mat',GenConsts.kDataHCP32URL, sid);
hcpsub = load_hcp_subject_32k(fn);
% cut subject
hcp_subj_cut = cut_on_sphere_32k(hcpsub);
hemi_cut = hcp_subj_cut.(lr);
F = hemi_cut.F;
pRF = hemi_cut.pRF;
ecc = pRF(:,2);
uv = hemi_cut.uv;  

%% show it for ecc type with uv coordinate
paramUV.face = F;
paramUV.vertex = uv;
paramUV.renderColor.valueType = 'ecc';
paramUV.renderColor.value = ecc;
paramUV.renderColor.ecc = ecc;
paramUV.renderColor.maxValue = 12;
genfigure01 = GenFigureRegHCP(paramUV);
genfigure01.draw().save('Fig_6_1');

% plot ang
renderColor.valueType = lr;
renderColor.value = pRF(:,1); 
paramUV.renderColor = renderColor;
genfigure02 = GenFigureRegHCP(paramUV);
genfigure02.draw().save('Fig_6_2');

%% Load template

load([GenConsts.kTemplatesURL,'HCPTemplate'],'template');
template_hemi = template.(lr);
feature = [template_hemi.ecc  template_hemi.ang]  ;

paramTMP.face = template_hemi.F;
paramTMP.vertex = template_hemi.uv;
paramTMP.renderColor.valueType = lr;
paramTMP.renderColor.value = feature(:,2);
genfigure03 = GenFigureRegHCP(paramTMP);
genfigure03.draw().setTitle('template angle');

% figure ecc
renderColor.valueType = 'ecc';
renderColor.value = feature(:,1); 
renderColor.maxValue = 8;
paramTMP.renderColor = renderColor;
genfigure04 = GenFigureRegHCP(paramTMP);
genfigure04.draw().setTitle('template ecc');


%% Register the subject to template 

moving.F = hemi_cut.F;
moving.uv = hemi_cut.uv;
moving.ecc =hemi_cut.pRF(:,2);
moving.ang =hemi_cut.pRF(:,1); 
moving.feature = [unit_cvt([moving.ecc moving.ang ],'p2c') moving.ecc ];
moving.feature = [moving.ecc, moving.ang/180*pi];
moving.R2 = hemi_cut.pRF(:,5); 

static.F = template_hemi.F;
static.uv = template_hemi.uv;
static.ecc =template_hemi.ecc;  
static.ang =template_hemi.ang;
% static.feature = [unit_cvt([static.ecc static.ang],'p2c') static.ecc];
  static.feature = [static.ecc, static.ang/180*pi];

global plotset
plotset =0;
[map,map_mu]  = disk_registration(static, moving,[],[]);

%% show new visual center
idvalid = ~isnan(template_hemi.ecc) & ~isnan(template_hemi.ang);
template_vis = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');

visxfun =  scatteredInterpolant(template_hemi.uv(idvalid,:), template_vis(idvalid,1),'natural'); % rampa
visyfun =  scatteredInterpolant(template_hemi.uv(idvalid,:), template_vis(idvalid,2),'natural');
eccang = unit_cvt([visxfun(map) visyfun(map)], 'c2p');
angchop = eccang(:,2); angchop(eccang(:,1)>8) =0;

paramVC.face = moving.F;
paramVC.vertex = map;
paramVC.renderColor.value = angchop;
paramVC.renderColor.valueType = lr;
genfigure04 = GenFigureRegHCP(paramVC);
genfigure04.draw().subPlot(template_hemi).save('Fig_6_3');

%%  generate figure 6.4 the eccentricity
wrapfeature =unit_cvt([visxfun(map), visyfun(map)],'c2p');

renderColor.value = wrapfeature(:,1);
renderColor.maxValue = 8;
renderColor.valueType = 'ecc';
paramVC.renderColor = renderColor;
genfigure05 = GenFigureRegHCP(paramVC);
genfigure05.draw().subPlot(template_hemi).save('Fig_6_4');



 

