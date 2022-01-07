%% Reproduce Figure 7
% Figure 7. Retinotopic map on the left-hemisphere of the first observer in 
% Studyforrest retinotopy dataset (Sengupta et al. 2016a). (a) eccentricity map,
%(b) polar angle map, (c) eccentricity map on the disk, and \
% (d) polar angle map on the disk. (e)-(h) shows the registered results correspondingly. 
clear;clc;close all;

%% load studyforrest subject file with file id
sid = 1;
sf_subj = load_studyforrest_subject_mat(sid);
sf_subj_cut = cut_on_sphere(sf_subj);

% show studyforrest figure of subject
ecc = sf_subj_cut.lh.pRF(:,2);
ang = sf_subj_cut.lh.pRF(:,1)+90; % they shifted the  axis
uv = sf_subj_cut.lh.uv;

paramSF.face = sf_subj_cut.lh.F;
paramSF.vertex = uv;
paramSF.renderColor.value = ecc;
paramSF.renderColor.valueType = 'ecc';
paramSF.renderColor.maxValue = 320;
paramSF.renderColor.ecc = ecc;
genFigureSF01 = GenFigureRegSF(paramSF);
genFigureSF01.draw();

% plot ang
paramSF.renderColor.value = ang;
paramSF.renderColor.valueType = 'lh';
paramSF.renderColor.ecc = ang;
genFigureSF02 = GenFigureRegSF(paramSF);
genFigureSF02.draw();

%% show figure with HCP template
% figure with 'lh' type.
load([GenConsts.kTemplatesURL,'HCPTemplate'],'template');
tem_hemi = template.lh;
paramHCP.face = tem_hemi.F;
paramHCP.vertex = tem_hemi.uv;
paramHCP.renderColor.value = tem_hemi.ang;
paramHCP.renderColor.valueType = 'lh';
genFigureSF03 = GenFigureRegSF(paramHCP);
genFigureSF03.draw().setTitle('template angle');

% figure with 'ecc' type.
paramHCP.renderColor.value = tem_hemi.ecc;
paramHCP.renderColor.maxValue = 16;
paramHCP.renderColor.valueType = 'ecc';
genFigureSF04 = GenFigureRegSF(paramHCP);
genFigureSF04.draw().setTitle('template ecc');

%% Register the subject to template 

moving.F = sf_subj_cut.lh.F;
moving.uv = sf_subj_cut.lh.uv;
moving.feature = [sf_subj_cut.lh.pRF(:,2)/320*16.5 (sf_subj_cut.lh.pRF(:,1)+90)/180*pi];  
moving.R2 = sf_subj_cut.lh.pRF(:,5);
 
static.F = tem_hemi.F;
static.uv = tem_hemi.uv;
static.feature = [tem_hemi.ecc tem_hemi.ang/180*pi]; 
 

global plotset
plotset =0;
[map,map_mu]  = disk_registration(static, moving,[],[]);
%%

%% show new visual center
color = sf_subj.lh.inflated.vertices*0 +[200 160 140]/255;

paramVS.face = sf_subj.lh.inflated.faces;
paramVS.vertex = sf_subj.lh.inflated.vertices;
paramVS.renderColor.color = color;
genFigureSF05 = GenFigureRegSF(paramVS);
genFigureSF05.draw()...
.subPlot(sf_subj_cut.lh.F, 1.01*sf_subj.lh.inflated.vertices(sf_subj_cut.lh.father,:), 'lh', sf_subj_cut.lh.pRF(:,1)+90)...
.setLightConfigurationExtra().setConfigurationExtra().save('Fig_7_2');

%% show inflated figure
genFigureSF06 = GenFigureRegSF(paramVS);
genFigureSF06.draw()...
.subPlot(sf_subj_cut.lh.F, 1.01*sf_subj.lh.inflated.vertices(sf_subj_cut.lh.father,:), 'ecc', sf_subj_cut.lh.pRF(:,2),350)...
.setLightConfigurationExtra().setConfigurationExtra().save('Fig_7_1');

%% show new visual center
feature = unit_cvt( [tem_hemi.ecc  tem_hemi.ang], 'p2c');
visxfun =  scatteredInterpolant(tem_hemi.uv, feature(:,1),'natural'); % rampa
visyfun =  scatteredInterpolant(tem_hemi.uv, feature(:,2),'natural');
wrapfeature =unit_cvt([visxfun(map), visyfun(map)],'c2p');

% figure of 'lh' type
paramNVS.face = sf_subj_cut.lh.F;
paramNVS.vertex = moving.uv;
paramNVS.renderColor.value = sf_subj_cut.lh.pRF(:,1)+90;
paramNVS.renderColor.valueType = 'lh';
genFigureSF07 = GenFigureRegSF(paramNVS);
genFigureSF07.draw().save('Fig_7_4');

% figure of 'ecc' type
paramNVS.renderColor.value = sf_subj_cut.lh.pRF(:,2);
paramNVS.renderColor.valueType = 'ecc';
paramNVS.renderColor.maxValue = 350;
genFigureSF08 = GenFigureRegSF(paramNVS);
genFigureSF08.draw().save('Fig_7_3');

%% show inflated figure with another vertices
genFigureSF09 = GenFigureRegSF(paramVS);
genFigureSF09.draw()...
.subPlot(sf_subj_cut.lh.F, 1.01*sf_subj.lh.inflated.vertices(sf_subj_cut.lh.father,:), 'lh', wrapfeature(:,2))...
.setLightConfigurationExtra().setConfigurationExtra().save('Fig_7_6');

% figure of 'ecc' type
genFigureSF10 = GenFigureRegSF(paramVS);
genFigureSF10.draw()...
.subPlot(sf_subj_cut.lh.F, 1.01*sf_subj.lh.inflated.vertices(sf_subj_cut.lh.father,:), 'ecc', wrapfeature(:,1), 16)...
.setLightConfigurationExtra().setConfigurationExtra().save('Fig_7_5');

%% The registered map
paramAB.face = moving.F;
paramAB.vertex = moving.uv;
paramAB.renderColor.value = wrapfeature(:,2).*(wrapfeature(:,1)<12);
paramAB.renderColor.valueType = 'lh';
genFigureSF11 = GenFigureRegSF(paramAB);
genFigureSF11.draw().subPlotBoundary(tem_hemi).save('Fig_7_8');

wrapfeature =unit_cvt([visxfun(map), visyfun(map)],'c2p');
paramAB.renderColor.value = wrapfeature(:,1);
paramAB.renderColor.valueType = 'ecc';
paramAB.renderColor.maxValue = 16;
genFigureSF12 = GenFigureRegSF(paramAB);
genFigureSF12.draw().subPlotBoundary(tem_hemi).save('Fig_7_7');


 
