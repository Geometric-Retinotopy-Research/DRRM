%% Reproduce the Figure 4
% Figure 4. Template and Subject Retinotopic Maps (synthetic data): (a) The predefined template; (b) noiseless retinotopic map of a subject; (c) retinotopic map of a subject with weak noise (PSNR = 20); (d) retinotopic map of a subject with strong noise (PSNR = 10); (e) ground truth displacement.  Red curves are eccentricity contours, black curves are polar angle contours, and some landmarks/anchors are marked in (a)-(d). 

clear;clc;close all;
%% load template file
[Ft,Vt,Et]=read_mfile([GenConsts.kGSLMURL, 'template.gslm']);

%% plot the register figure with the template and add the number description on each point.
param.face = Ft;
param.vertex = Vt(:,1:2);
param.EValue = Et;
param.renderColor.valueType = 'lh';
param.renderColor.value = Et.Vertex_vis(:,1)*180/pi;
figureRegSyn = GenFigureRegSyn(param);
figureRegSyn.draw().overlayPlot().subPlot(GenDataManager.getTemplateTxt()).save('Fig4_1');


%% plot the register figure with the subject and add the number description on each point.
% use the same way to plot as above 
[Fs,Vs,Es]=read_mfile([GenConsts.kGSLMURL, 'subject.gslm']);
 
param01.face = Fs;
param01.vertex = Vs(:,1:2);
param01.EValue = Es;
param01.renderColor.valueType = 'lh';
param01.renderColor.value = Es.Vertex_vis(:,1)*180/pi;
figureRegSyn01 = GenFigureRegSyn(param01);
figureRegSyn01.draw().overlayPlot().subPlot(GenDataManager.getTemplateTxt()).save('Fig4_2');

%% add small noise 
figureRegSynEn_S = GenFigureRegSyn(param01);
figureRegSynEn_S.draw().addSmallNoisy().overlayPlot().subPlot(GenDataManager.getTemplateTxt()).save('Fig4_3');

%%  add strong noise 
figureRegSynEn_L = GenFigureRegSyn(param01);
figureRegSynEn_L.draw().addStrongNoisy().overlayPlot().subPlot(GenDataManager.getTemplateTxt()).save('Fig4_4');

%% show the ground truth displacement 
roii =[];
for x = -6:0.3:2
    for y = -2.5:0.3:2.5
        d = sqrt((Vs(:,1)-x).^2 +(Vs(:,2)-y).^2 );
        [minval, id]=min(d);
        if (minval<0.05)
            roii =[roii; id];
        end
    end
end
figureRegSynEn_D = GenFigureRegSyn(param01);
figureRegSynEn_D.draw().addQuiver(roii).save('Fig4_5'); 