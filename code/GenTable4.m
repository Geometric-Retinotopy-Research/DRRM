%% Generate Table Table 4. 
% Performance of DRRM-registered retinotopic maps of all observers (N=15)
% of Studyforrest retinotopy dataset (Sengupta et al. 2016a)
% relative to structurally registered retinotopic maps (3T).

% The DataBus stores every evaulated metric; we display at the end
clear;clc;close all;
 
%%
global DataBus template_hemi hemi_cut
DataBus = [];
mkdir([GenConsts.kTablesURL, 'tmp_table4'])

for sid = [1:6 9 10 14:20]
    sfsub = load_studyforrest_subject_mat(sid);
	matfn = sprintf('%stmp_table4/running_sub_%s.mat', GenConsts.kTablesURL, sid);
    if exist(matfn,'file')
        DataBus = load(matfn).DataBus;        
        continue
    end
    
    sf_subj_cut = cut_on_sphere(sfsub);
    
    for lr = {'lh','rh'}
        lr = lr{1};
        close all
        %% ===============================LH ==================================
        % show it
        hemi_cut = sf_subj_cut.(lr);
        F = hemi_cut.F;
        
        uv = hemi_cut.uv;
        rendercolor = prf_value_2_color('ecc',hemi_cut.pRF(:,2),320); % 400 is the max pixel, corrosponding 16.5
        id =  isnan(hemi_cut.pRF(:,2));
        rendercolor(id,:) = rendercolor(id,:)*0+[200 160 140]/255;
        
        figure
        plot_surf(F, uv, rendercolor)
        axis off;
        
        % plot ang
        figure
        rendercolor = prf_value_2_color(lr, hemi_cut.pRF(:,1)+90);
        id =  isnan(hemi_cut.pRF(:,1));
        rendercolor(id,:) = rendercolor(id,:)*0+[200 160 140]/255;
        plot_surf(F, uv, rendercolor)
        axis off;
        %% Load template
        
		load([GenConsts.kTemplatesURL,'HCPTemplate'],'template');
        if lr=='lh'
            template_hemi = template.lh;
        else
            template_hemi = template.rh;
        end
        
        
        %% Register it
        
        moving.F = hemi_cut.F;
        moving.uv = hemi_cut.uv;
        moving.feature = [hemi_cut.pRF(:,2)  (hemi_cut.pRF(:,1) +90)/180*pi];
        moving.R2 = hemi_cut.pRF(:,5);
        moving.R2(moving.R2<0)=0;
        
        
        static.F = template_hemi.F;
        static.uv = template_hemi.uv;
        static.feature = [template_hemi.ecc template_hemi.ang/180*pi];
        
        evaulate(moving.uv, moving.feature);
 
        [map,map_mu]  = disk_registration(static, moving,[],[],5,2000);
        
        evaulate(map, moving.feature);
    end
    save(matfn,'DataBus')

end

save('Table4','DataBus')

%%
load('Table4','DataBus')

clc

i = 1:2:15;

mdata = mean([DataBus((i-1)*2+1,1)/320*16.5,     DataBus((i-1)*2+2,2),       DataBus((i-1)*2+1,4),     DataBus((i-1)*2+2,4),   DataBus((i-1)*2+1,6),  DataBus((i-1)*2+2,6),  DataBus((i-1)*2+1,8),     DataBus((i-1)*2+2,8)], 1);
 fprintf('%1.3f\t%d\t%1.5f\t%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',mdata)
i = 2:2:15;
mdata = mean([DataBus((i-1)*2+1,1)/320*16.5,     DataBus((i-1)*2+2,2),       DataBus((i-1)*2+1,4),     DataBus((i-1)*2+2,4),   DataBus((i-1)*2+1,6),  DataBus((i-1)*2+2,6),  DataBus((i-1)*2+1,8),     DataBus((i-1)*2+2,8)], 1);
 fprintf('%1.3f\t%d\t%1.5f\t%1.5f\t%1.3f\t%1.3f\t%1.3f\t%1.3f\n',mdata)



%% evaulate metric and add to DataBus
function evaulate(map, moving_feature)

global DataBus template_hemi hemi_cut


template_vis = unit_cvt([template_hemi.ecc template_hemi.ang],'p2c');
visxfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,1)); % rampa
visyfun =  scatteredInterpolant(template_hemi.uv, template_vis(:,2));
wrapfeature =unit_cvt([visxfun(map), visyfun(map)],'c2p');


if isfield(hemi_cut,'atlas_hcp')
    
    
    interested_id = find(hemi_cut.atlas_wang>=1 & hemi_cut.atlas_wang<=17  & ...
        wrapfeature(:,1)<8 & hemi_cut.atlas_wang~=11 &  hemi_cut.atlas_wang~=10 );
    
 else
    vareafun = scatteredInterpolant(template_hemi.uv, template_hemi.varea); % rampa
    atlas_hcp =  round(vareafun(map));
    
    interested_id = find(atlas_hcp>=1   );
    
    
end

 fMRI =  hemi_cut.fMRI(interested_id,:);


pRF0 = hemi_cut.pRF(interested_id,:);
pRF =  pRF0;

[pRF0(:,1),pRF0(:,2)] = vis2xy_pixel([pRF0(:,2), pRF0(:,1)]);
[pRF(:,1),pRF(:,2)] = vis2xy_pixel([wrapfeature(interested_id,1)  wrapfeature(interested_id,2)]);

[NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new]= compute_fMRI_metrics_sf(fMRI, pRF,pRF0);

% forget about very boundary 
Fnonbd = hemi_cut.F;
bid = find(vecnorm(hemi_cut.uv')>0.9);
delid = ismember(Fnonbd(:,1), bid) | ismember(Fnonbd(:,2), bid) | ismember(Fnonbd(:,3), bid);
Fnonbd(delid,:)=[];
flipnum = length(find(abs(compute_bc(Fnonbd, hemi_cut.uv, map))>1.005));


mov_vis = unit_cvt([moving_feature(:,1)  moving_feature(:,2)/pi*180],'p2c');
visualdiff = nanmean(vecnorm([visxfun(map(interested_id,:)), visyfun(map(interested_id,:))]'-mov_vis(interested_id,:)'));

DataBus = [DataBus; visualdiff, flipnum, NRMSE_raw, NRMSE_new, R2_raw,R2_new, AIC_raw, AIC_new  ];
end
