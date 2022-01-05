% generate Figure 5
clear ;clc;close all;
addpath('utilities')
addpath(genpath('C:\MATLABPackages'))

dirinfo = dir('data/*lh*.m');


% load all subjects and take mean visual prf as initial template
static.vis = 0;
for subi = 1:length(dirinfo)
    [Ft,Vt, Et]=read_mfile(['data/' dirinfo(subi).name]);
    static.vis = static.vis + prf2vis_uv(Et.Vertex_prf); % dont't use r, theta, use x,y, instead. The reason is add theta directly may cause problem
end
static.F = Ft;
static.V = Et.Vertex_uv;
static.E =  Et;
static.vis = static.vis/length(dirinfo);


%%
close all
for subi = 3
    
    [Fm,Vm, Em]=read_mfile(['data/' dirinfo(subi).name]);
    % Em=gf_smooth(Fm,Vm,Em);
    
    moving.F = Fm;
    moving.V = Em.Vertex_uv;
    moving.E = Em;
    moving.vis = prf2vis_uv(Em.Vertex_prf);
    
    clc
    landmark =[];
    target =[];
    %  Load Hand Draw LandMark
    landmarkfn = ['landmarks/' dirinfo(subi).name '.txt'];
    if exist(landmarkfn)
        landmarkinfo = load(landmarkfn);
        if(~isempty(landmarkinfo))
            
            landmark = [landmark ; landmarkinfo(:,1)];
            target = [target ; landmarkinfo(:,2:3)];
        end
        
    else
        mkdir('landmarks');
        figure('units','normalized','outerposition',[0 0 1 1])
        quick_plot(moving); hold on;
        if(~isempty(target))
            x = moving.V(landmark,1);
            y = moving.V(landmark,2);
            u = target(end,1)-x;
            v = target(end,2)-y;
            myquiver(x,y,u,v);
        end
        
        
        for i=1:15
            srcpt = ginput(1);
            if(norm(srcpt)>1)
                continue
            end
            [~, landmarkid]=min(vecnorm(moving.V' - srcpt'));
            landmark = [landmark; landmarkid];
            target =  [target;ginput(1)];
            x = moving.V(landmarkid,1);
            y = moving.V(landmarkid,2);
            u = target(end,1)-x;
            v = target(end,2)-y;
            myquiver(x,y,u,v);
        end
        
        landmarkinfo = [landmark target];
        save(landmarkfn, 'landmarkinfo','-ascii');
    end
    
    
    if(~isempty(target))
        quick_plot(moving);
        x = moving.V(landmark,1);
        y = moving.V(landmark,2);
        u = target(:,1)-x;
        v = target(:,2)-y;
        quiver(x,y,u,v,0,'LineWidth',1.5,'color','w');
    end
    
    
    [map,map_mu] = QCWIR(moving,static, landmark, target);
    % map= QCHR4RM(moving,static);
    
    maps{subi} = map;
    
    quick_plot(moving, map);
    
    
end


% convert prf to viusal uv coordinates in degree
function ptfnew = prf2vis_uv(prf)
ptfnew = prf;
ecc= prf(:,2);
theta =  prf(:,1)/180*pi;
ptfnew(:,1:2) = [ecc.*cos(theta) ecc.*sin(theta)];
end

function [ang,ecc] = vis2ang(vis)

ang = atan2(vis(:,2),vis(:,1) )*180/pi;

ecc = sqrt(vis(:,2).^2 + vis(:,1).^2);
end


% custom plot
function quick_plot(static, map)

figure
pltcolor =  prf_value_2_color('lh', vis2ang(static.vis));
pltcolor(static.vis(:,5)<5,:) = 0.75;

if nargin ==1
    plot_surf(static.F, static.E.Vertex_uv, 'FaceVertexCData',pltcolor);
else
    plot_surf(static.F, map, 'FaceVertexCData',pltcolor);
end
hold on;

landmarkfn = 'template_landmark.mat';
colorpool = ['y','m','c','w','k'];
if (exist(landmarkfn))
    load(landmarkfn)
    for i=1:length(landmark)
        plot(landmark{i}(:,1), landmark{i}(:,2), '-','color',colorpool(i), 'Linewidth',2);
    end
else
    
    atlaslabels =[1 5];
    for i =  1:length(atlaslabels)
        v_bd = get_boundary_vertex(static.F,static.E,atlaslabels(i));%
        plot(v_bd(:,1), v_bd(:,2), 'k');
        
    end
    
    for i=1:3
        x=[];    y=[];
        for j=1:40
            [x1,y1] = ginput(1);
            
            x =[x;x1];
            y =[y;y1];
            plot(x,y,'-o', 'color', colorpool(i));
            title(num2str(j));
        end
        landmark{i} =[x,y];
    end
    
    save(landmarkfn,'landmark');
end
axis off;

end


