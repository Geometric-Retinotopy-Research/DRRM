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


%% load data
movings=[];
for subi = 1:length(dirinfo)    
    [Fm,Vm, Em]=read_mfile(['data/' dirinfo(subi).name]);    
    moving.V = Em.Vertex_uv;
    moving.E = Em;
    moving.F = Fm; 
    moving.vis = prf2vis_uv(Em.Vertex_prf);
    
     
    
    %  Load Hand Draw LandMark
    lnadmark_fn = ['landmarks/' dirinfo(subi).name '.txt'];
    if exist(lnadmark_fn)
        landmarkinfo = load(lnadmark_fn);
        if(~isempty(landmarkinfo))
            
            landmark = landmarkinfo(:,1);
            target =landmarkinfo(:,2:3);
        end
        
    else
        landmark =[];
        target =[];
        mkdir('landmarks');
        figure('units','normalized','outerposition',[0 0 1 1])
        quick_plot(moving); hold on;              
        
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
            quiver(x,y,u,v,1,'color','w');
        end
        
        landmarkinfo = [landmark target];
        save(lnadmark_fn, 'landmarkinfo','-ascii');
    end
    
    
    movings{subi}.F = Fm;
    movings{subi}.V = Em.Vertex_uv;
    movings{subi}.E = Em;
    movings{subi}.vis = prf2vis_uv(Em.Vertex_prf);
    movings{subi}.landmark = landmark;
    movings{subi}.target = target;
    
    
end

%% load template boundary, if does not exist, create it
template_boundary_fn = 'template_boundary.mat';
if (exist(template_boundary_fn))
    load(template_boundary_fn)   
else    
    
    quick_plot(static);
    atlaslabels =[1 5];
    for i =  1:length(atlaslabels)
        v_bd = get_boundary_vertex(static.F,static.E,atlaslabels(i));%
        plot(v_bd(:,1), v_bd(:,2), 'k');        
    end
    colorpool = ['y','m','c','w','k'];
    for i=1:3
        x=[];    y=[];
        for j=1:40
            [x1,y1] = ginput(1);
            
            x =[x;x1];
            y =[y;y1];
            plot(x,y,'-o', 'color', colorpool(i));
            title(num2str(j));
        end
        template_boundary{i} =[x,y];
    end    
    save(template_boundary_fn,'template_boundary');
end


%%
engs = [];
for iter = 1:6
    
    maps =[];
    engs(iter) =0;
    for subi = 1:length(dirinfo)    
        moving = movings{subi};                
        [map,map_mu] = QCWIR(moving,static, moving.landmark, moving.target);        
        maps{subi} = map;         
        engs(iter) = engs(iter) + registrat_energy(moving, static, map);
    end
    
    
    % regenerate the new template, directly add mean map will make
    % overlapping, so we shall call registration again.     
    meanstatic = static;      meanstatic.vis=0;
    for i=1:length(dirinfo)
        meanstatic.vis = meanstatic.vis + mesh_interp( maps{i}, movings{i}.vis, static.V) ;
    end
    meanstatic.vis = meanstatic.vis/length(dirinfo);    
    [map2mean,~] = QCWIR(static, meanstatic, [], []); % reg template to mean value       
    

    
    % morphe the landmarks and boundaries
    for i=1:length(dirinfo)
        movings{subi}.target = mesh_interp(static.V, map2mean,  movings{subi}.V(movings{subi}.landmark,:)) ;    
    end

    for i = 1:length(template_boundary)
        template_boundary{i} = mesh_interp(static.V, map2mean, template_boundary{i}) ;   
    end
    
    
    % new template;
    static.V = map2mean;
    static.vis = meanstatic.vis;
    
         
    quick_plot(static, template_boundary, map2mean);
    title(sprintf('Iteration %d', iter))
    
end
%% calculate the energy before registration
eng0 =0;
for subi = 1:length(dirinfo)
    moving = movings{subi};
    eng0 = eng0 + registrat_energy(moving, static);
end
%%
figure;
h=plot([eng0 engs(1:6)],'-rs','LineWidth',4,...
    'MarkerSize',10,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0.5,0.5,0.5]);
set(gca,'FontSize',20)
xlabel('Iteration');
ylabel('Sum of Visual Mean Error');
grid on;


    



