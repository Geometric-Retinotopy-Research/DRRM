clear;clc;close all;

addpath(genpath('geometry-processing-package'));
addpath('utilities');
addpath('toolbox_fast_marching')

subjects = dir('../data/HCP/*lh.m');
rng(0)

R2th = 1;

subi = 184;
% Prepare data
fn = subjects(subi).name;
[Fm,Vm, Em]=read_mfile(sprintf('../data/HCP/%s.m',fn(1:end-2)));

% Load the full hemisphere
[Ffull,Vfull, Efull]=read_mfile(sprintf('../data/HCP/%s_ecc.m',fn(1:end-2)));

%

foveaid = 23589;

[D,S,Q] = perform_fast_marching_mesh(double(Vfull), double(Ffull), foveaid);

% Loop thru hemispheres
ind2del  = find(D > 90);

[Fout,Vout,father] = gf_remove_mesh_vertices(Ffull,Vfull, ind2del); % get a list of vertex id to delete

%%
figure(101)
try 
    Boundaries= load('boundary.mat').Boundaries;
catch
    Boundaries = {};
end


drawecc =0;
clf
draw_background(Ffull, Efull, Fout,  R2th, father, drawecc)
for i=1:length(Boundaries)
    if ~isempty(Boundaries{i}.u)
        plot(Boundaries{i}.u(:,1), Boundaries{i}.u(:,2),'w-', 'linewidth',1.5);
    end
end
 
bi = max(length(Boundaries),1);

while 1
    k = input('Select Menu by keyboard: \n 1-3. Draw line\n 5. Next Region \n 6. Select Region \n 7. Info&Switch to  Ecc\n 8. Clear current\n 9. Edit\n 0. Quit\n');
    try
        Boundaries{bi}.u 
    catch
        Boundaries{bi}.u =[];
    end
    switch k
        case 1 
            u = gui_get_boundary();            
            Boundaries{bi}.u = [Boundaries{bi}.u;u];
%         case 2
%             u = gui_get_boundary();            
%             Boundaries{bi}.u = [Boundaries{bi}.u;u];
%         case 3
%             u = gui_get_boundary();            
%             Boundaries{bi}.u = [Boundaries{bi}.u;u];
%             
%         case 4
%             u = gui_get_boundary();            
%             Boundaries{bi}.u = [Boundaries{bi}.u;u];
%             Boundaries{bi}.v = ones(size(Boundaries{bi}.u,1),2);
%             Boundaries{bi}.v(:,2) = nan;  % in format of [r \theta]
%             Boundaries{bi}.v(:,1) = 8;  % in format of [r \theta]
            
        case 5
            bi = bi+1;
        case 6
            s= sprintf('select the id of boundary, maximal is %d', length(Boundaries));
            bi = input(s);
        case 7
            disp(bi)
            disp(Boundaries{bi})
            drawecc = 1- drawecc;
            
        case 8 
            Boundaries{bi}.u =[];
        case 9
           s= sprintf('select the id of boundary, maximal is %d', length(Boundaries));
           bi = input(s);
           
           
           while 1
               plot(Boundaries{bi}.u(:,1), Boundaries{bi}.u(:,2),'g-o', 'linewidth',1.5);
               
               [x,y, button] = ginput(1);
               if button ==3
                   break
               end
               
               [~,mid]=min((Boundaries{bi}.u(:,1) - x).^2 +  (Boundaries{bi}.u(:,2) - y).^2);
               
               [x,y, button] = ginput(1);
               if(button==97)
                   Boundaries{bi}.u = [Boundaries{bi}.u(1:mid-1,:); x y; Boundaries{bi}.u(mid:end,:)];
               else
                   Boundaries{bi}.u(mid,:) =[x y];
               end 
              
           end
           
         case 0
            break
            
    end   
    
    save('boundary.mat','Boundaries');
    
    clf
    draw_background(Ffull, Efull, Fout,  R2th, father, drawecc)
    for i=1:length(Boundaries)    
        if ~isempty(Boundaries{i}.u)
            plot(Boundaries{i}.u(:,1), Boundaries{i}.u(:,2),'w-','linewidth',1.5);    
        end
    end
end


%%
clf
draw_background(Ffull, Efull, Fout,  R2th, father, drawecc)

for i=1:length(Boundaries)
    if ~isempty(Boundaries{i}.u)
        plot(Boundaries{i}.u(:,1), Boundaries{i}.u(:,2),'w-','linewidth',1.5);
        plot(Boundaries{i}.u(1:5,1), Boundaries{i}.u(1:5,2),'wo','linewidth',1.5);
    end
        
    Boundaries{i}.node =[];
    start = 1;
    Boundaries{i}.v =Boundaries{i}.u*0;
    
    for j=1:2
        [x,y, button] = ginput(1);
        [~,mid]=min((Boundaries{i}.u(:,1) - x).^2 +  (Boundaries{i}.u(:,2) - y).^2);       
        
        k = input('r=');
        Boundaries{i}.v(start:mid,1)=k;
        k = input('theta=');
        Boundaries{i}.v(start:mid,2)=k;
        
        start = mid+1;
        
    end
    
    k = input('r=');
    Boundaries{i}.v(start:end,1)=k;
    k = input('theta=');
    Boundaries{i}.v(start:end,2)=k;
    
    
    
    
end

save('boundary.mat','Boundaries');

%%
Boundaries= load('boundary.mat').Boundaries;
for i=1:length(Boundaries)
    theta= Boundaries{i}.v(:,2);
    id = find(abs(diff(theta))>0.5);
    if isempty(id)
        Boundaries{i}.v([1;end],1) =0;
    else
        Boundaries{i}.v(id:id+1,1) =0;
    end
    
end
save('boundary.mat','Boundaries');
%% 
function draw_background(Ffull, Efull, Fout,  R2th, father, drawecc)
 
R2thmean = 40;

% figure
if drawecc
    rendercolor = prf_value_2_color('ecc', Efull.Vertex_prf(:,2),8);
else
    rendercolor = prf_value_2_color('lh', Efull.Vertex_prf(:,1));
end
id = find(Efull.Vertex_prf(:,5)<R2th);
rendercolor(id,:)= repmat([189 140 124 ]/255,length(id),1);
colorout = rendercolor(father,:);

uv = disk_conformal_map(Fout, Efull.Vertex_Vpial(father,:));

plot_surf(Fout,uv,'FaceVertexCData',colorout, 'Edgecolor','none','FaceLighting','gouraud',...
    'AmbientStrength',0.5);

[Bv, Bw,~,Nb]=read_borderfile('../data/Q1-Q6_RelatedParcellation210.L.CorticalAreas_dil.32k_fs_LR.border');

for bi=1:length(father)
    father2child(father(bi)) = bi;
end

hold on;
for i=1:length(Bv)
    bi = Bv{i}+1;
    bid = father2child(bi);
    w = Bw{i};
    w(any(bid'==0),:) =[];
    bid(any(bid'==0),:)=[];
    
    
    
    
    uvi= uv(bid(:,1),:).*w(:,1) + uv(bid(:,2),:).*w(:,2) + uv(bid(:,3),:).*w(:,3);
    
    
    inid = inpolygon(uv(:,1),uv(:,2),uvi(:,1),uvi(:,2));
    R2cut = Efull.Vertex_prf(father,5);
    if mean(R2cut(inid))>R2thmean
        plot(uvi(:,1), uvi(:,2),'k-','Linewidth',0.5)       
        mc = mean(uvi)-[0.02 0];
%         text(mc(1),mc(2), Nb{i}(3:end-4),'Color',[0,0,0],'FontSize',12, 'BackgroundColor', [1,1,1])
  
    end
end

  
[Bv, Bw,~,Nb]=read_borderfile('../data/wang2015.L.32k_fs_LR.border');

for i=1:length(father)
    father2child(father(i)) = i;
end

hold on;
for i=1:length(Bv)
    bi = Bv{i}+1;
    bid = father2child(bi);
    w = Bw{i};
    w(any(bid'==0),:) =[];
    bid(any(bid'==0),:)=[];
    
    
    
    
    uvi= uv(bid(:,1),:).*w(:,1) + uv(bid(:,2),:).*w(:,2) + uv(bid(:,3),:).*w(:,3);
    
    
    inid = inpolygon(uv(:,1),uv(:,2),uvi(:,1),uvi(:,2));
    R2cut = Efull.Vertex_prf(father,5);
    if mean(R2cut(inid))>R2thmean
        plot(uvi(:,1), uvi(:,2),'r-','Linewidth',0.5) 
        mc = mean(uvi)-[0.02 0];
%          text(mc(1),mc(2), Nb{i},'Color',[0,0,0],'FontSize',12, 'BackgroundColor', [1,1,1])
        
    end
end
  
set(gca,'Position',[0.130000 0.110000 0.775000 0.815000])
set(gca,'Xlim',[-0.194912 0.622923])
set(gca,'Ylim',[-0.425429 0.392377])
set(gca,'Zlim',[-1.000000 1.000000])
set(gcf,'position',[1.000000 41.000000 1920.000000 963.000000])
hLegend = findobj(gcf, 'Type', 'Legend');
end



%% 
function xy = gui_get_boundary()

xy =[];

while 1
     
    [x,y, button] = ginput(1);
   
    
    if(button==3 || button==27)
        break
    end
    
  
    
    xy = [xy; x y];
    
    if(button==122)
        xy(end:-1:end-2,:) =[];
    end 
    try
    delete(h);       
    catch
    end
        
    h = plot(xy(:,1),xy(:,2),'w-','linewidth',1.5);
    
end 




end

