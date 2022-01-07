% generate the synthetic data based on the schira package
% The model has a lot parameters. Setting two sets of parameters, can
% generate two meshes. We consider on of them is tempalte in teh synthetic
% data, and another as the subject. 
%

clear; close all;clc;
addpath('schira_package')

%% Generating the template
% basic common Model Parameters:
a=1.5; % Foveal pole
b=90;  % Peripheral pole
% shear parameters alpha 1 to 3
V1linShear=1;
V2linShear=0.5;
V3linShear=0.4;

%  a few more parameters determening the precrision and time for computation
minEcc=0.05;
maxEcc=8;
% now open a small gui to set these Parameters

% now a few more settings that we provide no GUI support for.
isoEccRings=10;
isoPolarRays=50;

resolution=100; % resolution of the dots along the grid 
  
complexGrid=makeVisualGrid(minEcc,maxEcc,isoEccRings,isoPolarRays,resolution); 
%---------------------------------------------------------------------

[V1Grid,V2Grid,V3Grid]=assembleV1V3Complex(complexGrid,[V1linShear,V2linShear,V3linShear],0);


%executing the model
[V1cartx,V1carty] = bandedDoubleSech(V1Grid,a,b);
[V2cartx,V2carty] = bandedDoubleSech(V2Grid,a,b);
[V3cartx,V3carty] = bandedDoubleSech(V3Grid,a,b);


vx_t = [V1cartx,V2cartx,V3cartx];
vy_t = [V1carty,V2carty,V3carty];
label = [ones(1,length(V1carty)) 2*ones(1,length(V2carty)) 3*ones(1,length(V3carty)) ];

ecc = [complexGrid(1,:),complexGrid(1,:),complexGrid(1,:)];
ang = [complexGrid(2,:),complexGrid(2,:),complexGrid(2,:)];
% delete duplicated
[V,II,~]=unique([vx_t;vy_t]','rows');

ecc=ecc(II)';
ang = ang(II)';

ang = aroundpi(ang);


label = label(II)';


t=delaunay(V);
flag = islegalface(t,V, -3.838);
F = t(flag==1,:);



write_mfile('template.gslm','Face',F,'Vertex %d %f %f %f {vis=(%f %f %f) targe=(%f %f)}\n',[V zeros(size(V,1),1), ang, ecc, label, zeros(size(V,1),2)]);

improve_mesh_quality('template.gslm')
% copyfile('template.gslm','template_e.gslm')
% expand_mesh('template_e.gslm')

[F,V,E]=read_mfile('template.gslm');
%  V(E.Vertex_vis(:,2)<0,2) = NaN;
figure
plot_mesh(F, V, E.Vertex_vis(:,2));
shading interp;
colormap hsv;
view(0,90);
axis equal;
% axis off


 %% subject
 

% basic common Model Parameters:
a=0.5; % Foveal pole
b=90;  % Peripheral pole 
% shear parameters alpha 1 to 3
V1linShear=1.2;
V2linShear=0.5;
V3linShear=0.2; 
 
complexGrid=makeVisualGrid(minEcc,maxEcc,isoEccRings,isoPolarRays,resolution); 
[V1Grid,V2Grid,V3Grid]=assembleV1V3Complex(complexGrid,[V1linShear,V2linShear,V3linShear],0);
%executing the model
[V1cartx,V1carty] = bandedDoubleSech(V1Grid,a,b);
[V2cartx,V2carty] = bandedDoubleSech(V2Grid,a,b);
[V3cartx,V3carty] = bandedDoubleSech(V3Grid,a,b);


vx_s = [V1cartx,V2cartx,V3cartx];
vy_s = [V1carty,V2carty,V3carty];
label = [ones(1,length(V1carty)) 2*ones(1,length(V2carty)) 3*ones(1,length(V3carty)) ];
ang = [complexGrid(2,:),complexGrid(2,:),complexGrid(2,:)];
ecc = [complexGrid(1,:),complexGrid(1,:),complexGrid(1,:)];

% delete duplicated
[V,II,~]=unique([vx_s;vy_s]','rows');

ecc=ecc(II)';
ang = ang(II)';
label = label(II)';
targe = [vx_t(II); vy_t(II)]'; 

ang = aroundpi(ang);

t=delaunay(V);
flag = islegalface(t,V, -4.56);
F = t(flag==1,:);
 

write_mfile('subject.gslm','Face',F,'Vertex %d %f %f %f {vis=(%f %f %f) targe=(%f %f)}\n',[V,zeros(size(V,1),1), ang, ecc, label, targe]);

improve_mesh_quality('subject.gslm')
% copyfile('subject.gslm','subject_e.gslm')
% expand_mesh('subject_e.gslm')

 [F,V,E]=read_mfile('subject.gslm');
 
 V(E.Vertex_vis(:,2)<0,:) = NaN;
 figure
plot_mesh(F, V, E.Vertex_vis(:,2));
shading interp;
colormap hsv;
view(0,90);


function  ang = aroundpi(ang)
th = 0.018;
ang(abs(ang - pi/2)<th)=pi/2;
ang(abs(ang + pi/2)<th)=-pi/2;
ang(abs(ang)<th)=0;

end



function expand_mesh(fn)
% uniform mesh on disk

[F,V,E]=read_mfile(fn); 
[~, pt_disk]=unit_disk_mesh(5);
pt_ext = pt_disk*3 + [-4 0 ];

Vn = [V(:,1:2); pt_ext];
Fn = delaunay(Vn);

% 
bd = compute_bd(F);

V = V(:,1:2);
bdV = V([bd; bd(1)],:);
En = E;
for i=1:size(E.Vertex_vis,2)
    Fi = scatteredInterpolant(V, E.Vertex_vis(:,i));    
    En.vis(:,i) = Fi(Vn);
end


for i=1:size(E.Vertex_targe,2)
    Fi = scatteredInterpolant(V, E.Vertex_targe(:,i));
    En.targe(:,i) = Fi(Vn);
end


for i=1:2
    Fi = scatteredInterpolant(V, V(:,i));
    Vn(:,i) = Fi(Vn);
end


 in = inpolygon(Vn(:,1),Vn(:,2), bdV(:,1),bdV(:,2));
 En.vis(~in,:) = -1;
 En.targe(~in,:) = -1;
  
 
 Fn = delaunay(Vn);
 
write_mfile(fn,'Face',Fn,'Vertex %d %f %f %f {vis=(%f %f %f) targe=(%f %f)}\n',[Vn,zeros(size(Vn,1),1), En.vis, En.targe]);


end







% give n, to generate mesh on unit disk
function [Fn, Vn]=unit_disk_mesh(n)
    fn = sprintf('udisk_%d.mat',n);
    if(exist(fn))

        load(fn);
    else


        V=[0 0 0];
        for i =1:6
            V =[V; cos(2*pi/6*i) sin(2*pi/6*i) 0 ];
        end
        F=[1 2 3; 1 3 4; 1 4 5; 1 5 6; 1 6 7; 1 7 2];

        for i=1:n
            [V, F] = LoopSubdivisionLimited( V, F, 1e-2);
        end

        Fn = delaunay(V(:,1:2));
        Vn = disk_omt(Fn,V, disk_conformal_map(Fn,V));

         save(fn, 'Fn','Vn');
    end
end



function improve_mesh_quality(fn)
    % uniform mesh on disk
    [Fn, pt_disk]=unit_disk_mesh(5);


    %
    [F,V,E]=read_mfile(fn);
    uv = disk_harmonic_map(F,V);
    


    En = E; 
    for i=1:size(E.Vertex_vis,2)
        Fi = scatteredInterpolant(uv, E.Vertex_vis(:,i));
        En.vis(:,i) = Fi(pt_disk);
    end


    for i=1:size(E.Vertex_targe,2)
        Fi = scatteredInterpolant(uv, E.Vertex_targe(:,i));
        En.targe(:,i) = Fi(pt_disk);
    end

    Vn=[];
    for i=1:2
        Fi = scatteredInterpolant(uv, V(:,i));
        Vn(:,i) = Fi(pt_disk);
    end

    mu = compute_bc(Fn, [Vn,zeros(size(pt_disk,1),1)], disk_harmonic_map(Fn, [Vn,zeros(size(pt_disk,1),1)]));
    id = find(abs(mu>1));


    write_mfile(fn,'Face',Fn,'Vertex %d %f %f %f {vis=(%f %f %f) targe=(%f %f)}\n',[Vn,zeros(size(pt_disk,1),1), En.vis, En.targe]);

end