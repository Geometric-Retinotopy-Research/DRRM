clear;clc;close all;
% gen_syntehtic_Data()
[Ft,Vt,Et]=read_mfile('template.gslm');  % expanded mesh
[Fs,Vs,Es]=read_mfile('subject.gslm');
   

Vt = Vt(:,1:2); Vs = Vs(:,1:2);
disp_true = Es.Vertex_targe - Vs; 

Table1 = zeros(6,7);
 Nlevels = [20 10] ;
%% TPS: It only use landmark, as well as we find boundary

for i = 1:length(Nlevels)
    
% load  Human  landmark.
ptsrc = load(sprintf('subject_landmark_noise%d.txt',i));
ptdst = load('template_landmark.txt');
    
Vis_s = [Es.Vertex_vis(:,2).*cos(Es.Vertex_vis(:,1)) Es.Vertex_vis(:,2).*sin(Es.Vertex_vis(:,1))];
Vis_t = [Et.Vertex_vis(:,2).*cos(Et.Vertex_vis(:,1)) Et.Vertex_vis(:,2).*sin(Et.Vertex_vis(:,1))];

bd_s = compute_bd(Fs);
bd_t = compute_bd(Ft);

% find boundary target
bd_target = Vt(bd_t(find_nearest_id(Vis_t(bd_t,:), Vis_s(bd_s,:))),:);
bd_src = Vs(bd_s,:);


t0 = tic;
disp_tps = tpswarp2d(Vs, [ptsrc; bd_src], [ptdst; bd_target]); 
t1 =toc(t0);


[meanerr, maxerr, overlapping]=evaulate_metric(Fs, Vs, disp_tps, disp_true);
fprintf('TPS(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', Nlevels(i), meanerr, maxerr, overlapping);
Table1(1,i) = meanerr;
Table1(1,i+2) = maxerr;
Table1(1,i+4) = overlapping;
Table1(1,7) = Table1(1,7) + t1;

end
 


%%   Bayessian

addpath('../../../BensonRegister')
rng(0)
NoiseLevels = [0.173 0.5473];   
 
for i = 1:length(NoiseLevels)
    
    noiseI = NoiseLevels(i);
    
    %  add noise    
    Esnoise = Es.Vertex_vis;
    Esnoise(:,1:2) = Esnoise(:,1:2) + noiseI*rand(size(Esnoise(:,1:2)));    
    psnr(Es.Vertex_vis(:,1:2), Esnoise(:,1:2))
    
    ptsrc = load(sprintf('subject_landmark_noise%d.txt',i));
    ptdst = load('template_landmark.txt');
    
    
    
    Vis_s = [Es.Vertex_vis(:,2).*cos(Es.Vertex_vis(:,1)) Es.Vertex_vis(:,2).*sin(Es.Vertex_vis(:,1))];
    Vis_t = [Et.Vertex_vis(:,2).*cos(Et.Vertex_vis(:,1)) Et.Vertex_vis(:,2).*sin(Et.Vertex_vis(:,1))];
    
    bd_s = compute_bd(Fs);
    bd_t = compute_bd(Ft);
    
    % find boundary target
    bd_target = Vt(bd_t(find_nearest_id(Vis_t(bd_t,:), Vis_s(bd_s,:))),:);
    bd_src = Vs(bd_s,:);

    % Bayessian run registration outsize
     t0 = tic;
    disp_bayessian = bayessian_reg(Fs, Vs, Esnoise, Vt, Et.Vertex_vis, [ptsrc; bd_src], [ptdst ; bd_target]);
     t1 = toc(t0);
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, Vs, disp_bayessian, disp_true);
    fprintf('Bayessian(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', Nlevels(i), meanerr, maxerr, overlapping);   
     
    
    Table1(2,i) = meanerr;
    Table1(2,i+2) = maxerr;
    Table1(2,i+4) = overlapping;
    Table1(2,7) = Table1(2,7) + t1;


end


%% Gen Images

addpath(genpath('../../../lddmm'));
addpath(genpath('../../../DDemons'));
addpath(genpath('../../../QCHR'));
rng(0)
NoiseLevels = [0.173 0.5473];  
 
rect = [-6 -1 -2.5 2.5];
 N = 100;
for i = 1:length(NoiseLevels)
    
    
    % load  Human  landmark.
    ptsrc = load(sprintf('subject_landmark_noise%d.txt',i));
    ptdst = load('template_landmark.txt');


    %  add noise        
    noiseI = NoiseLevels(i);    
    Esnoise = Es.Vertex_vis;
    Esnoise(:,1:2) = Esnoise(:,1:2) + noiseI*rand(size(Esnoise(:,1:2)));    
    psnr(Es.Vertex_vis(:,1:2), Esnoise(:,1:2))
    
    
    ang_t_fn = 'ang_t.png';
    ecc_t_fn = 'ecc_t.png';
    
    ang_s_fn = sprintf('ang_s%d.png',i);
    ecc_s_fn = sprintf('ecc_s%d.png',i);
    
    value2image(Ft, Vt, Et.Vertex_vis(:,1), rect, ang_t_fn);
    value2image(Ft, Vt, Et.Vertex_vis(:,2), rect, ecc_t_fn);
    
    value2image(Fs, Vs, Esnoise(:,1), rect, ang_s_fn);
    value2image(Fs, Vs, Esnoise(:,2), rect, ecc_s_fn);
    
    ecc_t = im2double(imread(ecc_t_fn));
    ang_t = im2double(imread(ang_t_fn));
    ecc_s = im2double(imread(ecc_s_fn));
    ang_s = im2double(imread(ang_s_fn));         
    
    
    %  QCHR    
    
    Vis_s = [Es.Vertex_vis(:,2).*cos(Es.Vertex_vis(:,1)) Es.Vertex_vis(:,2).*sin(Es.Vertex_vis(:,1))];
    Vis_t = [Et.Vertex_vis(:,2).*cos(Et.Vertex_vis(:,1)) Et.Vertex_vis(:,2).*sin(Et.Vertex_vis(:,1))];
    
    bd_s = compute_bd(Fs);
    bd_t = compute_bd(Ft);
    
    % find boundary target
    bd_target = Vt(bd_t(find_nearest_id(Vis_t(bd_t,:), Vis_s(bd_s,:))),:);
    bd_src = Vs(bd_s,:);
    
     

    t0 = tic;     
    src_rc = xy2rc([ptsrc], N, rect);
    dst_rc = xy2rc([ptdst], N, rect);
    qchr_map = QCHR_reg(ang_s, ang_t, src_rc, dst_rc);
     
    for iter = 1:5                
        qchr_map = QCHR_reg(ecc_s, ecc_t, src_rc, dst_rc, qchr_map);     
        qchr_map = QCHR_reg(ang_s, ang_t, src_rc, dst_rc, qchr_map);
    end  
    t1= toc(t0);
    
   %%
    
    [tempface, tempvertex] = image_meshgen(size(ang_s,1),size(ang_s,2));        
    gridxy = rc2xy(tempvertex,N,rect);
    vxy =  rc2xy(qchr_map,N,rect) - gridxy;     
     
    Fx = scatteredInterpolant(gridxy(:,1),gridxy(:,2), vxy(:,1)); 
    Fy = scatteredInterpolant(gridxy(:,1),gridxy(:,2), vxy(:,2)); 
    disp_qchr = [Fx(Vs) Fy(Vs)];
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, Vs, disp_qchr, disp_true);
    fprintf('QCHR(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', Nlevels(i), meanerr, maxerr, overlapping);
    % load result and evaulate 
    
    Table1(5,i) = meanerr;
    Table1(5,i+2) = maxerr;
    Table1(5,i+4) = overlapping;
    Table1(5,7) = Table1(5,7) + t1;   
    
    
    
    %% LDDMM
    t0 = tic;
    res = lddmm(ang_s, ang_t);
    for iter = 1:2
        res = lddmm(ang_s, ang_t,  res.vx, res.vy);
        res = lddmm( ecc_s, ecc_t, res.vx, res.vy);
        close all;
    end  
    t1 = toc(t0);
    
    [x,y]=meshgrid(1:N,1:N);
    vx = res.vx(:,:, end);
    vy = res.vy(:,:, end);
    Fx = scatteredInterpolant(x(:),y(:), vx(:)); 
    Fy = scatteredInterpolant(x(:),y(:), vy(:)); 
    rc = xy2rc(Vs, N, rect);
    disp_lddmm = [Fx(rc) Fy(rc)];
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, Vs, disp_lddmm, disp_true);
    fprintf('LDDMM(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', Nlevels(i), meanerr, maxerr, overlapping);
    % load result and evaulate 
    
    Table1(4,i) = meanerr;
    Table1(4,i+2) = maxerr;
    Table1(4,i+4) = overlapping;
    Table1(4,7) = Table1(4,7) + t1;

    
    
    
     
    t0 = tic;
    res = logdemo(ang_s, ang_t);
    for iter = 1:5
        res = logdemo(ang_s, ang_t,  res.vx, res.vy);
        res = logdemo( ecc_s, ecc_t, res.vx, res.vy);
        close all;
    end  
    t1= toc(t0);
    
    [x,y]=meshgrid(1:N,1:N);
    vx = res.vx(:,:, end);
    vy = res.vy(:,:, end);
    Fx = scatteredInterpolant(x(:),y(:), vx(:)); 
    Fy = scatteredInterpolant(x(:),y(:), vy(:)); 
    rc = xy2rc(Vs, N, rect);
    disp_logdemo = [Fx(rc) Fy(rc)];
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, Vs, disp_logdemo, disp_true);
    fprintf('LogDemo(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', Nlevels(i), meanerr, maxerr, overlapping);
    % load result and evaulate 
    
    Table1(3,i) = meanerr;
    Table1(3,i+2) = maxerr;
    Table1(3,i+4) = overlapping;
    Table1(3,7) = Table1(3,7) + t1;
   
      

end

 


%% Proposed
clc
addpath('../../../Proposed')
rng(0)
NoiseLevels = [0.173 0.5473];


[Fte,Vte,Ete]=expand_mesh(Ft,Vt,Et); % expand the mesh to a disk domain, this prepare some space for registrtion
                                     %  just like image generation,which makes the domain rect
                                     
[Fse,Vse,Ese]=expand_mesh(Fs,Vs,Es);                                    


for i = 1:length(NoiseLevels)
    noiseI = NoiseLevels(i);
    %  add noise
    Esnoise = Ese;
    
    Esnoise.Vertex_vis(:,1:2) = Ese.Vertex_vis(:,1:2) + noiseI*rand(size(Ese.Vertex_vis(:,1:2)));
    
      
    Vis_s = [Ese.Vertex_vis(:,2).*cos(Ese.Vertex_vis(:,1)) Ese.Vertex_vis(:,2).*sin(Ese.Vertex_vis(:,1))];
    Vis_snoisy = [Esnoise.Vertex_vis(:,2).*cos(Esnoise.Vertex_vis(:,1)) Esnoise.Vertex_vis(:,2).*sin(Esnoise.Vertex_vis(:,1))];
    % weight
    Esnoise.Weight = 1./vecnorm(Vis_s' - Vis_snoisy')'+ 10; 
    Esnoise.Weight(Ese.n+1:end) = 0; 
    
    ptsrc = load(sprintf('subject_landmark_noise%d.txt',i));
    ptdst = load('template_landmark.txt');
    
    % qchr run registration outsize
    uv_s = disk_harmonic_map(Fse, Vse);
    uv_t = disk_harmonic_map(Fte, Vte);
    t0 = tic;
%     moving.F = Fse;
%     moving.V = uv_s;
%     moving.E = Esnoise;
%     static.F = Fte; 
%     static.V = uv_t; 
%     static.E = Ete; 
%     
%     
%     landmark = vertex_search(ptsrc, moving.V);
%     dst_id = vertex_search(ptdst,static.V); 
%     target = uv_t(dst_id,:);
%     [map,map_mu] = QCWIR_4synthetic(moving,static, landmark, target);
%      
   

    disp_proposed = proposed_reg(Fse, Vse, Esnoise, Fte, Vte, Ete, ptsrc, ptdst);
    t1= toc(t0);
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, Vs, disp_proposed(1:Ese.n,:), disp_true);
    fprintf('Proposed_reg(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', Nlevels(i), meanerr, maxerr, overlapping);
    
    Table1(6,i) = meanerr;
    Table1(6,i+2) = maxerr;
    Table1(6,i+4) = overlapping;
    Table1(6,7) = Table1(6,7) + t1;
end



%% display table


disp_table(Table1)




%%
% 
function value2image(F,V, vis, rect, fn)
N=100;
minvis = min(vis);
maxvis = max(vis);
% bd = compute_bd(F);
% bd =[ bd; bd(1)];

vis = (vis -minvis)/ (maxvis - minvis);
[xx,yy]=meshgrid(linspace(rect(1), rect(2),N), linspace(rect(3), rect(4),N));
F = scatteredInterpolant(V(:,1), V(:,2), vis);
zz = F(xx(:),yy(:));
% xa= xx(:);
% ya =yy(:);
% inp = inpolygon(xa,ya, V(bd,1), V(bd,2));
% zz(~inp) = 0;
img=reshape(zz, 100,100);
imwrite(img,fn);
end


function rc = xy2rc(xy, N, rect)
spanx = rect(2)-rect(1);
spany = rect(4)-rect(3);
rc = [(xy(:,1) - rect(1))/spanx*N (xy(:,2)- rect(3))/spany*N];
end

function xy = rc2xy(rc, N, rect)
spanx = rect(2)-rect(1);
spany = rect(4)-rect(3);
xy = [rc(:,1)/N*spanx + rect(1)  rc(:,2)/N*spany + rect(3)];
end


% Evaulate the registration metric
% Input: F is the subject face list, Nfx3
%        V: subject vertex list Nvx2 
%        disp: The registration result Nvx2
%        disp_true: The ground truth result   Nvx2
% Output: error of displacement 
% overlapping is calculated by counting the overlapping triangles of mesh F and V' = V+disp 
% There is no overlapping on ground truth.
function [meanerr, maxerr, overlapping]=evaulate_metric(F, V, disp, disp_true)
err = vectorNorm(disp-disp_true);
meanerr = nanmean(err(:));
maxerr = nanmax(err(:));
mu =compute_bc(F, V+disp_true, V+disp);
mu(isnan(mu))=[];
overlapping = length(find(abs(mu)>1.001));
end





function disp_table(table)

for i = 1:6    
   fprintf('%.2f/%.2f\t%.2f/%.2f\t%d/%d\t%.1f\n', table(i,1),  table(i,2),  table(i,3), table(i,4), table(i,5), table(i,6), table(i,7)) 
end
end


 