%% Reproduce Table 1
% Table 1. Comparing registration performance relative to the ground truth. 
% Each cell has two values, for the low and high fMRI noise conditions, respectively. 
% Landmarks (the circled positions in Fig. 4) were used if the method accepts them (marked with "*").

% The DataBus stores every evaulated metric; we display at the end


clear;clc;close all;
rng(0)

% gen_syntehtic_Data()
[Ft,Vt,Et]=read_mfile([GenConsts.kGSLMURL, 'template.gslm']);
[Fs,Vs,Es]=read_mfile([GenConsts.kGSLMURL, 'subject.gslm']);

Vt = Vt(:,1:2); Vs = Vs(:,1:2);
disp_true = Es.Vertex_targe - Vs;

Table1 = zeros(4,7); 
NoiseLevels = [0.173 0.5473];
    

rect = [-1 1 -1 1];

global N
N = 100; % image resolution NxN

%% expand to a disk
[Fse,Vse,Ese]=expand_mesh(Fs,Vs,Es); % subject
[Fte,Vte,Ete]=expand_mesh(Ft,Vt,Et); % template

% ensure to be a disk
uv_s = disk_harmonic_map(Fse, Vse);
uv_t = disk_harmonic_map(Fte, Vte);

n_before_expand = Ese.n;
%% TPS: It only use landmark, as well as we find boundary
for noisei = 1:length(NoiseLevels)
    
    %% prepare landmarks    % load  Human  landmark.
    ptsrc_0 = load(sprintf('%ssubject_landmark_noise%d.txt', GenConsts.kTXTURL, noisei));
    ptdst_0 = load([GenConsts.kTXTURL, 'template_landmark.txt']);   
    
    landmark = vertex_search(ptsrc_0, Vs);
    dst_id = vertex_search(ptdst_0,Vt);
    target = uv_t(dst_id,:);
    
    ptsrc = uv_s(landmark, :);
    ptdst = uv_t(dst_id, :);

 
    %% TPS
    
    bd_s = compute_bd(Fse); 
     
    bd_src = uv_s(bd_s,:);
    bd_target = bd_src;
     Esenoise = Ese.Vertex_vis;
     noiseI = NoiseLevels(noisei);
    Esenoise(:,1:2) = Esenoise(:,1:2) + noiseI*rand(size(Esenoise(:,1:2)));
    
    t0 = tic;
    
    disp_tps = tpswarp2d(uv_s, [ptsrc; bd_src], [ptdst; bd_target]);
    t1 =toc(t0);
    
    
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, uv_s(1:n_before_expand,:), disp_tps(1:n_before_expand,:), disp_true);
    signalpsnr = psnr(Ese.Vertex_vis(:,1:2), Esenoise(:,1:2));  
    
    fprintf('TPS(PSNR=%f): mean err = %f max error =%f  overlap = %d\n',signalpsnr, meanerr, maxerr, overlapping);
    Table1(1,noisei) = meanerr;
    Table1(1,noisei+2) = maxerr;
    Table1(1,noisei+4) = overlapping;
    Table1(1,7) = Table1(1,7) + t1;
    
     mu = compute_bc(Fs, uv_s(1:n_before_expand,:), uv_s(1:n_before_expand,:)+disp_tps(1:n_before_expand,:));
    Table1(1,noisei+7) = mean(abs(mu));
    
    %% Bayessian
    noiseI = NoiseLevels(noisei);
    
    %  add noise
    Esenoise = Ese.Vertex_vis;
    Esenoise(:,1:2) = Esenoise(:,1:2) + noiseI*rand(size(Esenoise(:,1:2)));
    signalpsnr = psnr(Ese.Vertex_vis(:,1:2), Esenoise(:,1:2));  
    
 
    Vis_t = [Ete.Vertex_vis(:,2).*cos(Ete.Vertex_vis(:,1)) Ete.Vertex_vis(:,2).*sin(Ete.Vertex_vis(:,1))];
    
    bd_s = compute_bd(Fse); 
     
    bd_src = Vse(bd_s,:);
    bd_target = bd_src;
    
    
    %% Bayessian run registration outsize
    t0 = tic;
    disp_bayessian = bayessian_reg(Fse, uv_s, Esenoise, uv_t, Ete.Vertex_vis, [ptsrc; bd_src], [ptdst; bd_target]);
    t1 = toc(t0);
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, uv_s(1:n_before_expand,:), disp_bayessian(1:n_before_expand,:), disp_true);
    fprintf('Bayessian(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', signalpsnr, meanerr, maxerr, overlapping);
    
    
    Table1(2,noisei) = meanerr;
    Table1(2,noisei+2) = maxerr;
    Table1(2,noisei+4) = overlapping;
    Table1(2,7) = Table1(2,7) + t1; 
    
        
    mu = compute_bc(Fs, uv_s(1:n_before_expand,:), uv_s(1:n_before_expand,:)+disp_bayessian(1:n_before_expand,:));
    Table1(2,noisei+7) = mean(abs(mu));
    %% DDemons
    
    t0 = tic;     
    
    %  add noise
    noiseI = NoiseLevels(noisei);
    Esenoise = Ese.Vertex_vis;
    Esenoise(:,1:2) = Esenoise(:,1:2) + noiseI*rand(size(Esenoise(:,1:2)));
    
    
    
    ang_t_fn = [GenConsts.kImageURL,'ang_t.png'];
    ecc_t_fn = [GenConsts.kImageURL,'ecc_t.png'];
    
    ang_s_fn = sprintf('%sang_s%d.png',GenConsts.kImageURL,noisei);
    ecc_s_fn = sprintf('%secc_s%d.png',GenConsts.kImageURL,noisei);
    
    value2image(Fte, uv_t, Ete.Vertex_vis(:,1), rect, ang_t_fn);
    value2image(Fte, uv_t, Ete.Vertex_vis(:,2), rect, ecc_t_fn);
    
    value2image(Fse, uv_s, Esenoise(:,1), rect, ang_s_fn);
    value2image(Fse, uv_s, Esenoise(:,2), rect, ecc_s_fn);
    
    ecc_t = im2double(imread(ecc_t_fn));
    ang_t = im2double(imread(ang_t_fn));
    ecc_s = im2double(imread(ecc_s_fn));
    ang_s = im2double(imread(ang_s_fn));
    
    
    res = logdemo(ang_s, ang_t);
    for iter = 1:200
        res = logdemo(ang_s, ang_t, res.vx, res.vy);
        res = logdemo(ecc_s, ecc_t, res.vx, res.vy);
    end
    t1= toc(t0);
    
    [x,y]=meshgrid(1:N,1:N);
    vx = res.vx(:,:, end);
    vy = res.vy(:,:, end);
    Fx = scatteredInterpolant(x(:)/N*2-1,y(:)/N*2-1, vx(:), 'natural');
    Fy = scatteredInterpolant(x(:)/N*2-1,y(:)/N*2-1, vy(:), 'natural');
 
    disp_logdemo = [Fx(uv_s) Fy(uv_s)]/N*2-1;
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, uv_s(1:n_before_expand,:), disp_logdemo(1:n_before_expand,:), disp_true);
    fprintf('LogDemo(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', signalpsnr, meanerr, maxerr, overlapping);
    % load result and evaulate
    
    Table1(3,noisei) = meanerr;
    Table1(3,noisei+2) = maxerr;
    Table1(3,noisei+4) = overlapping;
    Table1(3,7) = Table1(3,7) + t1;
    
    mu = compute_bc(Fs, uv_s(1:n_before_expand,:), uv_s(1:n_before_expand,:)+disp_logdemo(1:n_before_expand,:));
    Table1(3,noisei+7) = mean(abs(mu));
    
    %% Proposed
    
    noiseI = NoiseLevels(noisei);
    %  add noise
    Esenoise = Ese;
    
    Esenoise.Vertex_vis(:,1:2) = Ese.Vertex_vis(:,1:2) + noiseI*rand(size(Ese.Vertex_vis(:,1:2)));
    
    
    Vis_s = [Ese.Vertex_vis(:,2).*cos(Ese.Vertex_vis(:,1)) Ese.Vertex_vis(:,2).*sin(Ese.Vertex_vis(:,1))];
    Vis_snoisy = [Esenoise.Vertex_vis(:,2).*cos(Esenoise.Vertex_vis(:,1)) Esenoise.Vertex_vis(:,2).*sin(Esenoise.Vertex_vis(:,1))];
    % weight
    Esenoise.Weight = 1./vecnorm(Vis_s' - Vis_snoisy')'+ 10;
    Esenoise.Weight(n_before_expand+1:end) = 0;
    
    ptsrc = load(sprintf('%ssubject_landmark_noise%d.txt', GenConsts.kTXTURL, noisei));
    ptdst = load([GenConsts.kTXTURL, 'template_landmark.txt']);   
    

    
    
    moving.F = Fse;
    moving.uv = uv_s;
    moving.feature = Esenoise.Vertex_vis(:,1:2);
    moving.R2 = Esenoise.Weight;
    
    static.F = Fte;
    static.uv = uv_t;
    static.feature = Ete.Vertex_vis(:,1:2);

    
    t0 = tic;
    map =  disk_registration(static, moving, landmark, target, 2);
    disp_proposed = map - moving.uv;
    t1= toc(t0);
    [meanerr, maxerr, overlapping]=evaulate_metric(Fs, uv_s(1:n_before_expand,:), disp_proposed(1:n_before_expand,:), disp_true); % the first n is the interior
    fprintf('Proposed_reg(PSNR=%f): mean err = %f max error =%f  overlap = %d\n', signalpsnr, meanerr, maxerr, overlapping);
    
    Table1(4,noisei) = meanerr;
    Table1(4,noisei+2) = maxerr;
    Table1(4,noisei+4) = overlapping;
    Table1(4,7) = Table1(4,7) + t1;
    
    mu = compute_bc(Fs, uv_s(1:n_before_expand,:), uv_s(1:n_before_expand,:)+disp_proposed(1:n_before_expand,:));
    Table1(4,noisei+7) = mean(abs(mu));
    
end


%% display table

disp_table(Table1)
function disp_table(table)

    for i = 1:size(table,1)
        fprintf('%.3f/%.3f\t%.3f/%.3f\t%d/%d\t%.1f\n', table(i,1), ...
      table(i,2),  table(i,3), table(i,4), table(i,5), table(i,6), table(i,7)) 
        % mean(20db/10db) max(20db/10db) flip(20db/10db) time(20db/10db)
    end
end


%% Evaulate the registration metric
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
    mu =compute_bc(F, V, V+disp);
    mu(isnan(mu))=[];
    overlapping = length(find(abs(mu)>1.001));
end


