%% Reproduce the Figure 1g

clear;clc;close all;
subjects = dir([GenConsts.kDataMeshURL,'*lh.m']);
%%
subi =1;
% Prepare data
fn = subjects(subi).name;

% Load the full hemisphere from mat file 
[Ffull,Vfull, Efull]=read_mfile([GenConsts.kDataMeshURL fn]);

% find the distance & final state & index of cloest points.
% filter points which distance great than 60.
foveaid = 43052;
[D,S,Q] = perform_fast_marching_mesh(double(Vfull), double(Ffull), foveaid);
ind2del  = find(D > 60);

% get a list of vertex id to delete
[Fout,Vout,father] = gf_remove_mesh_vertices(Ffull,Vfull, ind2del); 

%% plot figure
% color with ecc type
param.face = Ffull;
param.vertex = Vfull;
param.EValue = Efull;

param.renderColor.valueType = 'ecc';
param.renderColor.value = Efull.Vertex_prf(:,2);
param.renderColor.maxValue = 8;
figure1 = GenFigure_g(param);
figure1.draw();

% color with lh type
param.renderColor.valueType = 'lh';
param.renderColor.value = Efull.Vertex_prf(:,1);
param.renderColor.maxValue = -1;
figure2 = GenFigure_g(param);
figure2.draw();



