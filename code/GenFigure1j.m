%% Reproduce the Figure 1j

clear;clc;close all;
subjects = dir([GenConsts.kDataMeshURL,'*lh.m']);
%%
subi =1;
% Prepare data
fn = subjects(subi).name;
% Load the full hemisphere
[Ffull,Vfull, Efull]=read_mfile([GenConsts.kDataMeshURL fn]);

% find the distance & final state & index of cloest points.
% filter points which distance great than 60.
foveaid = 43052;
[D,S,Q] = perform_fast_marching_mesh(double(Vfull), double(Ffull), foveaid);
ind2del  = find(D > 60);

% from conformal map get the coordinate -uv
[Fout,Vout,father] = gf_remove_mesh_vertices(Ffull,Vfull, ind2del); % get a list of vertex id to delete
uv = disk_conformal_map(Fout, Vout);


%% plot figur3
areaids = [1 4];

param.face = Fout;
param.vertex = uv;
param.EValue = Efull;
param.father = father;
param.renderColor.valueType = 'lh';
param.renderColor.value = Efull.Vertex_prf(:,1);
figure_j = GenFigure_j(param);
figure_j.draw().sub_plot(areaids).save('Fig1_k'); 
 
 
