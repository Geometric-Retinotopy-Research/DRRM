close all
figure;plot_mesh(Fs, uv_s, Es.Vertex_vis(:,1));title('start')
figure;plot_mesh(Fs,map, Es.Vertex_vis(:,1)); title('moving')
figure; plot_mesh(Ft,uv_t, Et.Vertex_vis(:,1)); title('target')


%%
figure; plot_mesh(face,map); title('map')

%%
figure
subplot(1,2,1)
plot_mesh(Fs, uv_s);
hold on;
quiver(uv_s(src_id,1), uv_s(src_id,2), dst_uv(:,1)- uv_s(src_id,1) , dst_uv(:,2)- uv_s(src_id,2),1)

subplot(1,2,2) 
plot_mesh(Ft, uv_t);

%%

figure;plot_mesh(Fs,map0); title('moving'); hold on
quiver(map0(:,1),map0(:,2), map(:,1)-map0(:,1),map(:,2)-map0(:,2),1)