clear;clc;close all;
% gen_syntehtic_Data()
[Ft,Vt,Et]=read_mfile('template.gslm');

%%
figure(101)
subplot(1,5,1)
h = gca;

plot_surf(Ft,Vt(:,1:2), Et.Vertex_vis(:,1))
colormap('hsv')
alpha(0.5);
hold on;
overlay_levelset(h,Ft,Vt(:,1:2),Et.Vertex_vis)
xlim([-6,-2]);
ylim([-2.5 2.5])
axis off



% click landmark 
if ~exist('template_landmark.txt')
    xy = ginput(21);
    save('template_landmark.txt','xy', '-ascii'); 
end
xy = load('template_landmark.txt');
plot(xy(:,1), xy(:,2), 'bo','MarkerSize',6,'Linewidth',2)
for i=1:size(xy,1)
    text(xy(i,1)+0.05, xy(i,2), num2str(i),'color','b','FontSize',16) ;
end





%%
figure(101)
subplot(1,5,2)
[Fs,Vs,Es]=read_mfile('subject.gslm');
 
h = gca;
plot_surf(Fs,Vs(:,1:2), Es.Vertex_vis(:,1))
colormap('hsv')
alpha(0.5);
hold on;
overlay_levelset(h,Fs,Vs(:,1:2),Es.Vertex_vis)

xlim([-6,-2]);
ylim([-2.5 2.5])
axis off



% click landmark 
if ~exist('subject_landmark.txt')
    xy = ginput(21);
    save('subject_landmark.txt','xy', '-ascii');
    
end
xy = load('subject_landmark.txt');
plot(xy(:,1), xy(:,2), 'bo','MarkerSize',6,'Linewidth',2)
for i=1:size(xy,1)
    text(xy(i,1)+0.05, xy(i,2), num2str(i),'color','b','FontSize',16) ;
end






%% add small noise 
rng(0)
Enoise = Es.Vertex_vis;
Enoise(:,1:2) = Enoise(:,1:2) + 0.173*rand(size(Enoise(:,1:2)));

psnr(Es.Vertex_vis(:,1:2), Enoise(:,1:2))
% Enoise(:,1:2) = Enoise(:,1:2) + 0.1*rand(size(Enoise(:,1:2))).*Es.Vertex_vis(:,2);

figure(101)
subplot(1,5,3)
h = gca;
plot_surf(Fs,Vs(:,1:2), Enoise(:,1))
colormap('hsv')
alpha(0.5);
hold on;
overlay_levelset(h,Fs,Vs(:,1:2),Enoise)

xlim([-6,-2]);
ylim([-2.5 2.5])
axis off

% click landmark 
if ~exist('subject_landmark_noise1.txt')
    xy = ginput(21);
    save('subject_landmark_noise1.txt','xy', '-ascii');
    
end
xy = load('subject_landmark_noise1.txt');
plot(xy(:,1), xy(:,2), 'bo','MarkerSize',6,'Linewidth',2)
for i=1:size(xy,1)
    text(xy(i,1)+0.05, xy(i,2), num2str(i),'color','b','FontSize',16) ;
end


%%  add strong noise 
rng(0)
Enoise = Es.Vertex_vis;
Enoise(:,1:2) = Enoise(:,1:2) + 0.5473*rand(size(Enoise(:,1:2)));


psnr(Es.Vertex_vis(:,1:2), Enoise(:,1:2))
% Enoise(:,1:2) = Enoise(:,1:2) + 0.1*rand(size(Enoise(:,1:2))).*Es.Vertex_vis(:,2);

figure(101)
subplot(1,5,4)
h = gca;
plot_surf(Fs,Vs(:,1:2), Enoise(:,1))
colormap('hsv')
alpha(0.5);
hold on;
overlay_levelset(h,Fs,Vs(:,1:2),Enoise)

xlim([-6,-2]);
ylim([-2.5 2.5])
axis off


% click landmark 
if ~exist('subject_landmark_noise2.txt')
    xy = ginput(21);
    save('subject_landmark_noise2.txt','xy', '-ascii');
    
end
xy = load('subject_landmark_noise2.txt');
plot(xy(:,1), xy(:,2), 'bo','MarkerSize',6,'Linewidth',2)
for i=1:size(xy,1)
    text(xy(i,1)+0.05, xy(i,2), num2str(i),'color','b','FontSize',16) ;
end



%% show the ground truth displacement 
figure(101)
subplot(1,5,5)
h = gca;
plot_surf(Fs,Vs(:,1:2), Es.Vertex_vis(:,1))
colormap('hsv')
alpha(0.5);
hold on;
% overlay_levelset(h,Fs,Vs(:,1:2),Es.Vertex_vis)

xlim([-6,-2]);
ylim([-2.5 2.5])
axis off
 

roii =[];
for x = -6:0.3:2
    for y = -2.5:0.3:2.5
         
        d = sqrt((Vs(:,1)-x).^2 +(Vs(:,2)-y).^2 );
        [minval, id]=min(d);
        if (minval<0.05)
            roii =[roii; id];
        end
    end
end
quiver(Vs(roii,1), Vs(roii,2), Es.Vertex_targe(roii,1)-Vs(roii,1) ,  Es.Vertex_targe(roii,2)-Vs(roii,2),1.0, 'Linewidth',1.2,'color',[0 0 0], 'MaxHeadSize',5)
 


 