%% Description  -- function overlay_levelset(h,F,V,vis)
%		add line on the figure overlay
%
%  Parameter(s): 
%		h[struct]         --  Current axes or chart
% 		F[double array]   --  connectivity of mesh
%		V[double array]   --  vertex of mesh
% 		vis[double array] --  Vertex visual coordiantes
%		
%%
function overlay_levelset(h,F,V,vis)
figure(123)
subplot(121)
[~,Hang]=tricontour(V,F,vis(:,1),5); 
subplot(122)
[~,Hecc]=tricontour(V,F,vis(:,2),[0:0.2:sqrt(8)].^2);
 
for i = 1:length(Hecc)
   Hi = Hecc(i) ;
   plot(h,Hi.Vertices(:,1), Hi.Vertices(:,2), 'r-','Linewidth',1);
end
 
for i = 1:length(Hang)
   Hi = Hang(i) ;
   plot(h,Hi.Vertices(:,1), Hi.Vertices(:,2), 'k-','Linewidth',1);
end
close(123);
end
