%%
% Description  -- function [map,map_mu]  = disk_registration(static, moving,landmark, target, maxiter, maxfixstep)
%		register the subject of disk.
%
% Parameter(s): 
%		static[struct]         -- the static mesh structure 
%		moving[struct]         -- the moving mesh structure
%		landmark[double array] -- landmark including boundary index
%		target[double array]   -- landmark target 
%		maxiter[int]           -- max value for loop
%       maxfixstep[int]        -- step for loop iterate
% Return: 
%		result[GenFigureBase]  -- GenFigureBase object
% 
%%
function [map,map_mu]  = disk_registration(static, moving,landmark, target, maxiter, maxfixstep)
if ~exist('maxiter','var')
    maxiter = 3;
end

if ~exist('maxfixstep','var')
    maxfixstep = 200;
end

bd = compute_bd(moving.F);
landmark_andbd = [landmark; bd];
target_andbd = [target; moving.uv(bd,:)];
[map,map_mu] = registration(moving,static, landmark_andbd, target_andbd, maxiter, maxfixstep);
end



function [map,map_mu] = registration(moving,static, landmark, target,maxiter, maxfixstep)
 
% Inputs
% moving : the moving mesh structure
% statuc : the static mesh structure
% landmark: Specify the landmark including boundary index
% target: landkmark target
%
% Output
% map : registration result
% map_mu : corresponding Beltrami coefficients
 
% initial map
 
face = moving.F;
vertex = moving.uv;
  
%% Use T-map to get an initial registration
[map,map_mu] = techimuller_map(face, vertex, landmark, target);
id = find(isnan(map(:,1)));
map(id,:) = vertex(id,:);
moving.uv = map;
%% refinement
stepsize =0.1;
for i = 1:maxiter
    [map, map_mu] = registration_rm(moving, static, landmark, target, stepsize, maxfixstep); % Intensity match registration
    moving.uv = map; % update moving.uv
    stepsize = stepsize*0.8;
%     moving.F = delaunay(moving.uv);
end

end

function [map,map_mu] = registration_rm(moving,static,landmark,target, stepsize, maxfixstep)
 
% registration mapping by both retinotopic coordinates and landmark constraints.
%
face = moving.F;
vertex = moving.uv;

% updated_map = nearby_search(moving,static);
 updated_map = simple_demon(moving,static, stepsize);
 
% smooth 
[smooth_disp,lambdas] = grid_laplacian_smooth(vertex, updated_map - vertex, moving.R2);
% smooth_disp(moving.R2<1,:)=0;
 
% adjust to ensure diffeomorphic
diff_map = vertex + smooth_disp; 
for t=1:maxfixstep
    mu = compute_bc(face,vertex, diff_map);
    if max(abs(mu))<1.01
        break
    else
        fprintf(' max(abs(mu))=%f\n', max(abs(mu)));
    end
%      uv_face = (diff_map(face(:,1),:)+diff_map(face(:,2),:)+diff_map(face(:,3),:))/3;
%      mu_smooth = laplacian_smooth(uv_face,real(mu)) + 1i*laplacian_smooth(uv_face,imag(mu));
    choped = mu_chop(mu,0.8, 0.8);        
    diff_map = linear_beltrami_solver(face,vertex,choped,landmark, target);
    figure(201)
    plot_mesh(face,diff_map);
    drawnow
    
end
map = diff_map;
id = isnan(map(:,1));
map(id,:) = vertex(id,:);

map_mu = mu;
end
 
   
function map = simple_demon(moving, static, stepsize)
 

alpha=2;
Tx=zeros(size(moving.uv,1),1); Ty=zeros(size(moving.uv,1),1);


Mvis = moving.feature; 
Svis = static.feature; 

Mvis(isnan(Mvis))=0;
Svis(isnan(Svis))=0;
dimfeature = size(moving.feature,2);
for i =1:dimfeature
    Sf{i} = scatteredInterpolant(static.uv(:,1), static.uv(:,2), Svis(:,i), 'natural');
end

delta = 0.01;  dx = [delta, 0]; dy = [0, delta];
 
 
map = moving.uv;
for itt=1 : 2
    %
    % update according to vis_x vis_y, accordingly
    % M is function of (u,v)
    for i=1:dimfeature
        xy = map + [Tx Ty];        

        Mf =  scatteredInterpolant(xy, Mvis(:,i)  , 'natural');
        
        M = Mf(xy) ;    S = Sf{i}(xy);     Idiff = (M - S);
        
        fprintf('demon mse %d, %f\n', itt,mse(Idiff));
        
        Mx = (Mf(xy +  dx) - M)/delta;     My = (Mf(xy + dy) - M)/delta;
        Sx = (Sf{i}(xy + dx) - S)/delta;    Sy = (Sf{i}(xy + dy) - S)/delta;
        
        Ux =  Idiff.*  ((Sx./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2)) +(Mx./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        Uy =  Idiff.*  ((Sy./((Sx.^2+Sy.^2)+alpha^2*Idiff.^2)) +(My./((Mx.^2+My.^2)+alpha^2*Idiff.^2)));
        
        Ux(isnan(Ux))=0; Uy(isnan(Uy))=0;
        Ux(isinf(Ux))=0; Uy(isinf(Uy))=0;
        % ignore the update beyond disk radio > 0.95
        id = find(vecnorm(xy')>0.95 );
        Ux(id)=0; Uy(id)=0;
        Tx=Tx+stepsize*Ux;    Ty=Ty+stepsize*Uy;
        Tx(isnan(Tx))=0; Ty(isnan(Ty))=0;

    end
     
end

map = xy;

end



function mu_new = mu_chop(mu,bound,constant)

mu_new = mu;
ind = abs(mu)>=bound;
mu_new(ind) = constant*(cos(angle(mu(ind)))+1i*sin(angle(mu(ind))));
mu_new(isnan(mu_new)) = 0;
mu_new(isinf(mu_new)) = 0;

end

      

function [uv_new, mu ]= techimuller_map(face, uv, landmark, target) 
uv_new = uv;
for i=1:20000
    % update target
    uv_new(landmark,:) = target;
    
    % compute bc
    mu = compute_bc(face,uv,uv_new);
    
    % chop mu
    id = find(abs(mu)>0.8);
    mu(id) = mu(id)./abs(mu(id))*0.8;
    
    % reduce mu
    meanmu =  nanmean(abs(mu));
    mu =meanmu*mu./abs(mu);
    mu(isnan(mu))=nanmean(mu);
    % update uv
    [uv_new,mu] = linear_beltrami_solver(face,uv,mu,landmark, target);      
    
    if max(abs(mu))<0.8
        break
    else
        fprintf('max(abs(mu))=%f\n',nanmax(abs(mu)));
    end
   
end

end
