function [F,V,E]=expand_mesh(Ft,Vt,Et)
 
Vt = Vt(:,1:2); 
bd = compute_bd(Ft);
[~, vdisk]=unit_disk_mesh(5);
vdiskbig = vdisk * 4 -[3.5 0];

% use only points outside boundary
[in, on] = inpolygon( vdiskbig(:,1), vdiskbig(:,2), Vt(bd,1),Vt(bd,2));
vdiskbig((in|on) ,:) =[];

for i = 1:length(bd)
   dis = vectorNorm( Vt(bd(i),:)  - vdiskbig);
   vdiskbig(dis<0.08,:) =[];
end
%


pt = [Vt;vdiskbig];

F = delaunay(pt); 
[in1,on1] = inpolygon(pt(F(:,1),1), pt(F(:,1),2), Vt(bd,1),Vt(bd,2));
[in2,on2] = inpolygon(pt(F(:,2),1), pt(F(:,2),2), Vt(bd,1),Vt(bd,2));
[in3,on3] = inpolygon(pt(F(:,3),1), pt(F(:,3),2), Vt(bd,1),Vt(bd,2));

pc = (pt(F(:,1),:) +  pt(F(:,2),:) +  pt(F(:,3),:))/3;
in = inpolygon(pc(:,1),pc(:,2), Vt(bd,1),Vt(bd,2));
 
inflag = (in1&~on1) | (in2&~on2) | (in3&~on3) | in;

F(inflag,:)=[];
F = [F; Ft];
V = pt;
E = Et;

N= size(V,1);
n = size(Vt,1);
E.Vertex_vis = zeros(N,3);
E.Vertex_vis(1:n,:) = Et.Vertex_vis;
E.Vertex_targe = zeros(N,2);
E.Vertex_targe(1:n,:) = Et.Vertex_targe;
E.n = n;
E.N =N;
  
 
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
    Vn = disk_harmonic_map(Fn,V);
    
     save(fn, 'Fn','Vn');
end
end
