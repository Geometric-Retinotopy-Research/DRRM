
%% 
% Description  --  function fq = sphere_interp(uvw, feature, uvw_q)
%		interpretate scalar features on surface
% Parameter(s): 
%     uvw[double array]        --  sites
%     feature[double array]    --  feature on site
%     uvw_q[double array]      --  query point on the spehre
% Return:
% 	  fq[double array]         -- the result of interpretation
% 
% uvw: sites
% feature : feature on site
% uvw_q: query point on the spehre


function fq = sphere_interp(uvw, feature, uvw_q)
dim = size(feature,2);

for i=1:dim
    idn = uvw(:,3)<1;   ids = uvw(:,3)> -1;

     u = north_prj(uvw(idn,:));
     Fi = scatteredInterpolant(u, feature(idn,i));    
     
    uq = north_prj(uvw_q);    
    bd = boundary(u,0);
    in = inpolygon(uq(:,1), uq(:,2), u(bd,1), u(bd,2));    
    fq(in,i) = Fi(uq(in,:)); 



    u = south_prj(uvw(ids,:));
    Fi = scatteredInterpolant(u, feature(ids,i));    
    uq = south_prj(uvw_q);
    bd = boundary(u,0);
    in = inpolygon(uq(:,1), uq(:,2),u(bd,1), u(bd,2));    
    fq(in,i) = Fi(uq(in,:));
end

end 
 


function u = north_prj(V)
u = double([V(:,1)./(100-V(:,3)) V(:,2)./(100-V(:,3)) ]);
end

function u = south_prj(V)
V(:,3) = -V(:,3);
u = double([V(:,1)./(100-V(:,3)) V(:,2)./(100-V(:,3)) ]);
end

