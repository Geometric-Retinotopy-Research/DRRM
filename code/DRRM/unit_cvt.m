%% 
% Description  --  function v=unit_cvt(v0,type)
%		utilities to convert units.			
% 
% Parameter(s): 
%     v0[double array]		-- raw data
%     type[string]          -- type of convertation, either c2p: x-y
%     to polar angle; or p2c polar angle to x-y
% 
% Return:
%     v[double array]		--  data from convert
% 
function v=unit_cvt(v0,type)
if strcmp(type, 'p2c')
    v = [v0(:,1).*cos(v0(:,2)/180*pi)  v0(:,1).*sin(v0(:,2)/180*pi)];
elseif strcmp(type, 'c2p')
    v = [sqrt(v0(:,1).^2 + v0(:,2).^2), atan2(v0(:,2), v0(:,1))*180/pi ];
%     v(v0(:,1),:)=NaN;
else
    error('type can only be p2c or c2p')
end

end