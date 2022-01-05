function displace = tpswarp2d(Vs, ptsrc, ptdst)

% Thin-plate spline mesh warping.
%  Dis = tpswarp2d(Vs, ps, pd)
%   input:
%       Vs: subject vertex mesh
%       ptsrc: 2d source landmark position n x 2
%       ptdst: 2d destin landmark position n x 2
%
%   output:
%       dis  : output displacement
%   Bookstein, F. L.
%   "Principal Warps: Thin Plate Splines and the Decomposition of Deformations."
%   IEEE Trans. Pattern Anal. Mach. Intell. 11, 567-585, 1989.
%
%   2019/9/26
%   Yanshuai
%   yanshuai@asu.edu.com

Ncorr = size(ptdst,1); %
% use a matrix K that evaluates the function U(rij)
K=zeros(Ncorr,Ncorr); % initilize the K matrix
for i=1:Ncorr
    for j = 1:Ncorr
        rij = norm(ptsrc(i,:) - ptsrc(j,:));
        K(i,j) = TPS_U(rij);
    end
end

P = [ones(Ncorr,1),ptsrc];  %The control point positions
L = [K,P;P',zeros(3,3)];
Y = [ptdst - ptsrc;zeros(3,2)]; %displacement of x, y, (it is independent)

% w = (W | a1 a2 a3)
w = L\Y; % will yield 2xm for x, y respectively

% Apply tranform: use w evaulate  displacement
displace = zeros(size(Vs,1), 2);
for dim = 1:2
    a1 = w(end-2,dim);
    a2 = w(end-1,dim);
    a3 = w(end,dim);    
    displace(:,dim) = a1 + a2* Vs(:,1) + a3*Vs(:,2);
    for i = 1:length(w) - 3
        displace(:,dim) =w(i, dim)*TPS_U( vecnorm(ptsrc(i,:)' - Vs') )';
    end
end
end

function k = TPS_U(ri)
rnonzero = ri;
rnonzero((ri==0))=realmin;
k = ri.*log(rnonzero);
end



