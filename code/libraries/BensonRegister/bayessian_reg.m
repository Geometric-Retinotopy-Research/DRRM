function displacement = bayessian_reg(Fs, Vs, ecc_ang_s, Vt, ecc_ang_t, ptsrc, ptdst)

% Register by calling nben library
% We use  weight=1 for non-landmark, and weight = 1e3 for landmark
% Please refer  N. C. Benson and J. Winawer, “Bayesian analysis of retinotopic maps,” Elife, vol. 7, Dec. 2018.

% landmark will use very big weight and  other vertices use small wegiht as anchors

% When searching anchors, we use landmark to interprete a initial wrap,
% then find the corrosponding anchor~

import nben.mesh.registration.* 

field = nben.mesh.registration.Fields();
faces = int32(Fs'-1); 
vertices = Vs(:,1:2)';

potential = field.newSum();
pot_edge = field.newHarmonicEdgePotential(faces, vertices);
pot_ang = field.newWellAnglePotential(faces, vertices);
pot_peri = field.newHarmonicPerimeterPotential(faces, vertices);


% anchor is the most important potential 
weight = 1e3*ones(size(Vs,1),1); % weigth =scale
ldmk_ids = find_nearest_id(Vs, ptsrc);
weight(ldmk_ids) = 1e6;
sig = 5*ones(size(Vs,1),1); 
[anchorids, anchorpos] = cal_anchors(Vs,ecc_ang_s, Vt, ecc_ang_t, ptsrc, ptdst);
pot_anchor = field.newGaussianAnchorPotential(weight, sig, 2, int32(anchorids-1), anchorpos', vertices);
potential.addField(pot_edge)
% potential.addField(pot_ang)
potential.addField(pot_peri)
potential.addField(pot_anchor)
 
minimizer = Minimizer(potential, vertices);
max_step_size = 0.001;
max_steps = int32(1000);
max_pe_change = 1;
minimizer.randomStep(max_pe_change, max_steps, max_step_size, logical(false));
displacement = minimizer.getX()' - Vs;


end



% find a empire registration by landmarks + value search
function [anchorids, anchorpos, target0] = cal_anchors(Vs,ecc_ang_s, Vt, ecc_ang_t, ptsrc, ptdst)

vis_s = [ecc_ang_s(:,2).*cos(ecc_ang_s(:,1)) ecc_ang_s(:,2).*sin(ecc_ang_s(:,1))  ];
vis_t = [ecc_ang_t(:,2).*cos(ecc_ang_t(:,1)) ecc_ang_t(:,2).*sin(ecc_ang_t(:,1))  ];

displace0 = tpswarp2d(Vs, ptsrc, ptdst);
target0 = Vs + displace0;
target = target0;

ldmk_ids = find_nearest_id(Vs, ptsrc);

% for each vertex, find the anchor target 
isanchor = zeros(size(Vs,1),1);
isanchor(ldmk_ids) =1;
for i=1:1000:size(Vs,1)
    rdis = vectorNorm(vis_t - vis_s(i,:)).* ((vectorNorm(Vt - Vs(i,:))<0.4)+1);    
    [minvis,j]=min(rdis);        
    if(minvis<0.05)
       isanchor(i) = 1;
    end
    target(i,:) = Vt(j,:);
end
% restore the landmark 
target(ldmk_ids,:) = target0(ldmk_ids,:);

% select anchor only 
anchorids = find(isanchor);
anchorpos = target(anchorids,:);
end


function ids = find_nearest_id( V, qv)

nq = size(qv,1);
ids= zeros(nq,1);
for i = 1:nq
    [~, id]=min(vectorNorm(V - qv(i,:)));
    ids(i) = id;
end

end





