function plot_subject(sub)
figure;  subplot(121)
hold on;
plot_surf(sub.lh.very_inflated.faces, sub.lh.very_inflated.vertices, prf_value_2_color('ecc', sub.lh.pRF(:,2),8))
plot_surf(sub.rh.very_inflated.faces, sub.rh.very_inflated.vertices, prf_value_2_color('ecc', sub.rh.pRF(:,2),8))
if isfield(sub.lh, 'fovid')
    plot3(sub.lh.very_inflated.vertices(sub.lh.fovid,1),sub.lh.very_inflated.vertices(sub.lh.fovid,2), sub.lh.very_inflated.vertices(sub.lh.fovid,3), 'ko','Markersize',4)
    plot3(sub.rh.very_inflated.vertices(sub.rh.fovid,1),sub.rh.very_inflated.vertices(sub.rh.fovid,2), sub.rh.very_inflated.vertices(sub.rh.fovid,3), 'ko','Markersize',4)
end

subplot(122)
hold on;
plot_surf(sub.lh.very_inflated.faces, sub.lh.very_inflated.vertices, prf_value_2_color('lh', sub.lh.pRF(:,1)))
plot_surf(sub.rh.very_inflated.faces, sub.rh.very_inflated.vertices, prf_value_2_color('rh', sub.rh.pRF(:,1)))
if isfield(sub.lh, 'fovid')
    plot3(sub.lh.very_inflated.vertices(sub.lh.fovid,1),sub.lh.very_inflated.vertices(sub.lh.fovid,2), sub.lh.very_inflated.vertices(sub.lh.fovid,3), 'ko','Markersize',4)
    plot3(sub.rh.very_inflated.vertices(sub.rh.fovid,1),sub.rh.very_inflated.vertices(sub.rh.fovid,2), sub.rh.very_inflated.vertices(sub.rh.fovid,3), 'ko','Markersize',4)
end

end