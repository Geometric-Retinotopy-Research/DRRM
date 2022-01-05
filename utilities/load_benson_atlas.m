function atlas = load_benson_atlas(lr)


fsavgfn = sprintf('bensondata/fsaverage/surf/%s.benson14_retinotopy.v4_0.sphere.reg',lr);
 
[V, faces] = freesurfer_read_surf(fsavgfn);

eccfn = sprintf('bensondata/fsaverage/surf/%s.benson14_eccen.v4_0.mgz',lr);
angfn = sprintf('bensondata/fsaverage/surf/%s.benson14_angle.v4_0.mgz',lr);
sigfn = sprintf('bensondata/fsaverage/surf/%s.benson14_sigma.v4_0.mgz',lr);
varfn = sprintf('bensondata/fsaverage/surf/%s.benson14_varea.v4_0.mgz',lr);

ecc   = squeeze(fs_load_mgh(eccfn));
ang   = squeeze(fs_load_mgh(angfn));
if lr=='lh'
    ang = 90-ang;
else
    ang = ang -270; 
end
ang(ang<0) = ang(ang<0) +360;
ang(ang>360) = ang(ang>360) -360;

sigma = squeeze(fs_load_mgh(sigfn));
varea = squeeze(fs_load_mgh(varfn));


atlas.F = faces;
atlas.V = V;
atlas.ecc = ecc;
atlas.ang = ang;
atlas.sigma = sigma;
atlas.varea = varea;

end
