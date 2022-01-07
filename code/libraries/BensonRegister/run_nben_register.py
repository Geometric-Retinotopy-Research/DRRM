import neuropythy as ny
import nibabel                      as     nib
import nibabel.freesurfer.io        as     fsio
import nibabel.freesurfer.mghformat as     fsmgh
import matplotlib.pyplot as plt

def _guess_surf_file(fl):
    # MGH/MGZ files
    try: return np.asarray(fsmgh.load(fl).dataobj).flatten()
    except Exception: pass
    # FreeSurfer Curv files
    try: return fsio.read_morph_data(fl)
    except Exception: pass
    # Nifti files
    try: return np.squeeze(nib.load(fl).dataobj)
    except Exception: pass
    raise ValueError('Could not determine filetype for: %s' % fl)
    
    
    

ny.config['data_cache_root']='data/'
sub = ny.data['benson_winawer_2018'].subjects['S1202']



polar_angle =  _guess_surf_file('data/benson_winawer_2018/retinotopy/S1202/lh_angle.mgz')

eccentricity =  _guess_surf_file('data/benson_winawer_2018/retinotopy/S1202/lh_eccen.mgz')

weight = _guess_surf_file('data/benson_winawer_2018/retinotopy/S1202/lh_vexpl.mgz')

m = ny.vision.register_retinotopy(hemi = sub.lh,  polar_angle = polar_angle, eccentricity = eccentricity, weight = weight)


#def register_retinotopy(hemi,
#                        model='benson17', model_hemi=Ellipsis,
#                        polar_angle=None, eccentricity=None, weight=None, pRF_radius=None,
#                        weight_min=0.1,
#                        eccentricity_range=None,
#                        partial_voluming_correction=False,
#                        radius_weight=1, field_sign_weight=1, invert_rh_field_sign=False,
#                        scale=20.0,
#                        sigma=Ellipsis,
#                        select='close',
#                        prior=None,
#                        resample=Ellipsis,
#                        radius=np.pi/3,
#                        max_steps=2000, max_step_size=0.05, method='random',
#                        yield_imap=False):
    
m.coordinates.T    
m.tess.faces.T