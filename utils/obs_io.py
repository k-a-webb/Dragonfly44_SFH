import numpy as np

def read_file_data( file_data=None, **extras ):
    import h5py
    data_dict = {}
    with h5py.File(file_data, 'r') as hfile:

        for key in ['wavelength', "Redshift", "mass", "sigma_smooth", "e_sigma_smooth","inres",  "maggies","spectrum","maggies_unc","unc"]:
            if key not in hfile.keys():
       	       	data_dict[key] = None
       	    else:
                v = np.copy( hfile[key] )
                if v.size>1:
                    data_dict[key] = v
                else:
                    data_dict[key] = v.item()

        for key in ['filternames']: # special case for string type data
            data_dict[key] = np.array( np.copy( hfile[key] ), dtype=str )

        key = "smoothtype"  # special case for string type data
        if key in hfile.keys(): data_dict[key] = np.copy( hfile[key] ).astype(str).item()
        else: data_dict[key] = None

        # copy mask for spectroscopy, if one is specified
        try: data_dict['mask'] = np.copy( hfile['mask'] )
        except: pass

    return data_dict

def mask_spectroscopic_regions( wavelength, mask,
                               wave_range=[0,int(1e10)],
                               mask_feature_Hbeta=False,
                               mask_feature_4960=False,
                               mask_feature_4733=False,
                               mask_feature_5010=False,
                               mask_feature_5040=False,
                               mask_feature_5080=False,
                               mask_feature_5170=False,
                               mask_feature_5230=False,
                               mask_feature_5270=False,
                               mask_feature_5330=False):

    # minimum and maxinimum wavelength range of spectrum to be included
    iwavemin = np.argmin(np.abs( wavelength - wave_range[0] ))
    iwavemax = np.argmin(np.abs( wavelength - wave_range[1] ))
    mask[:iwavemin] = False
    mask[iwavemax:] = False

    # mask_features
    if mask_feature_4733:
        il,ir = np.argmin(np.abs( wavelength- 4815 )), np.argmin(np.abs( wavelength-4855 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_4960:
        il,ir = np.argmin(np.abs( wavelength- 5063 )), np.argmin(np.abs( wavelength-5067 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_5010:
        il,ir = np.argmin(np.abs( wavelength-5112  )), np.argmin(np.abs( wavelength-5117 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_5040:
        il,ir = np.argmin(np.abs( wavelength-5145  )), np.argmin(np.abs( wavelength-5155 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_5080:
        il,ir = np.argmin(np.abs( wavelength-5185  )), np.argmin(np.abs( wavelength-5195 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_5170:
        il,ir = np.argmin(np.abs( wavelength-5273 )), np.argmin(np.abs( wavelength-5301 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_5230:
        il,ir = np.argmin(np.abs( wavelength-5335  )), np.argmin(np.abs( wavelength-5345 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_5270:
        il,ir = np.argmin(np.abs( wavelength-5372  )), np.argmin(np.abs( wavelength-5393 )) # wavelength in observed frame
        mask[il:ir] = False
    if mask_feature_5330:
        il,ir = np.argmin(np.abs( wavelength-5437  )), np.argmin(np.abs( wavelength-5447 )) # wavelength in observed frame
        mask[il:ir] = False

    if mask_feature_Hbeta:
        il,ir = np.argmin(np.abs( wavelength- 4957 )), np.argmin(np.abs( wavelength-4975 )) # wavelength in observed frame
        mask[il:ir] = False

    return mask
