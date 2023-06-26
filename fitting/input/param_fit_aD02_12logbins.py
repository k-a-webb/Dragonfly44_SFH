import os
import time, sys
import numpy as np

from prospect import prospect_args
from prospect.fitting import fit_model
from prospect.io import write_results as writer
from prospect.io import read_results as pread
from prospect.sources.constants import cosmo

from utils.transforms import build_agebins
from param_fit_setup import build_obs, build_sps_nonparametric, build_gp

from param_fit_setup import *

# --------------
# RUN_PARAMS
# --------------

# https://prospect.readthedocs.io/en/latest/usage.html
global_run_params = {
              # Obs data parameters
              'sfh': 3,
              'alphaD': 0.2,

              'nbins': 12,
              'tlims_first': [0.03, 0.1, 0.5, 1., 2., 3],
              'tlims_logspace': True,

              'file_data': path_data+"DF44/obs_phot_specKCWI_sigma110.h5",
              'outfile': 'aD1_phot_specKCWI_12logbins',

              'fit_sigma':    False, # whether to fit for the smoothing parameters (False usually)
              'fit_logzsol':  True, # Whether metallicity is a free or fixed parameter
              'fit_mass':  True, # Whether mass is a free or fixed parameter

              'fit_redshift':     False, # Whether redshift is a free or fixed parameter

              'fit_spectra':      True, # Fit spectroscopy?
              'fit_phot':         True, # Fit photometry?

              'npoly': 3, # Degress of Chebyshev polynomial for specphot calibration

              'opt_polynomial':   True, # whether to optimize a polynomial ratio between F_true and F_obs
              'fit_polynomial':   False, # whether to fit for a polynomial ratio between F_true and F_obs (ie includes prior)
              'rescale_spectrum': False, # If True, rescale the spectrum to have an average of 1 before doing anything.

              'fit_agn':     False, # fit agn emission
              'fit_duste':   True, # fit dust emission
              'fit_neb':     False, # fit nebular emission

              'fit_outlier_spec': True, # fit for outliers in spectroscopy
              'fit_noise_spec':   True, # Jitter parameter

              'wave_range': [3000,10000],
              'max_snr': 10000000,
              }

# --------------
# MODEL_PARAMS
# --------------

# def build_agebins( redshift, nbins, **extras ):
#     """
#     Define fixed age bins
#     redshift = 1.2 is the median redshift of the GOGREEN spectroscopic sample
#     Age bins are spaced:
#         0 < t < 30 Myr
#         30 < t < 100 Myr
#         100 < t < 500 Myr
#         500 Myr < t < 1 Gyr
#         nbins-4 linearly spaced bins
#         0.95*t_universe < t < t_universe
#     """
#     from prospect.sources.constants import cosmo
#     import astropy.units as u
#     tuniv = cosmo.age(redshift).value
#     tbinmax = (tuniv * 0.95)
#     #lim1, lim2, lim3, lim4 = 0.03, 0.1, 0.5, 1.
#     #agelims = [1e-9,lim1,lim2,lim3] + np.linspace(lim4,tbinmax,nbins-4).tolist() + [tuniv]
#     lims = [0.03, 0.1, 0.5, 1., 2., 3]
#     agelims = [1e-9] + lims[:-1] + np.logspace( np.log10(lims[-1]), np.log10(tbinmax), nbins-len(lims) ).tolist() + [tuniv]
#
#     agelims = np.log10( np.array(agelims) * 1e9)
#     agebins = np.array([agelims[:-1], agelims[1:]])
#     agebins = agebins.T
#
#     #agebins_Gyr = np.power(10., agebins) *1e-9 # Gyr
#     return agebins



# -----------
# Everything
# ------------

if __name__=='__main__':

    # - Parser with default arguments -
    parser = prospect_args.get_parser()

    # - Add custom arguments -
    parser.add_argument('--file_data', type=str, default=global_run_params['file_data'], help='Name of file with spectral information')

    parser.add_argument('--fit_spectra', type=int, default=1, help='Whether to fit spectra or not')
    parser.add_argument('--fit_phot', type=int, default=1, help='Whether to fit photometry or not')
    parser.add_argument('--fit_logzsol', type=int, default=1, help='Whether to fit logzsol')
    parser.add_argument('--fit_sigma', type=int, default=1, help='Whether to fit smoothing parameter')
    parser.add_argument('--wave_range', type=int, nargs="*", default=[0,10000])
    parser.add_argument('--max_snr', type=float, default=1000000, )
    parser.add_argument('--npoly', type=int, default=global_run_params['npoly'] )

    parser.add_argument('--mask_feature_Hbeta', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_4733', type=int, default=1, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_4960', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_5010', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_5040', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_5080', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_5170', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_5230', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_5270', type=int, default=0, help='Whether to mask weird feature')
    parser.add_argument('--mask_feature_5330', type=int, default=0, help='Whether to mask weird feature')

    args = parser.parse_args()
    run_params = vars(args)

    if not run_params['fit_spectra']:
        run_params['fit_redshift'] = False
        run_params['fit_sigma'] = False
        run_params['fit_noise_spec'] = False
        run_params['fit_outlier_spec'] = False
        run_params['opt_polynomial'] = False
        run_params['fit_polynomial'] = False
        run_params['rescale_spectrum'] = False

    if run_params['verbose']:
        print ('\nRun params:\n')
        for k,v in run_params.items():
            print ('    {}: {}'.format( k, v))
        print ('\n')
    for k,v in global_run_params.items():
        if k not in run_params.keys():
            run_params[k] =v
            if run_params['verbose']:
                print ('    {}: {}'.format( k, v))
    if run_params['verbose']: print( '\n')

    model = build_model( **run_params )
    obs = build_obs( **run_params )
    sps = build_sps_nonparametric(**run_params)
    noise = build_noise(**run_params)

    model_param_text = 'Model params\n'
    if run_params['verbose']: print ('\nModel params:')
    for key in model.init_config.keys() :
        txt = '    {}: {}'.format(key, model.init_config[key])
        if run_params['verbose']: print (txt)
        model_param_text += txt+'\n'
    if run_params['verbose']: print ('\n')

    run_params['modelparam_text'] = model_param_text
    run_params["sps_libraries"] = sps.ssp.libraries
    run_params["param_file"] = __file__

    if args.debug:
        print ('Stopping for debug')
        sys.exit()

    hfilename = "{0}_{1}_mcmc.h5".format(run_params['outfile'], int(time.time()))
    import h5py
    hfile = h5py.File(hfilename, 'a')
    writer.write_h5_header(hfile, run_params, model)
    writer.write_obs_to_h5(hfile, obs)
    print ('>>>> Writing to file {}'.format(hfilename))

    print ('>>>> Sampling ...')
    output = fit_model(obs, model, sps, noise, **run_params)

    print ('>>>> Writing to output ...')
    writer.write_hdf5(hfile, run_params, model, obs,
                      output["sampling"][0], output["optimization"][0],
                      tsample=output["sampling"][1],
                      toptimize=output["optimization"][1],
                      sps=sps)

    try:
        hfile.close()
    except(AttributeError):
        pass


    # Add more to output (that Prospect doesn't automatically save)
    with h5py.File(hfilename,'a') as hfile:

        # Save theta index because sometimes model informaion is not saved correctly
        thidx = hfile.create_group('theta_index')
        for key,val in model.theta_index.items():
            thidx[key] = np.arange( model.ndim )[ val ]
        hfile['agebins'] = model.params['agebins']

    print('>>>> Done <<<<')
