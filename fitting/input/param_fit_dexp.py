import os
import time, sys
import numpy as np

from prospect import prospect_args
from prospect.fitting import fit_model
from prospect.io import write_results as writer

from Dragonfly44_SFH.fitting.input.param_fit_setup_zsortafit import *

# --------------
# RUN_PARAMS
# --------------

# https://prospect.readthedocs.io/en/latest/usage.html
default_run_params = {
              # Obs data parameters
              'sfh': 4,
              'prior_tau_log': True,

              'file_data': path_data+"Dragonfly44/obs_phot_specKCWI_sigma110.h5",
              'path_outfile': path_fits,

              'fit_sigma':    False, # whether to fit for the smoothing parameters (False usually)
              'fit_logzsol':  True, # Whether metallicity is a free or fixed parameter
              'fit_mass':  True, # Whether mass is a free or fixed parameter

              'fit_redshift':     True, # Whether redshift is a free or fixed parameter

              'fit_spectra':      True, # Fit spectroscopy?
              'fit_phot':         True, # Fit photometry?

              'npoly': 3, # Degress of Chebyshev polynomial for specphot calibration
              'spec_norm':1e-8, # guess at the normalization of the spectrum (relative to photometry)

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

              'mask_feature_4733':True,
              'mask_feature_Hbeta':False,
              'mask_feature_4960':False,
              'mask_feature_5010':False,
              'mask_feature_5040':False,
              'mask_feature_5080':False,
              'mask_feature_5170':False,
              'mask_feature_5230':False,
              'mask_feature_5270':False,
              'mask_feature_5330':False,

              }


# -----------
# Everything
# ------------

if __name__=='__main__':

    # - Parser with default arguments -
    parser = prospect_args.get_parser()

    # default arguments
    # outfile
    # verbose
    # debug

    # - Add custom arguments -
    parser.add_argument('--file_data', type=str, default=default_run_params['file_data'], help='Name of file with spectral information')

    parser.add_argument('--fit_spectra', type=int, default=default_run_params['fit_spectra'], help='Whether to fit spectra or not')
    parser.add_argument('--fit_phot', type=int, default=default_run_params['fit_phot'], help='Whether to fit photometry or not')
    parser.add_argument('--fit_sigma', type=int, default=default_run_params['fit_sigma'], help='Whether to fit smoothing parameter')
    parser.add_argument('--fit_mass', type=int, default=default_run_params['fit_mass'], help='Whether to stellar mass')
    parser.add_argument('--fit_redshift', type=int, default=default_run_params['fit_redshift'], help='Whether to fit redshift')
    parser.add_argument('--prior_tau_log', type=int, default=default_run_params['prior_tau_log'], help='Whether the prior on tau is log scaled')

    parser.add_argument('--wave_range', type=int, nargs="*", default=[0,10000])
    parser.add_argument('--max_snr', type=float, default=1000000, )
    parser.add_argument('--npoly', type=int, default=default_run_params['npoly'] )

    args = parser.parse_args()
    run_params = vars(args)

    if not run_params['fit_spectra']:
        run_params['fit_sigma'] = False
        run_params['fit_noise_spec'] = False
        run_params['fit_outlier_spec'] = False
        run_params['opt_polynomial'] = False
        run_params['fit_polynomial'] = False
        run_params['rescale_spectrum'] = False

    if '/' not in run_params['outfile']:
        run_params['outfile'] = default_run_params['path_outfile']+run_params['outfile']

    # print to command line list of run_params
    if run_params['verbose']:
        print ('\nRun params:\n')
        for k,v in run_params.items():
            print ('    {}: {}'.format( k, v))
        print ('\n')
    for k,v in default_run_params.items():
        if k not in run_params.keys():
            run_params[k] =v
            if run_params['verbose']:
                print ('    {}: {}'.format( k, v))
    if run_params['verbose']: print( '\n')

    obs = build_obs( **run_params )
    model = build_model( **run_params )
    sps = build_sps(**run_params)
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
