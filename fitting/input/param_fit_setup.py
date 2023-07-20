import os
import time, sys
import numpy as np

from prospect import prospect_args
from prospect.fitting import fit_model
from prospect.io import write_results as writer
from prospect.io import read_results as pread
from prospect.sources.constants import cosmo

from Dragonfly44_SFH.utils import obs_io

import os
path_base = os.getcwd().split('Dragonfly44_SFH')[0] # gross, but works
path_fits = path_base+ "Dragonfly44_SFH/fitting/output/fits/"
path_data = path_base+ "Dragonfly44_SFH/data/"


# OBS
# -------------


def build_obs( data_dict=None,
               fit_spectra=True, fit_phot=True,
               max_snr_spec=100000., max_snr_phot=20., wave_range=[0,10000],
               norm_spec=True, file_data=None, \
               **extras ):

    if data_dict is None: data_dict = obs_io.read_file_data( file_data=file_data, **extras )

    obs = {} # Observation dictionary
    obs['Redshift'] = data_dict['Redshift']

    for key in ['phot_wave','filters','maggies','maggies_unc','spectrum','unc','mask','phot_mask','wavelength']:
        obs[key] = None

    ################################
    # photometric observations

    # copy over photometry, in units of maggies
    maggies, maggies_unc = [ np.copy( data_dict[key]) for key in ['maggies','maggies_unc']]

    # set minimum phot uncertainty (maxinimum S/N)
    if max_snr_phot is not None:
        where_unc_too_small = maggies_unc / maggies < 1./max_snr_phot
        maggies_unc[ where_unc_too_small ] = np.abs( maggies[ where_unc_too_small ] / max_snr_phot )   # inflate uncertainties

    # setup mask for photometry (all True by default)
    phot_mask = np.ones_like(maggies_unc, dtype=bool)

    # read in filters
    from sedpy.observate import load_filters
    phot_filters = load_filters( np.array( data_dict['filternames'], dtype=str) )
    phot_wave = np.array([ f.wave_effective for f in phot_filters ])

    if fit_phot:
        obs["maggies"]     = np.copy( maggies  )
        obs["maggies_unc"] = np.copy( maggies_unc )
        obs["phot_mask"]   = np.copy( phot_mask )
        obs["filters"]     = np.copy( phot_filters )
        obs["phot_wave"]   = np.copy( phot_wave )

    ################################
    # spectroscopic observations

    # copy over spectroscopy
    wavelength, spectrum, unc = [ np.copy( data_dict[key]) for key in ['wavelength','spectrum','unc'] ]

    # setup mask for spectrum
    if "mask" in data_dict.keys(): mask = np.copy( data_dict['mask'] ).astype(bool)
    else: mask = np.array(np.ones_like(wavelength), dtype=bool)
    mask = obs_io.mask_spectroscopic_regions( wavelength, mask, wave_range=wave_range, **extras)


    if fit_spectra:
        if max_snr_spec is not None:
            where_unc_too_small = spectrum / unc > max_snr_spec
            unc[ where_unc_too_small ] = np.abs( spectrum[ where_unc_too_small ] / max_snr_spec )   # inflate uncertainties

        # normalize spectra
        if norm_spec:
            norm = np.nanmedian( spectrum[ mask ] )
            spectrum /= norm
            unc /= norm

        obs['spectrum']   = np.copy( spectrum )
        obs['unc']        = np.copy( unc  )
        obs['mask']       = np.copy( mask )
        obs['wavelength'] = np.copy( wavelength  )

    from prospect.utils.obsutils import fix_obs
    obs = fix_obs(obs)  # This ensures all required keys are present and adds some extra useful info

    return obs

# --------------
# SPS Object
# --------------

def build_sps( sfh, zcontinuous=1, compute_vega_mags=False, **extras):
    if sfh == 3:
        return build_sps_nonparametric( zcontinuous, compute_vega_mags, **extras )
    else:
        raise Exception

def build_sps_nonparametric(zcontinuous=1, compute_vega_mags=False, **extras):
    from prospect.sources import FastStepBasis
    sps = FastStepBasis(zcontinuous=zcontinuous, compute_vega_mags=compute_vega_mags)
    return sps


# -----------------
# Noise Model
# ------------------

def build_noise(fit_noise_spec, **extras):
    if fit_noise_spec:
        from prospect.likelihood import NoiseModel
        from prospect.likelihood.kernels import Uncorrelated
        jitter = Uncorrelated(parnames = ['spec_jitter'])
        spec_noise = NoiseModel(kernels=[jitter],metric_name='unc',weight_by=['unc'])
        return spec_noise, None
    else:
        return None, None

# -----------------
# Gaussian Process
# ------------------

def build_gp(**extras):
    return None, None

# --------------
# MODEL_PARAMS
# --------------
# Define fixed age bins
def build_agebins( redshift, nbins, tuniv=None, tlims_first=[1e-9,0.03,0.1,0.5,1.], tlims_logspace=False, tbinmax=None, **extras ):
    """
    Define fixed age bins
    redshift = 1.2 is the median redshift of the GOGREEN spectroscopic sample
    Age bins are spaced:
        0 < t < 30 Myr
        30 < t < 100 Myr
        100 < t < 500 Myr
        500 Myr < t < 1 Gyr
        nbins-4 linearly spaced bins
        0.95*t_universe < t < t_universe
    """
    from prospect.sources.constants import cosmo # import cosmology assumed when fitting
    # if age of the Universe not provided, calculate based on redshift and cosmology
    if tuniv is None: tuniv = cosmo.age(redshift).value
    # if maximum time bin not provided, set based on 95\% of the age of the Universe
    if tbinmax is None: tbinmax = (tuniv * 0.95)

    # specify edges of the time bins
    # starts with fixed time bins
    agelims = tlims_first[:-1]
    # can edit as necessary, currently specifies linearlly spaced time bins
    if tlims_logspace:
        agelims += np.logspace( np.log10(tlims_first[-1]), np.log10(tbinmax), nbins-len(tlims_first)+1 ).tolist()
    else:
        agelims += np.linspace( tlims_first[-1], tbinmax, nbins-len(tlims_first)+1 ).tolist()
    # last time bin covers tbinmax to tuniv
    agelims += [tuniv]

    # convert to units of log(t/yr)
    agelims = np.log10( np.array(agelims) * 1e9)
    # convert from list of bin edges to to array of bins
    agebins = np.array([agelims[:-1], agelims[1:]])
    agebins = agebins.T

    # agebins_Gyr = np.power(10., agebins) *1e-9 # Gyr
    return agebins

def build_model_dirichlet( model_params, nbins, zred,
                           alphaD=1,
                           tlims_first=None, tlims_logspace=False, **extras ):
    """
    Add Dirichlet prior parameters
    """
    from prospect.models import priors
    from prospect.models.templates import TemplateLibrary

    model_params.update( TemplateLibrary["dirichlet_sfh"] )

    model_params["total_mass"] = {"N": 1,
                                   "init": 10**(8.48),
                                   "isfree": extras['fit_mass'],
                                   "prior": priors.LogUniform(mini=1e8, maxi=1e12),
                                   "units": "Solar masses formed",
                                  }

    agebins = build_agebins( zred, nbins, tlims_first=tlims_first, tlims_logspace=tlims_logspace, **extras )

    model_params['agebins']['N'] = nbins
    model_params['agebins']['init'] = agebins

    model_params['mass']['N'] = nbins
    model_params["z_fraction"]["N"] = nbins - 1
    model_params["z_fraction"]["init"] = np.zeros(nbins-1)

    beta = np.full( nbins-1, alphaD )
    alpha = np.cumsum( beta )[::-1]
    model_params["z_fraction"]["prior"] = priors.Beta( alpha=alpha, beta=beta, mini=0.0, maxi=1.0)

    return model_params

def build_model_continuity( model_params, nbins, zred, tlims_first=None, tlims_logspace=False,
                            stdT_mean=0., stdT_scale=0.3, stdT_df=2.,
                            **extras ):
    """
    Add continuity prior parameters
    """
    from prospect.models import priors
    from prospect.models.templates import TemplateLibrary

    model_params.update( TemplateLibrary["continuity_sfh"] )

    model_params["logmass"] = {"N": 1,
                               "init": 8.48,
                               "isfree": extras['fit_mass'],
                               "prior": priors.TopHat(mini=6, maxi=10),
                               "units": "log solar masses formed",
                              }

    agebins = build_agebins( zred, nbins, tlims_first=tlims_first, tlims_logspace=tlims_logspace, **extras )

    model_params['agebins']['N'] = nbins
    model_params['agebins']['init'] = agebins

    model_params['mass']['N'] = nbins

    model_params["logsfr_ratios"]["N"] = nbins - 1
    model_params["logsfr_ratios"]["init"] = np.zeros(nbins-1)

    if type(stdT_mean) is float: stdT_mean = np.full(nbins-1, stdT_mean)
    if type(stdT_scale) is float: stdT_scale = np.full(nbins-1, stdT_scale)
    if type(stdT_df) is float: stdT_df = np.full(nbins-1, stdT_df)

    model_params["logsfr_ratios"]["prior"] = priors.StudentT(mean=stdT_mean, scale=stdT_scale, df=stdT_df)


    return model_params

def build_model_UMtuned_continuity( model_params, nbins, zred, tlims_first=None, tlims_logspace=False,
                                    **extras ):
    """
    Add continuity prior parameters
    """
    from prospect.models import priors, transforms
    from prospect.models.templates import TemplateLibrary

    model_params.update( TemplateLibrary["beta"] )

    model_params['nzsfh'] = {'N': nbins+2,
                             'isfree': True,
                             'init': np.concatenate([[0.5,8,0.0], np.zeros(nbins-1)]),
                             'prior': priors_beta.NzSFH(zred_mini=1e-3, zred_maxi=15.0,
                                                    mass_mini=7.0, mass_maxi=12.5,
                                                    z_mini=-1.98, z_maxi=0.19,
                                                    logsfr_ratio_mini=-5.0, logsfr_ratio_maxi=5.0,
                                                    logsfr_ratio_tscale=0.3, nbins_sfh=nbins,
                                                    const_phi=True)}

    model_params["logmass"]["init"] = 8.48
    model_params["zred"]["init"] = zred

    model_params["logsfr_ratios"]["N"] = nbins-1
    model_params['mass']['N'] = nbins

    model_params['agebins']['N'] = nbins
    model_params['agebins']['init'] = transforms.zred_to_agebins_pbeta(np.atleast_1d(0.5), np.zeros(nbins))


    return model_params


def build_model_dust( model_params, dust_type=4, **extras ):
    """
    Add dust model parameters
    refer to FSPS manual for details
    """
    from prospect.models import priors, transforms

    # dust_type = 4 # Kriek
    # dust_type = 1 # MW like  # mwr = 3.1 by default, uvb = 1. by default
    # dust_type = 2 # Calzetti

    model_params['dust_type']['init'] = dust_type

    model_params["dust2"] = {"N": 1,
                             "init": 0.01,
                             "isfree": True,
                             "prior": priors.TopHat(mini=0, maxi=4.),
                             "units": "dust2",
                            }

    model_params["dust1"] = {"N": 1,
                             "isfree": False,
                             'depends_on': transforms.dustratio_to_dust1,
                             "init": 0.0,
                             "units": "optical depth towards young stars"
                            }
    if dust_type == 4:
        model_params["dust_ratio"] = {"N": 1,
                                      "isfree": True,
                                      "init": 1.0,
                                      "units": "ratio of birth-cloud to diffuse dust",
                                      "prior": priors.ClippedNormal(mini=0.0, maxi=2.0, mean=1.0, sigma=0.3)
                                     }

        model_params["dust_index"] = {"N": 1,
                                      "isfree": True,
                                      "init": 0.0,
                                      "units": "power-law multiplication of Calzetti",
                                      "prior": priors.TopHat(mini=-2.0, maxi=0.5)
                                     }
    return model_params


def build_model( data_dict=None, sfh=3,
                 alphaD=None, stdT_mean=None, spec_norm=1,
                 fit_duste=False, fit_agn=False, fit_neb=False, nbins=10, dust_type=4, **extras):
    """Construct a model.  This method defines a number of parameter
    specification dictionaries and uses them to initialize a
    `models.sedmodel.SedModel` object."""
    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors, sedmodel, transforms

    model_params = {}

    if data_dict is None: data_dict = obs_io.read_file_data( **extras )
    zred = data_dict['Redshift']

    if sfh==3: # nonparametric model
        if alphaD is not None:
            model_params = build_model_dirichlet( model_params, nbins, zred, alphaD=alphaD, **extras )
        elif stdT_mean is not None:
            model_params = build_model_continuity( model_params, nbins, zred,  stdT_mean=stdT_mean, **extras )
    elif sfh in [0,4]: # parametric model
        print('Error: parametric model not yet implemented. Exiting...')
        raise Exception

    # Non-parameteric SFH fitting for mass in flexible time bins with a smoothness prior
    model_params['imf_type']['init'] = 1 # Chabrier

    model_params["zred"] = {"N": 1,
                            "init": zred,
                            "isfree": extras["fit_redshift"],
#                            "prior": priors.TopHat(mini=np.max([0, zred-0.01]), maxi=zred+0.01),
                            "prior": priors.TopHat(mini=0, maxi=10),
                            "units": "spectroscopic redshift",
                            }

    model_params["logzsol"] = {"N": 1,
                                   "init": -1.26,
                                   "isfree": extras['fit_logzsol'],
                                   "units": "log Z/Zsol",
                                   "prior": priors.TopHat(mini=-2, maxi=0.19),
                                  }

    model_params = build_model_dust( model_params, dust_type, **extras )

    if extras["fit_sigma"]: # whether to smooth the model spectrum with a variable resolution (otherwise fixed)
        model_params.update(TemplateLibrary["spectral_smoothing"])

        model_params["smoothtype"]["init"] = data_dict["smoothtype"]

        model_params["sigma_smooth"]["isfree"]  = "True"
        model_params["sigma_smooth"]["init"]  = data_dict["sigma_smooth"]
        model_params["sigma_smooth"]["units"] = data_dict["smoothtype"]
        model_params["sigma_smooth"]["prior"] = priors.Normal(mean=data_dict['sigma_smooth'], sigma=5)

    assert not extras["fit_polynomial"], "Error: fit_polynomial not implemented"
    assert not (extras["fit_polynomial"] and extras["opt_polynomial"]), 'Error: Can either fit or opt polynomial, not both'
    if extras["fit_polynomial"]:

        model_params.update(TemplateLibrary["fit_speccal"])

        model_params["spec_norm"]["N"] = 1
        model_params["spec_norm"]["isfree"] = True
        model_params["spec_norm"]["prior"] = priors.Normal(mean=1, sigma=0.5)

        npoly = extras["npoly"]
        model_params["poly_coeffs"]["N"] = npoly
        model_params["poly_coeffs"]["init"] = np.zeros(npoly)
        polymax = 0.1 / (np.arange(npoly) + 1)
        model_params["poly_coeffs"]["prior"] = priors.TopHat(mini=-polymax, maxi=polymax)

    elif extras["opt_polynomial"]:
        model_params.update(TemplateLibrary["optimize_speccal"])

        npoly = extras["npoly"]
        model_params["spec_norm"]["N"] = spec_norm
        model_params["polyorder"]["init"] = npoly
        if "poly_regularization" in extras.keys():
            model_params["poly_regularization"]["init"] = extras['poly_regularization']

    if extras['fit_outlier_spec']:
        model_params["f_outlier_spec"] = {"N": 1,
                                          "isfree": True,
                                          "init": 0.01,
                                          "prior": priors.TopHat(mini=1e-5, maxi=0.5)}
        model_params["nsigma_outlier_spec"] = {"N": 1,
                                               "isfree": False,
                                               "init": 5.0}

    if extras['fit_noise_spec']:
        # add jitter
        model_params['spec_jitter'] = {"N": 1,
                                       "isfree": True,
                                       "init": 1.0,
                                       "prior": priors.TopHat(mini=1.0, maxi=3.0)}

    if fit_neb:
        model_params.update(TemplateLibrary["nebular"])

    if fit_duste:
        model_params.update(TemplateLibrary["dust_emission"])

        model_params['duste_umin']['isfree'] = True # MMP83 local MW intensity
        model_params['duste_qpah']['isfree'] = True # Percent mass fraction of PAHs in dust
        model_params['duste_gamma']['isfree'] = True # Mass fraction of dust in high radiation intensity

    if fit_agn:
        model_params.update(TemplateLibrary["agn"])

        model_params['fagn']['isfree'] = True # L_{AGN}/L_*
        model_params['add_agn_dust']['isfree'] = True # MMP83 local MW intensity
        model_params['agn_tau']['isfree'] = True # optical depth

    if extras["opt_polynomial"]:
        if (sfh==3) and (alphaD is not None): # Dirichlet SFH model
            # use edited version of the PolySpecModel prior
            # testing proved that the OG model fails when alphaD<1
            from Dragonfly44_SFH.fitting.prospect.models.sedmodel import PolySpecModel_dirichlet
            model = PolySpecModel_dirichlet(model_params)
        else:
            from Dragonfly44_SFH.fitting.prospect.models.sedmodel import PolySpecModel_diffspeccal
            model = PolySpecModel_diffspeccal(model_params)
            # model = sedmodel.PolySpecModel(model_params)
    elif extras["fit_polynomial"]:
        print("Error: Fitting coefficient of spec calibration not implemented")
        sys.exit()
    else:
    	model = sedmodel.SpecModel(model_params)

    if not extras["fit_sigma"]:
        model.params['sigma_smooth'] = data_dict["sigma_smooth"]
        model.params['smoothtype'] = data_dict["smoothtype"]
    model.params['inres'] = data_dict["inres"]

    return model
