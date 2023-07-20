
from Dragonfly44_SFH.fitting.input.param_fit_setup import *

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
    elif sfh in [1,4]: # parametric model
        model_params = build_model_dexp( model_params, **extras )
    else:
        print('Error: sfh model {} not yet implemented. Exiting...'.format(sfh))
        raise Exception

    # Non-parameteric SFH fitting for mass in flexible time bins with a smoothness prior
    model_params['imf_type']['init'] = 1 # Chabrier

    model_params["zred"] = {"N": 1,
                            "init": zred,
                            "isfree": extras["fit_redshift"],
                            "prior": priors.TopHat(mini=np.max([0, zred-0.01]), maxi=zred+0.01),
#                            "prior": priors.TopHat(mini=0, maxi=10),
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
            model = PolySpecModel_diffspeccal( model_params )
#            model = sedmodel.PolySpecModel(model_params)
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
