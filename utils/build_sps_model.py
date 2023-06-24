# --------------
# SPS Object
# --------------

def build_sps_sps(zcontinuous=1, compute_vega_mags=False, **extras):
    from prospect.sources import SSPBasis
    sps = SSPBasis(zcontinuous=zcontinuous, compute_vega_mags=compute_vega_mags)
    return sps


def build_model_sps( zred, opt_polynomial=True, fit_duste=False, fit_agn=False, fit_neb=False, npoly=4 ):


    from prospect.models.templates import TemplateLibrary
    from prospect.models import priors, sedmodel, transforms

    model_params = TemplateLibrary["ssp"]

    # Non-parameteric SFH fitting for mass in flexible time bins with a smoothness prior
    model_params['imf_type'] = {"N": 1, "isfree": False, "init": 1}        # Chabrier

    model_params["zred"]['init'] = zred
    model_params["dust2"]['init'] = 0.01

    model_params["logzsol"]['init'] = -1.3
    model_params["mass"]['init'] = 10**(8.5)
    model_params["tage"]['init'] = 12

    if opt_polynomial:
        model_params.update(TemplateLibrary["optimize_speccal"])
        model_params["polyorder"]["init"] = npoly

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

    if opt_polynomial:
        #model = sedmodel.PolySpecModel(model_params)
        model = sedmodel.PolySpecModel_v2(model_params)
    else:
        model = sedmodel.SpecModel(model_params)

    model.params['sigma_smooth'] = 100
    model.params['smoothtype'] = "vel"
    model.params['inres'] = ( 3e5 ) * (2.5/2.355 / 5000) # 64 km/s

    return model
