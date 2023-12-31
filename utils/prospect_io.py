import numpy as np
import h5py
import prospect.io.read_results as pread
from Dragonfly44_SFH.fitting.prospect.sources.constants import cosmo
import astropy.units as u
import os
from copy import deepcopy

from Dragonfly44_SFH.utils import transforms

def read_input(file_data=None, **extras):
    import h5py
    ddict = {}
    with h5py.File(file_data, 'r') as hfile:

        for key in ['wavelength',"Redshift","sigma_smooth","e_sigma_smooth","inres",  "maggies","spectrum","maggies_unc","unc"]:
            v = np.copy( hfile[key] )
            if v.size>1:
                ddict[key] = v
            else:
                ddict[key] = v.item()

        for key in ['filternames']:
            ddict[key] = np.array( np.copy( hfile[key] ), dtype=str )

        key = "smoothtype"
        ddict[key] = np.copy( hfile[key] ).astype(str).item()

        try: ddict['mask'] = np.copy( hfile['mask'] )
        except: pass

    from sedpy.observate import load_filters
    phot_filters = load_filters( np.array( ddict['filternames'], dtype=str) )
    ddict["filters"] = np.copy( phot_filters )
    ddict['phot_wave'] = np.array([ f.wave_effective for f in phot_filters ])

    return ddict

def read_results( result_file, more_groups=["draws"] , manual_agebins=None, sfh=None, **extras ):

    result = get_res( result_file, more_groups=more_groups, **extras )
    try:
        model = get_model( result, **extras )
    except Exception as e:
        print(e)
        model = None

    if sfh is None:
        if model is not None:
            sfh = int( model.params['sfh'] )

    if "agebins" not in result.keys():
        if (model is not None) and (manual_agebins is None):
            if 'agebins' in model.params:
                # if verbose: print('Warning: "agebins" not in result.')

                agebins = model.params['agebins']
                result['agebins'] = agebins

        elif manual_agebins is not None:
            result['agebins'] = manual_agebins


    from .prospect_postprocessing import add_missing_mass
    result = add_missing_mass( result, model )

    result = check_theta_index( result )

    from .prospect_postprocessing import add_mwa_to_chain
    result = add_mwa_to_chain( result, sfh=sfh )

    from .prospect_postprocessing import add_Av_to_chain
    result = add_Av_to_chain( result )

    from .prospect_postprocessing import add_logmass_stellar_to_chain
    result = add_logmass_stellar_to_chain( result )

    if 'agebins' in result.keys():
        from .prospect_postprocessing import add_sfr_quantiles_to_chain
        result = add_sfr_quantiles_to_chain( result, ctimescales=[0.5,0.9] )

    return result, result['obs'], model

def get_res( result_file, pickle=False, more_groups=[], manual_theta_labels=None, \
             verbose=False, start_index=0, thin=1, run_params_edit={}, **extras):
    """
    Building on prospect.io.read_results read Prospector h5py output and
    extra saved groups
    """
    if pickle: from prospect.io.read_results import unpick

    result= pread.read_hdf5(result_file, **extras)

    if "chain" not in result.keys():
        print("Error: fit {}, incomplete, no chain".format( result_file ))
        return None,None,None

    for key,v in run_params_edit.items():
        result['run_params'][key] = v

    try:
        mod = pread.get_model( result, **extras )
        result['model'] = mod
    except:
        result['model'] = None

    groups = ["theta_index", "fit", "posteriors",  "restframe"]+more_groups
    keys = ["agebins","stellar_mass_fraction"]

    with h5py.File(result_file, "r") as hf:

        if "sampling" not in hf.keys():
            print("File is incomplete, no sampling")
            print("rm ",result_file)
            return

        for group in groups:
            if group not in hf.keys():
                # print( group )
                continue
            result[group] = {}

            for k, v in hf[group].items():

                if type(hf[group][k])==h5py._hl.group.Group:
                    result[group][k] = {}
                    for k2,v2 in hf[group][k].items():
                        try:
                            result[group][k][k2] = json.loads(v2)
                        except:
                            result[group][k][k2] = np.array(v2)
                else:
                    try:
                        result[group][k] = json.loads(v)
                    except:
                        result[group][k] = np.array(v)

            for k, v in hf[group].attrs.items():
                if pickle:
                    try:
                        result[group][k] = unpick(v)
                    except:
                        result[group][k] = np.array(v)
                else:
                    try:
                        result[group][k] = json.loads(v)
                    except:
                        result[group][k] = np.array(v)

        for k in keys:
            if k not in hf.keys(): continue

            result[k] = np.array( hf[k] )


    if "chain" not in result.keys():
        print("WARNING: no chain in file {}".format(result_file))

    if result["chain"].ndim > 2:
        flatchain = result['chain'][:, start_index::thin, :]
        dims = flatchain.shape
        flatchain = flatchain.reshape(dims[0]*dims[1], dims[2])
        result["chain_unflattened"] = deepcopy( result['chain'] )
        result["chain"] = flatchain
    elif result["chain"].ndim == 2:
        flatchain = result["chain"][start_index::thin, :]
        result["chain"] = flatchain

    if 'theta_labels' not in result.keys():
        if verbose: print('Warning: "theta_labels" not in result.')
        try:
            model = get_model(res)
            result['theta_labels'] = model.theta_labels()
        except:
            if ( manual_theta_labels is None) :
                print('Error: provide "manual_theta_labels", default {}'.format( manual_theta_labels ))
            else:
                result['theta_labels'] = manual_theta_labels

    if "theta_index" not in result.keys():
        result= build_theta_index( result)

    if 'lnprobability' not in result.keys():
        if verbose: print('Warning: "lnprobability" not in result.')
        try:
            try: model
            except: model = get_model(res)
            agebins = model.params['agebins']
            result['lnprobability'] = result['lnlikelihood'] + model.prior_product( result['chain'] )
        except:
            pass

    return result

def build_theta_index( result):
    """
    Add to the resultdictionary and entry which specifies indices of the theta keys, "theta_index"
    """
    assert "theta_labels" in result.keys(), "Error: theta_labels not in res"

    theta_index = {}
    model = None
    model_nbins = 0

    # iterate through theta,
    for i,theta_label in enumerate(result['theta_labels']):
        if "logsfr_ratio" in theta_label:
            model = "continuity"
            model_nbins = i
            continue
        elif "z_fraction" in theta_label:
            model_nbins = i
            model = "dirichlet"
            continue
        else:
            theta_index[theta_label] = slice(i,i+1)

    assert model is not None, "Error: couldn't determine model type"
    if model=="continuity":
        for i,theta_label in enumerate(result['theta_labels']):
            if theta_label=='logsfr_ratios_1': break
        theta_index['logsfr_ratios'] = slice( i, model_nbins+1 )

    if model=="dirichlet":
        if   'z_fraction_0' in result['theta_labels']: first_key = 'z_fraction_0'
        elif 'z_fraction_1' in result['theta_labels']: first_key = 'z_fraction_1'
        else: raise Exception

        for i,theta_label in enumerate(result['theta_labels']):
            if theta_label==first_key: break
        theta_index['z_fraction'] = slice(i,model_nbins+1)

    result['theta_index'] = theta_index
    return result

def check_theta_index( result):
    """
    Go through "theta_index" entry in resultdictionary and check if all the keys are there,
    if not, add them
    """

    # if any entries are arrays (mistake in some save files), change to slices
    for th,v in result['theta_index'].items():
        if (type(v) is np.ndarray) or (type(v) is list):
            i0 = v[0]
            i1 = v[-1]
            result['theta_index'][th] = slice( i0, i1+1 )

    # check if log mass entry exists, add otherwise
    if 'logmass' not in result['theta_index'].keys():
        from .prospect_postprocessing import add_logmass_to_chain
        result = add_logmass_to_chain( result )

    if 'logzsol' not in result['theta_index'].keys():
        from .prospect_postprocessing import add_logzsol_to_chain
        result = add_logzsol_to_chain( result )


    return result

def get_model(result, file_data=None, file_data_path=None, param_file_path=None, **extras ):
    """
    Build Prospector model based on parameter file.
    Pretty janky setup.
    """

    if file_data is not None:
        result['run_params']['file_data'] = file_data
    elif file_data_path is not None:
        _, file_data = os.path.split( result['run_params']['file_data'] )
        result['run_params']['file_data'] = file_data_path+file_data

    param_file = result['run_params'].get('param_file', '')
    paramfile_text = result.get("paramfile_text", '')
    path, filename = os.path.split(param_file)
    modname = filename.replace('.py','')
    modname = modname.replace('.','')

    user_module = pread.import_module_from_string( paramfile_text, modname )
    try:
        model = user_module.load_model(**result['run_params'])
    except(AttributeError):
        try:
            model = user_module.build_model(**result['run_params'])
        except:
            try:
                dict_data = user_module.read_input( **result['run_params'] )
                model = user_module.build_model( dict_data, **result['run_params'] )
            except:
                model = None
                print('Error: failed to produce model')

    return model

def get_model_with_dict_data(res):
    """
    Build Prospector model based on parameter file.
    Use when load_model requiresultdict_data as input.
    Pretty janky setup.
    """
    param_file, paramfile_text = ( result['run_params'].get('param_file', ''), result.get("paramfile_text", '') )
    path, filename = os.path.split(param_file)
    modname = filename.replace('.py','')
    modname = modname.replace('.','')

    user_module = pread.import_module_from_string( paramfile_text, modname )
    dict_data = user_module.read_input( **result['run_params'] )
    try:
        model = user_module.load_model(dict_data, **result['run_params'])
    except(AttributeError):
        try:
            model = user_module.build_model(dict_data, **result['run_params'])
        except(AttributeError):

            model = user_module.build_model_continuity(dict_data,  **result['run_params'])

    return model


def get_models_from_posterior_draws( use_saved_draws, result, draws_quantiles=[0.16,0.84], **extras ):
    """
    Function to get saved SED models from posterior draws, saved in Prospector output files
    Note that the public version of Prospector does not automatically calculate this information
    The SED models are generated post-fitting, and appended to the output files.

    If use_saved_draws is False, return empty arrays
    """

    # default empty lists
    draws = { 'quantiles':{ 'qs':draws_quantiles } }
    for k in ["phot", "spec", "sed", "speccal", "sed2"]:
        draws[k] = []
        draws['quantiles'][k] = []

    # get SED models for draws from the posteriors (saved in result['draws'])
    # calculate quantiles

    if 'draws' not in result.keys():
        return None

    # spectrum
    if ('spectrum' in result['draws'].keys()):
        draws['spec'] = np.copy( result['draws']['spectrum'] )
        draws['quantiles']['spec'] = np.quantile( draws['spec'], q=draws_quantiles, axis=0)

    # photometry
    if 'photometry' in result['draws'].keys():
        draws['phot'] = np.copy( result['draws']['photometry'] )
        draws['quantiles']['phot'] = np.quantile( draws['phot'], q=draws_quantiles, axis=0)

    # full SED
    if 'sed' in result['draws'].keys():
        draws['sed'] = np.copy( result['draws']['sed'] )
        draws['quantiles']['sed'] = np.quantile( draws['sed'], q=draws_quantiles, axis=0)

    # full SED v2
    if 'sed2' in result['draws'].keys():
        draws['sed2'] = np.copy( result['draws']['sed2'] )
        draws['quantiles']['sed2'] = np.quantile( draws['sed2'], q=draws_quantiles, axis=0)

    # spectro-photometric calibration polynomial
    if 'speccal' in result['draws'].keys():
        x = np.copy( result['draws']['speccal'] )
        if x.size>0:
            draws['speccal'] = x
            draws['quantiles']['speccal'] = np.quantile( draws['speccal'], q=draws_quantiles, axis=0)

    return draws

def get_models_for_bestfits( result ):
    """
    Function to get saved SED models from bestfit solution, saved in Prospector output files

    Note that the public version of Prospector automatically calculate this information
    in result['bestfit'], however sometimes this is nonsense
    When calculating the SED models for draws from the posteriors (post-fitting), the bestfit model
    is also generated. In this case, can use use_saved_draws=True.

    """
    # default empty lists
    bestfits = {  }
    for k in ["phot", "spec", "sed", "speccal", "sed2", "wave_sed"]:
        bestfits[k] = []

    if 'draws' in result.keys():
        bestfits["spec"] = np.copy( result['draws']['bestfit']['spectrum'] )
        bestfits["phot"] = np.copy( result['draws']['bestfit']['photometry'] )
        bestfits["sed"] = np.copy( result['draws']['bestfit']['spectrum_div_speccal'] )
        bestfits["speccal"] = np.copy( result['draws']['bestfit']['speccal'] )

        try:
            bestfits["sed"] = np.copy( result['draws']['bestfit']['sed'] )
            bestfits["wave_sed"] = np.copy( result['draws']['bestfit']['wave_sed'] )
        except:
            pass

    else:
        bestfits["spec"] = np.copy( result['bestfit']['spectrum'] )
        bestfits["phot"] = np.copy( result['bestfit']['photometry'] )
        bestfits["speccal"] = np.copy( result['bestfit']['speccal'] )
        try:
            bestfits["sed"] = np.copy( result['bestfit']['spectrum_div_speccal'] )
        except:
            bestfits["sed"] = bestfits["spec"] / bestfits["speccal"]

    return bestfits

def get_prior_draws( model, size=int(1e5), params=None ):
    from prospect.plotting.utils import sample_prior
    prior_draws, labels = sample_prior( model, size )
    prior_theta_index = {}
    i = 0
    for p in labels:
        s = model.params[p].size
        prior_theta_index[p] = slice(i,i+s)
        i+=s

    prior_dict = {}
    if params is None: params = labels

    if int(model.params['sfh'])==3:
        agebins = model.params['agebins']
    else:
        agebins = None

    for p in params:

        x = None

        if p in labels:
            x = prior_draws[:, prior_theta_index[p] ]
            prior_dict[p] = np.squeeze(x)
        elif p=='sfr':
            x = chain_to_sfr( prior_draws, prior_theta_index, agebins=agebins, norm_by_mass=False )
        elif p=='ssfr':
            x = chain_to_sfr( prior_draws, prior_theta_index, agebins=agebins, norm_by_mass=True )
        elif p=='cmf':
            if 'sfr' in prior_dict.keys():
                prior_sfrs = prior_dict['sfr']
            else:
                prior_sfrs = chain_to_sfr( prior_draws, prior_theta_index, agebins=agebins, norm_by_mass=False )
            from prospect.plotting.sfh import sfh_to_cmf
            if model.params['sfh']==3:
                _, x = sfh_to_cmf( prior_sfrs, model.params['agebins'] )
            else:
                print('Warning: havent added function to calculate cmf for parametric models')
        elif p.lower()=='mwa':
            from Dragonfly44_SFH.utils.transforms import chain_to_mwa
            x = chain_to_mwa( prior_draws, prior_theta_index, agebins=agebins)

        elif p.startswith('t'):
            try:

                if agebins is None:
                    print('Warning: havent added function to calculate quantiles for parametric models')
                else:
                    from Dragonfly44_SFH.utils.transforms import chain_to_quantiles
                    qt = float( p[1:] )/100.
                    x = chain_to_quantiles( cmf_quantiles=[qt], chain=prior_draws,
                                            theta_index=prior_theta_index, agebins=agebins )


            except Exception as e:
                print(p,e)
                pass

        if x is None:
            print("Warning: havne't added function to calculate {}".format(p))

        prior_dict[p] = np.squeeze(x)

    return prior_dict


def get_sfh_priors( model, size=int(1e5) ):
    from prospect.plotting.utils import sample_prior

    prior_draws, label = sample_prior( model, size )

    prior_theta_index = {}
    i = 0
    for p in label:
        s = model.params[p].size
        prior_theta_index[p] = slice(i,i+s)
        i+=s

    from Dragonfly44_SFH.utils.transforms import chain_to_sfr
    agebins = model.params['agebins']

    prior_sfrs = chain_to_sfr( prior_draws, prior_theta_index, agebins=agebins, norm_by_mass=False )
    prior_ssfrs = chain_to_sfr( prior_draws, prior_theta_index, agebins=agebins, norm_by_mass=True )

    from prospect.plotting.sfh import sfh_to_cmf
    x, prior_cmfs = sfh_to_cmf( prior_sfrs, agebins )

    return dict( sfr=prior_sfrs, ssfr=prior_ssfrs, cmf=prior_cmfs, agebins=agebins )
