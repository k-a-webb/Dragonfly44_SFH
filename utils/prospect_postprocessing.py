import os
import numpy as np

def get_models( theta, obs, model, sps, **extras ):

    # mock obs to calculate full sed
    obs_mock = dict( wavelength=None, spectrum=None, unc=None, mask=None )
    for key in ['filters','filternames','maggies','maggies_unc']:
        obs_mock[key]=obs[key]

    # get models
    spectrum, _, _ = model.predict( theta, obs=obs, sps=sps)
    speccal = model._speccal.copy()
    sed1 = spectrum/speccal

    sed2, phot, mfrac = model.predict( theta, obs=obs_mock, sps=sps)
    sed_wavelength = sps.wavelengths*(1.+model.params['zred'])


    return dict( spectrum=spectrum, sed=sed1, sed2=sed2, speccal=speccal,
                  photometry=phot, mfrac=mfrac, sed_wavelength=sed_wavelength)

def save_bestfit_model( result_file=None, group=None, overwrite=False, **extras ):
    import h5py
    from .h5py_utils import add_data, add_group
    # from .prospect_io import read_results
    from prospect.io.read_results import get_sps, results_from

    # check if already saved to output
    with h5py.File(result_file,'r') as hfile:
        if group is None:
            hfile_group = hfile
        else:
            assert group in hfile.keys(), 'Error: "{}" not in {}'.format( group, result_file )
            hfile_group = hfile[group]

        if ( 'sed' in hfile_group['bestfit'].keys()) and not overwrite:
            print("Bestfit SED already saved, returning")
            return

    # result, obs, model = read_results( result_file, **extras )
    result, obs, model = results_from( result_file, **extras )
    sps = get_sps(result)

    params_bf = [ result['bestfit']['parameter'][ result['theta_index'][p] ][0]
                  for result in model.theta_labels() ]
    models_bf = get_models( params_bf, obs, model, sps )


    # save output
    with h5py.File( result_file, 'a' ) as hfile:
        if group is None:
            hfile_group = hfile
        else:
            if group in hfile.keys():
                hfile_group = hfile[group]
            else:
                hfile_group = add_group( hfile, group )
        hfile_group_bf = add_group( hfile_group, 'bestfit' )

        add_data( hfile_group_bf, 'sed', models_bf['sed'], overwrite=overwrite )
        add_data( hfile_group_bf, 'wave_sed', models_bf['sed_wavelength'], overwrite=overwrite )

    return

def draw_seds( model, chain, weights=None, n_seds=int(100),
               obs=None, sps=None, obs2=None, theta_samples=None, **extras ):

    if theta_samples is None:
        from prospect.plotting.utils import sample_posterior
        # theta_samples = sample_posterior( chain, weights, nsample=n_seds )
        chain_inds = np.arange( chain.shape[0] ).astype(int)
        theta_samples, theta_indices = sample_posterior( chain, weights, nsample=n_seds, extra=chain_inds )

    spectra, phots, speccals, seds, mfracs, = [], [], [], [], []
    seds2 = []

    fit_spec=False
    if 'spectrum' in obs.keys():
        if obs['spectrum'] is not None:
            fit_spec=True

    for theta in theta_samples:
        s, p, m = model.predict(theta, obs=obs, sps=sps)
        spectra.append(s)
        phots.append(p)
        mfracs.append(m)

        if fit_spec: speccals.append( model._speccal.copy() )
        seds.append( model._spec )

        if obs2 is not None:
            s, _,_ = model.predict(theta, obs=obs2, sps=sps)
            seds2.append(s)

    return dict( spectrum=np.array(spectra),
                 photometry=np.array(phots),
                 speccal=np.array(speccals),
                 sed=np.array(seds),
                 mfrac=np.array(mfracs),
                 sed2=np.array(seds2),
                 theta_samples=np.array(theta_samples),
                 theta_indices=np.array(theta_indices),
               )

def save_posterior_draws_models( result_file=None, overwrite=False, n_seds=100,
                                 save_bestfit=True, **extras ):
    import h5py
    from .h5py_utils import add_data, add_group
    # from .prospect_io import read_results
    from prospect.io.read_results import get_sps, results_from

    # check if already saved to output
    with h5py.File(result_file,'r') as hfile:

        if 'draws' in hfile.keys():
            ndraws = np.shape( hfile['draws']['draws'] )[0]
            if (ndraws >= n_seds) and not overwrite:
                print("SED draws already saved, returning")
                return

    # result, obs, model = read_results( result_file, **extras )
    result, obs, model = results_from( result_file, **extras )
    sps = get_sps(result)

    # mock obs to calculate full sed
    obs_mock = dict( wavelength=None, spectrum=None, unc=None, mask=None )
    for key in ['filters','filternames','maggies','maggies_unc']:
        obs_mock[key]=obs[key]
    sed_wavelength = sps.wavelengths*(1.+model.params['zred'])

    draws = draw_seds( model, result['chain'], weights=result['weights'], n_seds=n_seds,
                       obs=obs, sps=sps, obs2=obs_mock )

    if save_bestfit:
        params_bf = result['bestfit']['parameter']
        models_bf = get_models( params_bf, obs, model, sps )


    # save output
    with h5py.File( result_file, 'a' ) as hfile:

        if 'draws' in hfile.keys():
            hfile_draws = hfile['draws']
        else:
            hfile_draws = add_group( hfile, 'draws' )

        add_data( hfile_draws, 'draws', draws['theta_samples'], overwrite=overwrite )
        add_data( hfile_draws, 'drawn_inds', draws['theta_indices'], overwrite=overwrite )

        add_data( hfile_draws, 'spectrum', draws['spectrum'], overwrite=overwrite )
        add_data( hfile_draws, 'photometry', draws['photometry'], overwrite=overwrite )
        add_data( hfile_draws, 'sed', draws['sed'], overwrite=overwrite )
        add_data( hfile_draws, 'speccal', draws['speccal'], overwrite=overwrite )
        add_data( hfile_draws, 'mfrac', draws['mfrac'], overwrite=overwrite )
        add_data( hfile_draws, 'sed2', draws['sed2'], overwrite=overwrite )

        add_data( hfile_draws, 'wave_sed', sed_wavelength, overwrite=overwrite )

        if save_bestfit:
            if 'bestfit' in hfile_draws.keys():
                hfile_bf = hfile_draws['bestfit']
            else:
                hfile_bf = add_group( hfile_draws, 'bestfit' )

            add_data( hfile_bf, 'wave_sed', sed_wavelength, overwrite=overwrite )
            add_data( hfile_bf, 'photometry', models_bf['photometry'], overwrite=overwrite )
            add_data( hfile_bf, 'spectrum', models_bf['spectrum'], overwrite=overwrite )
            add_data( hfile_bf, 'spectrum_div_speccal', models_bf['sed'], overwrite=overwrite )
            add_data( hfile_bf, 'speccal', models_bf['speccal'], overwrite=overwrite )
            add_data( hfile_bf, 'mfrac', models_bf['mfrac'], overwrite=overwrite )
            add_data( hfile_bf, 'sed', models_bf['sed2'], overwrite=overwrite )

    return

def save_mfrac( result_file, overwrite=False, **extras ):
    import h5py
    from .h5py_utils import add_data, add_group
    # from .prospect_io import read_results
    from prospect.io.read_results import get_sps, results_from

    # check if already saved to output
    with h5py.File(result_file,'r') as hfile:

        if ( 'chain_mfrac' in hfile.keys()) and not overwrite:
            print("mfrac's already saved, returning")
            return

    # result, obs, model = read_results( result_file, **extras )
    result, obs, model = results_from( result_file, **extras )
    sps = get_sps(result)

    mfracs = []
    from prospect.utils.plotting import hist_samples
    flatchain, pnames = hist_samples(result)
    for theta in flatchain:
        s, p, m = model.predict(theta, obs=obs, sps=sps)
        mfracs.append(m)
    mfracs = np.squeeze( np.array( mfracs ) )

    # save output
    with h5py.File( result_file, 'a' ) as hfile:
        add_data( hfile, 'chain_mfrac', mfracs, overwrite=overwrite )

    return

def add_to_chain( result, x, label, label2=None, **extras ):
    chain = result['chain']
    x = np.squeeze( x )

    # pay attention to dimension of the chain when adding new entry
    if chain.ndim>2:
        result['chain'] = np.dstack([ chain, x ])
    else:
        result['chain'] = np.vstack([ chain.T, x ]).T

    N,M = np.shape( chain )
    result['theta_index'][label] = slice( M, M+1 )

    # if duplicate nameing convention, add second entry to theta_index
    if label2 is not None:
        result['theta_index'][label2] = slice( M, M+1 )

    return result

def add_to_bestfit( result, x, label ):
    # also copy bestfit values
    result['bestfit']['parameter'] = np.append( result['bestfit']['parameter'], x )
    # result['bestfit']['parameter'] = np.append( result['bestfit']['parameter'], x )

    return result

def add_Av_to_chain( result ):

    from .transforms import chain_to_Av
    x = chain_to_Av( **result )
    result = add_to_chain( result, x, 'Av' )

    return result

def add_mwa_to_chain( result ):

    from .transforms import chain_to_mwa
    x = chain_to_mwa( **result )
    result = add_to_chain( result, x, 'mwa', label2='MWA' )

    return result

def add_logmass_to_chain( result ):
    theta_index = result['theta_index']
    theta_index_keys = theta_index.keys()

    # if using mass-metallicity prior, copy entry with name "logmass"
    if 'massmet' in theta_index_keys:
        idx_logmass = theta_index['massmet'][0]

        result['theta_index']['logmass'] = slice( idx_logmass, idx_logmass+1 )
        result['theta_index']['massmet_1'] = slice( idx_logmass, idx_logmass+1 )

    # if using mass-metallicity prior, copy entry with name "logmass"
    elif 'massmet_1' in theta_index_keys:
        result['theta_index']['logmass'] = theta_index['massmet_1']

    # if total_mass exists, calculated log-mass and make new entry
    elif ('total_mass' in theta_index_keys):
        from .transforms import chain_to_param
        x = chain_to_param( result['chain'], result['theta_index'], 'total_mass' )
        x = np.log10( x )
        result = add_to_chain( result, x, 'logmass' )

    elif ('mass' in theta_index_keys):

        from .transforms import chain_to_param
        x = chain_to_param( result['chain'], result['theta_index'], 'mass' )
        x = np.log10( x )
        result = add_to_chain( result, x, 'logmass' )

    else:
        print('Caution: No logmass or massmet/massmet_1 or total_mass parameter!')

    return result

def add_logzsol_to_chain( result ):

    # follow same logic to add metallicity information

    if 'massmet' in result['theta_index'].keys():
        idx_logzsol = result['theta_index']['massmet'][1]
        result['theta_index']['logzsol'] = slice( idx_logzsol, idx_logzsol+1 )
        result['theta_index']['massmet_2'] = slice( idx_logzsol, idx_logzsol+1 )

    elif 'massmet_2' in result['theta_index'].keys():
        result['theta_index']['logzsol'] = result['theta_index']['massmet_2']

    else:
        print('Caution: No logzsol or massmet parameter!')

    return result


def add_missing_mass( result, model, **extras ):
    thidx_keys = result['theta_index']

    if ('total_mass' not in thidx_keys) and ('logmass' not in thidx_keys):
        if model is None:
            print("Warning: Mass is missing, but model is None, can't fill... ")
            return result

        chain = result['chain']

        if 'total_mass' in model.params:
            total_mass = np.full( chain.shape[:-1], model.params['total_mass'] )
            logmass = np.log10( total_mass )
        elif 'logmass' in model.params:
            logmass = np.full( chain.shape[:-1], model.params['logmass'] )
            total_mass =10**logmass
        elif 'mass' in model.params:
            total_mass = model.params['mass'].sum()
            total_mass = np.full( chain.shape[:-1], total_mass )
            logmass = np.log10( total_mass )

        total_mass = np.squeeze(total_mass)
        logmass = np.squeeze(logmass)

        npar = np.shape(chain)[-1]
        result['chain'] = np.vstack([ chain.T, total_mass ]).T
        result['chain'] = np.vstack([ result['chain'].T, logmass ]).T

        result['theta_index']['total_mass'] = slice( npar, npar+1 )
        result['theta_index']['logmass'] = slice( npar+1, npar+2 )

    return result

def add_mfrac_to_chain( result ):

    chain_shape = result['chain'].shape
    x = np.full( chain_shape[:-1], np.nan )

    if 'chain_mfrac' in result.keys():
        x = np.copy( result['chain_mfrac'] )

    elif ('draws' in result.keys()):
        if 'mfrac' in result['draws'].keys():
            x_test = result['draws']['mfrac']
            if x_test.shape == chain_shape[:-1]:
                x = x_test
            elif 'drawn_inds' in result['draws'].keys():
                inds = result['draws']['drawn_inds']
                x[inds] = result['draws']['mfrac']
            else:
                print('Warning: mfrac saved in draws, but not used')

    else:
        print('Warning: mfrac not saved to output')

    result = add_to_chain( result, x, 'mfrac' )
    return result

def add_logmass_stellar_to_chain( result ):

    from .transforms import chain_to_param, logmass_mfrac_to_logmass_stellar

    if 'logmass' not in result['theta_index'].keys():
        result  = add_logmass_to_chain( result )

    if 'logmass' not in result['theta_index'].keys():
        print('Warning: no logmass, cant compute logmass_stellar')
        return result

    logmass_total = chain_to_param( param='logmass', **result )

    if 'mfrac' not in result['theta_index'].keys():
        result  = add_mfrac_to_chain( result )
    mfrac = chain_to_param( param='mfrac', **result )

    x = logmass_mfrac_to_logmass_stellar( logmass_total, mfrac )
    result = add_to_chain( result, x, 'logmass_stellar' )
    return result

def add_sfr_quantiles_to_chain( result, ctimescales=[0.5,0.9], **extras):

    from .transforms import chain_to_sfr
    sfrs = chain_to_sfr( norm_by_mass=True, **result )

    from prospect.plotting.sfh import sfh_to_cmf
    x, cmfs = sfh_to_cmf( sfrs, result['agebins'] )

    Nx,M = np.shape(cmfs)
    Nt = len(ctimescales)

    txs_lookback = np.full((Nx,Nt),np.nan)
    from scipy.interpolate import interp1d
    for ii in range(Nx):
        txs_lookback[ii,:] = interp1d( cmfs[ii,:], x )(ctimescales)
    txs_cosmic = x[-1] - txs_lookback

    for ii,t in enumerate(ctimescales):
        # in units of cosmic time
        str_t = 't{:.0f}'.format(t*100)
        x = np.squeeze( txs_cosmic[:,ii] )

        if str_t not in result['theta_index'].keys():
            result = add_to_chain( result, x, str_t )

    if (0.5 in ctimescales) and (0.9 in ctimescales) and ('t90mt50' not in result['theta_index'].keys()):
        i50 = np.argmin(np.abs( np.array(ctimescales) - 0.5 ))
        i90 = np.argmin(np.abs( np.array(ctimescales) - 0.9 ))
        x = np.squeeze( txs_cosmic[:,i90] -  txs_cosmic[:,i50] )
        result  = add_to_chain( result, x, 't90mt50' )

    return result
