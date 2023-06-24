import numpy as np

"""
Some things copied from prospect.models.transforms + added functionality for 2D arrays
Some convenience functions added
"""

__all__ = ["sfr_to_mwa", \
           "logsfr_ratios_to_masses", \
           "logsfr_ratios_to_sfrs", \
           "build_agebins", \
           "zfrac_to_sfrac", \
           "zfrac_to_masses", \
           "zfrac_to_sfr", \
           "sfrac_to_logsfr_ratios", \
           "chain_to_sfr", \
           "chain_to_mwa", \
           "chain_to_masses", \
           "sfr_dt", \
           "cumulate_masses", \
           "chain_to_logmass", \
           "chain_to_fitparam", \
           "chain_to_Av", \
           "chain_to_ssfr_dt", \
           "chain_to_logsfr_ratios_ind", \
           "chain_to_z_fraction_ind", \
           "chain_to_sfr_ind", \
           "chain_to_ssfr_ind", \
           "chain_to_cummass_ind", \
           ]

####################

def chain_to_param( chain, theta_index, param, **extras ):
    """
    Copy entry for parameter from the chain, based on index specified in theta_index
    Pay attention to the dimension of the chain: MCMC and nested sampling chains have different dimensions
    """
    if chain.ndim>2: # MCMC
        x = chain[:, :, theta_index[param]]
    else: # nested sampling
        x = chain[:, theta_index[param]]
    return x

def chain_to_sfr( chain, theta_index, agebins=None, norm_by_mass=True,
                  mass=None, times=None, **extras ):
    """
    Extract SFRs from Prospector output chain without the bother of specifying the SFH prior
    input:
    - chain = results["chain"], where 'results' is the output of prospect.io.read_results
    - theta_index = a dictionary with parameter names (e.g., "logmass") with the index in the chain
    - agebins = SFH time bins
    - (optional) norm_by_mass = True, sets mass = 1 such that returns SFRs
    - (optional) mass = None, use this mass rather than the distribution in the chain
    Output:
    sfrs (or sSFR if mass provided, and norm_by_mass=False)
    """

    # nonparametric SFH model with continuity prior
    if "logsfr_ratios" in theta_index.keys():
        logsfr_ratios = chain_to_param( chain, theta_index, 'logsfr_ratios')
        if norm_by_mass: logmass = 0
        elif mass is None: chain_to_param( chain, theta_index, 'logmass').T
        else: logmass = mass
        sfrs = logsfr_ratios_to_sfrs( logmass=logmass, logsfr_ratios=logsfr_ratios, agebins=agebins )

    # nonparametric SFH model with Dirichlet prior
    elif "z_fraction" in theta_index.keys():
        z_fraction = chain_to_param( chain, theta_index, 'z_fraction')
        if norm_by_mass:
            total_mass = 1
        elif (mass is None) and ("total_mass" in theta_index.keys()):
            total_mass = chain_to_param( chain, theta_index, 'total_mass').T
        elif (mass is None) and ("logmass" in theta_index.keys()):
            total_mass = 10** chain_to_param( chain, theta_index, 'logmass').T
        elif mass is not None:
            total_mass = mass
        else:
            total_mass = 1
        sfrs = zfrac_to_sfr( total_mass=total_mass, z_fraction=z_fraction, agebins=agebins )

    # parametric SFH model
    elif 'tau' in theta_index.keys():
        from prospect.plotting.sfh import parametric_sfr

        tages = chain_to_param( chain, theta_index, 'tage')
        taus = chain_to_param( chain, theta_index, 'tau')
        masses = chain_to_param( chain, theta_index, 'mass')

        sfrs = []
        for tage,tau,mass in zip( tages, taus, masses ):
            sfr = parametric_sfr( tage=tage, \
                                   times=times, \
                                   **dict( tau=tau, mass=mass )
                                 )
            sfrs.append( sfr )
        sfrs = np.vstack( sfrs )

    else:
        return None

    return sfrs

def chain_to_mwa( chain, theta_index, agebins=None, sfh=None, **extras ):
    """
    Take output from 'chain_to_sfr' and convert it to mass-weighted age.
    """
    # if 'agebins' exists, nonparametric SFH model
    if agebins is not None:
        sfrs = chain_to_sfr( chain, theta_index, agebins, **extras )
        mwa = sfr_to_mwa( agebins, sfrs, **extras )

    # parametric model
    elif sfh is not None:
        power = int( sfh > 3 ) # 0 = declining exponential, 1 = delayed declining exponential

        from prospect.plotting.sfh import parametric_mwa
        mwa = parametric_mwa( tau=chain_to_param( chain, theta_index, 'tau'), \
                              tage=chain_to_param( chain, theta_index, 'tage'), \
                              power = power,
                            )
    return mwa

def chain_to_masses( chain, theta_index, agebins=None, times=None, **extras ):
    """
    Extract masses from Prospector output chain without the bother of specifying the SFH prior
    input:
    - chain = results["chain"], where 'results' is the output of prospect.io.read_results
    - theta_index = a dictionary with parameter names (e.g., "logmass") with the index in the chain
    - agebins = SFH time bins
    Output:
    masses
    """
    # nonparametric SFH model with continuity prior
    if "logsfr_ratios" in theta_index.keys():
        logsfr_ratios = chain_to_param( chain, theta_index, 'logsfr_ratios')
        masses = logsfr_ratios_to_masses(logmass=0, logsfr_ratios=logsfr_ratios, agebins=agebins).T

    # nonparametric SFH model with Dirichlet prior
    elif "z_fraction" in theta_index.keys():
        z_fraction = chain_to_param( chain, theta_index, 'z_fraction')
        masses = zfrac_to_masses(total_mass=1, z_fraction=z_fraction, agebins=agebins)

    # parametric SFH model
    elif 'tau' in theta_index.keys():
        from prospect.plotting.sfh import parametric_cmf

        tages = chain_to_param( chain, theta_index, 'tage')
        taus = chain_to_param( chain, theta_index, 'tau')
        masses = chain_to_param( chain, theta_index, 'mass')

        all_masses = []
        for tage,tau,mass in zip( tages, taus, masses ):
            m = parametric_cmf( tage=tage, \
                                times=times, \
                                **dict( tau=tau, mass=mass )
                              )
            all_masses.append( m )
        masses = np.vstack( all_masses )

    return masses

def chain_to_logmass( chain, theta_index, **extras ):
    mass = chain_to_param( chain, theta_index, 'mass')
    return np.log10( mass )

def chain_to_Av( chain, theta_index, **extras ):
    dust2 = chain_to_param( chain, theta_index, 'dust2')[:,0]
    dust_ratio, dust_index = None,None
    if "dust_ratio" in theta_index.keys():
        dust_ratio = chain_to_param( chain, theta_index, 'dust_ratio')[:,0]
    if "dust_index" in theta_index.keys():
        dust_index = chain_to_param( chain, theta_index, 'dust_index')[:,0]

    from utils.dust_transforms import get_attenuation_curve,tau_to_extinct
    tau = get_attenuation_curve( 5500., dust_ratio=dust_ratio,
                                       dust2=dust2,
                                       dust_index=dust_index,
                                )
    extinct = tau_to_extinct( tau )
    return extinct[0]

def sfr_to_mwa(agebins, sfrs, **extras):
    """
    t_mw = ( int_0^t_obs t SFR(t) dt ) / ( int_0^t_obs SFR(t) dt)
    """
    if sfrs is None: return None

    agebins_Gyr = np.power(10., agebins) *1e-9 # Gyr
    mts = np.median(agebins_Gyr, axis=1)
    dts = np.diff(agebins_Gyr)[:,0]
    if sfrs.ndim>1:
        t_mws = np.sum( sfrs * mts * dts, axis=1) / np.sum(  sfrs * dts, axis=1)
    else:
        t_mws = np.sum( sfrs * mts * dts ) / np.sum(  sfrs * dts )
    return t_mws

#################### the following are adapted Prospector scripts

def logsfr_ratios_to_masses(logmass=0, logsfr_ratios=None, agebins=None, **extras):
    """This converts from an array of log_10(SFR_j / SFR_{j+1}) and a value of
    log10(\Sum_i M_i) to values of M_i.  j=0 is the most recent bin in lookback
    time.
    """
    nbins = agebins.shape[0]
    sratios = 10**np.clip(logsfr_ratios, -100, 100) # numerical issues...
    dt = (10**agebins[:, 1] - 10**agebins[:, 0])

    if sratios.ndim>1:
        coeffs = np.array([ (1. / np.prod(sratios[:,:i], axis=1)) * (np.prod(dt[1: i+1]) / np.prod(dt[: i])) for i in range(nbins)])
        m1 = (10**logmass) / np.sum(coeffs, axis=0)
        try:
            m1 = ( m1.T * coeffs.T).T
        except:
            m1 = m1.T * coeffs
    else:
        coeffs = np.array([ (1. / np.prod(sratios[:i]))           * (np.prod(dt[1: i+1]) / np.prod(dt[: i])) for i in range(nbins)])
        m1 = (10**logmass) / coeffs.sum()
        m1 = m1 * coeffs

    return m1

def logsfr_ratios_to_sfrs(logmass=0, logsfr_ratios=None, agebins=None, **extras):
    """Convenience function
    """
    masses = logsfr_ratios_to_masses(logmass=logmass, logsfr_ratios=logsfr_ratios, agebins=agebins)
    dt = (10**agebins[:, 1] - 10**agebins[:, 0])

    if masses.ndim>1:
        nx,ny = np.shape(masses)[1], len(agebins)
        dt = np.tile(dt, (nx,1)).T
        sfrs = masses / dt
        if np.shape(logsfr_ratios)+(0,1) != np.shape(sfrs): sfrs = sfrs.T
        return sfrs
    else:
        return masses / dt

# Functions for the SFH prior
def zfrac_to_sfrac(z_fraction=None, **extras):
    z_fraction = np.atleast_2d( z_fraction )
    shape = list(z_fraction.shape)
    shape[-1] += 1
    sfr_fraction = np.zeros( shape )
    sfr_fraction[:, 0] = 1. - z_fraction[:, 0]
    for i in range(1, shape[-1]-1):
        sfr_fraction[:,i] = np.prod( z_fraction[:,:i], axis=-1 ) * ( 1. - z_fraction[:,i] )
    sfr_fraction[:,-1] = 1. - np.sum( sfr_fraction[:,:-1], axis=-1 )
    sfr_fraction = np.squeeze( sfr_fraction )
    return sfr_fraction

# functions copied from prospect.models.transforms (and adapeted for 2D arrays)
def zfrac_to_masses(z_fraction=None, total_mass=1, agebins=None, sfr_fraction=None, **extras ):
    if sfr_fraction is None: sfrac = zfrac_to_sfrac(z_fraction)
#     sfr_fraction = np.atleast_2d( np.copy( sfr_fraction ) )
    else: sfrac = np.copy(sfr_fraction)

    # convert to mass fractions
    time_per_bin = np.diff(np.power(10, agebins), axis=-1)[:,0]
    sfrac *= time_per_bin
    mtot = np.atleast_1d( sfrac.sum(axis=-1) )
    mass_fraction = sfrac / mtot[:,None]
    masses = np.atleast_2d(total_mass) * mass_fraction.T
    return masses.T

def zfrac_to_sfr(total_mass=1, z_fraction=None, agebins=None, sfr_fraction=None, **extras):
    time_per_bin = np.diff(np.power(10, agebins), axis=-1)[:,0]
    masses = zfrac_to_masses(total_mass=total_mass, z_fraction=z_fraction, agebins=agebins, sfr_fraction=sfr_fraction)
    return masses / time_per_bin

def sfrac_to_logsfr_ratios( sfrac, **extras ):
    N = len( sfrac )
    logsfr_ratios = np.full( N-1, np.nan )

    for i in range( N-1 ):
        logsfr_ratios[i] = np.log10( sfrac[i]/sfrac[i+1])
    return logsfr_ratios

#################### miscellaneous

def zred_to_age( zred=None, **extras ):
    from prospect.sources.constants import cosmo
    return cosmo.age(zred).value

def dust2_to_Av( dust2 ):
    """ 2.5 * log10( exp(1) ) * tau_dust """
    print("Caution: This function does not work as intended for dust_type=4")
    return 2.5*np.log10(np.exp(1)) * dust2

# Define fixed age bins
def build_agebins( redshift, ncomp, tuniv=None, tlims_first=[0.03,0.1,0.5,1.], tlims_logspace=False, tbinmax=None, **extras ):
    """
    Define fixed age bins
    redshift = 1.2 is the median redshift of the GOGREEN spectroscopic sample
    Age bins are spaced:
        0 < t < 30 Myr
        30 < t < 100 Myr
        100 < t < 500 Myr
        500 Myr < t < 1 Gyr
        ncomp-4 linearly spaced bins
        0.95*t_universe < t < t_universe
    """
    from prospect.sources.constants import cosmo # import cosmology assumed when fitting
    # if age of the Universe not provided, calculate based on redshift and cosmology
    if tuniv is None: tuniv = cosmo.age(redshift).value
    # if maximum time bin not provided, set based on 95\% of the age of the Universe
    if tbinmax is None: tbinmax = (tuniv * 0.95)

    # specify edges of the time bins
    agelims = [1e-9] # close enough to zero
    # starts with fixed time bins
    agelims += tlims_first[:-1]
    # can edit as necessary, currently specifies linearlly spaced time bins
    if tlims_logspace:
        agelims += np.logspace( np.log10(tlims_first[-1]), np.log10(tbinmax), nbins-len(tlims_first) ).tolist()
    else:
        agelims += np.linspace( tlims_first[-1], tbinmax, ncomp-len(tlims_first) ).tolist()
    # last time bin covers tbinmax to tuniv
    agelims += [tuniv]

    # convert to units of log(t/yr)
    agelims = np.log10( np.array(agelims) * 1e9)
    # convert from list of bin edges to to array of bins
    agebins = np.array([agelims[:-1], agelims[1:]])
    agebins = agebins.T

    # agebins_Gyr = np.power(10., agebins) *1e-9 # Gyr
    return agebins

def add_Av_to_chain( result ):
    if 'Av' in result['theta_index'].keys(): return result

    chain = result['chain']

    # calculate mass-weighted age,
    x = chain_to_Av( chain, result['theta_index'] )

    # pay attention to dimension of the chain when adding new entry
    if chain.ndim>2:
        result['chain'] = np.dstack([ chain, x ])
    else:
        result['chain'] = np.vstack([ chain.T, x ]).T

    N,M = np.shape( chain )
    result['theta_index']['Av'] = slice( M, M+1 )

    return result

# def cumulate_masses( masses, **extras ):
#     """ Cumulative sum of masses  """
#     if masses.ndim == 1:
#         return np.cumsum( masses[::-1] )[::-1]
#     else:
#         return np.cumsum( masses[:,::-1], axis=1 )[:,::-1]

def sfh_to_tquantiles_binnedSFH( cmf=None, agebins=None, cmf_quantiles=[0.5], **extras ):
    from scipy.interpolate import interp1d

    ages = np.unique( agebins )

    age_quantiles = np.array([ interp1d( c, ages  )(cmf_quantiles) for c in cmf ])

    return age_quantiles, ages[-1] - age_quantiles

def chain_to_cmf( chain, theta_index, agebins=None, **extras ):
    sfrs = chain_to_sfr( norm_by_mass=True, chain=chain, theta_index=theta_index, agebins=agebins )

    from prospect.plotting.sfh import sfh_to_cmf
    x, cmfs = sfh_to_cmf( sfrs, agebins )

    return cmfs

def chain_to_quantiles( chain, theta_index, agebins=None, cmf_quantiles=[0.5], **extras ):
    cmfs = chain_to_cmf( chain, theta_index, agebins=agebins )
    agebins_Gyr = 10**(agebins-9)

    _, ctime_quantiles = sfh_to_tquantiles_binnedSFH( cmf=cmfs, agebins=agebins_Gyr, cmf_quantiles=cmf_quantiles )

    return ctime_quantiles
