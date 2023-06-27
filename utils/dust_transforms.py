import numpy as np
__all__ = ["get_attenuation_curve_diffuse", \
           "get_attenuation_curve_birth", \
           "get_attenuation_curve", \
           "tau_to_extinct", \
          ]

############# Charlot+Fall
from sedpy.attenuation import noll

def get_attenuation_curve_diffuse(wave, dust2=0.0, dust_index=1.0, **extras):
    """ Calculate diffuse dust attenuation curve
    Parameters
    ----------
    wave: Float or 1-D Array
        Wavelengths (Angstroms) at which attenuation curve should be evaluated
    dust_index: Float
        Slope parameter of Noll+09 dust attenuation parametrization--difference between true slope and that of Calzetti+00 curve, with positive values signifying shallower curves
    dust2: Float
        Diffuse dust optical depth; also referred to as tau throughout this document

    Returns
    -------
    Diffuse dust attenuation curve at given wavelength(s)
    """
    if (len(dust2)>1) & (np.ndim(wave)<2):
        wave =  np.tile( wave, ( len(dust2),1)).T
    Eb = 0.85 - 1.9*dust_index
    return noll(wave, tau_v=dust2, delta=dust_index, c_r=0.0, Ebump=Eb)

def get_attenuation_curve_birth(wave, dust_ratio=1.0, dust2=0.0, dust_index=1.0, **extras):
    """ Calculate birth cloud dust attenuation curve
    Parameters
    ----------
    wave: Float or 1-D Array
        Wavelengths (Angstroms) at which attenuation curve should be evaluated
    d1: Float
        Birth cloud dust optical depth

    Returns
    -------
    Birth dust attenuation curve at given wavelength(s); inverse law from Charlot+Fall 00 assumed
    """
    if (len(dust2)>1) & (np.ndim(wave)<2):
        wave =  np.tile( wave, ( len(dust2),1)).T

    dust2 = get_attenuation_curve_diffuse(wave, dust2=dust2, dust_index=dust_index)
    dust1 = dust_ratio * dust2

    return dust1*(wave/5500.0)**(-1)

def tau_to_extinct( tau ):
    return 2.5*np.log10(np.exp(1)) * tau

def get_attenuation_curve( wave, dust_ratio=1.0, dust2=0.0, dust_index=1.0, **extras):
    tau_diffuse = get_attenuation_curve_diffuse(wave, dust2=dust2, dust_index=dust_index, **extras)
    tau_birth = get_attenuation_curve_birth(wave, dust_ratio=dust_ratio, dust2=dust2, dust_index=dust_index, **extras)
    return tau_diffuse+tau_birth

############# MW

import astropy.units as u

def attenuation_law_MW_infrared( x, mwr=3.1, threshs=[0.1,1.1], **extras ):
        x = np.array(x)
        if threshs is not None:
            x = x[ (threshs[0]<=x) & (x<threshs[1]) ]
            if len(x)<1: return []

        a_x =  0.574 * x**(1.61)
        b_x = -0.527 * x**(1.61)
        return a_x + b_x / mwr

def attenuation_law_MW_optical_NIR( x, mwr=3.1, threshs=[1.1,3.3], **extras ):
        x = np.array(x)
        if threshs is not None:
            x = x[ (threshs[0]<=x) & (x<threshs[1]) ]
            if len(x)<1: return []

        a_coeffs = [1, 0.17699, -0.50447, -0.02427,  0.72085,  0.01979, -0.77530, 0.32999]
        b_coeffs = [0, 1.41338,  2.28305,  1.07233, -5.38434, -0.62251,  5.3026, -2.09002]
        a_x = np.polyval( a_coeffs[::-1], x-1.82 )
        b_x = np.polyval( b_coeffs[::-1], x-1.82 )
        return a_x + b_x / mwr

def attenuation_law_MW_NUV( x, mwr=3.1, uvb=1.0, threshs=[3.3,5.9], **extras ):
        x = np.array(x)
        if threshs is not None:
            x = x[ (threshs[0]<=x) & (x<threshs[1]) ]
            if len(x)<1: return []

        a_x = lambda xx: 1.752 - 0.316*xx - 0.104 / ((xx-4.67)**2 + 0.341) * uvb
        b_x = lambda xx: -3.09 + 1.825*xx + 1.206 / ((xx-4.62)**2 + 0.263) * uvb
        # this hack parameter is not in the original CCM89
        # parameterization.  It is a hack designed to result in
        # a smooth profile in the presence of a variable UVB strength
        hack = ( 3.3/x )**6.*( attenuation_law_MW_optical_NIR( 3.3, mwr=mwr, threshs=None ) - \
                                ( a_x(3.3) + b_x(3.3) / mwr )
                              )
        return a_x(x) + b_x(x) / mwr + hack

def attenuation_law_MW_MUV( x, mwr=3.1, uvb=1.0, threshs=[5.9,8], **extras ):
        x = np.array(x)
        if threshs is not None:
            x = x[ (threshs[0]<=x) & (x<threshs[1]) ]
            if len(x)<1: return []

        F_ax = -0.04473 * (x-5.9)**2 - 0.009779*(x-5.9)**3
        a_x =  1.752 - 0.316*x - 0.104 / ((x-4.67)**2 + 0.341) * uvb + F_ax
        F_bx =  0.2130  * (x-5.9)**2 + 0.1207  *(x-5.9)**3
        b_x = -3.09  + 1.825*x + 1.206 / ((x-4.62)**2 + 0.263) * uvb + F_bx
        return a_x + b_x / mwr

def attenuation_law_MW_FUV( x, mwr=3.1, threshs=[8,12], **extras ):
        x = np.array(x)
        if threshs is not None:
            x = x[ (threshs[0]<=x) & (x<threshs[1]) ]
            if len(x)<1: return []
        y = x-8
        a_x = -1.073 - 0.628 * y + 0.137 * y**2 - 0.070*y**3
        b_x = 13.670 + 4.257 * y - 0.420 * y**2 + 0.374*y**3
        return a_x + b_x / mwr

def attenuation_law_MW_FUV_constant( x, mwr=3.1, thresh=[12], **extras ):
        x = np.array(x)
        if thresh is not None:
            x = x[ (thresh<=x) ]
            if len(x)<1: return []
        return attenuation_law_MW_FUV( np.full_like(x, 12), mwr=mwr, threshs=None, **extras )

def attenuation_law_MW( wave, dust2=1., mwr=3.1, uvb=1.0, **extras ):

    if ( type(wave) == int ):
        wave = float( wave )
    if ( type(wave) == float ):
        wave = np.array([wave])
    if (type(wave) == u.quantity.Quantity):
        wave_um = wave.to(u.um).value
    else:
        print("'wave' must be in units of Angstroms")
        wave_um = wave * 1e-4

    x = 1/wave_um

    Aratio = np.hstack([
                        attenuation_law_MW_infrared( x, mwr=mwr ),
                        attenuation_law_MW_optical_NIR( x, mwr=mwr ),
                        attenuation_law_MW_NUV( x, mwr=mwr, uvb=uvb ),
                        attenuation_law_MW_MUV( x, mwr=mwr, uvb=uvb ),
                        attenuation_law_MW_FUV( x, mwr=mwr ),
                        attenuation_law_MW_FUV_constant( x, mwr=mwr ),
                       ])

    if type(dust2) != float:
        if (len(dust2)>1) & (np.ndim(wave)<2):
            Aratio =  np.tile( Aratio, ( len(dust2),1)).T

    return wave[ (0.1<=x)  ], dust2 * Aratio


############## Calzetti 2000 #############

def get_calzetti2000_atten( wave, dust2, wave_units='A', dd63=6300.00, R = 4.05, **extras ):
    # https://iopscience.iop.org/article/10.1086/308692/pdf
    if wave is not np.array: wave = np.array([wave])
    klam = np.zeros_like( wave )
    if wave_units == 'A': wave = wave*1e-4 # convert to microns
    elif wave_units != 'um': print('Units can be "A" or "um", otherwise add conversion') ; raise Exception

    if np.any( wave > dd63 ):
        sel = wave > dd63
        klam[sel] = 2.659*( -1.857 + 1.04 / wave[sel]  ) + R

    if np.any( wave <= dd63 ):
        sel = wave <= dd63
        klam[sel]  = 2.659*( -2.156 + 1.509 / wave[sel] - 0.198 / wave[sel]**2 + 0.011 / wave[sel]**3 ) + R

    # klam = Alam / Es(B-V)
    # Es(B-V) = (0.44 +- 0.03) E(B-V)
    # Rv = Av / E(B-V)

    # klam = Alam / Es(B-V) = Alam / ( 0.44 E(B-V) ) =  Alam * Rv / ( 0.44 * Av )
    # Alam / Av = klam * 0.44 / Rv

    Alam_div_Av = klam * 0.44  / R
    return Alam_div_Av * dust2

############## Kriek 2013 #############


def get_kriek2013_atten( wave, dust2, wave_units='A', dust1=0, dust_index=1, R=4.05, lamv=5500., **extras ):
    """
    https://arxiv.org/pdf/1308.1099.pdf
    1) Alam = Av / 4.05 ( klam + Dlam ) ( lam / lam_V )**dust_index
    2) Dlam = Eb ( lam dlam )**2 / [ (lam**2 - lam0**2) + (lam dlam)**2 ]
    3) Eb = (0.85 +- 0.09 ) - (1.9 +- 0.4 ) * dust_index
    """
    if wave_units == 'um': lamv*= 1e-4 # convert from A to um
    elif wave_units != 'A': print('Units can be "A" or "um", otherwise add conversion') ; raise Exception


    Dlam = get_drude_bump( wave, dust_index )
    Alam_div_Av_C00 = get_calzetti2000_atten( wave, np.ones_like(dust2), wave_units=wave_units, R=R )
    Alam = dust2 * ( Alam_div_Av_C00 + Dlam/R ) * (wave/lamv)**dust_index

    return Alam

def get_drude_bump( wave, dust_index, dlam=350., lamuvb=2175.0, wave_units='A', **extras ):
    if wave_units == 'um': dlam*= 1e-4 ; lamuvb *= 1e-4 # convert from A to um
    elif wave_units != 'A': print('Units can be "A" or "um", otherwise add conversion') ; raise Exception

    Eb = 0.85 - 1.9 * dust_index
    Dlam = Eb * ( wave * dlam )**2 / ( ( wave**2 - lamuvb**2 ) + (wave * dlam)**2 )
    return Dlam
