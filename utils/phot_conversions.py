import numpy as np


def get_flux_Jy( mAB=None, maggie=None, **extras):
    """
    per unit frequency
    """
    if maggie is None:
        maggie  = get_maggie( mAB=mAB )

    flux_Jy = maggie / 3631.
    return flux_Jy

def get_flux_Jy_unc( maggie=None, maggie_unc=None, flux_Jy=None, mAB_unc=None, mAB=None, **extras):
    """
    mAB = -2.5 log maggie = -2.5 log f + ZP
    f = 10^(-0.4 (mAB-ZP) )
    sigma_f = ( ln 10 / -2.5 ) * f / mAB_unc
    """

    if maggie_unc is None:
        maggie_unc = get_maggie_unc( maggie=maggie, mAB=mAB, mAB_unc=mAB_unc, flux_Jy=flux_Jy, **extras )
    flux_Jy_unc = get_flux_Jy( maggie=maggie_unc )

#     if maggie_unc is not None:
#         flux_Jy_unc = get_flux_Jy( maggie=maggie_unc )

#     elif mAB_unc is not None:
#         if flux_Jy is None:
#             flux_Jy = get_flux_Jy( mAB=mAB )
#         flux_Jy_unc = ( flux_Jy / mAB_unc ) * ( np.log(10) / -2.5 )

#     else:
#         print("Error: provide input ")
#         raise Exception

    return np.abs( flux_Jy_unc )

def get_maggie( mAB=None, flux_Jy=None, **extras ):
    """
     mAB = -2.5 log maggie
     maggie = 10^( mAB / -2.5 ) = 10^( -0.4 * mAB )
     1 maggie = 3631 Jy
    """
    if mAB is not None:    maggie = np.power( 10., mAB / -2.5 )
    elif flux_Jy is not None: maggie = 3631. * flux_Jy
    else:
        print("Error: provide input ")
        raise Exception
    return maggie

def get_maggie_unc( maggie=None, mAB=None, mAB_unc=None, flux_Jy=None, flux_Jy_unc=None, snr=None, **extras ):
    """
    mAB = -2.5 log maggie
    mAB_unc = (-2.5/ln10) * ( maggie_unc/maggie ) = (-2.5/ln10) * 1./snr
    maggie_unc = mAB_unc * maggie / (-2.5/ln10)

    1 maggie = 3631 Jy
    maggie_unc = 3631 * flux_Jy_unc
    """
#     print( "mAB_unc: {}, flux_Jy_unc: {}".format( mAB_unc, flux_Jy_unc ))
    if (mAB_unc is not None):
        if maggie is None:
            if mAB is not None:
                maggie = get_maggie( mAB=mAB, **extras)
            elif flux_Jy is not None:
                maggie = get_maggie( flux_Jy=flux_Jy, **extras)
            else:
                print("Error: provide input ")
                raise Exception

        maggie_unc  = ( np.log(10) / -2.5 ) * mAB_unc * maggie

    elif (flux_Jy_unc is not None):
        maggie_unc = get_maggie( flux_Jy=flux_Jy_unc, **extras)

    else:
        print("Error: provide input ")
        raise Exception

    return np.abs( maggie_unc )

def get_mAB( flux_Jy=None, maggie=None, zeropoint=8.9, **extras):
    """
    m = -2.5 log f + ZP
    in AB magnitudes and flux in Jy, ZP = 8.9
    or
    mAB = -2.5 log maggie
    """
    if maggie is None:
        maggie = get_maggie( flux_Jy=flux_Jy )

    mAB = -2.5*np.log10( maggie )

    return mAB

def get_mAB_unc( maggie=None, maggie_unc=None, flux_Jy=None, flux_Jy_unc=None, **extras):
    """
    mAB = -2.5 log f + ZP = -2.5 log maggie
    sigma_m = ( -2.5 / ln 10 ) * sigma_ x / x = ( -2.5 / ln 10 ) * 1/snr
    """
    if flux_Jy is not None:
        if flux_Jy_unc is None:
            flux_Jy_unc = get_flux_Jy_unc( maggie_unc=maggie_unc )
        snr = flux_Jy / flux_Jy_unc
    elif maggie is not None:
        if maggie_unc is None:
            maggie_unc = get_maggie_unc( maggie_unc=maggie_unc )
        snr = maggie / maggie_unc
    else:
        print("Error: provide input ")
        raise Exception

    mag_unc = (-2.5/np.log(10)) / snr

    return np.abs( mag_unc )

# def get_sigma_m_from_snr( snr ):
#     return 2.5 * np.log10( 1 + 1./snr )
#
# def get_snr_from_sigma_m( sigma_m ):
#     return 1./( 10**(0.4*sigma_m) - 1 )

# def get_sigma_x_with_snr( x, snr):
# #     return np.log10( 1 + 1./snr ) * np.log(10) * x
#     return np.log( 1 + 1./snr ) * np.log(10) * x
