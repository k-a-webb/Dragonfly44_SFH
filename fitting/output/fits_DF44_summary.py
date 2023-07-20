path_fits_DF44 = 'Dragonfly44/'

#################################################################
# Nonparametric SFH model
# Dirichlet prior, alpha_D = x, or aDx in label
# aD = 1 prefers to distribute SF equally among time bins
# aD < 1 increasingly prefers to distribute SF un-equally among time bins
#################################################################

# alphaD = 1, "smooth", "dispersed", or "exteneded" version of the SF prior
# specKCWI = KCWI spectroscopy
# phot = UV--NIR
# specMANGA = specMANGA

# phot_specKCWI = specKCWI + phot fit simultaneously
# phot_specMANGA = specMANGA + phot fit simultaneously

fit_aD1_phot_specKCWI = path_fits_DF44+"fit_DF44_aD1_phot_specKCWI_zsortafit_npoly8_1649733960_mcmc.h5"
fit_aD1_phot_specKCWI_v2 = path_fits_DF44+"fit_DF44_aD1_phot_specKCWI_zfixed_npoly3_1649733960_mcmc.h5"

fit_aD1_specKCWI = path_fits_DF44+ "fit_DF44_aD1_specKCWI_massfixed_zsortafit_npoly8_1634082827_mcmc.h5" # fixed mass
# fit_aD1_specKCWI = path_fits_DF44+ "fit_DF44_aD1_specKCWI_massfixed_zsortafit_npoly8_1650322035_mcmc.h5" # fixed mass
fit_aD1_specKCWI_massfit = path_fits_DF44+"fit_DF44_aD1_specKCWI_zsortafit_npoly8_1634071305_mcmc.h5" # free mass

fit_aD1_phot    = path_fits_DF44+"fit_DF44_aD1_phot_1649889938_mcmc.h5"
fit_aD1_phot_v2 = path_fits_DF44+"fit_DF44_aD1_phot_1687841281_mcmc.h5"

# minium age of time bins is 1 Gyr (i.e., suppress late burst of SF)
fit_aD1_older1Gyr_phot_specKCWI    = path_fits_DF44+"fit_DF44_aD1_older1Gyr_phot_specKCWI_zsortafit_npoly8_1649879742_mcmc.h5"
fit_aD1_older1Gyr_phot_specKCWI_v2 = path_fits_DF44+"fit_DF44_aD1_older1Gyr_phot_specKCWI_zfixed_npoly3_1687841282_mcmc.h5"
fit_aD1_older1Gyr_specKCWI = None
fit_aD1_older1Gyr_phot = None

# degrade S/N of spectroscopy (i.e., increase uncertainties) to test S/N dependence of the results
fit_aD1_phot_specKCWI_snr5 = path_fits_DF44+"fit_DF44_aD1_phot_specKCWI_zsortafit_npoly8_maxsnr5_1649879779_mcmc.h5"
fit_aD1_specKCWI_snr5 = path_fits_DF44+"fit_DF44_aD1_specKCWI_massfixed_maxsnr5_zsortafit_npoly8_1635282103_mcmc.h5"
fit_aD1_phot_snr5 = None

fit_aD1_phot_specKCWI_snr10 = path_fits_DF44+"fit_DF44_aD1_phot_specKCWI_zsortafit_npoly8_maxsnr10_1649879768_mcmc.h5"
fit_aD1_specKCWI_snr10 = path_fits_DF44+"fit_DF44_aD1_specKCWI_massfixed_maxsnr10_zsortafit_npoly8_1635282112_mcmc.h5"
fit_aD1_phot_snr10 = None

fit_aD1_phot_specKCWI_snr15 = path_fits_DF44+"fit_DF44_aD1_phot_specKCWI_zsortafit_npoly8_maxsnr15_1649708931_mcmc.h5"
fit_aD1_specKCWI_snr15 = path_fits_DF44+"fit_DF44_aD1_specKCWI_massfixed_maxsnr15_zsortafit_npoly8_1634850512_mcmc.h5"
fit_aD1_phot_snr15 = path_fits_DF44+"fit_DF44_aD1_phot_1649708925_mcmc.h5"
fit_aD1_phot_snr15 = None

fit_aD1_phot_specKCWI_snr20 = path_fits_DF44+"fit_DF44_aD1_phot_specKCWI_zsortafit_npoly8_maxsnr20_1649879779_mcmc.h5"
fit_aD1_specKCWI_snr20 = path_fits_DF44+"fit_DF44_aD1_phot_specKCWI_massfixed_zsortafit_npoly8_maxsnr20_1634689586_mcmc.h5"
fit_aD1_phot_snr20 = None

# MANGA data instead of KCWI
fit_aD1_specMANGA_photgi = None


# alphaD = 0.2, "bursty" or "concentrated" version of the SF prior

fit_aD02_phot_specKCWI    = path_fits_DF44+"fit_DF44_aD02_phot_specKCWI_zsortafit_npoly8_1649878702_mcmc.h5"
fit_aD02_phot_specKCWI_v2 = path_fits_DF44+"fit_DF44_aD02_phot_specKCWI_zfixed_npoly3_1687841283_mcmc.h5"
fit_aD02_specKCWI = None
fit_aD02_phot = path_fits_DF44+'fit_DF44_aD02_phot_1687841281_mcmc.h5'

# minimum age of time bins is 1 Gyr (i.e., suppress late burst of SF)
fit_aD02_min1Gyr_phot_specKCWI    = path_fits_DF44+"fit_DF44_aD02_older1Gyr_phot_specKCWI_zsortafit_npoly8_1649979206_mcmc.h5"
fit_aD02_min1Gyr_phot_specKCWI_v2 = path_fits_DF44+"fit_DF44_aD02_older1Gyr_phot_specKCWI_zfixed_npoly3_1687843283_mcmc.h5"
fit_aD02_min1Gyr_specKCWI = None
fit_aD02_min1Gyr_phot = None

# alphaD = 0.3

fit_aD03_phot_specKCWI = None
fit_aD03_specKCWI = None
fit_aD03_phot = path_fits_DF44+"fit_DF44_aD03_phot_1649712254_mcmc.h5"

# alphaD = 0.5

fit_aD05_phot_specKCWI    = path_fits_DF44+"fit_DF44_aD05_phot_specKCWI_zsortafit_npoly8_1649707131_mcmc.h5"
fit_aD05_phot_specKCWI_v2 = path_fits_DF44+"fit_DF44_aD05_phot_specKCWI_zfixed_npoly3_1687841282_mcmc.h5"
fit_aD05_specKCWI = None
fit_aD05_phot = path_fits_DF44+"fit_DF44_aD05_phot_1687841280_mcmc.h5"

#################################################################
# Nonparametric SFH model
# constant SFR Continuity prior ( StudentT(0,0.3,2) prior on log-ratio of SFR of adjacent time bins)
#################################################################

fit_csfrcont_phot_specKCWI = None
fit_csfrcont_phot_specKCWI = path_fits_DF44+"fit_DF44_csfrcont_phot_specKCWI_zsortafit_npoly8_1689720284_mcmc.h5"
fit_csfrcont_phot_specKCWI_v2 = path_fits_DF44+"fit_DF44_csfrcont_phot_specKCWI_zfixed_npoly3_1687841282_mcmc.h5"
fit_csfrcont_specKCWI = None
fit_csfrcont_phot = path_fits_DF44+"fit_DF44_csfrcont_phot_1649697330_mcmc.h5"

#################################################################
# Parametric SFH model, delayed declining tau model
#################################################################

# log-uniform prior on tau
fit_dexp_phot_specKCWI = path_fits_DF44+"fit_DF44_dexp_phot_specKCWI_zsortafit_npoly8_1649700118_mcmc.h5"
# fit_dexp_phot_specKCWI = path_fits_DF44+"fit_DF44_dexp_phot_specKCWI_zsortafit_npoly8_1649962577_mcmc.h5"
# fit_dexp_phot = path_fits_DF44+"fit_DF44_dexp_phot_1649697850_mcmc.h5"
fit_dexp_phot = path_fits_DF44+"fit_DF44_dexp_phot_1649971330_mcmc.h5"

# uniform prior on tau
fit_dexp_utau_phot_specKCWI = None
fit_dexp_utau_phot = path_fits_DF44+"fit_DF44_dexp_utau_phot_1649697872_mcmc.h5"


# variants
# fit_DF44_aD02_phot_specKCWI_fixeddust001_zfit_npoly8_1650297773_mcmc
