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

fit_aD1_phot_specKCWI = "fit_aD1_1649733960_mcmc.h5"
fit_aD1_specKCWI = "fit_aD1_specKCWI_phot_specKCWI_fixmass_1634082827_mcmc.h5" # fixed mass
fit_aD1_specKCWI_freemass = "fit_aD1_specKCWI_phot_specKCWI_1634071305_mcmc.h5" # free mass
fit_aD1_phot = "fit_aD1_phot_1649889938_mcmc.h5"

# minium age of time bins is 1 Gyr (i.e., suppress late burst of SF)
fit_aD1_min1Gyr_phot_specKCWI = "fit_aD1_old_1649879742_mcmc.h5"
fit_aD1_min1Gyr_specKCWI = None
fit_aD1_min1Gyr_phot = None

# degrade S/N of spectroscopy (i.e., increase uncertainties) to test S/N dependence of the results
fit_aD1_phot_specKCWI_snr5 = "fit_aD1_maxsnr5_1649879779_mcmc.h5"
fit_aD1_specKCWI_snr5 = "fit_aD1_specKCWI_phot_specKCWI_fixmass_maxsnr5_1635282103_mcmc.h5"

fit_aD1_phot_specKCWI_snr10 = "fit_aD1_maxsnr10_1649879768_mcmc.h5"
fit_aD1_specKCWI_snr10 = "fit_aD1_specKCWI_phot_specKCWI_fixmass_maxsnr10_1635282112_mcmc.h5"

fit_aD1_phot_specKCWI_snr15 = "fit_aD1_maxsnr15_1649708931_mcmc.h5"
fit_aD1_specKCWI_snr15 = "fit_aD1_specKCWI_phot_specKCWI_fixmass_maxsnr15_1634850512_mcmc.h5"
fit_aD1_phot_snr15 = "fit_aD1_maxsnr15_phot_1649708925_mcmc.h5"

fit_aD1_phot_specKCWI_snr20 = "fit_aD1_maxsnr20_1649879779_mcmc.h5"
fit_aD1_specKCWI_snr20 = "fit_aD1_specKCWI_phot_specKCWI_fixmass_maxsnr20_1634689586_mcmc.h5"

# MANGA data instead of KCWI
fit_aD1_specMANGA_photgi = None


# alphaD = 0.2, "bursty" or "concentrated" version of the SF prior

fit_aD02_phot_specKCWI = "fit_aD02_1649878702_mcmc.h5"
fit_aD02_specKCWI = None
fit_aD02_phot = None

# minium age of time bins is 1 Gyr (i.e., suppress late burst of SF)
fit_aD02_min1Gyr_phot_specKCWI = "fit_aD02_old_1649979206_mcmc.h5"
fit_aD02_min1Gyr_specKCWI = None
fit_aD02_min1Gyr_phot = None

# alphaD = 0.3

fit_aD03_phot_specKCWI = None
fit_aD03_specKCWI = None
fit_aD03_phot = "fit_aD03_phot_1649712254_mcmc.h5"

# alphaD = 0.5

fit_aD05_phot_specKCWI = "fit_aD05_1649707131_mcmc.h5"
fit_aD05_specKCWI = None
fit_aD05_phot = "fit_aD05_phot_1649696855_mcmc.h5"

#################################################################
# Nonparametric SFH model
# constant SFR Continuity prior ( StudentT(0,0.3,2) prior on log-ratio of SFR of adjacent time bins)
#################################################################

fit_csfrcont_phot_specKCWI = "fit_csfrcont_1657666517_mcmc.h5"
# fit_csfrcont_phot_specKCWI_temp = "fit_csfrcont_mask4733_mask5170_flogzsol_fixsmooth_1627320346_mcmc.h5"
fit_csfrcont_specKCWI = ""
fit_csfrcont_phot = "fit_csfrcont_phot_1649697330_mcmc.h5"

#################################################################
# Parametric SFH model, delayed declining tau model
#################################################################

# log-uniform prior on tau
fit_dexp_phot_specKCWI = "fit_dexp_1649700118_mcmc.h5"
# fit_dexp_phot_specKCWI = "fit_dexp_1649962577_mcmc.h5"
# fit_dexp_phot = "fit_dexp_phot_1649697850_mcmc.h5"
fit_dexp_phot = "fit_dexp_phot_1649971330_mcmc.h5"

# uniform prior on tau
fit_dexp_utau_phot_specKCWI = None
fit_dexp_utau_phot = "fit_dexp_utau_phot_1649697872_mcmc.h5"
