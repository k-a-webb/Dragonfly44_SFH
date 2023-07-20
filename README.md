# Dragonfly44_SFH

This repository contains all the data used to analyze the star formation history (SFH) of Dragonfly 44 (DF44) as discussed in Webb+2022 (https://ui.adsabs.harvard.edu/abs/2022MNRAS.516.3318W/abstract) based on high-S/N and high-resolution rest-frame optical KCWI spectroscopy and NUV--NIR photometry.

#### Analysis presented in Webb+2022

The scripts used to generate the figures in Webb+2022 are provided in
> analysis/Webb2022_figures/*

Note that some figures have slight errors, as noted in
>  analysis/Webb2022_figures/errata.txt

These errors have been corrected in the notebooks shared here.

#### Additional analysis for Dragonfly 44

as shown in various notebooks in
> analysis/

Additional fits with different SFH models/priors:
- nonparametic model with constant SFR "continuity" prior
- nonparametic model with Dirichlet prior with other $\alpha_\mathrm{D}$ values
- delayed declining exponential SFR model


#### Requirements:

- numpy (v1.25.1)
- matplotlib (v3.7.2)
- scipy (v1.11.1)
- astropy (v5.3.1)
- pandas (v2.0.3)
- seaborn (v0.12.2)
- [corner](https://corner.readthedocs.io/en/latest/) (v2.2.2)
- [ipython](https://ipython.org/install.html), [Jupyter notebook](https://jupyter.org/install) (v4.0.3)
- [Prospector](https://prospect.readthedocs.io/en/latest/) (v1.2.0):
  - [and dependent packages for Prospector](https://prospect.readthedocs.io/en/latest/installation.html#requirements)
  - h5py (v3.9.0)
  - dynesty (v2.1.2)



Setup details:
- In order to construct the filter objects used by Prospector, via sedpy, copy filter files for Dragonfly44 from
  > data/Dragonfly44/sedpy_filters/

  to the appropriate sedpy folder

  > sedpy/data/filters/

#### Acknowledgments
