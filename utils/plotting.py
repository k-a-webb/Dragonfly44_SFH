import numpy as np
import matplotlib.pyplot as plt
from copy import deepcopy
from matplotlib.legend_handler import HandlerLine2D, HandlerTuple

from prospect.plotting.figuremaker import *
from prospect.plotting import pretty
from prospect.plotting.corner import  get_spans, allcorner

saveparams = {'bbox_inches':'tight', 'pad_inches':0.05, 'dpi':200, 'format':'pdf'}

fig_width_one = 242.26653/72.*2. *0.9
fig_width_two = 513.11743/72.*2. *0.9
textheight = 657.3189/72.*2. *0.9

import matplotlib as mpl
# https://matplotlib.org/stable/tutorials/introductory/customizing.html
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['font.size'] = 12
mpl.rcParams['axes.labelsize'] = 12
mpl.rcParams['legend.fontsize'] = 12
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['xtick.labelsize'] = 10

from Dragonfly44_SFH.utils.prospect_io import get_models_from_posterior_draws, get_models_for_bestfits

def plot_photometry( ax, obs, bestfits, draws=None,
                     obs_params={'marker':'D', 'color':'c', 'ms':12},
                     posts_params={'color':'goldenrod', },
                     bestfit_params={'lw':2.4, 'color':'k', 'marker':'o', 's':160},
                     legend_items={},
                     plot_draws=True,
                     **extras ):

    if obs['maggies'] is not None:

        # copy information from obs
        x, y, ey, mask = [ np.copy( obs[key] )
                           for key in [ "phot_wave","maggies","maggies_unc","phot_mask" ]]

        # plot masked observations, with low alpha
        l = ax.errorbar( x[~mask], y[~mask], yerr=ey[~mask],
                         mec=obs_params['color'], mfc='None', fmt=' ',
                         capsize=2, zorder=2, mew=0.7, alpha=0.2, **obs_params)
        legend_items['obs']['phot_masked'] = dict( handle=l, label='Observed photometry, masked' )

        # plot unmasked observations
        l = ax.errorbar( x[mask], y[mask], yerr=ey[mask],
                         mec=obs_params['color'], mfc='None', fmt=' ',
                         capsize=2, zorder=2, mew=0.7, alpha=1,
                         **obs_params)
        legend_items['obs']['phot'] = dict( handle=l, label='Observed photometry' )

        # plot bestfit photometry model for masked points
        l = ax.scatter( x[~mask], bestfits['phot'][~mask],
                        edgecolor=bestfit_params['color'], facecolor='None',
                        zorder=3, alpha=0.2, **bestfit_params)
        legend_items['bestfit']['phot_masked'] = dict( handle=l, label='Bestfit photometry, masked' )

        # plot bestfit photometry model for unmasked points
        l = ax.scatter( x[mask], bestfits['phot'][mask],
                        edgecolor=bestfit_params['color'], facecolor='None',
                        zorder=3, alpha=1, **bestfit_params)
        legend_items['bestfit']['phot'] = dict( handle=l, label='Bestfit photometry' )

    try:
        ax.plot( bestfits['wave_sed'], bestfits['sed'],
                  color=bestfit_params['color'], lw=0.5, zorder=0 )
    except:
        # sometimes the 'sed' entry is nonsense, requires post-processing step... not always done
        pass

    # plot SED model draws from posterior, quantiles
    if plot_draws and (draws is not None):
        l = ax.fill_between( bestfits['wave_sed'], draws['sed2'][0], draws['sed2'][1],
                             color=posts_params['color'], alpha=0.4, lw=0, label=None )
        legend_items['draws'] = {'handle':l, 'label':'68% CRs'}


    ax.set( xscale='log',
            yscale='log',
            ylabel=r'$\nu f_\nu$ (maggies)',
          )

    return ax, legend_items

def plot_photometry_residuals( ax, obs, bestfits,
                     draws=None, plot_draws=True,
                     obs_params={'marker':'D', 'color':'c', 'ms':12},
                     posts_params={'color':'goldenrod'},
                     bestfit_params={'lw':2.4, 'color':'k', 'marker':'o', 's':160},
                     legend_items={},
                     draws_quantiles=[0.16,0.84],
                     **extras ):

    if obs['maggies'] is None: return ax, legend_items

    # copy information from obs
    x, y, ey, mask = [ np.copy( obs[key] )
                       for key in [ "phot_wave","maggies","maggies_unc","phot_mask" ]]

    # draw horizontal line on residual plot at zero (ideal residuals)
    ax.axhline( 0, color='k', ls='--', zorder=0.1)

    # calculate residual of obsevations and bestfit
    chi = ( y - bestfits['phot'] ) /ey
    chi2_phot = np.sum(chi[mask]**2)/np.sum(mask).astype(float)

    # plot residuals
    bestfit_params_dupe = deepcopy( bestfit_params )
    if 's' in bestfit_params_dupe.keys():
        bestfit_params_dupe['s'] = bestfit_params_dupe['s']/2.
    ax.scatter( x[mask], chi[mask],
                 edgecolor=bestfit_params_dupe['color'], facecolor='None', zorder=2,
                 alpha=1, **bestfit_params_dupe)

    # plot SED model draws from posterior, quantiles
    if plot_draws and (draws is not None):
        chis = (y - draws['phot']) /ey
        qs = np.quantile( chis, np.sort([0.5]+draws_quantiles), axis=0 )
        dqs = np.diff( qs, axis=0 )

        ax.errorbar( x[mask], qs[1], yerr=dqs, \
                      color=bestfit_params_dupe['color'], mec=bestfit_params_dupe['color'],
                      elinewidth=3, marker='.',
                      ms=0, fmt=' ', capsize=7, zorder=-10, lw=1.5, alpha=0.7 )


    ymax = np.max( np.abs( chi[mask] )) # calculate axis ylim based on residuals, make symmetric
    ax.set( xscale='log',
            xlabel=u'Observed wavelength (\u00c5)',
            ylabel=r'$\chi_\mathrm{bestfit}$',
            ylim=(-ymax/0.8,ymax/0.8),
           )

    return ax, chi2_phot

def plot_spectroscopy( ax, obs, bestfits, draws=None, plot_draws=True, zobs=0,
                     obs_params={'marker':'D', 'color':'c', 'ms':12},
                     posts_params={'color':'goldenrod', },
                     bestfit_params={'lw':2.4, 'color':'k', 'marker':'o', 's':160},
                     legend_items={},
                     **extras ):

    if obs['spectrum'] is None: return ax, legend_items
    if 'lw' not in bestfit_params.keys(): bestfit_params['lw'] = 1
    if 'color' not in bestfit_params.keys(): bestfit_params['color'] = 'r'

    x, y, ey, mask = [ np.copy( obs[key] ) for key in [ "wavelength","spectrum","unc","mask" ]]
    x /= (1.+zobs) # observed-frame to rest-frame

    # make a second mask to use when plotting bestfit spectrum, spanning min- max-wavelength
    mask_bf = np.ones_like(x, dtype=bool)
    mask_bf[ (x < x[mask][0])|(x > x[mask][-1]) ] = False

    # y[~mask] = np.nan

    # plot uncertainty region of observed spectrum
    l = ax.fill_between( x, y-ey, y+ey,
                         step='post', color='lightgrey', lw=0, zorder=0)
    legend_items['obs']['spec_masked'] = dict( handle=l, label='Observed spectroscopy, masked' )

    # plot obseved spectrum
    l, = ax.step( x, y,
               where='post', lw=1, color=obs_params['color'], zorder=3,)
    legend_items['obs']['spec'] = dict( handle=l, label='Observed spectroscopy' )

    # plot (normalized) bestfit spectrum
    l, = ax.step( x[mask_bf], bestfits['spec'][mask_bf],
               where='post', color=bestfit_params['color'], lw=bestfit_params['lw'], zorder=3 )
    legend_items['bestfit']['spec'] = dict( handle=l, label='Observed spectroscopy' )

    if plot_draws and (draws is not None):
     l = ax.fill_between( x[mask_bf],
                            draws['quantiles']['spec'][0,:][mask_bf], draws['quantiles']['spec'][1,:][mask_bf],
                            color=posts_params['color'], alpha=0.2, zorder=2, lw=0)
     legend_items['draws']['spec'] = dict( handle=l, label='Spectroscopy' )

    ax.set( xscale='linear',
         yscale='linear',
         ylabel=r'$\nu f_\nu$ (normalized)',
         xticklabels=[], \
         xlim=(x[mask][0]*0.99, x[mask][-1]/0.99),
         ylim=( np.min((y-ey)[mask])*0.99, np.max((y+ey)[mask])/0.99),
       )

    return ax, legend_items

def plot_spectroscopy_residuals( ax, obs, bestfits, draws=None, plot_draws=True, zobs=0,
                                 obs_params={'marker':'D', 'color':'c', 'ms':12},
                                 posts_params={'color':'goldenrod', },
                                 bestfit_params={'lw':2.4, 'color':'k', 'marker':'o', 's':160},
                                 legend_items={},
                                 draws_quantiles=[0.16,0.84],
                                 **extras ):

    if obs['spectrum'] is None: return ax, legend_items

    x, y, ey, mask = [ np.copy( obs[key] ) for key in [ "wavelength","spectrum","unc","mask" ]]
    x /= (1.+zobs) # observed-frame to rest-frame

    # make a second mask to use when plotting bestfit spectrum, spanning min- max-wavelength
    mask_bf = np.ones_like(x, dtype=bool)
    mask_bf[ (x < x[mask][0])|(x > x[mask][-1]) ] = False

    y[~mask] = np.nan

    ax.axhline( 0, color='k', ls='--', lw=0.5, zorder=0)

    chi = (y-bestfits['spec'])/ey
    chi[~mask] = np.nan
    chi2_spec = np.nansum( chi[mask]**2 ) / np.sum(mask).astype(float)

    ax.step( x[mask], chi[mask],
           where='post', color=bestfit_params['color'], lw=bestfit_params['lw'], zorder=3)

    if plot_draws and (draws is not None):
     chis = (y-draws['spec'])/ey
     qs_chi = np.quantile( chis, q=draws_quantiles, axis=0)
     ax.fill_between( x[mask], qs_chi[0,:][mask], qs_chi[1,:][mask],
                      color=posts_params['color'], alpha=0.2, zorder=2, lw=0)

    ax.set( ylabel=r'$\chi_\mathrm{bestfit}$',
        xticklabels=[],
      )

    return ax, chi2_spec

def plot_speccal( ax, obs, bestfits, draws=None, plot_draws=True, zobs=0, norm=1e8,
                 obs_params={'marker':'D', 'color':'c', 'ms':12},
                 posts_params={'color':'goldenrod'},
                 bestfit_params={'lw':2.4, 'color':'k', 'marker':'o', 's':160},
                 legend_items={},
                 **extras ):

    if obs['spectrum'] is None: return ax, legend_items

    x, y, ey, mask = [ np.copy( obs[key] ) for key in [ "wavelength","spectrum","unc","mask" ]]
    x /= (1.+zobs) # observed-frame to rest-frame

    # make a second mask to use when plotting bestfit spectrum, spanning min- max-wavelength
    mask_bf = np.ones_like(x, dtype=bool)
    mask_bf[ (x < x[mask][0])|(x > x[mask][-1]) ] = False

    ax.step( x[mask_bf], bestfits['speccal'][mask_bf]/norm,
          where='post', color=bestfit_params['color'], lw=1.5, zorder=3)

    # if plot_draws and (draws is not None):
    #     ax.fill_between( x[mask_bf],
    #                      draws['quantiles']['speccal'][0,:][mask_bf]/norm,
    #                      draws['quantiles']['spec'][1,:][mask_bf]/norm,
    #                      color=posts_params['color'], alpha=0.2, zorder=0, lw=0 )


    ax.set(
        ylabel='Calibration\nPolynomial\n'+r'($\times 10^{{{:.0f}}}$)'.format( np.log10(norm) ),
      )

    return ax

def plot_relative_spectroscopy( ax, obs1, obs2, spec1, spec2, zobs=0,
                                params={'color':'k', 'lw':1.5}, reverse_order=False, **extras ):

    x, mask1 = [ np.copy( obs1[k]) for k in ['wavelength','mask'] ]
    mask2 = np.copy( obs2['mask'] )
    x /= (1.+zobs) # observed-frame to rest-frame

    x[~mask1] = np.nan
    x[~mask2] = np.nan

    ax.axhline( 1, color='k', ls='--', lw=0.5, zorder=0)

    y = spec1 / spec2
    if reverse_order: y = 1./y

    ax.plot( x, y, **params )

    ax.set_yticks( np.arange(0.9,1.1,0.005), minor=True )
    ax.set_ylim(0.992,1.013)
    ax.set_ylabel( "Relative\nchange\nof bestfits" )

    return ax


def plot_obs_and_fits(result, obs, zobs,
             # col_models="c", col_bf='limegreen', col_obs='k',
             fig=None, axes=None, figsize=(10,.5*10), return_fig=False,
             outfname=None,
             xlim_phot = (1.5e3,7e4),
             draws_quantiles=[0.16,0.84],
             obs_params={'marker':'D', 'color':'c', 'ms':12},
             posts_params={'color':'goldenrod'},
             bestfit_params={'lw':2.4, 'color':'k', 'marker':'o', 's':160},
             plot_legend=True, show_chi2s=True, label=None,
             **extras,
            ):
    """
    Plot the observations and model SEDs
    If model SEDs are saved in the output in result['draws'], plot (plot_draws = True)
    If bestfit model SEDs is saved in result['draws']['bestfit'], plot
        Note that there also exists result['bestfit'], however this is sometimes buggy

    Lots of hardcoded things in this function.
    e.g., number of axes determined by which observations are included in fit (phot, spec, or both)

    res: result dictionary
    obs: observation dictionary

    col_models: colour of models
    col_bf: colour of bestfit models
    col_obs: colour of observations

    fig: figure element, if not specified with use default
    figsize: dimensions of default figure

    return_fig: flag to indicate whether to close plot (plt.show) or return

    outfname: if not None, save plot with this name

    xlim_phot: hardcoded xlim of photometry axis, used when photometry not included in fit
    draws_quantiles: when plotting SED models drawn from the posterior, show region between these quantiles
    """

    if True: # draws and bestfits, check if exist, make flags based on checks
        plot_draws = True
        if "draws" not in result.keys():
            # print("  draws not in file")
            plot_draws = False
        elif "bestfit" not in result['draws'].keys():
            # print("  draws/bestfit not in file")
            plot_draws = False

        draws = get_models_from_posterior_draws( plot_draws, result, draws_quantiles, **extras )
        bestfits = get_models_for_bestfits( result )

    if True: # setup figure and axes
        # if figure not already defined, and spectrum included in fit, make figure with six rows
        if (obs['spectrum'] is not None) and (fig is None):
            fig, axs = plt.subplots(6,1, figsize=figsize,
                                    gridspec_kw={'height_ratios':[2.5,1,1.1,2.5,1,1]}, sharex=False)
            ax_phot,ax_rphot,axx,ax_spec,ax_rspec,ax_pspec = axs.flatten()
            axx.axis('off')
        # else, figure MUST have six rows, and we now label each
        elif (obs['spectrum'] is not None) and (fig is not None):
            if axes is None:
                ax_phot,ax_rphot,_,ax_spec,ax_rspec,ax_pspec = fig.get_axes()[:6]
            else:
                ax_phot,ax_rphot,_,ax_spec,ax_rspec,ax_pspec = axes[:6]
        # if figure not already defined, and spectrum NOT included in fit, make figure with two rows
        elif (obs['spectrum'] is None) and (fig is None):
            if axes is None:
                fig, axs = plt.subplots(2,1, figsize=(figsize[0],figsize[1]/2.),
                                        gridspec_kw={'height_ratios':[3.5,1]}, sharex=True)
                ax_phot,ax_rphot = axs.flatten()
            else:
                ax_phot,ax_rphot = axes
        # else, figure MUST have two rows, and we now label each
        elif (obs['spectrum'] is None) and (fig is not None):
            if axes is None:
                ax_phot,ax_rphot = fig.get_axes()
            else:
                ax_phot,ax_rphot = axes
            # label all unused axes as None
            axx,ax_spec,ax_rspec,ax_pspec = [ None for i in range(4) ]

    # set parameter for annotations
    params_ann = {"xycoords":"axes fraction", "xy":(0,0), "fontsize":12, "va":"top"}

    legend_items = { 'obs':{}, 'bestfit':{}, 'draws':{} }

    if True: # plot photometry

        # if photometry included in the observations, plot it and bestfit
        if (obs['maggies'] is not None):

            ax_phot, legend_items = plot_photometry( ax=ax_phot, obs=obs, bestfits=bestfits,
                                                 plot_draws=plot_draws, draws=draws,
                                                  obs_params=obs_params,
                                                  posts_params=posts_params,
                                                  bestfit_params=bestfit_params,
                                                  legend_items=legend_items,
                                                   )
            ax_rphot, chi2_phot = plot_photometry_residuals( ax=ax_rphot, obs=obs, bestfits=bestfits,
                                                  obs_params=obs_params,
                                                  posts_params=posts_params,
                                                  bestfit_params=bestfit_params,
                                                  legend_items=legend_items,
                                                   )

            # label chi2
            label_chi2_phot = r"Bestfit $\chi^2$/N$_\mathrm{phot}$="+"{:.2f}".format( chi2_phot )
            if show_chi2s:
                ax_phot.annotate( label_chi2_phot, xytext=(0.01,0.97), **params_ann)
            else: print( label, label_chi2_phot )

            ax_phot.set( xticklabels=[] )
            ax_rphot.set_xlim( ax_phot.get_xlim())

        # if photometry NOT included in the observations, just plot bestfit
        else:
            ax_rphot.axis('off') # turn off residual axis, since won't include anything
            # find ylim based on limits of bestfit in xlim_phot
            y_in_xlim = bestfits['sed'][ ( xlim_phot[0] <= bestfits['wave_sed']) & (bestfits['wave_sed'] <= xlim_phot[1])]
            ylim = ( np.min(y_in_xlim), np.max(y_in_xlim) )
            ax_phot.set( xlabel=u'Observed wavelength (\u00c5)' )

    if True: # spectrum

        # if spectrum included in observations, plot
        if (obs['spectrum'] is not None):


            ax_spec, legend_items = plot_spectroscopy( ax=ax_spec, obs=obs, bestfits=bestfits,
                                                  draws=draws, zobs=zobs,
                                                  plot_draws=plot_draws,
                                                  obs_params=obs_params,
                                                  posts_params=posts_params,
                                                  bestfit_params=bestfit_params,
                                                  legend_items=legend_items,
                                                   )
            ax_rspec, chi2_spec = plot_spectroscopy_residuals( ax=ax_rspec, obs=obs, bestfits=bestfits, draws=draws, zobs=zobs,
                                                  plot_draws=plot_draws,
                                                  obs_params=obs_params,
                                                  posts_params=posts_params,
                                                  bestfit_params=bestfit_params,
                                                  legend_items=legend_items,
                                                  draws_quantiles=draws_quantiles,
                                                   )
            ax_pspec = plot_speccal( ax=ax_pspec, obs=obs, bestfits=bestfits, draws=draws, zobs=zobs,
                                                  plot_draws=plot_draws,
                                                  obs_params=obs_params,
                                                  posts_params=posts_params,
                                                  bestfit_params=bestfit_params,
                                                  legend_items=legend_items,
                                                   )


            label_chi2_spec = r"Bestfit $\chi^2$/N$_\mathrm{spec}$="+"{:.2f}".format( chi2_spec )
            if show_chi2s: # label chi2
                ax_spec.annotate( label_chi2_spec, xytext=(0.01,0.97), **params_ann)
            else: print( label, label_chi2_spec )

            ax_spec.set( xticklabels=[] )
            ax_rspec.set( xticklabels=[] )
            ax_rspec.set_xlim(ax_spec.get_xlim())
            ax_pspec.set_xlim(ax_spec.get_xlim())

            if zobs>0: xlabel=u'Rest wavelength (\u00c5)'
            else: xlabel=u'Observed wavelength (\u00c5)'
            ax_pspec.set_xlabel( xlabel )

    if zobs > 0: # add second x-axis with observed and rest wavelengths
        if ax_spec is not None:
            xlim = ax_spec.get_xlim()
            xlim_rf = np.array(xlim)*(1+zobs)
            ax_spect = ax_spec.twiny()
            ax_spect.set_xlim( xlim_rf )

            ax_pspec.set_yscale('log')

            x = np.copy( obs['wavelength'] )
            mask = np.copy( obs['mask'] )

            if (x[mask][0] > 4500):
                ax_spect.set_xticks( np.arange( 4700,5400+1, 100) )
                [ ax.set_xticks( np.arange( 4600,5300+1, 100) ) for ax in [ax_spec,ax_rspec,ax_pspec] ]
                [ ax.set_xticks( np.arange( 4600,5300+1, 50), minor=True ) for ax in [ax_spec,ax_rspec,ax_pspec] ]
            else:
                [ ax.xaxis.set_major_locator(MultipleLocator(100)) for ax in [ax_spec,ax_spect,ax_rspec,ax_pspec] ]
                [ ax.xaxis.set_minor_locator(MultipleLocator(10)) for ax in [ax_spec,ax_spect,ax_rspec,ax_pspec] ]

    if plot_legend: # make legend
        if  (obs['spectrum'] is not None):
            leg_params = dict( numpoints=1, loc=4, handler_map={tuple: HandlerTuple(ndivide=None)}, markerscale=0.7 )
            if plot_draws and (obs['maggies'] is not None) and (obs['spectrum'] is not None):
                legend_handles = [ ( legend_items['obs']['phot']['handle'], legend_items['obs']['spec']['handle'] ),
                                   ( legend_items['bestfit']['phot']['handle'], legend_items['bestfit']['spec']['handle'] ),
                                   legend_items['draws']['handle'],
                                 ]
                legend_labels = ['Observations',"Bestfit model","68% CR of randomly drawn models"]
                ax_phot.legend( legend_handles, legend_labels, **leg_params)

            elif plot_draws and (obs['maggies'] is None):
                legend_handles = [ legend_items['obs']['spec']['handle'],
                                   legend_items['bestfit']['spec']['handle'],
                                   legend_items['draws']['handle'],
                                 ]
                legend_labels = ['Observations',"Bestfit model","68% CR of randomly drawn models"]
                ax_phot.legend( legend_handles, legend_labels, **leg_params)
            elif not plot_draws and (obs['maggies'] is not None):
                legend_handles = [ ( legend_items['obs']['phot']['handle'], legend_items['obs']['spec']['handle'] ),
                                   ( legend_items['bestfit']['phot']['handle'], legend_items['bestfit']['spec']['handle'] ),
                                 ]
                legend_labels = ['Observations',"Bestfit model"]
                ax_phot.legend( legend_handles, legend_labels, **leg_params)
            else:
                legend_handles = [ legend_items['obs']['phot']['handle'],
                                   legend_items['bestfit']['phot']['handle'],
                                   legend_items['draws']['handle'],
                                 ]
                legend_labels = ['Observations',"Bestfit model"]
                ax_phot.legend( legend_handles, legend_labels, **leg_params)

    fig.subplots_adjust(wspace=0.05, hspace=0.05)
    if outfname is not None:
        plt.savefig(outfname, **saveparams)
    if return_fig:
        return fig
    else:
        plt.show()


def compare_two_fits(plot_params1, plot_params2, draws_quantiles=[0.16,0.84], zobs=0,
                     figsize=(12,12), reverse_order=False, **extras):
    """
    Plot and compare two fits (to the same data, including both photometry and spectroscopy)
    """

    for plot_params in [plot_params1, plot_params2]:
        if 'draws' not in plot_params:
            plot_params['draws'] = get_models_from_posterior_draws( True, plot_params['result'], draws_quantiles=draws_quantiles )
        if 'bestfits' not in plot_params:
            plot_params['bestfits'] = get_models_for_bestfits( plot_params['result'] )


    fig, axes = plt.subplots(2+1+3+1,1, figsize=figsize, \
                        gridspec_kw={'height_ratios':[2.5,1, 0.4, 3,1,1,1]}, sharex=False)

    ax_phot,ax_rphot, aoff, ax_spec,ax_rspec,ax_pspec, ax_cspec = axes

    aoff.axis('off')

    # set parameter for annotations
    params_ann = {"xycoords":"axes fraction", "xy":(0,0), "fontsize":12, "va":"top"}

    legend_items1 = { 'obs':{}, 'bestfit':{}, 'draws':{} }
    legend_items2 = { 'obs':{}, 'bestfit':{}, 'draws':{} }

    for legend_items, plot_params in [[legend_items1,plot_params1],[legend_items2,plot_params2]]:

        ax_phot, legend_items1 = plot_photometry( ax=ax_phot,
                                                   legend_items=legend_items,
                                                   **plot_params,
                                                  )
        ax_rphot, chi2_phot = plot_photometry_residuals( ax=ax_rphot, **plot_params, )

        ax_spec, legend_items = plot_spectroscopy( ax=ax_spec,
                                                   legend_items=legend_items,
                                                   **plot_params,
                                                 )
        ax_rspec, chi2_spec = plot_spectroscopy_residuals( ax=ax_rspec,
                                                           draws_quantiles=draws_quantiles,
                                                           **plot_params,
                                                         )
        ax_pspec = plot_speccal( ax=ax_pspec, **plot_params, )

    if True: # axis limits and labels
        ax_phot.set_ylim(1e-12,1e-7)
        ax_phot.set_xlim(1.5e3,1e5)

        ax_phot.set( xticklabels=[] )
        ax_rphot.set_xlim(ax_phot.get_xlim())

        [ ax.set( xticklabels=[] ) for ax in [ax_spec, ax_rspec, ax_pspec ] ]
        [ ax.set_xlim( ax_spec.get_xlim() ) for ax in [ax_rspec, ax_pspec, ax_cspec ] ]

        if zobs>0: xlabel=u'Rest wavelength (\u00c5)'
        else: xlabel=u'Observed wavelength (\u00c5)'
        ax_cspec.set_xlabel( xlabel )

        for ax in axes:
            ax.tick_params(direction='out', length=5, which='major')
            ax.tick_params(direction='out', length=3, which='minor')

        from matplotlib.ticker import MultipleLocator
        for ax in [ax_spec,ax_rspec,ax_pspec,ax_cspec]:
            ax.xaxis.set_major_locator(MultipleLocator( 100 ))
            ax.xaxis.set_minor_locator(MultipleLocator( 20 ))


    if True: # compare bestfit spectra

        ax_cspec = plot_relative_spectroscopy( ax=ax_cspec,
                                               obs1=plot_params1['obs'],
                                               obs2=plot_params2['obs'],
                                               spec1=plot_params1['bestfits']['spec'],
                                               spec2=plot_params2['bestfits']['spec'],
                                               zobs=zobs, reverse_order=reverse_order,
                                            )

    if True: # top legend
        ax_legend = ax_phot
        from matplotlib.legend_handler import HandlerLine2D, HandlerTuple


        l1a = ax_legend.errorbar([-1],[-1], yerr=[1], color=plot_params1['obs_params']['color'], mec=plot_params1['obs_params']['color'], mfc='None', fmt=' ', capsize=2, ms=12, marker='o', mew=0.7)
        l1b, = ax_legend.plot([],[], color=plot_params1['obs_params']['color'])
        l2a = ax_legend.scatter([],[], edgecolor=plot_params1['bestfit_params']['color'], facecolor='None', marker=plot_params1['bestfit_params']['marker'], zorder=2, s=160, lw=1.5, )
        l2b, = ax_legend.plot([],[], color=plot_params1['bestfit_params']['color'])
        l3a = ax_legend.scatter([],[], edgecolor=plot_params2['bestfit_params']['color'], facecolor='None', marker=plot_params2['bestfit_params']['marker'], zorder=2, s=160, lw=1.5, )
        l3b, = ax_legend.plot([],[], color=plot_params2['bestfit_params']['color'])
        l4 = ax_legend.fill_between([],[], lw=0, color=plot_params1['bestfit_params']['color'], alpha=0.3)
        l5 = ax_legend.fill_between([],[], lw=0, color=plot_params2['bestfit_params']['color'], alpha=0.3)

        handles = [(l1a,l1b),(l2a,l2b),(l3a,l3b),(l4,l5)]
        labels = ['Observations','','',"68% CRs"]
        for i,pp in enumerate( [plot_params1,plot_params2] ):
            if ('label' in pp['bestfit_params'].keys()):
                labels[i+1] = pp['bestfit_params']['label']
            elif ('label' in pp.keys()):
                labels[i+1] = pp['label']
            else:
                labels[i+1] = 'Bestfit {}'.format(i+1)


        ax_legend.legend( handles, labels,
                          numpoints=1, loc=4, markerscale=0.5, bbox_to_anchor=[0.8,0.02],
                          handler_map={tuple: HandlerTuple(ndivide=None)},
                           )
    return fig

def plot_sfh(ax, result, prior_draws=None, show_bestfit=True, show_priors=True,
             posts_params={'color':'c', 'lw':1.5},
             label='',
             priors_params={'facecolor':'None', 'edgecolor':'goldenrod', 'hatch':'//', 'lw':0 },
             bestfit_params={'marker':'X', 'edgecolor':'k', 'facecolor':'None', 's':80, 'lw':0.9},
             quantiles=[0.16,0.84],
             norm_by_mass=False,
             **extras,
            ):
    from Dragonfly44_SFH.utils.misc_utils import weighted_quantile
    from Dragonfly44_SFH.utils.transforms import chain_to_sfr

    agebins = np.copy( result['agebins'] )
    agebins_Gyr = np.power(10, agebins-9)
    ncomp = np.shape( agebins )[0]

    x = np.unique( agebins_Gyr )
    mx = x[:-1]+np.diff(x)/2. # midpoint of bins

    sfrs_post = chain_to_sfr( norm_by_mass=norm_by_mass, **result )

    ibest = np.argmax( result["lnprobability"] )
    sfrs_best = sfrs_post[ibest,:]

    if show_bestfit:
        ax.scatter( mx, sfrs_best,
                    zorder=3, label=' '.join(["Bestfit",label]), **bestfit_params)

    quantiles = [0.5]+quantiles
    w = result['weights']
    qs50, qs_1, qs_2 = np.array([ weighted_quantile( x, quantiles, w )
                            for x in sfrs_post.T ]).T
    ax.step( x, np.append( qs50, qs50[-1] ),
             where='post', zorder=2, label=' '.join(["Posterior median",label]), **posts_params)
    ax.fill_between( x, np.append( qs_1, qs_1[-1] ), np.append( qs_2, qs_2[-1] ),
                     alpha=0.3, color=posts_params['color'], step='post', lw=0, zorder=1,
                     label=' '.join(["68% CRs of posterior",label]))

    if show_priors and (prior_draws is not None):
        if norm_by_mass: prior_key = 'ssfr'
        else: prior_key = 'sfr'
        qs50, qs_1, qs_2 = np.array([ np.quantile( x, q=quantiles )
                                for x in prior_draws[prior_key].T ]).T
        ax.step( x, np.append(qs50, qs50[-1]),
                 where='post', zorder=1, color=priors_params['edgecolor'], ls='--',
                 label=' '.join(['Prior median',label]) )
        ax.fill_between( x, np.append( qs_1, qs_1[-1] ), np.append( qs_2, qs_2[-1] ),
                         step='post', zorder=1,
                         label=' '.join(['68% CRs of prior',label]), **priors_params)

    ax.set(xlim=[-0.1,x[-1]+0.1], yscale='log')
    if norm_by_mass: ylabel = r'sSFR  (yr$^{-1}$)'
    else: ylabel = r'SFR  (M$_\odot$/yr)'
    ax.set_ylabel( ylabel )
    ax.set_xlabel(r'Lookback time  (Gyr)' )

    return ax

def plot_cmf(ax, result, prior_draws=None, show_bestfit=True, show_priors=True,
             style='violin',
             posts_params={'color':'c', 'lw':0},
             priors_params={'facecolor':'None', 'edgecolor':'goldenrod', 'hatch':'//', 'lw':0 },
             bestfit_params={'marker':'X', 'edgecolor':'k', 'facecolor':'None', 's':80, 'lw':0.9},
             quantiles=[0.16,0.84],
             xscale='log', vwidths=0.2,
             **extras,
            ):

    assert style in ['violin','step'], "Error: style must be one of: violin, step"

    agebins = np.copy( result['agebins'] )
    agebins_Gyr = np.power(10, agebins-9)
    xs = np.unique( agebins_Gyr )

    # get posteriors for the sfr
    from Dragonfly44_SFH.utils.transforms import chain_to_sfr
    sfrs_post = chain_to_sfr( norm_by_mass=True, **result )

    # covert to cmf
    from prospect.plotting.sfh import sfh_to_cmf
    x, cmfs_post = sfh_to_cmf( sfrs_post, agebins )
    if xscale=='log': x = np.log10(x)

    # best bestfit solution
    ibest = np.argmax( result["lnprobability"] )
    cmfs_best = cmfs_post[ibest,:]

    # reweight, if weights exist (i.e., used nested sampling)
    if 'weights' in result.keys():
        from dynesty.utils import resample_equal
        cmfs_post = resample_equal( cmfs_post, result['weights'])

    if style=='violin':

        if show_priors and (prior_draws is not None):
            parts = ax.violinplot( prior_draws['cmf'],
                                   positions=x,
                                   widths=vwidths*2,
                                   showmeans=False, showmedians=True, showextrema=False,
                                 )
            for pc in parts['bodies']:
                pc.set_facecolor( 'None'  )
                pc.set_edgecolor( priors_params['edgecolor']  )
                pc.set_linewidth( 1  )
                pc.set_alpha( 0.5 )
                pc.set_zorder( 1 )

            for p in ['cmedians']:
                if p not in parts.keys(): continue
                parts[p].set_color( priors_params['edgecolor'] )
                parts[p].set_zorder(1)
            parts['cmedians'].set_linewidth(1.5)

        parts = ax.violinplot( cmfs_post,
                               positions=x,
                               widths=vwidths,
                               showmeans=False, showmedians=True, showextrema=False,
                             )

        for pc in parts['bodies']:
            pc.set_facecolor( posts_params['color'] )
            pc.set_edgecolor( 'k'  )
            pc.set_linewidth( 0.5  )
            pc.set_alpha( 0.3 )
            pc.set_zorder( 2 )

        for p in ['cmedians']:
            if p not in parts.keys(): continue
            parts[p].set_color( posts_params['color'] )
            parts[p].set_zorder(3)
        parts['cmedians'].set_linewidth(2.5)


        q50 = np.median( cmfs_post, axis=0 )
        ax.plot( x, q50, color=posts_params['color'], alpha=0.2, zorder=0, lw=4 )

        if show_bestfit:
            ax.scatter( x, cmfs_best,
                        zorder=3, **bestfit_params)


    elif style=='step':

        if show_priors and (prior_draws is not None):
            qs50, qs_1, qs_2 = np.array([ np.quantile( x, q=[0.5]+quantiles )
                                    for x in prior_draws['cmf'].T ]).T
            ax.step( x, qs50, where='post', zorder=1, color=priors_params['edgecolor'], ls='--' )
            ax.fill_between( x, qs_1, qs_2,
                             step='post', zorder=1, **priors_params)


        qs50, qs_1, qs_2 = np.array([ np.quantile( x, q=[0.5]+quantiles )
                                for x in cmfs_post.T ]).T
        ax.step( x, qs50, where='post', zorder=1, color=posts_params['color'], ls='-' )
        ax.fill_between( x, qs_1, qs_2,
                         alpha=0.3, step='post', zorder=1, **posts_params)

        if show_bestfit:
            mx = np.median( agebins_Gyr , axis=1 )
            ax.scatter( mx, cmfs_best[:-1],
                        zorder=3, **bestfit_params)


    if xscale=='log':
        xticks = ([0.1,0.5,1,2,3,5,8,13])
        log_xticks = np.log10(xticks)
        ax.set_xticks( log_xticks )
        ax.set_xticklabels( xticks )

        ax.set_xlim(-2, np.log10(xs[-1]) )

    ax.set_ylabel( r'Cumulative  $M_\ast/\mathrm{M}_\odot$' )
    ax.set_xlabel(r'Lookback time  (Gyr)' )

    return ax

from corner import corner as triangle
def plot_corner( result, showpars=None,  return_fig=False, xparams = {}, fig=None, color='c', **extras ):
    from Dragonfly44_SFH.utils.misc_utils import flatten

    if showpars is None: showpars = ['logmass','logzsol','dust2','MWA','spec_norm']
    Npar = len(showpars)

    fchain, weights, brange = flatten(result, showpars=showpars, q=[0.05,0.95])


    xparams2 = {"plot_datapoints": True, "plot_density": True, "fill_contours": False, "bins":30,\
                'show_titles':False, 'quantiles':[0.5], 'color':color, 'smooth':1, 'range':brange }
    xparams2.update( xparams )
    if "labels" not in xparams2.keys(): xparams2['labels'] = showpars


    sel = weights > 0
    cornerfig = triangle(fchain[sel,:], weights=weights[sel], fig=fig, **xparams2)
    if return_fig: return cornerfig
    else: plt.show()

def chain_to_struct_kw( chain, model, theta_index, names ):
    """Given a (flat)chain (or parameter dictionary) and a model, convert the
    chain to a structured array

    :param chain:
        A chain, ndarry of shape (nsamples, ndim) or a dictionary of
        parameters, values of which are numpy datatypes.

    :param model:
        A ProspectorParams instance

    :returns struct:
        A structured ndarray of parameter values.

    KW: edited from prospect.io.write_resutls.chain_to_struct to accomodate
    chains which include derived parameters. Instead of using the model,
    use the theta_index key in the results dictionary
    """
    from Dragonfly44_SFH.utils.transforms import chain_to_param

    indict = type(chain) == dict
    if indict:
        return dict_to_struct(chain)
    else:
        n = np.prod(chain.shape[:-1])

    dt = []
    for p in names:
        if p in model.params.keys():
            s = model.params[p].shape
            x = (p, model.params[p].dtype, s )
        else:
            xx = chain_to_param( chain, theta_index, p )
            s = xx.shape[-1]
            x = (p, type(xx), tuple([s])  )
        dt.append( x )

    struct = np.zeros(n, dtype=np.dtype(dt))
    for i, p in enumerate(names):
        inds = theta_index[p]

        x = chain[..., inds]

        if p in model.params.keys():
            s = model.params[p].shape[0]
        else:
            s = x.shape[-1]

        struct[p] = x.reshape( -1, s )

    return struct

class FigureMaker_kw( FigureMaker ):
    """
    Edited version of prospect.plotting.figuremaker allowing for input result which chain
    includes derived parameters (so dimension of chain different from model dimension).
    Result dictionary must include 'theta_index' entry which specifies the indices of
    each entry of theta (including derived parameters)
    """

    def __init__(self, results_file="", result=None, show=None, nufnu=False, microns=True,
                 n_seds=-1, prior_samples=10000, **extras):

        self.results_file = results_file
        self.prior_samples = prior_samples
        self.n_seds = n_seds
        self.nufnu = nufnu
        self.microns = microns
        if results_file or result:
            self.read_in(results_file, result, **extras)
        self.spec_best = self.phot_best = None
        if show is not None:
            self.show = show

    def read_in(self, results_file=None, result=None, model=None ):
        """Read a prospector results file, cache important components,
        and do any parameter transformations.  The cached attributes are:

        * `obs` - The `obs` dictionary ised for the fit
        * `model` - The model used for the fit, if it could be reconstructed.
        * `chain` - Structured array of parameter vector samples
        * `weights` - Corresponding weights for each sample
        * `ind_best` - Index of the sample with the highest posterior probability
        * `parchain` - The chain transformed to the desired derived parameters via
                       the `convert` method.

        :param results_file: string
            full path of the file with the prospector results.
        """
        if results_file:
            self.result, self.obs, self.model = reader.results_from(results_file)
        elif result:
            self.result = result
            self.obs = result['obs']

        if model:
            self.model = model
        if self.model is None:
            self.model = reader.get_model(self.results)
        self.sps = None

        # if fit with emcee, need to reshape the chain
        # convert dimensions of the chain to (Nwalkers*Nsamples, Npar)
        chain = self.result["chain"]
        if chain.ndim>2:
            dims = chain.shape
            chain = chain.reshape( (dims[0]*dims[1], dims[2]) )
        self.result["chain"] = chain

        self.chain = chain_to_struct_kw( self.result["chain"], self.model,
                                         self.result["theta_index"], self.result['theta_index'].keys(),
                                       )

        self.weights = self.result.get("weights", None)
        self.ind_best = np.argmax(self.result["lnprobability"])
        self.parchain = self.convert(self.chain)

    # bug fixes
    def plot_corner(self, caxes, **extras):
        """Example to make a corner plot of the posterior PDFs for the
        parameters listed in `show`.

        :param caxes: ndarray of axes of shape (nshow, nshow)
        """
        xx = np.squeeze([ np.array(self.parchain[p]) for p in self.show ])
        labels = [pretty.get(p, p) for p in self.show]
        spans = get_spans(None, xx)
        caxes = allcorner(xx, labels, caxes, weights=self.weights, span=spans,
                          color=self.pkwargs["color"], hist_kwargs=self.hkwargs,
                          label_kwargs=self.label_kwargs,
                          tick_kwargs=self.tick_kwargs, max_n_ticks=4, **extras)
        # plot priors
        if self.prior_samples > 0:
            self.show_priors(np.diag(caxes), spans, smooth=0.05, **self.rkwargs)

    def make_axes(self, figsize=None):
        """Make a set of axes and assign them to the object.
        """
        self.caxes = plt.subplots(len(self.show), len(self.show), figsize=figsize)

    def make_seds(self, full=False):
        """Generate and cache the best fit model spectrum and photometry.
        Optionally generate the spectrum and photometry for a number of
        posterior samples.

        Populates the attributes `*_best` and `*_samples` where `*` is:
        * spec
        * phot
        * sed
        * cal

        :param full: bool, optional
            If true, generate the intrinsic spextrum (`sed_*`) over the entire wavelength
            range.  The (restframe) wavelength vector will be given by
            `self.sps.wavelengths`
        """
        if self.sps is None:
            self.build_sps()

        # --- best sample ---
        xbest = np.hstack([ self.result["chain"][self.ind_best, self.result["theta_index"][p]]
                            for p in self.model.free_params ]).reshape(-1)
        blob = self.model.predict(xbest, obs=self.obs, sps=self.sps)
        self.spec_best, self.phot_best, self.mfrac_best = blob
        self.cal_best = self.model._speccal.copy()
        if full:
            from copy import deepcopy
            dummy = deepcopy(self.obs)
            dummy["wavelength"] = None
            dummy["spectrum"] = None
            s, p, _ = self.model.predict(xbest, obs=dummy, sps=self.sps)
        else:
            dummy = None
        self.sed_best = self.model._sed.copy()

        # --- get SED samples ---
        if self.n_seds > 0:
            blob = self.draw_seds(self.n_seds, dummy=dummy)
            self.spec_samples, self.phot_samples, self.sed_samples, self.cal_samples = blob
        else:
            self.spec_samples = np.atleast_2d(self.spec_best)
            self.phot_samples = np.atleast_2d(self.phot_best)
            self.sed_samples = np.atleast_2d(self.sed_best)
            self.cal_samples = np.atleast_2d(self.cal_best)

from matplotlib import transforms as mpl_transforms

def multicolor_text(ax, fig, x0,y0, string_list, color_list, **kw):
    """
    Take a list of strings ``ls`` and colors ``lc`` and place them next to each
    other, with text ls[i] being shown in color lc[i].

    This example shows how to do both vertical and horizontal text, and will
    pass all keyword arguments to plt.text, so you can set the font size,
    family, etc.
    """
    t = ax.transData

    #horizontal version
    for s,c in zip(string_list, color_list ):
        text = ax.text(x0,y0, s, color=c, transform=t, **kw)
        text.draw(fig.canvas.get_renderer())
        ex = text.get_window_extent()
        t = mpl_transforms.offset_copy(text._transform, x=ex.width, units='dots')

def plot_columns( table, ax=None, col_x=None, col_y=None, params={}, ebparams={}, x_off=0., y_off=0., add_err=False, x_err=0., y_err=0., **extras  ):
    from Dragonfly44_SFH.utils.misc_utils import get_data_x3_from_data, add_error

    x, ex = get_data_x3_from_data( table, [col_x,'em_'+col_x,'ep_'+col_x])
    y, ey = get_data_x3_from_data( table, [col_y,'em_'+col_y,'ep_'+col_y])

    x = x+x_off
    y = y+y_off

    if add_err:
        ex[0,:] = add_error( ex[0,:], -x_err[0] )
        ex[1,:] = add_error( ex[1,:], x_err[1] )
        ey[0,:] = add_error( ey[0,:], -y_err[0] )
        ey[1,:] = add_error( ey[1,:], y_err[1] )

    cb = ax.scatter(  x, y, **params )
    ax.errorbar( x, y, xerr=ex, yerr=ey, **ebparams )

    return cb

def plot_cols_wcolor( table, axes=None, col_x=None, cols_y=None, col_c=None, params=None, ebparams=None, return_cb=False, params_colors=None, **extras  ):
    from Dragonfly44_SFH.utils.misc_utils import get_data_x3_from_data

    x, ex = get_data_x3_from_data( table, [col_x,'em_'+col_x,'ep_'+col_x])
    if col_c is not None:
        c, ec = get_data_x3_from_data( table, [col_c,'em_'+col_c,'ep_'+col_c])
        if np.all(np.isnan(c)): col_c = None

    for j,col_y in enumerate( cols_y ):
        y,ey = get_data_x3_from_data( table, [col_y,'em_'+col_y,'ep_'+col_y])
        try:
            params.update( extras )
            params.update( params_colors )

            if col_c is not None:
                cb = axes[j].scatter(  x, y, c=c, **params )
            else:
                axes[j].scatter(  x, y, **params )
                cb = None

            if ebparams is not None:
                axes[j].errorbar( x, y, xerr=ex, yerr=ey, **ebparams )
        except:
            pass

    if return_cb: return cb


def draw_connection( tab1, tab2, params_c=None, cols_x=None, cols_y=None, axes=None, xoff=0, yoff=0, **extras ):
    from Dragonfly44_SFH.utils.misc_utils import get_data_x3_from_data

    for i,col_x in enumerate( cols_x ):
        axes_i = axes[i]

        x1, _ = get_data_x3_from_data( tab1, [col_x,'em_'+col_x,'ep_'+col_x])
        x2, _ = get_data_x3_from_data( tab2, [col_x,'em_'+col_x,'ep_'+col_x])

        for j,col_y in enumerate( cols_y ):
            ax = axes_i[j]

            y1,_ = get_data_x3_from_data( tab1, [col_y,'em_'+col_y,'ep_'+col_y])
            y2,_ = get_data_x3_from_data( tab2, [col_y,'em_'+col_y,'ep_'+col_y])
            y1 += yoff
            y2 += yoff
            x1 += xoff
            x2 += xoff

            if type(x1)==np.ndarray: x1=float(x1)
            if type(x2)==np.ndarray: x2=float(x2)
            if type(y1)==np.ndarray: y1=float(y1)
            if type(y2)==np.ndarray: y2=float(y2)

            ax.plot(  [x1,x2], [y1,y2], lw=2, color='w', zorder=params_c['zorder']-0.001  )
            ax.plot(  [x1,x2], [y1,y2], **params_c )
