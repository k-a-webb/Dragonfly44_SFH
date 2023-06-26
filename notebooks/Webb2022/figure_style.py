
from matplotlib import rcParams
params = {'backend': 'ps', 'axes.labelsize': 12, 'axes.titlesize': 12, 'font.size': 12, \
          'legend.fontsize': 12, 'xtick.labelsize': 10, 'ytick.labelsize': 10, \
           'font.family': 'sans serif' }
rcParams.update(params)

saveparams = {'bbox_inches':'tight', 'pad_inches':0.05, 'dpi':250}

fig_width_one = 244./72.*2.
fig_width_two = 508./72.*2.
textheight = 682./72.*2.

import matplotlib.ticker as ticker
from matplotlib.legend_handler import HandlerTuple

import matplotlib as mpl
mpl.rcParams['hatch.linewidth'] = 0.5

color_aD1 = '#d71746ff'
color_aD02 = '#42d4f4ff'

marker_params_aD1 =dict( marker='D', s=160, edgecolor=color_aD1, facecolor='None' )
marker_params_aD02 =dict( marker='s', s=120, edgecolor=color_aD02, facecolor='None' )
