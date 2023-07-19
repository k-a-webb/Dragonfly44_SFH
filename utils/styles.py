
parameter_labels = dict( dust2=r"$\hat{\tau}_\mathrm{dust,~diffuse}$",
                         Av=r"$A_\mathrm{V}$",
                         logzsol=r"$\log(~Z_\ast/\mathrm{Z}_\odot~)$",
                         logmass=r"$\log(~M_{\ast,~\mathrm{total}}/\mathrm{M}_\odot~)$",
                         logmass_stellar=r"$\log(~M_{\ast}/\mathrm{M}_\odot~)$",
                         mwa=r"$t_\mathrm{mass}$",
                         t50=r'$t_\mathrm{50}$',
                         t90=r'$t_\mathrm{90}$',
                         t50mt90=r'$t_\mathrm{5s0}-t_\mathrm{90}$',
                         sfr=r'$\mathrm{SFR}$',
                         log_sfr=r'$\log \mathrm{SFR}$',
                         ssfr=r'$\mathrm{SFR}~/~M_\mathrm{formed}$',
                         log_ssfr=r'$\log \left( \mathrm{SFR}~/~M_\mathrm{formed} \right)$',
                        )

parameter_units = dict( dust2=None,
                        Av=None,
                         logzsol=None,
                         logmass=None,
                         logmass_stellar=None,
                         mwa="Gyr",
                         t50="Gyr",
                         t90="Gyr",
                         t50mt90="Gyr",
                         sfr=r"$\mathrm{M}_\odot~\mathrm{yr}^{-1}$",
                         log_sfr=r"$\mathrm{M}_\odot~\mathrm{yr}^{-1}$",
                         ssfr=r"$\mathrm{yr}^{-1}$",
                         log_ssfr=r"$\mathrm{yr}^{-1}$",
                        )

parameter_labels_with_units = {}
for k,v in parameter_labels.items():
    if parameter_units[k] is None:
        parameter_labels_with_units[k] = v
    else:
        parameter_labels_with_units[k] = '{}  ({})'.format(v, parameter_units[k])
parameter_labels_with_units['log_ssfr'] = r'$\log\left( \mathrm{SFR}~/~M_\mathrm{formed} ~~ (\mathrm{yr}^{-1})~\right) $'
parameter_labels_with_units['log_sfr'] = r'$\log\left( \mathrm{SFR} ~~ (\mathrm{M}_\odot~\mathrm{yr}^{-1})~\right) $'


parameter_bounds = dict( dust2=[0,2],
                         Av=[0,3],
                         logzsol=[-2,0.2],
                         logmass=[7,10],
                         logmass_stellar=[7,10],
                         mwa=[0,14],
                         t50=[0,14],
                         t90=[0,14],
                         t50mt90=[0,14],
                        )


sfh_labels = dict( aD02=r"Dirichlet($\alpha=0.2$)",
                       aD03=r"Dirichlet($\alpha=0.3$)",
                       aD05=r"Dirichlet($\alpha=0.5$)",
                       aD1=r"Dirichlet($\alpha=1$)",
                       csfrcont=r"Constant SFR Continuity",
                       dexp=r"Delayed-$\tau$",
                        )

# integer used by Prospector/FSPS to identify type of SFH model
sfh_sfhs = dict( aD02=3,
                       aD03=3,
                       aD05=3,
                       aD1=3,
                       csfrcont=3,
                       dexp=4,
                        )
