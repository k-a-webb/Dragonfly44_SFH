import numpy as np
import pandas as pd
from astropy.cosmology import FlatLambdaCDM

def convert_entry( x ):
    try:
        return np.float(x)
    except:
        return np.nan

def read_yagi2016_table4( file_path="", **extras ):
    # yagi+2016
    # https://iopscience.iop.org/article/10.3847/0067-0049/225/1/11#apjsaa2413t5
    # https://iopscience.iop.org/article/10.3847/1538-4365/aa6e4e

    #    1-  3  A3     ---           ID        Subaru UDG ID
    #  103-106  F4.1   arcsec        ReS       Effective radius from SExtractor
    #  108-111  F4.1   mag/arcsec2   SuBrmS    Mean surface brightness in ReS from model (2)
    #  139-142  F4.2   kpc           Re1       ? Effective radius from single Sersic fitting
    #  144-147  F4.2   kpc         e_Re1       ? Uncertainty in Re1
    #  197-200  F4.2   kpc           Re2       ? Effective radius of Sersic component
    #                                            from double component fitting
    #  202-205  F4.2   kpc         e_Re2       ? Uncertainty in Re2

    ID,ID2,ReS,SuBrmS,SuBrS,SuBr0,Re1,e_Re1,Re2,e_Re2 = [ np.array([]) for i in range(10)]
    with open( file_path+'yagi2016_table4.csv','r') as infile:
        for ii,lin in enumerate(infile.readlines()):
            if ii<70: continue
            ID     = np.append( ID,     int(lin[:4]))
            ID2    = np.append( ID2,    (lin[43-1:89]).strip(' ') )
            ReS    = np.append( ReS,    convert_entry(lin[103-1:106] ))
            SuBrmS = np.append( SuBrmS, convert_entry(lin[108-1:111] ))
            SuBrS = np.append( SuBrS, convert_entry(lin[108-1:111] ))
            SuBr0 = np.append( SuBr0, convert_entry(lin[176-1:179] ))
            Re1    = np.append( Re1,    convert_entry(lin[139-1:142] ))
            e_Re1  = np.append( e_Re1,  convert_entry(lin[144-1:147] ))
            Re2    = np.append( Re2,    convert_entry(lin[197-1:200] ))
            e_Re2  = np.append( e_Re2,  convert_entry(lin[202-1:205] ))

    df_yagi2016 = pd.DataFrame( dict( ID=ID, ID2=ID2, ReS=ReS, SuBrmS=SuBrmS, SuBrS=SuBrS, Re1=Re1, e_Re1=e_Re1, SuBr0=SuBr0, Re2=Re2, e_Re2=e_Re2 ))

    # Re2 = Effective radius of Sersic component from double component fitting
    # Re1 = Effective radius from single Sersic fitting (used by FM18)
    # ReS = Effective radius from SExtractor
    # SuBr0 = Central surface brightness of single Sersic model (2)
    # SuBrmS = Mean surface brightness in ReS from model (2) (Sextractor)
    # SuBrS = Surface brightness at ReS from model (2)
    df_yagi2016['r_eff'] = df_yagi2016['Re1'].values
    df_yagi2016['e_r_eff'] = df_yagi2016['e_Re1'].values

    df_yagi2016['mu'] = df_yagi2016['SuBrmS'].values
    df_yagi2016['e_mu'] = np.full( len(df_yagi2016), np.nan)

    df_yagi2016['Yagi'] = ['Yagi{:.0f}'.format(idn) for idn in df_yagi2016['ID'].values ]
    df_yagi2016['IBG'] = np.full( len(df_yagi2016), np.nan )
    df_yagi2016['DF'] = np.full( len(df_yagi2016), np.nan )
    df_yagi2016['GMP'] = np.full( len(df_yagi2016), np.nan )
    for idx in df_yagi2016.index.values:
        idn2 = df_yagi2016.loc[idx,'ID2']

        if ("DF" in idn2) or ('GMP' in idn2) or ('IBG' in idn2):
            idn2 = idn2.split()
            for idn3 in idn2:
                if ('DF' in idn3): df_yagi2016.loc[idx,'DF'] = idn3
                elif ('GMP' in idn3): df_yagi2016.loc[idx,'GMP'] = idn3
                elif ('IBG' in idn3): df_yagi2016.loc[idx,'IBG'] = idn3.replace('IBG-','')

    return df_yagi2016

def read_alibi2020( file_path="", df_yagi2016=None, **extras ):
    # alabi + 2020
    # https://academic.oup.com/mnras/article/496/3/3182/5859958#supplementary-data

    #V:  total AB V-band magnitude (K-corrected and corrected for Galactic extinction)
    #Verr: uncertainty on V-band magnitude (mag)
    #Re_V:  V-band circularized effective radius (kpc)
    #ReV_err: uncertainty on V-band circularized effective radius (kpc)
    #SB_V:  V-band mean surface brightness within the effective radius (mag/arcsec^2)
    #SB_V_err: uncertainty on V-band mean surface brightness (mag/arcsec^2)

    #R:  total AB R-band magnitude (K-corrected and corrected for Galactic extinction)
    #Rerr: uncertainty on R-band magnitude (mag)
    #Re_R:  R-band circularized effective radius (kpc)
    #ReR_err: uncertainty on R-band circularized effective radius (kpc)
    #SB_R:  R-band mean surface brightness within the effective radius (mag/arcsec^2)
    #SB_R_err: uncertainty on R-band mean surface brightness (mag/arcsec^2)

    if df_yagi2016 is None:
        print('Warning: Need "df_yagi2016" to add extra information. ')
        if make_input_if_not_given:
            print('\tReading df_yagi2016...')
            df_yagi2016 = read_yagi2016_table4( file_path )

    headers = ['ID','RA','Dec','V','Verr','re_V','ReV_err','SB_V','SB_V_err','n_V','nV_err','q_V','qV_err','PA_V','R','Rerr','re_R','ReR_err','SB_R','SB_R_err','n_R','nR_err','q_R','qR_err','PA_R','V-R','comment']
    df_alabi2020 = pd.read_csv(file_path+'alabi2020.csv', delim_whitespace=True, skiprows=32, names=headers)
    df_alabi2020 = df_alabi2020.loc[:,['ID','re_V','ReV_err','SB_V','SB_V_err','re_R','ReR_err','SB_R','SB_R_err','comment']]

    filt = "R" # "V"
    df_alabi2020['r_eff'] = df_alabi2020['re_{}'.format(filt)].values
    df_alabi2020['e_r_eff'] = df_alabi2020['Re{}_err'.format(filt)].values

    df_alabi2020['mu'] = df_alabi2020['SB_{}'.format(filt)].values
    df_alabi2020['e_mu'] = df_alabi2020['SB_{}_err'.format(filt)].values

    for col in ['Yagi','DF','GMP','IBG']:
        df_alabi2020[col] = np.full( len(df_alabi2020), np.nan )

    for idx in df_alabi2020.index.values:
        idn2 = df_alabi2020.loc[idx,'comment']

        if ("DF" in idn2) or ('GMP' in idn2) or ('Y' in idn2):
            idn2 = idn2.split(',')
            for idn3 in idn2:
                if ('DF' in idn3): df_alabi2020.loc[idx,'DF'] = idn3
                elif ('GMP' in idn3): df_alabi2020.loc[idx,'GMP'] = idn3
                elif idn3.startswith('Y'):
                    idn_yagi = 'Yagi{}'.format( idn3[1:]  )
                    df_alabi2020.loc[idx,'Yagi'] = idn_yagi
                    if df_yagi2016 is not None:
                        df_alabi2020.loc[idx,'DF']  = df_yagi2016.query("Yagi=='{}'".format( idn_yagi ))['DF'].values[0]
                        df_alabi2020.loc[idx,'IBG'] = df_yagi2016.query("Yagi=='{}'".format( idn_yagi ))['IBG'].values[0]



    return df_alabi2020

def read_alabi2018_table1( file_path="", **extras ):
    df_alabi2018 = pd.read_csv( file_path+'alabi2018_table1.csv', skiprows=5, delim_whitespace=True)
    return df_alabi2018

def read_ruizlara2018_table2( file_path="", df_yagi2016=None, df_alabi2020=None, make_input_if_not_given=False, **extras ):
    # ruiz-lara+ 2018

    if df_yagi2016 is None:
        print('Warning: Need "df_yagi2016" to add extra information. ')
        if make_input_if_not_given:
            print('\tReading df_yagi2016...')
            df_yagi2016 = read_yagi2016_table4( file_path )

    if df_alabi2020 is None:
        print('Warning: Need "df_alabi2020" to add extra information')
        if make_input_if_not_given:
            print('\tReading df_alabi2020...')
            df_alabi2020 = read_alibi2020( file_path, df_yagi2016=df_yagi2016 )


    cosmo_ruizlara2018 = FlatLambdaCDM(H0=69.7, Om0=0.281)
    df_ruizlara2018 = pd.read_csv( file_path+'ruizlara2018_table2.csv')

    if df_yagi2016 is not None:
        #  add other names
        for idx in df_ruizlara2018.index.values:
            idn = df_ruizlara2018.loc[idx,'Galaxy']
            if "Yagi" in idn: idn = idn.replace('Yagi 0','Yagi').replace('Yagi ','Yagi')

            cols = ['DF','Yagi','GMP']
            for col1 in cols:
                if col1 in idn:
                    df_ruizlara2018.loc[idx,col1] = idn

                    if idn in df_yagi2016[col1].values:
                        idx2 = df_yagi2016.query('{}=="{}"'.format(col1,idn)).index.values[0]

                        for col2 in cols:
                            if col1==col2: continue
                            df_ruizlara2018.loc[idx,col2] = df_yagi2016.loc[idx2,col2]

    x,ex_m,ex_p = df_ruizlara2018[['M*/1E8 (Msun)', 'M*_m', 'M*_p']].values.T
    (x,ex_m,ex_p) = np.array([x,ex_m,ex_p])*1e8
    x2 = np.log10(x)
    ex2_m = ex_m / x / np.log(10)
    ex2_p = ex_m / x / np.log(10)

    df_ruizlara2018['logmass'] = x2
    df_ruizlara2018['em_logmass'] = ex2_m
    df_ruizlara2018['ep_logmass'] = ex2_p

    x,ex_m,ex_p = df_ruizlara2018[['log(AgeLW[yr]) (dex)', 'lAge_m', 'lAge_p']].values.T
    x2 = 10**x
    ex2_m = x2 * ex_m * np.log(10)
    ex2_p = x2 * ex_p * np.log(10)
    (x2,ex2_m,ex2_p) = np.array([x2,ex2_m,ex2_p]) * 1e-9 # yr to Gyr
    df_ruizlara2018['age_LW'] = x2
    df_ruizlara2018['em_age_LW'] = ex2_m
    df_ruizlara2018['ep_age_LW'] = ex2_p

    df_ruizlara2018['lbt50'] = df_ruizlara2018['t50 (Gyr)']
    df_ruizlara2018['lbt90'] = df_ruizlara2018['t90 (Gyr)']

    df_ruizlara2018['tH'] = cosmo_ruizlara2018.age(df_ruizlara2018['z']).value

    df_ruizlara2018['age_MW'] = df_ruizlara2018['t50 (Gyr)']

    df_ruizlara2018['logzsol'] = df_ruizlara2018['[M/H]LW (dex)']
    df_ruizlara2018['em_logzsol'] = df_ruizlara2018['[M/H]_m'].values
    df_ruizlara2018['ep_logzsol'] = df_ruizlara2018['[M/H]_p'].values

    df_ruizlara2018['t90mt50'] = df_ruizlara2018['t50 (Gyr)'] - df_ruizlara2018['t90 (Gyr)']

    df_ruizlara2018['t50'] = df_ruizlara2018['tH'] - df_ruizlara2018['lbt50']
    df_ruizlara2018['t90'] = df_ruizlara2018['tH'] - df_ruizlara2018['lbt90']


    if True: # surface brightness

        df_ruizlara2018.loc[:,'mu_self']       = df_ruizlara2018.loc[:,'muB (mag arcsec-2)']
        df_ruizlara2018.loc[:,'r_eff_self']    = df_ruizlara2018.loc[:,'Reff (kpc)']


        for i in df_ruizlara2018.index.values:
            idn = df_ruizlara2018.loc[i,'Galaxy']

            if idn.startswith('OGS'):
                df_ruizlara2018.loc[i,'mu_self_corr'] = df_ruizlara2018.loc[i,'muB (mag arcsec-2)'] - 0.9
                continue

            if df_yagi2016 is not None: # yagi

                found = None

                if idn.startswith('Yagi') and (df_yagi2016 is not None):
                    idx = int( idn.replace('Yagi ','') )
                    found = df_yagi2016.query( 'ID == {}'.format( idx ))

                elif idn.startswith('DF') and (df_yagi2016 is not None):
                    for idx in df_yagi2016.index.values:
                        if idn in df_yagi2016.loc[idx,'ID2']:
                            found = df_yagi2016.loc[[idx],:]
                            break

                if found is not None:
                    df_ruizlara2018.loc[i,'mu_y16']       = found['mu'].values[0]
                    df_ruizlara2018.loc[i,'em_mu_y16']    = found['e_mu'].values[0]
                    df_ruizlara2018.loc[i,'ep_mu_y16']    = found['e_mu'].values[0]

                    df_ruizlara2018.loc[i,'r_eff_y16']    = found['r_eff'].values[0]
                    df_ruizlara2018.loc[i,'em_r_eff_y16'] = found['e_r_eff'].values[0]
                    df_ruizlara2018.loc[i,'ep_r_eff_y16'] = found['e_r_eff'].values[0]
                else:
                    print( idn, 'not found in Yagi+2016')

            if df_alabi2020 is not None:

                found = None
                if idn.startswith('Yagi'):
                    idn2 = idn.replace('Yagi ','Y').replace('Y0','Y')

                elif idn.startswith('DF') and (df_yagi2016 is not None):
                    for idx in df_yagi2016.index.values:
                        if idn in df_yagi2016.loc[idx,'ID2']:
                            idn2 = "Y{:.0f}".format( df_yagi2016.loc[idx,'ID'] )
                            break

                if (df_alabi2020 is not None):
                    for idx in df_alabi2020.index.values:

                        if idn2 in df_alabi2020.loc[idx,'comment']:
                            found = df_alabi2020.loc[idx]
                            break

                if found is not None:
                    df_ruizlara2018.loc[i,'mu_a20']       = found['mu']
                    df_ruizlara2018.loc[i,'em_mu_a20']    = found['e_mu']
                    df_ruizlara2018.loc[i,'ep_mu_a20']    = found['e_mu']

                    df_ruizlara2018.loc[i,'r_eff_a20']    = found['r_eff']
                    df_ruizlara2018.loc[i,'em_r_eff_a20'] = found['e_r_eff']
                    df_ruizlara2018.loc[i,'ep_r_eff_a20'] = found['e_r_eff']
                else:
                    print( idn, 'not found in Alabi+2020')


#     df_ruizlara2018[['Galaxy','r_eff_self','r_eff_y16','r_eff_a20','mu_a20','age_MW','age_LW','FeH','logzsol','t90mt50',  'DF','Yagi','GMP']]
    return df_ruizlara2018


def read_ferremateu2018_table1_and_2( file_path="", df_yagi2016=None, df_alabi2020=None, df_alabi2018=None,
                                make_input_if_not_given=False, **extras ):
    # ferre-mateu + 2018
    cosmo_ferremateu2018 = FlatLambdaCDM(H0=70, Om0=0.27)
    tH = cosmo_ferremateu2018.age(0.023156).value

    if df_yagi2016 is None:
        print('Warning: Need "df_yagi2016" to add extra information. ')
        if make_input_if_not_given:
            print('\tReading df_yagi2016...')
            df_yagi2016 = read_yagi2016_table4( file_path )

    if df_alabi2020 is None:
        print('Warning: Need "df_alabi2020" to add extra information')
        if make_input_if_not_given:
            print('\tReading df_alabi2020...')
            df_alabi2020 = read_alibi2020( file_path, df_yagi2016=df_yagi2016 )

    if df_alabi2018 is None:
        print('Warning: Need "df_alabi2018" to add extra information')
        if make_input_if_not_given:
            print('\tReading df_alabi2018...')
            df_alabi2018 = read_alabi2018_table1( file_path )

    df_ferremateu2018 = pd.read_csv( file_path+'ferremateu2018_table2.csv')

    # for col in all_cols:
    #     if col not in df_ferremateu2018.columns:
    #         df_ferremateu2018[col] = np.full( len(df_ferremateu2018), np.nan )

    if df_yagi2016 is not None: # add other names
        for idx in df_ferremateu2018.index.values:
            idn = df_ferremateu2018.loc[idx,'Galaxy']
            if "Yagi" in idn: idn = idn.replace('Yagi0','Yagi')

            cols = ['DF','Yagi','GMP',"IBG"]
            for col1 in cols:
                if col1 in idn:
                    df_ferremateu2018.loc[idx,col1] = idn

                    if (idn in df_yagi2016[col1].values) and (df_yagi2016 is not None):
                        idx2 = df_yagi2016.query('{}=="{}"'.format(col1,idn)).index.values[0]

                        for col2 in cols:
                            if col1==col2: continue
                            df_ferremateu2018.loc[idx,col2] = df_yagi2016.loc[idx2,col2]

                elif (idn.startswith("J")) and (col1=="IBG"):

                    idn = idn.replace("J","")
                    df_ferremateu2018.loc[idx,col1] = idn

                    if idn in df_yagi2016[col1].values:
                        idx2 = df_yagi2016.query('{}=="{}"'.format(col1,idn)).index.values[0]

                        for col2 in cols:
                            if col1==col2: continue
                            df_ferremateu2018.loc[idx,col2] = df_yagi2016.loc[idx2,col2]

    df1_ferremateu2018 = pd.read_csv( file_path+'ferremateu2018_table1.csv')
    df_ferremateu2018['r_eff_self'] = df1_ferremateu2018['Re (kpc)'].values

    df_ferremateu2018['logmass'] = np.log10( df_ferremateu2018['M* (M⊙)'])

    df_ferremateu2018['age_MW'] = df_ferremateu2018['t50 (Gyr)']
    df_ferremateu2018['tH'] = np.full( len(df_ferremateu2018), tH )

    df_ferremateu2018['age_LW'] = df_ferremateu2018['Age_LW (Gyr)']
    df_ferremateu2018['em_age_LW'] = df_ferremateu2018['unc(Age)'].values
    df_ferremateu2018['ep_age_LW'] = df_ferremateu2018['unc(Age)'].values

    df_ferremateu2018['FeH'] = df_ferremateu2018['[Fe/H] (dex)']
    df_ferremateu2018['em_FeH'] = df_ferremateu2018['unc([Fe/H])'].values
    df_ferremateu2018['ep_FeH'] = df_ferremateu2018['unc([Fe/H])'].values

    df_ferremateu2018['logzsol'] = df_ferremateu2018['[Z/H] (dex)']
    df_ferremateu2018['em_logzsol'] = df_ferremateu2018['unc([Z/H])'].values
    df_ferremateu2018['ep_logzsol'] = df_ferremateu2018['unc([Z/H])'].values

    df_ferremateu2018['t90mt50'] = df_ferremateu2018['Δt90 (Gyr)']

    df_ferremateu2018['lbt50'] = df_ferremateu2018['t50 (Gyr)']
    df_ferremateu2018['lbt90'] = df_ferremateu2018['t90 (Gyr)']
    df_ferremateu2018['t50'] = df_ferremateu2018['tH'] - df_ferremateu2018['lbt50']
    df_ferremateu2018['t90'] = df_ferremateu2018['tH'] - df_ferremateu2018['lbt90']

    if True: # surface brightness

        # self == yagi
        for i in df_ferremateu2018.index.values:
            idn = df_ferremateu2018.loc[i,'Galaxy']

            if df_yagi2016 is not None: # Yagi + self
                if idn.startswith('Yagi'):
                    idx = int( idn.replace('Yagi','') )
                    found = df_yagi2016.query( 'ID == {}'.format( idx ))

                    for suff in ['_y16']:
                        df_ferremateu2018.loc[i,'mu'+suff]       = found['mu'].values[0]
                        df_ferremateu2018.loc[i,'em_mu'+suff]    = found['e_mu'].values[0]
                        df_ferremateu2018.loc[i,'ep_mu'+suff]    = found['e_mu'].values[0]

                        df_ferremateu2018.loc[i,'r_eff'+suff]    = found['r_eff'].values[0]
                        df_ferremateu2018.loc[i,'em_r_eff'+suff] = found['e_r_eff'].values[0]
                        df_ferremateu2018.loc[i,'ep_r_eff'+suff] = found['e_r_eff'].values[0]

            if (df_alabi2020 is not None) and (df_alabi2018 is not None): # Alabi+2020 then 2018

                found = None
                for col in ['DF','Yagi','GMP','IBG']:
                    idn = df_ferremateu2018.loc[i,col]
                    if type( idn ) == str:
                        if idn in df_alabi2020[col].values:
                            found = df_alabi2020.query('{}=="{}"'.format(col, idn))
                        break

                if found is not None:
                    df_ferremateu2018.loc[i,'mu_a20']       = found['mu'].values[0]
                    df_ferremateu2018.loc[i,'em_mu_a20']    = found['e_mu'].values[0]
                    df_ferremateu2018.loc[i,'ep_mu_a20']    = found['e_mu'].values[0]

                    df_ferremateu2018.loc[i,'r_eff_a20']    = found['r_eff'].values[0]
                    df_ferremateu2018.loc[i,'em_r_eff_a20'] = found['e_r_eff'].values[0]
                    df_ferremateu2018.loc[i,'ep_r_eff_a20'] = found['e_r_eff'].values[0]
                else:
                    print(idn,"Not found in Alabi+2020")

                    if "+" in idn:
                        idn1, idn2 = idn.split('+')

                        for i2 in df_alabi2018.index.values:
                            if (idn1 in df_alabi2018.loc[i2,'Galaxy']) & (idn2 in df_alabi2018.loc[i2,'Galaxy']):

                                found = df_alabi2018.loc[i2]
                                break
                    else:
                        i2 = df_alabi2018.query('Galaxy=="{}"'.format(idn)).index.values[0]
                        found = df_alabi2018.loc[i2]

                    if found is not None:
                        try:
                            df_ferremateu2018.loc[i,'mu_a20']       = found['mu0'] - 1.18 # to be consistent with Alabi+2020
                            df_ferremateu2018.loc[i,'em_mu_a20']    = 0.47
                            df_ferremateu2018.loc[i,'ep_mu_a20']    = 0.47
                        except:
                            pass
                        df_ferremateu2018.loc[i,'r_eff_a20']    = found['Re'] - 0.08
                        df_ferremateu2018.loc[i,'em_r_eff_a20']    = 0.47
                        df_ferremateu2018.loc[i,'ep_r_eff_a20']    = 0.47
                    else:
                        print( idn, 'not found in Alabi+2018')

#     df_ferremateu2018[['Galaxy','r_eff_self','r_eff_y16','r_eff_a20','mu_a20','age_MW','age_LW','FeH','logzsol','t90mt50',  'DF','Yagi','GMP']]
    return df_ferremateu2018

def read_gu2018_table1( file_path="", df_yagi2016=None, df_alabi2020=None, df_alabi2018=None,
                        make_input_if_not_given=False, **extras ):
    if df_yagi2016 is None:
        print('Warning: Need "df_yagi2016" to add extra information. ')
        if make_input_if_not_given:
            print('Reading df_yagi2016...')
            df_yagi2016 = read_yagi2016_table4( file_path )

    if df_alabi2020 is None:
        print('Warning: Need "df_alabi2020" to add extra information')
        if make_input_if_not_given:
            print('Reading df_alabi2020...')
            df_alabi2020 = read_alibi2020( file_path, df_yagi2016=df_yagi2016 )

    if df_alabi2018 is None:
        print('Warning: Need "df_alabi2018" to add extra information')
        if make_input_if_not_given:
            print('Reading df_alabi2018...')
            df_alabi2018 = read_alabi2018_table1( file_path )

    # gu + 2018

    # r_eff from vandokkum+2015a
    df_gu2018 = pd.read_csv( file_path+'gu2018_table1.csv')
    df_gu2018.rename( columns={'Target':'Galaxy', "[Fe/H]":"FeH", "em_[Fe/H]":"em_FeH", "ep_[Fe/H]":"ep_FeH", \
                                "r_eff":"r_eff_self", "em_r_eff":"em_r_eff_self", "ep_[r_eff":"ep_[r_eff_self"}, inplace=True)

#     for col in all_cols:
#         if col not in df_gu2018.columns:
#             df_gu2018[col] = np.full( len(df_gu2018), np.nan )

    # add other names
    if df_yagi2016 is not None:
        for idx in df_gu2018.index.values:
            idn = df_gu2018.loc[idx,'Galaxy']

            cols = ['DF','Yagi','GMP',"IBG"]
            for col1 in cols:
                if col1 in idn:
                    df_gu2018.loc[idx,col1] = idn

                    if df_yagi2016 is not None:
                        if idn in df_yagi2016[col1].values:
                            idx2 = df_yagi2016.query('{}=="{}"'.format(col1,idn)).index.values[0]

                            for col2 in cols:
                                if col1==col2: continue
                                df_gu2018.loc[idx,col2] = df_yagi2016.loc[idx2,col2]

    x,ex_m,ex_p = df_gu2018[['log(age/Gyr)','em_log(age/Gyr)','ep_log(age/Gyr)'] ].values.T
    x2 = 10**x
    ex2_m = x2 * ex_m * np.log(10)
    ex2_p = x2 * ex_p * np.log(10)
    (x2,ex2_m,ex2_p) = np.array([x2,ex2_m,ex2_p])
    df_gu2018['age_LW'] = x2
    df_gu2018['em_age_LW'] = ex2_m
    df_gu2018['ep_age_LW'] = ex2_p

    if True: # surface brightnesses etc
        for idx in df_gu2018.index.values:
            idn = df_gu2018.loc[idx,'Galaxy']
            if idn=="DF7": idn='DF07'

            if True and (df_yagi2016 is not None): # yagi

                try:
                    idx_yagi2016 = df_yagi2016.query('ID2=="{}"'.format(idn)).index.values[0]
                except IndexError:
                    print(idn, 'not found to reference in yagi 2016')
                    continue

                idx_yagi2016 = df_yagi2016.query('ID2=="{}"'.format(idn)).index.values[0]
                df_gu2018.loc[idx,'r_eff_y16']    = df_yagi2016.loc[idx_yagi2016,'r_eff']
                df_gu2018.loc[idx,'em_r_eff_y16'] = df_yagi2016.loc[idx_yagi2016,'e_r_eff']
                df_gu2018.loc[idx,'ep_r_eff_y16'] = df_yagi2016.loc[idx_yagi2016,'e_r_eff']

                df_gu2018.loc[idx,'mu_y16']    = df_yagi2016.loc[idx_yagi2016,'mu']
                df_gu2018.loc[idx,'em_mu_y16'] = df_yagi2016.loc[idx_yagi2016,'e_mu']
                df_gu2018.loc[idx,'ep_mu_y16'] = df_yagi2016.loc[idx_yagi2016,'e_mu']

            if True and (df_yagi2016 is not None) and (df_alabi2020 is not None): # alabi

                idn2 = 'Y{:.0f}'.format( df_yagi2016.loc[idx_yagi2016,'ID'] )

                found = None
                for idx2 in df_alabi2020.index.values:
                    if idn2 in df_alabi2020.loc[idx2,'comment']:
                        found = df_alabi2020.loc[idx2]
                        break

                if found is not None:
                    df_gu2018.loc[idx,'r_eff_a20']    = found['r_eff']
                    df_gu2018.loc[idx,'em_r_eff_a20'] = found['e_r_eff']
                    df_gu2018.loc[idx,'ep_r_eff_a20'] = found['e_r_eff']

                    df_gu2018.loc[idx,'mu_a20']    = found['mu']
                    df_gu2018.loc[idx,'em_mu_a20'] = found['e_mu']
                    df_gu2018.loc[idx,'ep_mu_a20'] = found['e_mu']
                else:
                    print(idn, 'not found in alabi 2020')
    return df_gu2018

def read_gu2018_extras( file_path="", df_gu2018=None, quantiles=[0.16,0.5,0.84], **extras ):

    if df_gu2018 is None: print("ERROr: I was lazy, and did not write this to not provide 'df_gu2018' as input")
    logzsol = np.loadtxt( file_path+'df44_zh_posterior_gu2018_for_kristi_01_24_2022.dat')
    qs = np.quantile( logzsol, quantiles )
    dqs = np.diff(qs)
    idx = df_gu2018.query('(Galaxy=="DF44") & (constraints=="spectra")').index.values[0]
    df_gu2018.loc[idx,'logzsol'] = qs[1]
    df_gu2018.loc[idx,'em_logzsol'] = dqs[0]
    df_gu2018.loc[idx,'ep_logzsol'] = dqs[1]


    logzsol = np.loadtxt( file_path+'UDG1_MPL11.csv', unpack=True, skiprows=1, delimiter=',')[2]
    qs = np.quantile( logzsol, quantiles )
    dqs = np.diff(qs)
    idx = df_gu2018.query('(Galaxy=="DF44") & (constraints=="both")').index.values[0]
    df_gu2018.loc[idx,'logzsol'] = qs[1]
    df_gu2018.loc[idx,'em_logzsol'] = dqs[0]
    df_gu2018.loc[idx,'ep_logzsol'] = dqs[1]

#     df_gu2018[['Galaxy','constraints','r_eff_self','r_eff_y16','r_eff_a20','mu_a20','age_MW','age_LW','FeH','logzsol','t90mt50',  'DF','Yagi','GMP']]
    return df_gu2018


def read_martinnavarro2019(  ):
    # Martin-Navarro + 2019 DGSAT1

    # Martin-Navarro+2019 only quotes t_50 for DGSAT 1 (5 Gyr) but once we know what kind of metrics we want to compare I’ll email him to get the fits/values

    # STECKMAP
    # The best-fitting model and residuals are shown in Fig. 2, where the most prominent absorption features are clearly noticeable. DGSAT I shows an extended SFH, with a mass-weighted age of 8.1 ± 0.4 Gyr, and luminosity-weighted age of 2.8 ± 0.5 Gyr. It took ∼5 Gyr for DGSAT I to form 50 per cent of its stellar mass. This result is in agreement with the photometric study of Pandya et al. (2018), where they estimated a luminosity-weighted age of ∼3 Gyr.

    # pPXF
    # allows for multiple metallicities per age bin. Thus, a metallicity spread at a given age should be captured by pPXF. With this addi- tional method, we recover a luminosity-weighted age of 3.5 Gyr, and a mass-weighted value of 7.9 Gyr, in excellent agreement with the STECKMAP analysis.

    # STARLIGHT

    #  α-enhanced ([Mg/Fe] = + 0.4), metal-poor ([M/H] = −1.7)
    #  [M/H] = −1.8 ± 0.4 dex assuming a mass-weighted age of 8.1 Gyr, as expected for a low-σ galaxy. However, metallicity measurements based on magnesium spectral features ([M/H]Mg = −1.0 ± 0.5) are higher than th
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    dict_martinnavarro2019 = dict( Galaxy="DGSAT I", z=0.018179,
                                   r_eff_self_I=4.7, em_r_eff_self_I=0.2, ep_r_eff_self_I=0.2, # Martínez-Delgado 2016 https://iopscience.iop.org/article/10.3847/0004-6256/151/4/96
                                   mu_self_V=24.8,  em_mu_self_V=0.2, ep_mu_self_V=0.2,
                                   mu_self_I=24.0,  em_mu_self_I=0.2, ep_mu_self_I=0.2,

                                   age_MW=8.1, em_age_MW=0.4, ep_age_MW=0.4,
                                   age_LW=2.8, em_age_LW=0.5, ep_age_LW=0.5,
                                   lbt95=0.4, lbt75=3.2, lbt50=8.1, e_lbt50=0.4, # t50 confusion between fig3 and text t50=5.0,

                                   age_MW_pPXF=7.9, age_LW_pPXF=3.5,
                                  )

#     for col in all_cols:
#         if col not in dict_martinnavarro2019.keys():
#             dict_martinnavarro2019[col] = np.nan

    dict_martinnavarro2019['FeH'] = -1.8 # from text [M/H]
    dict_martinnavarro2019['em_FeH'] = 0.4 # from text [M/H]
    dict_martinnavarro2019['ep_FeH'] = 0.4 # from text [M/H]

    # dict_martinnavarro2019['[Fe/H]'] = -2.7 # from fig 7
    # dict_martinnavarro2019['em_[Fe/H]'] = 0.3
    # dict_martinnavarro2019['ep_[Fe/H]'] = 0.3

    dict_martinnavarro2019['logmass'] = np.log10(4.8e8)

    dict_martinnavarro2019['tH'] = cosmo.age(dict_martinnavarro2019['z']).value
    dict_martinnavarro2019['lbt90'] = 0.9 # from fig 3

    dict_martinnavarro2019['t50'] = dict_martinnavarro2019['tH'] - dict_martinnavarro2019['lbt50']
    dict_martinnavarro2019['t90'] = dict_martinnavarro2019['tH'] - dict_martinnavarro2019['lbt90']
    dict_martinnavarro2019['t90mt50'] = dict_martinnavarro2019['t90'] - dict_martinnavarro2019['t50']

    # assume uncertainties are similar ....
    for col in ['em_lbt50','ep_lbt50','em_t50','ep_t50', 'em_lbt90','ep_lbt90','em_t90','ep_t90', 't90mt50']:
        dict_martinnavarro2019[col] = dict_martinnavarro2019['e_lbt50']

    df_martinnavarro2019 = pd.Series(dict_martinnavarro2019).to_frame().T
#     df_martinnavarro2019[['Galaxy', 'r_eff_self_I','mu_self_I','mu_self_V','age_MW','age_LW','FeH','logzsol','t90mt50', 'DF','Yagi','GMP']]
    return df_martinnavarro2019

def read_chan2018(  file_path="" ):
    cosmo = FlatLambdaCDM(H0=70, Om0=0.3)

    # chan + 2018 (FIRE simulation)
    df_chan2018 = pd.read_csv( file_path+'chan2018_table2_inpart.csv')
    df_chan2018.rename( columns={'Observed':'Galaxy', ' [Fe/H]':'FeH', ' age_MW (Gyr)':'age_MW', ' r_eff':'r_eff_self', ' mu_c':'mu_self'}, inplace=True)

    df_chan2018['logmass'] = np.log10( df_chan2018[' M (x1e8 Msun)']*1e8 )

#     for col in all_cols:
#         if col not in df_chan2018.columns:
#             df_chan2018[col] = np.full( len(df_chan2018), np.nan )

    df_chan2018['t90'] =  df_chan2018[' tq (Gyr)']
    df_chan2018['lbt90'] = cosmo.age(0).value - df_chan2018[' tq (Gyr)']
    df_chan2018['tH'] = np.full( len(df_chan2018), cosmo.age(0).value )
    df_chan2018['lbt50'] = df_chan2018['age_MW']
    df_chan2018['t50'] =  cosmo.age(0).value - df_chan2018['lbt50']
    df_chan2018['t90mt50'] = df_chan2018['t50'] - df_chan2018['t90']

    # df_chan2018

#     df_chan2018[['Galaxy','r_eff_self','mu_self','age_MW','age_LW','FeH','logzsol','t90','t90mt50']]
    return df_chan2018

def read_villaume2022( path_data, quantiles=[0.16,0.5,0.84], **extras ):
    #  results from Villaume+2022
    from copy import deepcopy

    dict_v22_1compSFH = dict( Galaxy='DF44', z = 0.02142 )
    dict_v22_1compSFH['source'] = '1_comp'

    # initialize with null entries
#     for col in all_cols:
#         if col not in dict_v22_1compSFH.keys():
#             dict_v22_1compSFH[col] = np.nan

    dict_v22_2compSFH = deepcopy(dict_v22_1compSFH)

    dict_v22_1compSFH['FeH'] = -1.33
    dict_v22_1compSFH['em_FeH'] = 0.04
    dict_v22_1compSFH['ep_FeH'] = 0.05

    dict_v22_1compSFH['age_LW'] = 10.2
    dict_v22_1compSFH['em_age_LW'] = 0.9
    dict_v22_1compSFH['ep_age_LW'] = 0.7

    file_data_V22_Z = path_data+'df44_zh_posterior_for_kristi_10_12_2021.dat'
    villaume2021_logzsol = np.loadtxt( file_data_V22_Z, unpack=1)
    qs = np.quantile(villaume2021_logzsol, quantiles)
    dqs = np.diff(qs)
    for key in ['logzsol']:
        dict_v22_1compSFH[key] = qs[1]
        dict_v22_1compSFH['em_'+key] = dqs[0]
        dict_v22_1compSFH['ep_'+key] = dqs[1]

    dict_v22_2compSFH['source'] = '2_comp'

    dict_v22_2compSFH['FeH'] = -1.29
    dict_v22_2compSFH['em_FeH'] = 0.02
    dict_v22_2compSFH['ep_FeH'] = 0.028

    dict_v22_2compSFH['age_LW'] = 9.8
    dict_v22_2compSFH['em_age_LW'] = 0.7
    dict_v22_2compSFH['ep_age_LW'] = 0.9

    dict_v22_2compSFH['age_MW'] = 9.7
    dict_v22_2compSFH['em_age_MW'] = 0.9
    dict_v22_2compSFH['ep_age_MW'] = 1.1

    df_villaume2022 = pd.concat([ pd.Series(dict_v22_1compSFH).to_frame().T, pd.Series(dict_v22_2compSFH).to_frame().T ])
    return df_villaume2022

def read_tremmel2022( path_data='', **extras  ):
    return pd.read_csv( path_data+'tremmel_data_combined.csv')
