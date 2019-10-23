import numpy as np
import pandas as pd
import astropy.units as u

def getTESSlum(TIC):
    '''
    Calculate the luminosity of a TESS target using Gaia colors
    interpolated onto main sequence isochrone
    Parameters
    ----------
    TIC : numpy array of integers or single int
        List of TICs or single TIC to crossmatch with Gaia
    
    Returns
    -------
    lum - log10 of the TESS luminosity. If a list of TICs was
    passed, lum is also a list. Any TICs that are not in the
    Gaia crossmatch return a luminosity of NaN
    '''
    tics_are_list = isinstance(TIC, list)
    iso = pd.read_csv('isochrones.csv', comment='#')
    tic_gaia = pd.read_csv('TESSgaia1to15.csv', comment='#')
    TIC_in_gaia = np.isin(tic_gaia['ticid'], TIC)
    TIC_in_match = np.isin(TIC, tic_gaia['ticid'])

    n = 1
    if tics_are_list:
     n = len(TIC)
    else:
        if len(tic_gaia[TIC_in_gaia]) == 0:
            print("TIC not found in Gaia crossmatch")
            exit()
    dist = np.full(n, np.nan)
    GBp = np.full(n, np.nan)
    GRp = np.full(n, np.nan)
    dist[TIC_in_match] = tic_gaia['r_est'][TIC_in_gaia].values*u.pc
    GBp[TIC_in_match] = tic_gaia['phot_bp_mean_mag'][TIC_in_gaia]
    GRp[TIC_in_match] = tic_gaia['phot_rp_mean_mag'][TIC_in_gaia]

    Bp_min_Rp = GBp - GRp
    Bp_iso = iso['G_BPbrmag'].values[::-1]
    Rp_iso = iso['G_RPmag'].values[::-1]
    G_iso = iso['Gmag'].values[::-1]
    T_iso = iso['TESSmag'].values[::-1]
    Bp_min_Rp_iso = Bp_iso - Rp_iso

    # Only use main sequence stars
    ms_mask = (G_iso > 4) & (Bp_min_Rp_iso > 0) & (Bp_min_Rp_iso < 4.5)
    T_int = np.interp(Bp_min_Rp, Bp_min_Rp_iso[ms_mask], T_iso[ms_mask])

    # Zero point TESS flux (from Sullivan 2017)
    Tf0 = 4.03e-6*u.erg/u.s/u.cm**2

    # TESS apparent magnitude
    t = T_int + 5*np.log10(dist) - 5
    flux = 10**(-t/2.5)*Tf0
    lum = (4*np.pi*((dist*u.pc).to(u.cm))**2*flux).value

    # Return a single value if we were only passed a single TIC
    if not tics_are_list:
        lum = lum[0]

    return np.log10(lum)
