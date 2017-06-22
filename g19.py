#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Light curve analysis of G1.9+0.3, using my light curve model. Will add
# more description here in future.
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


import snrlightcurve_freqparam as snr_ind
import numpy as np
import matplotlib.pyplot as plt
reload(snr_ind)

def flux_to_lum(flux, fluxerr, dist, e_dist):
    lum = 1.0e24*1.2*flux*dist*dist
    lumerr = 1.0e24*np.sqrt((2.4*flux*dist*e_dist)**2.0 + (1.2*dist*dist*fluxerr)**2.0)
    return (lum, lumerr)


#~~~~~~~~~PROPERTIES OF G1.9~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
lum1, e_lum1 = flux_to_lum(0.74*1.0e3, 0.038*1.0e3, 8.5e-3, 0.0) #Condon (1998) 1.425 GHz assumed to be taken in 1993, taken from Green (2008)
lum2, e_lum2 = flux_to_lum(0.935*1.0e3, 0.047*1.0e3, 8.5e-3, 0.0) #1.4 GHz Green (2008) measurement
rad = 2.0 #Reynolds (2008) said mean RADIUS of X-ray profile = 2 pc

raderr = 0.0
freq = 1.425e9 #GHz
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


def Ekin_Mej_ParameterSpace_G19(density=0.01, epse=0.0001, pp=2.2, age=140., freq=1.425e9):
    e51array = np.logspace(-1, 1, 800)
    mejarray = np.linspace(0.1,2.0,800)
    mejBest_0, e51Best_0 = [],[]
    age1, age2 = 126., 140. #These are kinematic age upper limits based on uniform expansion.
    for i, e51 in enumerate(e51array):
        for j, mej in enumerate(mejarray):
            model_R_1, model_L_1, vel = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp,\
                                                       freq=freq, tsnap=age1, sntype='ia')
            model_R_2, model_L_2, vel = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp,\
                                                       freq=freq, tsnap=age2, sntype='ia')
            if model_L_1 > lum1-e_lum1 and model_L_2 < lum2+e_lum2:
                mejBest_0.append(mej), e51Best_0.append(e51)
        
    return (mejBest_0, e51Best_0)

def radius_ED(t, E, n0, mej):
    return 1.29*((t/100.)**0.7)*(E**0.35)*(n0**(-0.1))*(mej**(-0.25))

def Ekin_Mej_RadiusConstraint_G19(rad=1.0, density=0.01, age=140., epse=1.0e-4, pp=2.2, freq=1.4e9):

    e51array = np.logspace(-1,1,800)
    mejarray = np.linspace(0.1, 2.0, 800)
    y_e51 = np.zeros(mejarray.size)
    for i, mej in enumerate(mejarray):
        rad_mej = np.zeros(e51array.size)
        for j, e51 in enumerate(e51array):
            rad_mej[j], lum, vel = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=age, sntype='ia')
        bestval = e51array[np.argmin(np.absolute(rad_mej - rad))]
        if (rad_mej[np.argmin(np.absolute(rad_mej-rad))] - rad)<0.01:
            y_e51[i]=bestval

    return (mejarray, y_e51)

def Ekin_Mej_VelocityConstraint_G19(vel=1.0, density=0.01, age=140., epse=1.0e-4, pp=2.2, freq=1.4e9):

    e51array = np.logspace(-1,1,800)
    mejarray = np.linspace(0.1, 2.0, 800)
    y_e51 = np.zeros(mejarray.size)
    for i, mej in enumerate(mejarray):
        vel_mej = np.zeros(e51array.size)
        for j, e51 in enumerate(e51array):
            rad, lum, vel_mej[j] = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=age, sntype='ia')
        bestval = e51array[np.argmin(np.absolute(vel_mej - vel))]
        y_e51[i]=bestval

    return (mejarray, y_e51)

def Ekin_Mej_LumConstraint_G19(lum_obs=1.0e24, density=0.01, age=140., epse=1.0e-4, pp=2.2, freq=1.4e9):

    e51array = np.logspace(-1,1,800)
    mejarray = np.linspace(0.1, 2.0, 800)
    y_e51 = np.zeros(mejarray.size)
    for i, mej in enumerate(mejarray):
        lum = np.zeros(e51array.size)
        for j, e51 in enumerate(e51array):
            rad, lum[j], vel = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=age, sntype='ia')
        bestval = e51array[np.argmin(np.absolute(np.log10(lum) - np.log10(lum_obs)))]
      #  if (np.log10(lum[np.argmin(np.absolute(np.log10(lum)-np.log10(lum_obs)))]) - np.log10(lum_obs))<0.01:
        y_e51[i]=bestval

    return (mejarray, y_e51)

def Ekin_Mej_RadiusConstraint_G19_twoages(density=0.01, epse=0.0001, pp=2.2, age=140., freq=1.425e9):
    e51array = np.logspace(-1, 1, 800)
    mejarray = np.linspace(0.1,2.0,800)
    mejBest_0, e51Best_0 = [],[]
    age1, age2 = 118., 140. #Uniform expansion, spanning observations between 1985 and 2007, Reynolds (2008)
    rad1, rad2 = 1.72, 2.0

    for i, e51 in enumerate(e51array):
        for j, mej in enumerate(mejarray):
            model_R_1, model_L_1, vel = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp,\
                                                       freq=freq, tsnap=age1, sntype='ia')
            model_R_2, model_L_2, vel = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp,\
                                                       freq=freq, tsnap=age2, sntype='ia')
            if model_R_1 >= rad1 and model_R_2 <= rad2:
                mejBest_0.append(mej), e51Best_0.append(e51)
        
    return (mejBest_0, e51Best_0)
