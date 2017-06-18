#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Program to determine the minimum age by ploting the parameter space 
# for G1.9+0.3 for changes in fixed parameters.
#
# Date: 06/13/2017
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


import snrlightcurve_freqparam as snr_ind
reload(snr_ind)
import g19
reload(g19)
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
fileext = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/SN1885/Diagnostics/'
params = {'axes.linewidth':1.5,'lines.linewidth':1.3,'xtick.labelsize':20,'ytick.labelsize':20,\
          'xtick.major.size':7,'xtick.major.width':2,'ytick.major.size':7,'ytick.major.width':2,\
          'xtick.minor.size':4,'xtick.minor.width':1.5,'ytick.minor.size':4,'ytick.minor.width':1.5}
plt.rcParams.update(params)
plt.rcParams.update({'figure.autolayout': True})

def flux_to_lum(flux, fluxerr, dist, e_dist):
    lum = 1.0e24*1.2*flux*dist*dist
    lumerr = 1.0e24*np.sqrt((2.4*flux*dist*e_dist)**2.0 + (1.2*dist*dist*fluxerr)**2.0)
    return (lum, lumerr)

def ageConstraintPlots(velmin, velmax):
    pdf_obj = PdfPages(fileext+'G19Age_LowLim_velmax_{}.pdf'.format(str(velmax)))
    print '######## Age Constraint for Max Velocity = ', velmax, ' km/s ########\n\n'
    for age_now in ages:
        for density in dens:

            print 'age = ', age_now, 'density = ', density
            age_bef = age_now-16
            
            total_expan = expan_rate*(age_now - age_bef)
            rad_bef = rad_now*(1.0 - (total_expan/100.))
            vel_now = ((rad_now - rad_bef)*pc_to_km)/((age_now-age_bef)*yr_to_sec)
            
            mejarray_vel_0, y_e51_vel_high_0 = g19.Ekin_Mej_VelocityConstraint_G19(vel=velmin, density=density, epse=epse, \
                                                                                       pp=pp, age=age_now)
            mejarray_vel_0, y_e51_vel_low_0 = g19.Ekin_Mej_VelocityConstraint_G19(vel=velmax, density=density, epse=epse, \
                                                                                      pp=pp, age=age_now)

            mejarray_rad_0, y_e51_rad_high_0 = g19.Ekin_Mej_RadiusConstraint_G19(rad=rad_now, density=density, epse=epse, \
                                                                                     pp=pp, age=age_now)
            mejarray_rad_0, y_e51_rad_low_0 = g19.Ekin_Mej_RadiusConstraint_G19(rad=rad_bef, density=density, epse=epse, \
                                                                                    pp=pp, age=age_bef)

            plt.figure(figsize=(9, 9))
            axes_1 = plt.gca()

            if 1:
                axes_1.fill_between(x=mejarray_vel_0, y1=y_e51_vel_low_0, y2=y_e51_vel_high_0, \
                     color='#1b9e77', alpha=0.7, label='Velocity constraint; {}'.format(density)+units[0])

            if 1:
                axes_1.fill_between(x=mejarray_rad_0, y1=y_e51_rad_low_0, y2=y_e51_rad_high_0, \
                     color='#d95f02', alpha=0.9, label='Radius constraint; {}'.format(density)+units[0])


            axes_1.set_title('Age = {} yr'.format(age_now), fontsize=20)
            axes_1.set_xlim(0.2,1.7)
            axes_1.set_ylim(0.1, 10.)
            axes_1.set_xlabel(r'Ejecta Mass [$\rm{M_{\odot}}$]', fontsize=24)
            axes_1.set_ylabel(r'Kinetic Energy $E_{51}$ [$\rm{\times\ 10^{51}\ ergs}$]', fontsize=24)
            axes_1.set_yscale('log')
            axes_1.legend(loc=2, fontsize=20)
            pdf_obj.savefig()
    pdf_obj.close()

#PROPERTIES OF G1.9
lum1, e_lum1 = flux_to_lum(0.74*1.0e3, 0.038*1.0e3, 8.5e-3, 0.0)
#Condon (1998) 1.425 GHz assumed to be taken in 1993, taken from Green (2008)
lum2, e_lum2 = flux_to_lum(0.935*1.0e3, 0.047*1.0e3, 8.5e-3, 0.0)
#1.4 GHz Green (2008) measurement

expan_rate = 0.65 # % per year, based on Green 2009
rad_now = 2.0
pc_to_km = 3.086e13 #km
yr_to_sec = 3.154e7 #s
freq=1.425e9
epse=1.0e-4
pp=2.2
units = [r'cm$^{-3}$']

ages = [80, 90, 100, 110, 120]
dens = [0.03, 0.1]


ageConstraintPlots(5300, 14000)

