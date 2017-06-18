import snrlightcurve_freqparam as snr_ind
reload(snr_ind)
import g19
reload(g19)
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt
import numpy as np
fileext = '/Users/sumits2k/Desktop/Research/SNResearch2/RadioSNRs/SN1885/Diagnostics/ParameterSpaces/'
params = {'axes.linewidth':1.5,'lines.linewidth':1.3,'xtick.labelsize':20,'ytick.labelsize':20,\
          'xtick.major.size':7,'xtick.major.width':2,'ytick.major.size':7,'ytick.major.width':2,\
          'xtick.minor.size':4,'xtick.minor.width':1.5,'ytick.minor.size':4,'ytick.minor.width':1.5}
plt.rcParams.update(params)
plt.rcParams.update({'figure.autolayout': True})
pdf_obj = PdfPages(fileext+'SN1885_Changes_in_pp.pdf')

epse = 0.0001

fluxsn1885 = 0.0038*3.0 #mJy, Chomiuk's upper limit from document
dist = 0.785 #McConnachie 2005 - http://adsabs.harvard.edu/abs/2005MNRAS.356..979M
e_dist = 0.025
lumsn1885 = 1.0e24*1.2*fluxsn1885*dist*dist
yerr = 1.0e24*np.sqrt((2.4*fluxsn1885*dist*e_dist)**2.0)
rad1885 = 1.52 #Fesen07
raderr = 0.15
freq = 6.23e9 #Hz, observation from Laura
pp_array = [2.2, 2.3, 2.4, 2.5]
density = [0.01]
print '\n ####### Checking changes in pp  ######### \n'

for pp in pp_array:
    print 'pp = ', pp
    e51array = np.logspace(-1, 1, 800)
    mejarray = np.linspace(0.1,2.0,800)
    M, e = np.meshgrid(mejarray, e51array)
    stat_1 = np.zeros((e51array.size, mejarray.size))
    stat_2 = np.zeros((e51array.size, mejarray.size))
    stat_3 = np.zeros((e51array.size, mejarray.size))
    
    for i, e51 in enumerate(e51array):
        for j, mej in enumerate(mejarray):
            rad, lum, vel = snr_ind.lightcurve(n0=density[0], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
            stat_1[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[1], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_2[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[2], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_3[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.


    plt.figure(figsize=(10,11))
    plt.rc('font', family='serif')
    cs = plt.contourf(M,e,stat_1,40, cmap='gray_r')
#plt.contourf(M,e,stat_2,40,cmap='gray_r',alpha=0.6)
#plt.contourf(M,e,stat_3,40,cmap='gray_r',alpha=0.4)
#Change 'alpha' with respect to stat_1. That means alpha=0.8 will give you the lighter color, not alpha=0.2

#plt.annotate(s='', xy=(1.4,1.4), xytext=(1.4,0.9), arrowprops={'arrowstyle':'<|-|>', \
 #                                                      'lw':2, 'color':'red'})
    if 1:
        plt.text(0.8, 1.5, 'Normal Ia\n(Badenes 07)', color='#CC3300', fontsize=18) 
        plt.axhspan(0.9, 1.4, color='#CC3300', alpha=0.5)
#plt.vlines(x=1.4, ymin=0.9, ymax=1.4, color='r', lw=3.0)
    plt.title('pp = {}'.format(str(pp))+r', $\epsilon_e$ = '+str(epse)+', density = '+str(density[0]), fontsize=16)
    plt.yscale('log')
    plt.ylim(0.1,5.0)
    plt.xlim(0.1, 1.7)
    plt.xlabel(r'Ejecta Mass [$\rm{M_{\odot}}$]', fontsize=24)
    plt.ylabel(r'Kinetic Energy [$\rm{\times\ 10^{51}\ ergs}$]', fontsize=24)
    plt.legend(numpoints=1, loc=2, fontsize=20, borderpad = 0.6, handletextpad=0.6, labelspacing=0.5)
    pdf_obj.savefig()

print '\n ####### Checking changes in density for pp = 2.5 ######### \n'

pp=2.5
density = [0.03, 0.04, 0.05, 0.06, 0.07]

for dens in density:
    print 'density = ', dens
    e51array = np.logspace(-1, 1, 800)
    mejarray = np.linspace(0.1,2.0,800)
    M, e = np.meshgrid(mejarray, e51array)
    stat_1 = np.zeros((e51array.size, mejarray.size))
    stat_2 = np.zeros((e51array.size, mejarray.size))
    stat_3 = np.zeros((e51array.size, mejarray.size))
    
    for i, e51 in enumerate(e51array):
        for j, mej in enumerate(mejarray):
            rad, lum, vel = snr_ind.lightcurve(n0=dens, mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
            stat_1[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[1], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_2[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[2], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_3[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.


    plt.figure(figsize=(10,11))
    plt.rc('font', family='serif')
    cs = plt.contourf(M,e,stat_1,40, cmap='gray_r')
#plt.contourf(M,e,stat_2,40,cmap='gray_r',alpha=0.6)
#plt.contourf(M,e,stat_3,40,cmap='gray_r',alpha=0.4)
#Change 'alpha' with respect to stat_1. That means alpha=0.8 will give you the lighter color, not alpha=0.2

#plt.annotate(s='', xy=(1.4,1.4), xytext=(1.4,0.9), arrowprops={'arrowstyle':'<|-|>', \
 #                                                      'lw':2, 'color':'red'})
    if 1:
        plt.text(0.8, 1.5, 'Normal Ia\n(Badenes 07)', color='#CC3300', fontsize=18) 
        plt.axhspan(0.9, 1.4, color='#CC3300', alpha=0.5)
#plt.vlines(x=1.4, ymin=0.9, ymax=1.4, color='r', lw=3.0)
    plt.title('pp = {}'.format(str(pp))+r', $\epsilon_e$ = '+str(epse)+', density = '+str(dens), fontsize=16)
    plt.yscale('log')
    plt.ylim(0.1,5.0)
    plt.xlim(0.1, 1.7)
    plt.xlabel(r'Ejecta Mass [$\rm{M_{\odot}}$]', fontsize=24)
    plt.ylabel(r'Kinetic Energy [$\rm{\times\ 10^{51}\ ergs}$]', fontsize=24)
    plt.legend(numpoints=1, loc=2, fontsize=20, borderpad = 0.6, handletextpad=0.6, labelspacing=0.5)
    pdf_obj.savefig()

print '\n ####### Checking changes in epse ######### \n'

pp=2.2
density = 0.005
epse_array = [1.0e-4, 4.0e-4, 7.0e-4, 1.0e-3]

for epse in epse_array:
    print 'epse = ', epse
    e51array = np.logspace(-1, 1, 800)
    mejarray = np.linspace(0.1,2.0,800)
    M, e = np.meshgrid(mejarray, e51array)
    stat_1 = np.zeros((e51array.size, mejarray.size))
    stat_2 = np.zeros((e51array.size, mejarray.size))
    stat_3 = np.zeros((e51array.size, mejarray.size))
    
    for i, e51 in enumerate(e51array):
        for j, mej in enumerate(mejarray):
            rad, lum, vel = snr_ind.lightcurve(n0=density, mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
            stat_1[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[1], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_2[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[2], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_3[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.


    plt.figure(figsize=(10,11))
    plt.rc('font', family='serif')
    cs = plt.contourf(M,e,stat_1,40, cmap='gray_r')
#plt.contourf(M,e,stat_2,40,cmap='gray_r',alpha=0.6)
#plt.contourf(M,e,stat_3,40,cmap='gray_r',alpha=0.4)
#Change 'alpha' with respect to stat_1. That means alpha=0.8 will give you the lighter color, not alpha=0.2

#plt.annotate(s='', xy=(1.4,1.4), xytext=(1.4,0.9), arrowprops={'arrowstyle':'<|-|>', \
 #                                                      'lw':2, 'color':'red'})
    if 1:
        plt.text(0.8, 1.5, 'Normal Ia\n(Badenes 07)', color='#CC3300', fontsize=18) 
        plt.axhspan(0.9, 1.4, color='#CC3300', alpha=0.5)
#plt.vlines(x=1.4, ymin=0.9, ymax=1.4, color='r', lw=3.0)
    plt.title('pp = {}'.format(str(pp))+r', $\epsilon_e$ = '+str(epse)+', density = '+str(density), fontsize=16)
    plt.yscale('log')
    plt.ylim(0.1,5.0)
    plt.xlim(0.1, 1.7)
    plt.xlabel(r'Ejecta Mass [$\rm{M_{\odot}}$]', fontsize=24)
    plt.ylabel(r'Kinetic Energy [$\rm{\times\ 10^{51}\ ergs}$]', fontsize=24)
    plt.legend(numpoints=1, loc=2, fontsize=20, borderpad = 0.6, handletextpad=0.6, labelspacing=0.5)
    pdf_obj.savefig()

print '\n ####### Checking changes in densities for epse = 1.0e-3 ######### \n'

pp=2.2
density = [0.001, 0.003, 0.005, 0.007, 0.009]
epse = 1.0e-3

for dens in density:
    print 'dens = ', dens
    e51array = np.logspace(-1, 1, 800)
    mejarray = np.linspace(0.1,2.0,800)
    M, e = np.meshgrid(mejarray, e51array)
    stat_1 = np.zeros((e51array.size, mejarray.size))
    stat_2 = np.zeros((e51array.size, mejarray.size))
    stat_3 = np.zeros((e51array.size, mejarray.size))
    
    for i, e51 in enumerate(e51array):
        for j, mej in enumerate(mejarray):
            rad, lum, vel = snr_ind.lightcurve(n0=dens, mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
            stat_1[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[1], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_2[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.
        #rad, lum, vel = snr_ind.lightcurve(n0=density[2], mej=mej, e51=e51, epse=epse, pp=pp, freq=freq, tsnap=127, sntype='ia')
        #stat_3[i,j] = 1. if (lum<=lumsn1885) & (rad>=rad1885) & (vel>=12500.) else 0.


    plt.figure(figsize=(10,11))
    plt.rc('font', family='serif')
    cs = plt.contourf(M,e,stat_1,40, cmap='gray_r')
#plt.contourf(M,e,stat_2,40,cmap='gray_r',alpha=0.6)
#plt.contourf(M,e,stat_3,40,cmap='gray_r',alpha=0.4)
#Change 'alpha' with respect to stat_1. That means alpha=0.8 will give you the lighter color, not alpha=0.2

#plt.annotate(s='', xy=(1.4,1.4), xytext=(1.4,0.9), arrowprops={'arrowstyle':'<|-|>', \
 #                                                      'lw':2, 'color':'red'})
    if 1:
        plt.text(0.8, 1.5, 'Normal Ia\n(Badenes 07)', color='#CC3300', fontsize=18) 
        plt.axhspan(0.9, 1.4, color='#CC3300', alpha=0.5)
#plt.vlines(x=1.4, ymin=0.9, ymax=1.4, color='r', lw=3.0)
    plt.title('pp = {}'.format(str(pp))+r', $\epsilon_e$ = '+str(epse)+', density = '+str(dens), fontsize=16)
    plt.yscale('log')
    plt.ylim(0.1,5.0)
    plt.xlim(0.1, 1.7)
    plt.xlabel(r'Ejecta Mass [$\rm{M_{\odot}}$]', fontsize=24)
    plt.ylabel(r'Kinetic Energy [$\rm{\times\ 10^{51}\ ergs}$]', fontsize=24)
    plt.legend(numpoints=1, loc=2, fontsize=20, borderpad = 0.6, handletextpad=0.6, labelspacing=0.5)
    pdf_obj.savefig()

pdf_obj.close()


