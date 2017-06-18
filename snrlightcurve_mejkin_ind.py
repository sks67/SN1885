
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Light Curve Analysis of SN 1885, keeping Ejecta energy and mass
# independent, and varying them as parameters. 
# 
# For M31
# 1) Distance = 785 Mpc McConnachie (2005), MNRAS, 356, 979                                                                                                  # 2) flux limit = 2.7 muJy/beam
# 3) beam size = 1 arcsec
# 4) frequency = 4.86 GHz 
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import numpy as np
import matplotlib.pyplot as py
import math as mt
        
def isDetect(fluxLim=0.0038, beamSize=0.04, diamSNR=0., lumSNR=0.):
    dist = 8.5e-3 #Mpc
    lum = lambda x: 1.0e24*1.2*x*dist*dist  #x is radio flux in mJy
    numBeams = 1.0 if diamSNR<=beamSize else np.pi*((diamSNR/(2.0*beamSize))**2)
    lumLimit = lum(fluxLim*np.sqrt(numBeams))
    return lumSNR<lumLimit

def lightcurve(n0 = 1.0, e51 = 1.0, mej = 2.0, pp=2.5,  epse=0.0, tsnap=500.0,sntype='cc',stop_vel200=False,stop_sblimit=False):
    #Declaring constants for Luminosity Calculation
    if n0 == 0.:
        return ( 0., 0., 0.)

    #pp = 2.5
    compf = 4.0
    ff = 0.38
    nei = 1.14
    cl = 6.27e18
    c5 = 9.68e-24
    c6 = 8.10e-41
    dist = 8.5*1.0e-3*1000*3.086e18  #m33 distance, for now
    nu = 1.425e9 
    mu=1.4
    mp = 1.67e-24
    ismrho = n0*mu*mp #in units of cm^-3
    me = 9.11e-28  #in grams
    c = 3.0e10
 
    #Characteristic variables - TM99
    tch = 423*(e51**(-0.5))*(mej**(5.0/6.0))*(n0**(-1.0/3.0)) #years
    rch = 3.07*(mej**(1.0/3.0))*(n0**(-1.0/3.0)) #pcs
    vch = 7090*(e51**0.5)*(mej**(-0.5)) #km/s
    
    #Check whether to apply Ia or CC physics of TM99
    if(sntype=='ia'):
        tstar_st = 0.481
        t_st0 = tstar_st*tch
        tstar_ed = tsnap/tch
        vstar_ed = 0.805*tstar_ed**(-(3.0/10.0)) if tsnap<t_st0 else 0.569*((1.42*tstar_ed - 0.2935)**(-3.0/5.0))
        rstar_ed = 1.15*tstar_ed**(7.0/10.0) if tsnap<t_st0 else (1.42*tstar_ed - 0.2935)**(2.0/5.0)
        v_ed = vstar_ed*vch  #in km/s                                                                                                                                         
        r_ed = rstar_ed*rch  #in par
        
    elif(sntype=='cc'):
        tstar_st = 0.424
        t_st0 = tstar_st*tch
        tstar_ed = tsnap/tch
        vstar_ed = 0.906*tstar_ed**(-(1.0/4.0)) if tsnap<t_st0 else 0.569*((1.42*tstar_ed - 0.28)**(-3.0/5.0))
        rstar_ed = 1.20906*tstar_ed**(3.0/4.0) if tsnap<t_st0 else (1.42*tstar_ed - 0.28)**(2.0/5.0)
        v_ed = vstar_ed*vch  #in km/s                                                                                                                                         
        r_ed = rstar_ed*rch  #in par
    #Check if SNR is past Radiative phase
    vrad = 200.
    if stop_vel200:
        if (v_ed<=vrad):
            return ( 0.,0.,0.)
#            return (0., 0.)

    
   
    #Luminosity Calculation
    rho0 = ismrho/((1.0e-24)*mu*1.67)
    bism = 9e-6*rho0**0.47    #magnetic field for galaxies scales as 0.47 (Krutcher et al. 1999)                                                        
    rr1 = r_ed*(3.086e18) #LAURA - cm                                                                                                                   
    vv1 = v_ed*1.0e5 #LAURA - cm/s                                                                                                                      
#    bb0 = np.sqrt(8.0*mt.pi*epsb*ismrho*vv1*vv1)
    
####### EPSILON_B AND EPSILON_E ###############################
    #epse = 0.01#4.0e-3
    f = 0.5
    epsCR = 0.06
    v_A = bism/np.sqrt(4.0*np.pi*ismrho)
    AlfMach0 = vv1/v_A

    if (vv1/c) >= 1/AlfMach0:
        epsb_up = f*epsCR*(vv1/c)
    elif epsCR <= (1.5e-3*AlfMach0 + 0.06):
        epsb_up = (f*epsCR)/AlfMach0
    else:
        epsCR = 1.5e-3*AlfMach0 + 0.06
        epsb_up = (f*epsCR)/AlfMach0

    alpha = (1.0 + 2.0*compf*compf)/3.0
    epsb_down = alpha*epsb_up


    B_up = np.sqrt(8*np.pi*epsb_up*ismrho*vv1*vv1)
    bb1 = np.sqrt(alpha)*B_up
###############################################################

    if bb1<(4.0*bism):
       bb1=4.0*bism
    

    gammam = (mu*epse*(pp-2)*((vv1/c)**2)*mp)/((pp-1)*compf*nei*me)
    if gammam<1.0:
        gammam=1.0

    if 0: #Check epsilon_B
        return (epsb_up, gammam)

    elow = gammam*me*c**2  #Eq. 10 Chevalier, 98                                                                                                                       
    #            n0 = (alpha*(bb0**2)*(pp-2)*(elow**(pp-2.0)))/(8.0*mt.pi)
    n_0 = (pp-2)*epse*ismrho*(vv1**2)*(elow**(pp-2))
    ss11 = (4.0/3.0)*ff*rr1
    nu11 = 2.0*cl*((ss11*c6*n_0)**(2.0/(pp+4.0)))*(bb1**((pp+2.0)/(pp+4.0)))
    ss1 = (c5/c6)*(bb1**(-0.5))*((nu11/(2.0*cl))**2.5)
    #This test should be fine, since nu/nu1 becomes 1 very early for SN 1a (10^-4 years), while CCSN even earlier. t_st ~ 10^2-10^3 years. 
    if tsnap<t_st0:
        jj1 = ((nu/nu11)**2.5)*(1.0-np.exp(-1.0*(nu/nu11)**(-1.0*(pp+4.0)/2.0)))
    else:
        jj1 = ((nu/nu11)**2.5)*((nu/nu11)**(-1.0*(pp+4.0)/2.0))
    f_nu1 = (ss1*jj1*mt.pi*rr1**2)/(1.0e-26*dist**2)
    lluni = f_nu1*1.0e-26*4.0*mt.pi*dist**2

    if stop_sblimit:
        if isDetect(diamSNR=2.0*r_ed, lumSNR=lluni):
            r_ed, lluni, v_ed = 0., 0., 0.
    
    if lluni==0.0 or r_ed==0.0 or v_ed==0.0:
        lluni,r_ed,v_ed = 0.0, 0.0, 0.0
          
    return (r_ed, lluni, v_ed)
 
def lightcurve_Full(n0 = 1.0, e51 = 1.0, mej = 2.0, epse=0.0, pp=2.5, sntype='cc', plot_LightCurve=False, plot_Radius = False, plot_EjectaVelocity = False, return_vals = True, stop_vel200 = True, stop_sblimit = False):

    tim = np.logspace(-6,6,5000)
    lum = np.zeros_like(tim)
    rad = np.zeros_like(tim)
    vel = np.zeros_like(tim)
    gammam = np.zeros_like(tim)
    for ind_t,t in enumerate(tim):
         rad[ind_t], lum[ind_t], vel[ind_t] = lightcurve(n0=n0,e51=e51,mej=mej,epse=epse,pp=pp,tsnap=t,sntype=sntype,stop_vel200=stop_vel200, stop_sblimit=stop_sblimit)
    
    #Label Strings    
    nh_string = r'$n_H = {}$'.format(n0)+r' $\rm{g/cm^{-3}}$'
    e51_string = r'$E = {}$'.format(e51)+r' $\times \rm{10^{51} ergs}$'
    mej_string = r'$M = {}$'.format(mej)+r' $\rm{M_{\odot}}$'
   # epsb_string = r'$\epsilon_B = {}$'.format(epsb)
   # epse_string = r'$\epsilon_E = {}$'.format(epse)
    
    #Plots
    
    if plot_LightCurve:
        fig = py.figure(1)
        ax = fig.add_subplot(111)
        ax.plot(tim,lum,'.',lw=1.5)
        #ax.text(0.7,0.5,nh_string+'\n'+e51_string+'\n'+mej_string+'\n'+epsb_string+'\n'+epse_string,fontsize=15,transform=ax.transAxes)
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_xlabel(r't [years]',fontsize='18')
        ax.set_ylabel(r'1.4 GHz Luminosity (ergs/s/Hz)',fontsize='18')
        ax.set_xlim(100,1.0e6)
        ax.set_ylim(1.0e22,1.0e26)                                                                                                                
        ax.legend()
  
    if plot_Radius:
        fig2 = py.figure(2)
        ax2 = fig2.add_subplot(111)
        #ax2.axvline(t_st0,linestyle='--',lw=1.5)
        ax2.plot(tim,rad,'r.',lw=1.5)#,label=r'$E = {}$'.format(str(delt_array[ind])))                                                                      
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel(r't [years]',fontsize='18')
        ax2.set_ylabel(r'Radius [pc]',fontsize='18')
        ax2.set_xlim(100,1.0e6)
        ax2.legend(loc=4)

    if plot_EjectaVelocity:
        fig3 = py.figure(3)
        ax3 = fig3.add_subplot(111)
       # ax3.axvline(t_st0,linestyle='--',lw=1.5)
        ax3.plot(tim,vel,'r-')
        ax3.set_xscale('log')
        ax3.set_yscale('log')
        ax3.set_xlabel(r't [years]',fontsize='18')
        ax3.set_ylabel(r'Velocity [km/s]',fontsize='18')
        ax3.set_xlim(100,1.0e6)
        ax3.legend(loc=3)
        py.tick_params(labelsize='16')

        py.show()

    if return_vals: 
        return (tim,rad,lum,vel)
