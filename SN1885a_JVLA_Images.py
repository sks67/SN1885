import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from astropy.io import fits
import sys

mpl.rcParams['figure.figsize'] = [7.0, 6.0]
mpl.rcParams['figure.dpi'] = 80
mpl.rcParams['savefig.dpi'] = 100

mpl.rcParams['font.size'] = 12
mpl.rcParams['legend.fontsize'] = 'large'
mpl.rcParams['figure.titlesize'] = 'medium'

import aplpy
from astropy.units import cds
from astropy.coordinates import Angle


path = '/Users/sumits2k/Dropbox/SN1885a/Images/'

#~~~~~~~~ Recenter units ~~~~~~~~~#
imageCen_RA = Angle('00h42m43.7s').degree
imageCen_Dec = Angle('41d16m06s').degree
image_width = (18.0*cds.arcs).to('deg').value
image_height = (14.0*cds.arcs).to('deg').value
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~ SN1885A center ~~~~~~~~~~#
sn1885Cen_RA = Angle('00h42m43.03s').degree
sn1885Cen_Dec = Angle('41d16m04.5s').degree
sn1885_Rad = (0.4*cds.arcs).to('deg').value
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

fig = plt.figure(figsize=(18, 7))

image_Hubble_file = 'sand_ca2_koords2.fits'
image_Hubble = aplpy.FITSFigure(path+image_Hubble_file, figure=fig, subplot=[0.1, 0.1, 0.4, 0.8])
image_Hubble.show_grayscale(stretch='linear')
image_Hubble.recenter(imageCen_RA, imageCen_Dec, width=image_width, height=image_height)
image_Hubble.show_circles(sn1885Cen_RA, sn1885Cen_Dec, radius=sn1885_Rad, edgecolor='b', facecolor='', lw=3.0)
image_Hubble.show_contour(levels=[0.1, 0.2, 0.3, 0.4], colors='b', stretch='log')
image_Hubble.add_label(0.2, 0.8, 'M31*\nDouble Nucleus', relative=True)
image_Hubble.add_label(0.9, 0.35, 'SN 1885A', relative=True)
image_Hubble.add_label(0.5, 1.05, 'Optical (Hubble)', relative=True, size=17)
#image_Hubble.add_colorbar('top')
image_Hubble.tick_labels.set_xformat('hh:mm:ss')
image_Hubble.set_theme('publication')

if 1:    
    image_6GHz_file = 'M31_6GHZ_COMB_R0.FITS'
    image_6GHz = aplpy.FITSFigure(path+image_6GHz_file, figure=fig, subplot=[0.45, 0.1, 0.4, 0.8])
    image_6GHz.show_grayscale(stretch='linear', vmin=0.0)
    image_6GHz.recenter(imageCen_RA, imageCen_Dec, width=image_width, height=image_height)
    image_6GHz.show_circles(sn1885Cen_RA, sn1885Cen_Dec, radius=sn1885_Rad, edgecolor='b', facecolor='', lw=3.0)
    image_6GHz.show_contour(path+image_Hubble_file, levels=[0.1, 0.2, 0.3, 0.4], colors='b')
    image_6GHz.tick_labels.set_xformat('hh:mm:ss')
    image_6GHz.add_label(0.5, 1.05, 'JVLA, 6 GHz, $\sigma=1.27$ $\mu$Jy/beam', relative=True, size=17)
#image_6GHz.add_colorbar('top')
#image_6GHz.colorbar.set_axis_label_text(r'$\mu$Jy')
#image_6GHz.colorbar.set_labels(['1','2','3','4','5'])
    image_6GHz.tick_labels.hide_y()
    image_6GHz.axis_labels.hide_y()
    image_6GHz.set_theme('publication')
    
path = '/Users/sumits2k/Dropbox/SN1885a/Images/'

#~~~~~~~~ Recenter units ~~~~~~~~~#
imageCen_RA = Angle('00h42m43.03s').degree
imageCen_Dec = Angle('41d16m04.5s').degree
image_width = (1.5*cds.arcs).to('deg').value
image_height = (1.0*cds.arcs).to('deg').value
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#~~~~~~~ SN1885A center ~~~~~~~~~~#
sn1885Cen_RA = Angle('00h42m43.03s').degree
sn1885Cen_Dec = Angle('41d16m04.5s').degree
sn1885_Rad = (0.4*cds.arcs).to('deg').value
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


image_6GHz_file = 'M31_6GHZ_COMB_R0.FITS'
image_6GHz_zoom = aplpy.FITSFigure(path+image_6GHz_file, figure=fig, subplot=[0.8, 0.1, 0.4, 0.8])
image_6GHz_zoom.show_grayscale(stretch='linear', vmin=0.0)
image_6GHz_zoom.recenter(imageCen_RA, imageCen_Dec, width=image_width, height=image_height)
image_6GHz_zoom.show_circles(sn1885Cen_RA, sn1885Cen_Dec, radius=sn1885_Rad, edgecolor='b', facecolor='', lw=3.0)
sigma = 1.9e-6
image_6GHz_zoom.show_contour(levels=[1.4*sigma, 2.8*sigma, 4.6*sigma], colors='g')
image_6GHz_zoom.tick_labels.set_xformat('hh:mm:ss')
image_6GHz_zoom.add_label(0.5, 1.05, 'JVLA, zoomed-in', relative=True, size=17)
#image_6GHz.add_colorbar('top')
#image_6GHz.colorbar.set_axis_label_text(r'$\mu$Jy')
#image_6GHz.colorbar.set_labels(['1','2','3','4','5'])
image_6GHz_zoom.tick_labels.hide_y()
image_6GHz_zoom.axis_labels.hide_y()
image_6GHz_zoom.set_theme('publication')

fig.canvas.draw()
fig.savefig('../Diagnostics/SN1885_3Panel_Hubble_6GHz.pdf', dpi=150)




