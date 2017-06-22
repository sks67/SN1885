#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# A short calculation of what are the expected
# energy range for SN Iax
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

import numpy as np
import astropy.units as u
from astropy.units import cds
mejlow = 0.2*cds.Msun.to(u.cgs.g)  #Solar masses in grams
mejhigh = 0.8*cds.Msun.to(u.cgs.g)
mejmed = 0.6*cds.Msun.to(u.cgs.g)

vellow = 5000*u.kilometer.to(u.cgs.cm)/u.s
velhigh = 6000*u.kilometer.to(u.cgs.cm)/u.s

print 'Ejecta Mass = ', mejmed 
print 'Ejecta Energy lower limit = ', 0.5*mejmed*vellow*vellow
print 'Ejecta Energy upper limit = ', 0.5*mejmed*velhigh*velhigh
print 'Ejecta Mass = ', mejlow 
print 'Ejecta Energy lower limit = ', 0.5*mejlow*vellow*vellow
print 'Ejecta Energy upper limit = ', 0.5*mejlow*velhigh*velhigh
print 'Ejecta Mass = ', mejhigh
print 'Ejecta Energy lower limit = ', 0.5*mejhigh*vellow*vellow
print 'Ejecta Energy upper limit = ', 0.5*mejhigh*velhigh*velhigh
