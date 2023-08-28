### Dimensional Conversions and Physical Constants ###################################################
import numpy as np

hbarc   = 0.197326968		# eV micro-meter  --> 1(eV)^-1 = 0.19732 micro-meter
c		= 0.299792458		# 1 micro-meter/femto-second


micron		= 1./hbarc		# 1 micro-meter in inverse eV : 1micron	== 5.06773 eV^-1
fsec		= c/hbarc		# 1 femto-second in inverse eV: 1fsec   == 1.51926 eV^-1
psec		= 1e3*fsec

finestruct	= 1./137.035999
elec_charge	= np.sqrt(4*np.pi*finestruct)	# electron charge in natural units
elec_mass	= 511.e3		# electron rest mass in eV



joule_per_eV	= 1.60217653e-19	# electron charge in SI units, conversion from eV to Joule
joule		    = 1/joule_per_eV

