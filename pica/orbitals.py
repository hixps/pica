import numpy as np
from scipy.interpolate import interp1d
from scipy import integrate
from random import random

from .__init__ import __version__



"""
Routines for loading atomic orbitals from file

"""


element_names = {1: "H", 2: "HE", 3: "LI", 4: "BE", 5: "B", 6: "C", 7: "N",
                 8: "O", 9: "F", 10: "NE", 11: "NA", 12: "MG", 13: "AL",
                 14: "SI", 15: "P", 16: "S", 17: "CL", 18: "AR", 19: "K",
                 20: "CA", 21: "SC", 22: "TI", 23: "V", 24: "CR", 25: "MN",
                 26: "FE", 27: "CO", 28: "NI", 29: "CU", 30: "ZN", 31: "GA",
                 32: "GE", 33: "AS", 34: "SE", 35: "BR", 36: "KR", 37: "RB",
                 38: "SR", 39: "Y", 40: "ZR", 41: "NB", 42: "MO", 43: "TC",
                 44: "RU", 45: "RH", 46: "PD", 47: "AG", 48: "CD", 49: "IN",
                 50: "SN", 51: "SB", 52: "TE", 53: "I", 54: "XE", 55: "CS",
                 56: "BA", 57: "LA", 58: "CE", 59: "PR", 60: "ND", 61: "PM",
                 62: "SM", 63: "EU", 64: "GD", 65: "TB", 66: "DY", 67: "HO",
                 68: "ER", 69: "TM", 70: "YB", 71: "LU", 72: "HF", 73: "TA",
                 74: "W", 75: "RE", 76: "OS", 77: "IR", 78: "PT", 79: "AU",
                 80: "HG", 81: "TL", 82: "PB", 83: "BI", 84: "PO", 85: "AT",
                 86: "RN", 87: "FR", 88: "RA", 89: "AC", 90: "TH", 91: "PA",
                 92: "U", 93: "NP", 94: "PU", 95: "AM", 96: "CM", 97: "BK",
                 98: "CF", 99: "ES", 100: "FM", 101: "MD", 102: "NO", 103: "LR"
                 }



# def SampleFromCDF( InvCDF, nsamples=1 ):
#     r        = np.random.uniform( 0 , 1, nsamples )
#     p_sample = InvCDF( r )
#     return p_sample

def InterpolateCDF( p, prob_density ):

    # p, PDF = CalculatePDF( *args )
    norm    = np.trapz(prob_density*p,np.log(p))
    PDF     = prob_density / norm

    print (PDF.shape)
    # Forward integration with subtraction
    cdf = integrate.cumtrapz( PDF*p ,  np.log( p ), initial=0 )
    CDF = cdf + (1.-cdf[-1])
    return p, CDF


def InterpolateInvCDF(p, prob_density ):
    p,CDF  = InterpolateCDF(p, prob_density )

    InvCDF = interp1d( CDF, p, bounds_error=False, assume_sorted=True, fill_value=(0, 1) ) 
    # Inverse CDF as callable function
    return InvCDF
    

class OrbitalSampler():

    def __init__(self, input_dict):

        self.input_dict = input_dict

        self.read_orbital_parameters()
        self.source_filename = f'{self.path}/{self.Z:03d}{element_names[self.Z]}.GAM'

        # meta data lines in the source file have to be commented out
        _density = np.loadtxt(self.source_filename).T
        p        = _density[0]
        density  = _density[1:] 
        prob_density = density[self.subshells.index(self.shell)] * p**2

        self.InvCDF = InterpolateInvCDF(p, prob_density )



    def read_orbital_parameters(self):

        # extract detector parameters
        self.Z           = int(self.input_dict['orbital']['Z'])
        self.path        = str(self.input_dict['orbital']['path'])
        self.shell       = str(self.input_dict['orbital']['shell']).upper()


    @property
    def subshells(self):
        with open(self.source_filename, 'r') as file:
            lines = file.readlines()
            for line in lines:
                if line.find("---") != -1:
                    x = lines.index(line)
            subshells = lines[x+1].strip().split()
        return subshells[2:]

    def SampleFromCDF( self, nsamples ):
        r        = np.random.uniform( 0 , 1, nsamples )
        p_sample = self.InvCDF( r )
        return p_sample

