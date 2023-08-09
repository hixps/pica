import numpy as np


"""
Auxiliary routines
"""



def beam_covariance_matrix(beam_size,rms_angle,L):


    cov_at_IP_Focus = np.array([[beam_size**2, 0], [0, rms_angle**2]])
     


    #Inverse of beam transport matrix
    M = np.array([[1, -L], [0, 1]])
    #Calculating the Covariance matrix at the Interaction Point
    cov=(M@cov_at_IP_Focus@M.T)
    return (cov)


def gaussian_sum_normalized( mean , rms , nbeam ):
    """
    Help generate normalized Gaussian distribtions
    """

    values        = mean + np.linspace(-3*rms, 3*rms, nbeam)
    dvalues       = values[1] - values[0]
    beam_weights  = np.exp(-(values-mean)**2/2/rms**2) / np.sqrt(2*pi) / rms * dvalues

    return values, beam_weights



