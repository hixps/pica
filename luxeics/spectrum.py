import numpy as np
# from numpy import *
# from scipy import *

from .constants import elec_mass,finestruct

"""
Routines for the radiation integral
"""


def energy_spectum(L,sigma):
    # the Fourier transform of the laser envelope, replaces the delta-function in the usual cross section


    # specify spectrum for cos^2 pulse
    laser_spec  = np.pi**4 * np.sin( (L-1) * sigma )**2 / ( (L-1)*np.pi**2 - (L-1)**3*sigma**2)**2

    return laser_spec



def F0_building_block(x,y):
    return 0.25*(x/y+y/x) + 1./x - 1./y + (1./x-1./y)**2

def F3_building_block(x,y):
    return - (1./x - 1./y) - (1./x-1./y)**2

def F11_building_block(x,y):
    return 1./x-1./y + 0.5

def F22_building_block(x,y):
    return 0.25*(x/y+y/x)*(1+2/x-2/y)

def F33_building_block(x,y):
    return (1/x-1/y)**2 + (1/x-1/y) + 0.5


def L_parameter_linear(y,s,b0):
    """
    Determine the value of L for energy spectrum
    use energy momentum conservation and that 
        1) in linear Compton scattering: y = 2b0*L(1-s)
    """
    L           = 0.5*y / b0/(1-s)
    return L


def L_parameter_a0correction(y,s,b0,a0):
    
    # approximately take into account the intensity-dependent red-shift 
    L           = 0.5* (y + s*a0**2/2.*b0) / b0/(1-s)
    return L


def Compton_spectrum_linear( U_in , a0 , omega0 , sigma , omega, theta, phi , S3 , a0_freq_correction=False ):
    # define the Baier-Katkov-type recoil factor E/(E-omega) = 1 /(1-omega/E)
    energy_in    = U_in[0]*elec_mass
    s            = omega/energy_in   # fractional energy loss


    print (s)
    print ((U_in[0] + U_in[3])/(omega*(1+np.cos(theta))))

    # direction of the outgoing photon
    nx = np.sin(theta)*np.cos(phi)
    ny = np.sin(theta)*np.sin(phi)
    nz = np.cos(theta)

    # quantum energy paramter == laser frequency in electron rest frame in units of electron rest mass
    b0          =     omega0 * ( U_in[0] + U_in[3]) / elec_mass

    # LL y-parameter (86-15), == 2* outgoing kappa-factor of Jauch-Rohrlich (up to mass)
    # y = (m^2-u)/m^2 = 2 K_out*U_in / m 
    y          = 2 * omega*( U_in[0]  - U_in[1] *nx - U_in[2] *ny - U_in[3] *nz ) / elec_mass

    if a0_freq_correction:
        L = L_parameter_a0correction(y,s,b0,a0)
    else:
        L = L_parameter_linear(y,s,b0)

    # LL x-parameter (86-15), == 2* incoming kappa-factor of Jauch-Rohrlich (up to mass)
    # x = (s-m^2)/m^2 = 2 L K_in*U_in / m
    # the L is due to the laser spectrum, == L * b0


    x          = 2 * L * omega0 * ( U_in[0] + U_in[3]) / elec_mass


    print (L)
    print (x)
    print (y)
    print (x/y+y/x)
    print ('----------')

    #X           = kap_in/kap_out + kap_out/kap_in + 2*(1/kap_in-1/kap_out) + (1/kap_in-1/kap_out)**2
    # -->  this needs to be modified for linear laser polarization!!!!

    # polarization degree with regard to scattering plane
    xi3         = S3 * np.cos(2*phi)


    # print (S3)
    # print (xi3)

    F           = F0_building_block(x,y) + xi3 * F3_building_block(x,y)

    print (all(F>0))

    prefactor   = finestruct * a0**2 * omega / (32*np.pi**2) / elec_mass**2 / b0**2 / (1-s)


    laser_spec  = energy_spectum(L,sigma)



    return prefactor * laser_spec * F * 8


def StokesParameters( U_in , a0 , omega0 , sigma , omega, theta, phi , S3 , a0_freq_correction=False ):
    energy_in    = U_in[0]*elec_mass
    s            = omega/energy_in   # fractional energy loss


    # direction of the outgoing photon
    nx = np.sin(theta)*np.cos(phi)
    ny = np.sin(theta)*np.sin(phi)
    nz = np.cos(theta)

    # quantum energy paramter == laser frequency in electron rest frame in units of electron rest mass
    b0          =     omega0 * ( U_in[0] + U_in[3]) / elec_mass

    print (b0)

    # LL y-parameter (86-15), == 2* outgoing kappa-factor of Jauch-Rohrlich (up to mass)
    # y = (m^2-u)/m^2 = 2 K_out*U_in / m 
    y          = 2 * omega*( U_in[0]  - U_in[1] *nx - U_in[2] *ny - U_in[3] *nz ) / elec_mass

    if a0_freq_correction:
        L = L_parameter_a0correction(y,s,b0,a0)
    else:
        L = L_parameter_linear(y,s,b0)


    # x          = 2 * L * omega0 * ( U_in[0] + U_in[3]) / elec_mass
    x          = 2 * omega0 * ( U_in[0] + U_in[3]) / elec_mass

    print ('vvvvvvvvvv')
    print (np.amin(L),np.amax(L))
    print (x)
    print (y)
    print (x/y+y/x)
    print ('^^^^^^^^^^')


    #X           = kap_in/kap_out + kap_out/kap_in + 2*(1/kap_in-1/kap_out) + (1/kap_in-1/kap_out)**2
    # -->  this needs to be modified for linear laser polarization!!!!

    # polarization degree with regard to scattering plane
    xi3         = S3 * np.cos(2*phi)
    xi1         = S3 * np.sin(2*phi)

    F           = F0_building_block(x,y) + xi3 * F3_building_block(x,y)

    Stokes1     = xi1 * F11_building_block(x,y) / F
    Stokes2     = 0*Stokes1
    Stokes3     = (F3_building_block(x,y) + xi3*F33_building_block(x,y)) / F 

    return Stokes1,Stokes2,Stokes3





def radiation_spectrum_linear( U_in , a0 , omega0 , sigma , omega, theta, phi , a0_freq_correction=True ):
    """
    Use a linear Compton based model 

    U_in        : normalized four-momentum of the incoming electron
    a0          : normalized laser vector potential at peak of pulse
    omega0      : laser frequency
    sigma       : laser pulse duration
    omega       : outgoing photon frequency
    theta       : outgoing photon polar angle
    phi         : outgoing photon azimuthal angle
    a0_freq_correction : switch on the intensity-dependent red-shift  

    """

    # define the Baier-Katkov-type recoil factor E/(E-omega) = 1 /(1-omega/E)
    energy_in    = U_in[0]*elec_mass
    s            = omega/energy_in



    # direction of the outgoing photon
    nx = np.sin(theta)*np.cos(phi)
    ny = np.sin(theta)*np.sin(phi)
    nz = np.cos(theta)


    # quantum energy paramter == laser frequency in electron rest frame in units of electron rest mass
    b0          =     omega0 * ( U_in[0] + U_in[3]) / elec_mass


    # outgoing kappa-factor
    kap_out     = omega*( U_in[0]  - U_in[1] *nx - U_in[2] *ny - U_in[3] *nz ) / elec_mass
 
    if a0_freq_correction:
        # approximately take into account the intensity-dependent red-shift 
        L           = (kap_out + s*a0**2/2.) / b0 / (1-s)
    else:
        L           = kap_out / b0 / (1-s)


    # incoming kappa-factor, the L is due to the laser spectrum, == L * b0
    kap_in      = L * omega0 * ( U_in[0] + U_in[3]) / elec_mass


    X           = kap_in/kap_out + kap_out/kap_in + 2*(1/kap_in-1/kap_out) + (1/kap_in-1/kap_out)**2
    # -->  this needs to be modified for linear laser polarization!!!!




    prefactor   = finestruct * a0**2 * omega / (32*np.pi**2) / elec_mass**2 / b0**2 / (1-s)


    laser_spec  = energy_spectum(L,sigma)



    return prefactor * laser_spec * X * 2


