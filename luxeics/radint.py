import numpy as np
# from numpy import *
# from scipy import *
from scipy.integrate import *

# from .filon import filon_volkov_1,filon_volkov_2
from .constants import elec_mass,finestruct

"""
Routines for the radiation integral
"""


# def radint( t , func , phase , phasederiv_start , phasederiv_end ):
    
#     """
#     Evaluate the radiation integral by explicitly treating the boundary terms according to 

#     \int_{-\infty}^\infty dt func(t) * exp( 1j*phase(t) ) 
#         = 
#         \int_{t_i}^{t_f} dt func(t) * exp( 1j*phase(t) ) 
#         + 1j * func(t) * exp(1j*phase(t))/ phase'(t) |_{t_i}^{t_f}

#     t:     1d array of timesteps
#     func:  1d array of pre-exponential fucntion values
#     phase: nd array of phase function values, time must be first axis. Frequency/angle axes to follow
#     phasederiv_start , phasederiv_end: values of the derivative of the phase at the start and endpoints of t

#     The time sequence needs to have an odd number of nodes for the filon quadrature
#     """

#     ts , fs , gs = t.shape , func.shape , phase.shape

#     if not ts[0] == fs[0] == gs[0]:
#         raise IndexError('t,func,phase must have same number of timesteps')
#     if ts[0] % 2 == 0:
#         raise IndexError('There must be an odd number of timesteps')

#     pdim = len(phase.shape)-1

#     # filon quadrature of the bulk of the integral from t_i to t_f
#     if pdim==1:
#         bulk          = filon_volkov_1( t , func , phase )
#     elif pdim==2:
#         bulk          = filon_volkov_2( t , func , phase )


#     surface_start = 1j * func[0]  * np.exp(1j*phase[0] ) / phasederiv_start
#     surface_end   = 1j * func[-1] * np.exp(1j*phase[-1]) / phasederiv_end

#     return bulk + surface_end - surface_start


# def calc_phase_1( track , baier_omega , theta , phi ):
#     t,x,y,z,ut,ux,uy,uz = track.T



#     n_times_x =  t -  x*nx -  y*ny -  z*nz
#     phase     = baier_omega[np.newaxis,:] * n_times_x[:,np.newaxis]

    
#     k_times_u_start = baier_omega*( ut[0]  - ux[0] *nx - uy[0] *ny - uz[0] *nz )
#     k_times_u_end   = baier_omega*( ut[-1] - ux[-1]*nx - uy[-1]*ny - uz[-1]*nz )

#     return phase , k_times_u_start , k_times_u_end


# def calc_phase_2theta(track , baier_omega , theta , phi ):
#     t,x,y,z,ut,ux,uy,uz = track.T


#     nx = np.sin(theta)*np.cos(phi)
#     ny = np.sin(theta)*np.sin(phi)
#     nz = np.cos(theta)


#     n_times_x =  t[:,np.newaxis] -  x[:,np.newaxis]*nx -  y[:,np.newaxis]*ny -  z[:,np.newaxis]*nz

#     # print (n_times_x.shape)

#     phase     = baier_omega[np.newaxis,:] * n_times_x[:,np.newaxis]

#     # print (phase.shape)
    
#     k_times_u_start = baier_omega*( ut[0]  - ux[0] *nx - uy[0] *ny - uz[0] *nz )
#     k_times_u_end   = baier_omega*( ut[-1] - ux[-1]*nx - uy[-1]*ny - uz[-1]*nz )

#     print (k_times_u_start.shape,k_times_u_end.shape)


#     return phase , k_times_u_start , k_times_u_end



# def quasiclass_radiation_vector( track , tau , omega , theta , phi ):
#     """
#     Evaluate the four vector integral

#     I_mu = \int d tau u_mu(tau) exp( i k_nu * x^nu(tau) )

#     tau ... proper time
#     u   ... four-velocity dx/dtau
#     k   ... four-momentum of the emitted photon, parametrized by (omega,theta,phi)

#     parameters:
#     ===========

#     track ... 2d array representing the orbit of a single particle, shape (ntau,8)
#               the 2nd dimension represents the phase space point (t,x,y,z,ut=gamma,ux,uy,uz)
#     tau   ... 1d array of proper time value, len = ntau

#     omega ... 1d array of frequencies
#     theta ... scalar, polar emission angle
#     phi   ... scalar, azimuthal emission angle


#     """

#     # for the filon quadrature we have to make sure that ntau is odd, if it is even remove the last timestep from tau and tracks
#     ntau=len(tau)
#     if ntau % 2 == 0:
#         tau   = tau[:-1]
#         ntau  = len(tau)
#         track = track[:-1,:]


#     # phase space checks

#     o      = omega*theta*phi
#     pshape = o.shape
#     pdim   = len(pshape)




#     # unpack the track
#     t,x,y,z,ut,ux,uy,uz = track.T


#     # define the Baier-Katkov recoil factor E/(E-omega) = 1 /(1-omega/E)
#     energy_in    = ut[0]*elec_mass
#     s            = omega/energy_in
#     baier_factor = 1./(1.-s)
#     baier_omega  = baier_factor*omega


#     # here it becomes interesting...need to allow multidimensional theta and phi arrays. maye just define the appropriate cases ?
#     if pdim==1:
#         # only omega has a dimension
#         phase , dphase_start , dphase_end = calc_phase_1(track,baier_omega,theta,phi)

#     elif pdim==2:
#         # assume that theta has the 2nd dimension
#         if pdim != len(theta.shape):
#             raise IndexError
#         print (pshape,pdim)
#         phase , dphase_start , dphase_end = calc_phase_2theta(track,baier_omega,theta,phi)





#     # radvec = [ radint( tau , ui , phase , dphase_start , dphase_end  ) for ui in (ut,ux,uy,uz) ] # as list
#     radvec = np.asarray( [ radint( tau , ui , phase , dphase_start , dphase_end  ) for ui in (ut,ux,uy,uz) ] ) # as array

#     # print (radvec.shape)
#     # returns shape of (4,pshape)
#     return radvec





# def radiation_spectrum(track, tau, omega, theta, phi ):

#     gamma0            = track[0,4]

#     # cartesian components
#     radvec            = quasiclass_radiation_vector( track , tau , omega , theta=theta , phi=phi )

#     squared_radvec    = abs(radvec[3])**2 + abs(radvec[2])**2 + abs(radvec[1])**2 - abs(radvec[0])**2 
#     radvec0           = abs(radvec[0])**2 

#     energy_in         = elec_mass * gamma0
#     energy_out        = energy_in - omega

#     prefactor         = omega * finestruct / (4*np.pi**2)


#     c1_BK             = 0.5*(energy_out**2 + energy_in**2)/energy_out**2
#     c2_BK             = elec_mass**2 * omega**2 / energy_in**2 / energy_out**2


#     quasiclass_spec   = prefactor * (c1_BK * squared_radvec + c2_BK*radvec0 )

#     print (quasiclass_spec.shape)
#     return quasiclass_spec



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


    # incoming kappa-factor, the L is due to the laser spectrum
    kap_in      = L * omega0 * ( U_in[0] + U_in[3]) / elec_mass


    #
    X           = kap_in/kap_out + kap_out/kap_in + 2*(1/kap_in-1/kap_out) + (1/kap_in-1/kap_out)**2


    prefactor   = finestruct *a0**2 * omega / (32*np.pi**2) / elec_mass**2 / b0**2 / (1-s)


    # specify spectrum for cos^2 pulse
    laser_spec  = np.pi**4 * np.sin( (L-1) * sigma )**2 / ( (L-1)*np.pi**2 - (L-1)**3*sigma**2)**2



    return prefactor * laser_spec * X * 2



if __name__ == '__main__':

    from constants import elec_mass,finestruct


    # # config parameters
    # omega0 = 1.55*3
    # a0     = 0.5
    # pplus  = 2*16.5e9/elec_mass

    # sigma = 2*pi*25/0.65

    # # define photon phase space
    # omega_ref = omega0*pplus**2 / ( 1  + 2*omega0*pplus/elec_mass)

    # omega = linspace(0.8*omega_ref , 1.05*omega_ref ,500)
    # theta = linspace(0,5e-6,6)

    # # laser phase
    # phi   = linspace(-sigma,sigma,5000)

    # # proper time
    # tau   = phi / omega0 / pplus

    # # laser envelope
    # gL    = cos( pi/2. * phi/sigma)**2 * (abs(phi)<= sigma)

    # # trajectory
    # px = a0 * sin(phi) * gL
    # py = a0 * cos(phi) * gL

    # # plot(phi,px)
    # # input()

    # pminus = ( 1 + px**2 + py**2 ) / pplus

    # pz    = (pplus - pminus)/2.
    # gamma = (pplus + pminus)/2. 


    # t     = cumtrapz( gamma, x=tau, initial=0 )
    # x     = cumtrapz( px   , x=tau, initial=0 )
    # y     = cumtrapz( py   , x=tau, initial=0 )
    # z     = cumtrapz( pz   , x=tau, initial=0 )


    # # particle orbit into single array, shape: (ntau, 8)
    # track = np.stack( [t,x,y,z,gamma,px,py,pz] , axis = -1)


    # spectrum = radiation_spectrum(track , tau , omega[:,newaxis] , theta[newaxis,:] , phi , elec_mass)
    # # print(spectrum)
    # import pylab
    # pylab.ion()
    # pylab.plot(omega,spectrum)
    # input()


