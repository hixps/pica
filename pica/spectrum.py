import numpy as np

from .constants import elec_mass,finestruct

"""
Routines for the radiation integral
"""


def laser_spectum(L,sigma):
    # the Fourier transform of the laser envelope, replaces the delta-function in the usual cross section

    # specify spectrum for cos^2 pulse
    laser_spec  = np.pi**4 * np.sin( (L-1) * sigma )**2 / ( (L-1)*np.pi**2 - (L-1)**3*sigma**2 )**2

    return laser_spec


def L_parameter_linear(y,s,b0):
    """
    Determine the value of L for energy spectrum
    use energy momentum conservation and that 
        1) in linear Compton scattering: y = 2b0*L(1-s)
    """
    # print ('no a0 correction')
    _L           = 0.5*y / b0/(1-s)
    return _L


def L_parameter_a0correction(y,s,b0,a0):
    # approximately take into account the intensity-dependent red-shift 
    # print ('a0 correction active')
    _L           = 0.5* (y + s*a0**2/2.*b0) / b0/(1-s)
    return _L



class _Compton_Spectrum():

    def __init__(self, U_in , a0 , omega, theta, phi , omega0 , sigma, laser_poldegree , laser_polangle , a0_freq_correction ):

        self.U_in = U_in
        self.a0   = a0

        self.omega  = omega
        self.theta  = theta
        self.phi    = phi


        self.omega0 = omega0
        self.sigma  = sigma
        self.laser_poldegree   = laser_poldegree
        self.laser_polangle    = laser_polangle
        self.a0_freq_correction = a0_freq_correction


        self.nx = np.sin(theta)*np.cos(phi)
        self.ny = np.sin(theta)*np.sin(phi)
        self.nz = np.cos(theta)

        
        self.energy_in    = self.U_in[0]*elec_mass
        # self.s            = self.omega/self.energy_in   # fractional energy loss

        self.s           = self.omega * ( 1 + self.nz) /  (self.U_in[0] + self.U_in[3])/elec_mass

        # quantum energy paramter == laser frequency in electron rest frame in units of electron rest mass
        self.b0           =     omega0 * ( U_in[0] + U_in[3]) / elec_mass



    @property
    def y(self):
        self._y = 2 * self.omega*( self.U_in[0]  - self.U_in[1] *self.nx 
                                                 - self.U_in[2] *self.ny 
                                                 - self.U_in[3] *self.nz ) / elec_mass
        return self._y

    
    @property
    def x(self):
        # LL x-parameter (86-15), == 2* incoming kappa-factor of Jauch-Rohrlich (up to mass)
        # x = (s-m^2)/m^2 = 2 L K_in*U_in / m
        # the L is due to the laser spectrum, == L * b0
        self._x = 2 * self.L * self.omega0 * ( self.U_in[0] + self.U_in[3]) / elec_mass
        return self._x
    
    @property
    def L(self):
        if self.a0_freq_correction:
            self._L = L_parameter_a0correction(self.y,self.s,self.b0,self.a0)
        else:
            self._L = L_parameter_linear(self.y,self.s,self.b0)
        return self._L


    @property
    def prefactor(self):
        self._prefactor   = finestruct * self.a0**2 * self.omega / (4*np.pi**2) / elec_mass**2 / self.b0**2 / (1-self.s)
        return self._prefactor

    @property
    def prefactor_ptarmigan_benchmark(self):
        self._prefactor   = 0.5 * finestruct * self.a0**2 * self.omega / (4*np.pi**2) / elec_mass**2 / self.b0**2 / (1-self.s)
        return self._prefactor




class Compton_Spectrum_Full(_Compton_Spectrum):

    def __init__(self, *args, **kwargs):
        _Compton_Spectrum.__init__(self, *args, **kwargs)

    @property
    def Ein(self):
        # incident polarization vector
        Ein_0 = 0
        Ein_1 = np.cos(self.laser_polangle)
        Ein_2 = np.sin(self.laser_polangle)
        Ein_3 = 0
        return (Ein_0, Ein_1, Ein_2, Ein_3)

    @property
    def Ein_ortho(self):
        # needed if incident laser is not fully polarized
        Ein_0 = 0
        Ein_1 = np.sin(self.laser_polangle)
        Ein_2 = -np.cos(self.laser_polangle)
        Ein_3 = 0
        return (Ein_0, Ein_1, Ein_2, Ein_3)
    


    @property
    def Eout_1(self):
        # outgoing polarization vector 1
        E1_0   = 0
        E1_1   = np.cos(self.theta) * np.cos(self.phi)
        E1_2   = np.cos(self.theta) * np.sin(self.phi)
        E1_3   = -np.sin(self.theta) 
        return (E1_0, E1_1, E1_2, E1_3)

    @property
    def Eout_2(self):
        # outgoing polarization vector 2
        E2_0   = 0
        E2_1   = -np.sin(self.phi)
        E2_2   = np.cos(self.phi)
        E2_3   = 0  
        return (E2_0, E2_1, E2_2, E2_3)  

    @property
    def U_out(self):
        # outgoing normalized momentum
        U_out_0 = self.U_in[0] + ( self.L*self.omega0 - self.omega        )/elec_mass
        U_out_1 = self.U_in[1] + (                    - self.omega*self.nx)/elec_mass
        U_out_2 = self.U_in[2] + (                    - self.omega*self.ny)/elec_mass
        U_out_3 = self.U_in[3] + (-self.L*self.omega0 - self.omega*self.nz)/elec_mass
        return (U_out_0,U_out_1,U_out_2,U_out_3)
    


    def dot4(self, X, Y ):
        return X[0]*Y[0] - X[1]*Y[1] - X[2]*Y[2] - X[3]*Y[3]

    @property
    def M0(self):
        return 0.5*(self.x/self.y+self.y/self.x) - 1

    # @property
    def M1(self,E_initial):
        Ein_E1   = self.dot4( E_initial, self.Eout_1)
        Ein_Uin  = self.dot4( E_initial, self.U_in )
        Ein_Uout = self.dot4( E_initial, self.U_out )
        E1_Uin   = self.dot4( self.Eout_1, self.U_in )
        E1_Uout  = self.dot4( self.Eout_1, self.U_out )
        return Ein_E1 - 2*Ein_Uin*E1_Uout/self.x + 2*Ein_Uout*E1_Uin/self.y

    # @property
    def M2(self,E_initial):
        Ein_E2   = self.dot4( E_initial, self.Eout_2)
        Ein_Uin  = self.dot4( E_initial, self.U_in )
        Ein_Uout = self.dot4( E_initial, self.U_out )
        E2_Uin   = self.dot4( self.Eout_2, self.U_in )
        E2_Uout  = self.dot4( self.Eout_2, self.U_out )
        return Ein_E2 - 2*Ein_Uin*E2_Uout/self.x + 2*Ein_Uout*E2_Uin/self.y



    def cross_section(self):

        m0     = self.M0
        m1     = self.M1(self.Ein)
        m2     = self.M2(self.Ein)
        m1o    = self.M1(self.Ein_ortho)
        m2o    = self.M2(self.Ein_ortho)

        F_in          = 0.5*(m0 +  m1**2  + m2**2 )
        F_in_ortho    = 0.5*(m0 +  m1o**2 + m2o**2)

        F             = 0.5*(1+self.laser_poldegree)*F_in + 0.5*(1-self.laser_poldegree)*F_in_ortho

        laser_spec  = laser_spectum(self.L,self.sigma)
        # return self.prefactor * laser_spec * F
        return self.prefactor_ptarmigan_benchmark * laser_spec * F


    def StokesParameters( self ):

        m0     = self.M0
        m1     = self.M1(self.Ein)
        m2     = self.M2(self.Ein)
        m1o    = self.M1(self.Ein_ortho)
        m2o    = self.M2(self.Ein_ortho)

        rho_11 = m0 + 2*m1**2
        rho_12 =      2*m1*m2
        rho_22 = m0 + 2*m2**2


        rho_11o = m0 + 2*m1o**2
        rho_12o =      2*m1o*m2o
        rho_22o = m0 + 2*m2o**2

        w       = 0.5*(1+self.laser_poldegree)
        wo      = 0.5*(1-self.laser_poldegree)


        Stokes1     = (w*(rho_11-rho_22)+wo*(rho_11o-rho_22o))/(w*(rho_11+rho_22) + wo*(rho_11o+rho_22o) ) # 0-90 degree polarization
        Stokes2     = 2*(w*rho_12+wo*rho_12o)/(w*(rho_11+rho_22) + wo*(rho_11o+rho_22o) )                  # diagonal polarization
        Stokes3     = 0*Stokes1                     # circular polarization degree, identically 0 for LP laser, not supported yet

        return Stokes1,Stokes2,Stokes3



