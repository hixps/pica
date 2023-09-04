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



class Compton_Spectrum():

    def __init__(self, U_in , a0 , omega0 , sigma , omega, theta, phi , laser_poldegree , laser_polangle , a0_freq_correction ):

        self.U_in = U_in
        self.a0   = a0
        self.omega0 = omega0
        self.sigma  = sigma

        self.omega  = omega
        self.theta  = theta
        self.phi    = phi

        
        self.laser_poldegree   = laser_poldegree
        self.laser_polangle    = laser_polangle


        self.a0_freq_correction = a0_freq_correction


        self.nx = np.sin(theta)*np.cos(phi)
        self.ny = np.sin(theta)*np.sin(phi)
        self.nz = np.cos(theta)

        
        self.energy_in    = self.U_in[0]*elec_mass
        self.s            = self.omega/self.energy_in   # fractional energy loss

        # self.s2           = self.omega * ( 1 + self.nz) /  (self.U_in[0] + self.U_in[3])/elec_mass

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
    def prefactor2(self):
        self._prefactor   = finestruct * self.a0**2 * self.omega / (4*np.pi**2) / elec_mass**2 / self.b0**2 / (1-self.s)
        return self._prefactor




class Compton_Spectrum_Full(Compton_Spectrum):

    def __init__(self, *args, **kwargs):
        Compton_Spectrum.__init__(self, *args, **kwargs)

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
        #needed if incident laser is not fully polarized
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
        return (E2_0,E2_1,E2_2,E2_3)  

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
        return self.prefactor * laser_spec * F


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



class Compton_Spectrum_Landau(Compton_Spectrum):

    def __init__(self, *args, **kwargs):
        Compton_Spectrum.__init__(self, *args, **kwargs)

        # only x-polarization should be somewhat ok ... others need to be checked
        self.S1 =   self.laser_poldegree * np.sin(2*self.laser_polangle)
        self.S3 = - self.laser_poldegree * np.cos(2*self.laser_polangle)

    def cross_section(self):

        F           = self.F0_building_block(self.x,self.y) + self.xi3 * self.F3_building_block(self.x,self.y)
        laser_spec  = laser_spectum(self.L, self.sigma)
        # prefactor   = finestruct * self.a0**2 * self.omega / (4*np.pi**2) / elec_mass**2 / self.b0**2 / (1-self.s)

        return self.prefactor * laser_spec * F


    def StokesParameters( self ):

        F           = self.F0_building_block(self.x,self.y) + self.xi3 * self.F3_building_block(self.x,self.y)

        Stokes1     = self.xi1 * self.F11_building_block(self.x,self.y) / F
        Stokes2     = 0*Stokes1
        Stokes3     = (self.F3_building_block(self.x,self.y) + self.xi3*self.F33_building_block(self.x,self.y)) / F 

        return Stokes1,Stokes2,Stokes3

    @property
    def xi1(self):
        return self.S3 * np.sin(2*self.phi)

    @property
    def xi3(self):
        # polarization degree with regard to scattering plane
        return self.S3 * np.cos(2*self.phi)


    def F0_building_block(self,x,y):
        return 0.25*(x/y + y/x) + (1./x - 1./y) + (1./x-1./y)**2

    def F3_building_block(self,x,y):
        return - (1./x - 1./y) - (1./x-1./y)**2

    def F11_building_block(self,x,y):
        return 1./x-1./y + 0.5

    def F22_building_block(self,x,y):
        return 0.25*(x/y+y/x)*(1+2/x-2/y)

    def F33_building_block(self,x,y):
        return (1/x-1/y)**2 + (1/x-1/y) + 0.5




class Compton_Spectrum_Greiner(Compton_Spectrum):

    def __init__(self, *args, **kwargs):
        Compton_Spectrum.__init__(self, *args, **kwargs)


    def cross_section(self):

        Epsilon_pi = - self.U_in[1]
        Epsilon_pf = - self.U_in[1] + self.omega/elec_mass * self.nx

        # print (Epsilon_pi)
        # print (Epsilon_pf)
        # print (Epsilon_pi/Epsilon_pf)

        Polarizationbracket = 2*(Epsilon_pi)/self.x - 2*Epsilon_pf/self.y

        # print (list(zip(self.phi/np.pi,Polarizationbracket))[:100])
        # print (0.25*(self.x/self.y + self.y/self.x))
        # print (0.5 * Polarizationbracket**2)

        F           = 0.25*(self.x/self.y+self.y/self.x) - 0.5 * Polarizationbracket**2


        # prefactor   = finestruct * self.a0**2 * self.omega / (4*np.pi**2) / elec_mass**2 / self.b0**2 / (1-self.s)
        laser_spec  = laser_spectum(self.L,self.sigma)
        return self.prefactor * laser_spec * F


    def StokesParameters( self ):

        # incident polarization vector
        Ein_0 = 0
        Ein_1 = 1
        Ein_2 = 0
        Ein_3 = 0

        # outgoing polarization vector 1
        E1_0   = 0
        E1_1   = np.cos(self.theta) * np.cos(self.phi)
        E1_2   = np.cos(self.theta) * np.sin(self.phi)
        E1_3   = -np.sin(self.theta) 

        # outgoing polarization vector 2
        E2_0   = 0
        E2_1   = -np.sin(self.phi)
        E2_2   = np.cos(self.phi)
        E2_3   = 0

        # outgoing normalized momentum
        U_out_0 = self.U_in[0] + ( self.L*self.omega0 - self.omega        )/elec_mass
        U_out_1 = self.U_in[1] + (                    - self.omega*self.nx)/elec_mass
        U_out_2 = self.U_in[2] + (                    - self.omega*self.ny)/elec_mass
        U_out_3 = self.U_in[3] + (-self.L*self.omega0 - self.omega*self.nz)/elec_mass


        # now all the scalar products

        Ein_E1   = Ein_0*E1_0 - Ein_1*E1_1 - Ein_2*E1_2 - Ein_3*E1_3
        Ein_E2   = Ein_0*E2_0 - Ein_1*E2_1 - Ein_2*E2_2 - Ein_3*E2_3

        Ein_Uin  = Ein_0*self.U_in[0] - Ein_1*self.U_in[1] - Ein_2*self.U_in[2] - Ein_3*self.U_in[3]
        Ein_Uout = Ein_0*U_out_0 - Ein_1*U_out_1 - Ein_2*U_out_2 - Ein_3*U_out_3

        E1_Uin   = E1_0*self.U_in[0] - E1_1*self.U_in[1] - E1_2*self.U_in[2] - E1_3*self.U_in[3]
        E1_Uout  = E1_0*U_out_0 - E1_1*U_out_1 - E1_2*U_out_2 - E1_3*U_out_3

        E2_Uin   = E2_0*self.U_in[0] - E2_1*self.U_in[1] - E2_2*self.U_in[2] - E2_3*self.U_in[3]
        E2_Uout  = E2_0*U_out_0 - E2_1*U_out_1 - E2_2*U_out_2 - E2_3*U_out_3

        M0       = 0.5*(self.x/self.y+self.y/self.x) - 1
        M1       = Ein_E1 - 2*Ein_Uin*E1_Uout/self.x + 2*Ein_Uout*E1_Uin/self.y
        M2       = Ein_E2 - 2*Ein_Uin*E2_Uout/self.x + 2*Ein_Uout*E2_Uin/self.y



        Polarizationbracket = 2*Ein_Uin/self.x - 2*Ein_Uout/self.y
        F           = 0.25*(self.x/self.y+self.y/self.x) - 0.5 * Polarizationbracket**2
        print (M0+M1**2+M2**2)
        print (2*F)


        rho_11 = M0 + 2 * M1*M1
        rho_12 =      2 * M1*M2
        rho_22 = M0 + 2 * M2*M2


        Stokes2     = 2*rho_12/(rho_11+rho_22)
        Stokes1     = (rho_11-rho_22)/(rho_11+rho_22)
        Stokes3     = 0*Stokes1


        return Stokes1,Stokes2,Stokes3



def Compton_spectrum_linear( U_in , a0 , omega0 , sigma , omega, theta, phi , S3 , a0_freq_correction=False ):
    energy_in    = U_in[0]*elec_mass
    s            = omega/energy_in   # fractional energy loss


    # print (s)
    # print ((U_in[0] + U_in[3])/(omega*(1+np.cos(theta))))

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


    # print (L)
    # print (x)
    # print (y)
    # print (x/y+y/x)
    # print ('----------')

    #X           = kap_in/kap_out + kap_out/kap_in + 2*(1/kap_in-1/kap_out) + (1/kap_in-1/kap_out)**2
    # -->  this needs to be modified for linear laser polarization!!!!

    # polarization degree with regard to scattering plane
    xi3         = S3 * np.cos(2*phi)


    # print (S3)
    # print (xi3)

    F           = F0_building_block(x,y) + xi3 * F3_building_block(x,y)

    # print (all(F>0))

    prefactor   = finestruct * a0**2 * omega / (32*np.pi**2) / elec_mass**2 / b0**2 / (1-s)


    laser_spec  = laser_spectum(L, sigma)



    return 8 * prefactor * laser_spec * F




def Compton_spectrum_linear_Greiner( U_in , a0 , omega0 , sigma , omega, theta, phi , S3 , a0_freq_correction=False ):
    energy_in    = U_in[0]*elec_mass
    s            = omega/energy_in   # fractional energy loss

    print ('xsection=Greiner')
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

    x          = 2 * L * omega0 * ( U_in[0] + U_in[3]) / elec_mass


    Epsilon_pi = - U_in[1]
    Epsilon_pf = - U_in[1] + omega/elec_mass * nx

    Polarizationbracket = 2*(Epsilon_pi)/x - 2*Epsilon_pf/y

    F          = 0.25*(x/y + y/x) - 0.5 * Polarizationbracket**2

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

    # print (b0)

    # LL y-parameter (86-15), == 2* outgoing kappa-factor of Jauch-Rohrlich (up to mass)
    # y = (m^2-u)/m^2 = 2 K_out*U_in / m 
    y          = 2 * omega*( U_in[0]  - U_in[1] *nx - U_in[2] *ny - U_in[3] *nz ) / elec_mass

    if a0_freq_correction:
        L = L_parameter_a0correction(y,s,b0,a0)
    else:
        L = L_parameter_linear(y,s,b0)


    # x          = 2 * L * omega0 * ( U_in[0] + U_in[3]) / elec_mass
    x          = 2 * omega0 * ( U_in[0] + U_in[3]) / elec_mass

    # print ('vvvvvvvvvv')
    # print (np.amin(L),np.amax(L))
    # print (x)
    # print (y)
    # print (x/y+y/x)
    # print ('^^^^^^^^^^')


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


