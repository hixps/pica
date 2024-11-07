import numpy as np
import yaml
import h5py

import dacite


from .constants import hbarc, c, elec_mass, finestruct, joule_per_eV, elementary_charge_si
from .spectrum  import Compton_Spectrum_Full as Compton_Spectrum
from .__init__  import __version__
from .auxiliary import beam_covariance_matrix, gaussian_sum_normalized
from .input     import PICA_Config, H5Reader
from .output    import H5Writer


"""
Routines for a collision of a beam with a laser pulse
"""


class ICSAnalysis():

    def __init__(self):
        pass

    @property
    def phi(self):
        _,K1,K2,_ = self.K_photon
        phi =  np.arctan2(K2,K1)
        return phi


    @property
    def S1_scatteringplane(self):
        _,Stokes1,Stokes2,_ = self.S_photon
        phi = self.phi
        Stokes_rot1 = np.cos(2*phi) * Stokes1 + np.sin(2*phi) * Stokes2
        return Stokes_rot1

    @property
    def S2_scatteringplane(self):
        _,Stokes1,Stokes2,_ = self.S_photon
        phi = self.phi
        Stokes_rot2 = - np.sin(2*phi) * Stokes1 + np.cos(2*phi) * Stokes2
        return Stokes_rot2

    
    @property
    def LinearPolarizationDegree(self):
        _,Stokes1,Stokes2,_ = self.S_photon
        PolDegree = np.sqrt(Stokes1**2 + Stokes2**2)
        return PolDegree

    @property
    def thetax(self):
        _,K1,_,K3 =  self.K_photon
        _thetax = np.arctan2(K1,K3)
        return _thetax

    @property
    def thetay(self):
        _,_,K2,K3 =  self.K_photon
        _thetay = np.arctan2(K2,K3)
        return _thetay  

    @property
    def theta(self):
        return np.sqrt( self.thetax**2 + self.thetay**2 )
    

    @property
    def xoffset(self):
        return self.thetax * self.config.beam.baseline

    @property
    def yoffset(self):
        return self.thetay * self.config.beam.baseline


class ICSSimulation(H5Writer, H5Reader, ICSAnalysis):

    def __init__(self, input_filename ):

        self.filename = input_filename

        with open( self.filename + '.yml', 'r' ) as stream:
            self.input_dict = yaml.load(stream, Loader=yaml.SafeLoader)

        self.config    = dacite.from_dict(data_class=PICA_Config, data=self.input_dict)



        self.Compton_Spectrum = Compton_Spectrum

        self.number_electrons = self.config.beam.charge / elementary_charge_si
        self.electron_weight  = self.number_electrons   / self.config.control.beam.sample_electrons


        # energy specific for cos^2 envelope
        self.total_energy    = 3/32 * self.config.laser.omega0 * ( 0.5*self.config.laser.a0**2 ) * elec_mass**2 / finestruct * self.config.laser.w0**2 * self.config.laser.sigma / hbarc**2
        self.total_energy_J  = self.total_energy * joule_per_eV




        print (f' >> pulse rescaling bias: {self.config.control.laser.pulse_rescale_bias:02.2f}:'
               f' sigma = {self.config.laser.sigma:.2f}    -> {self.config.laser.sigma/self.config.control.laser.pulse_rescale_bias:.2f}')
        print (f'                                :'
               f' TFWHM = {self.config.laser.TFWHM:.2f} fs -> {self.config.laser.TFWHM/self.config.control.laser.pulse_rescale_bias:.2f} fs')


        # initialize ouput arrays
        self.X_photon       = np.empty((4,0), dtype=float)
        self.K_photon       = np.empty((4,0), dtype=float)
        self.S_photon       = np.empty((4,0), dtype=float)
        self.W_photon       = np.empty(0, dtype=float)
        self.xi_peak        = np.empty(0, dtype=float)


        self.X_electron     = np.empty((4,0), dtype=float)
        self.P_electron     = np.empty((4,0), dtype=float)
        self.W_electron     = np.empty(0, dtype=float)


    def run(self):

        sampling_batch_sizes = self._get_batch_sizes()


        for j,n in enumerate( sampling_batch_sizes ):
            # iterate through batches
            self.MC_sampling_one_batch(j,n)

        self.write_file()




    def _get_batch_sizes(self):

        n_samp = self.config.control.beam.sample_electrons // self.config.control.beam.sample_batch_size
        r_samp = self.config.control.beam.sample_electrons % self.config.control.beam.sample_batch_size

        if r_samp:
            sampling_batch_sizes = [self.config.control.beam.sample_batch_size]*int(n_samp) + [r_samp]
        else:
            sampling_batch_sizes = [self.config.control.beam.sample_batch_size]*int(n_samp)


        print (f'n_batches={n_samp}, remainder={r_samp}, number_of_particles={sum(sampling_batch_sizes):.4g}={self.config.control.beam.sample_electrons:.4g}')
        print (f' >> number_electrons   = {self.number_electrons}')
        print (f' >> sample_electrons   = {self.config.control.beam.sample_electrons}')
        print (f' >> electron weight    = {self.electron_weight:.4g}')

        return sampling_batch_sizes



    def generate_electron_beam(self, number_sample_electrons):
        
        
        # rms angles for x and y space
        rms_angle_X   = self.config.beam.emittanceX / self.config.beam.sigmaX / self.config.beam.gamma
        rms_angle_Y   = self.config.beam.emittanceY / self.config.beam.sigmaY / self.config.beam.gamma
        
        
        #Mean and Covariant matrices for bivariate gaussian
        
        mean_x = [0,0]  # x-offset and x' offset for focus 
        cov_x = beam_covariance_matrix(self.config.beam.sigmaX, rms_angle_X, self.config.beam.focus_z )

        mean_y = [0,0]  # y-offset and y' offset for focus
        cov_y = beam_covariance_matrix(self.config.beam.sigmaY, rms_angle_Y, self.config.beam.focus_z )

        
        #Sampling (x,x') and (y,y')
        x,theta_x = np.random.multivariate_normal(mean_x, cov_x, number_sample_electrons).T
        y,theta_y = np.random.multivariate_normal(mean_y, cov_y, number_sample_electrons).T
        
        
        
        gamma        = np.random.normal( self.config.beam.gamma , self.config.beam.gamma*self.config.beam.energyspread , number_sample_electrons )
        beta         = np.sqrt(1-1./gamma**2)

        
        pz0          = gamma*beta*np.cos(theta_x)*np.cos(theta_y)
        px0          = gamma*beta*np.sin(theta_x)
        py0          = gamma*beta*np.sin(theta_y)
        pt0          = np.sqrt( 1 + px0**2 + py0**2 + pz0**2 )

        return x,theta_x,y,theta_y,gamma,pt0,px0,py0,pz0



    def MC_sampling_one_batch(self, jj, number_sample_electrons):



        print (f'  > batch {jj:03d}: {number_sample_electrons:.2g} macroelectrons of weight {self.electron_weight:.5g}')
        
        # define electron beam with correlated phase space
        x,theta_x,y,theta_y,gamma,pt0,px0,py0,pz0 = self.generate_electron_beam( number_sample_electrons )


        # radial position of electron at the interaction point        
        r            = np.sqrt( x**2 + y**2 )

        # peak value of laser intensity at electron radial location
        xi_peak      = self.config.laser.a0 * np.exp( -r**2/self.config.laser.w0**2 )


        U_in         = np.asarray([pt0,px0,py0,pz0])

        
        # random photon kinematics in detector range
        omega        = np.random.uniform( *self.config.detector.omega, number_sample_electrons ) 
        theta        = np.random.uniform( *self.config.detector.theta, number_sample_electrons ) 
        phi          = np.random.uniform( *self.config.detector.phi  , number_sample_electrons )
        # omega        = np.random.uniform(*self.detector.omega[0], self.omega_detector[1], number_sample_electrons) 
        # theta        = np.random.uniform(self.theta_detector[0], self.theta_detector[1], number_sample_electrons) 
        # phi          = np.random.uniform(self.phi_detector[0]  , self.phi_detector[1]  , number_sample_electrons)
        

        Spectrum_object = self.Compton_Spectrum( U_in , xi_peak , 
                                                    self.config.laser.omega0 , self.config.laser.sigma / self.config.control.laser.pulse_rescale_bias , 
                                                    omega, theta, phi, 
                                                    self.config.laser.poldegree, self.config.laser.polangle, 
                                                    self.config.control.detector.a0_freq_correction )

        rad_int         = Spectrum_object.cross_section()


        spectrum = self.config.control.laser.pulse_rescale_bias * rad_int * np.sin(theta)
        spec_max = np.amax(spectrum)



        s_val        = np.random.uniform(0             , spec_max,       number_sample_electrons )


        samplingvolume     = spec_max * (self.config.detector.omega[1]-self.config.detector.omega[0]) \
                                      * (self.config.detector.theta[1]-self.config.detector.theta[0]) \
                                      * (self.config.detector.phi[1]  -self.config.detector.phi[0])

        base_photon_weight = samplingvolume * self.electron_weight


        # selector1 are the particles that have been successfully sampled
        selector1          = s_val < spectrum
        number_photons     = sum(selector1)
        print ('   base photon weight :' , base_photon_weight)
        print ('   number photons     :' , number_photons)



        # sampled_gamma   = gamma[selector1]
        # sampled_theta_x = theta_x[selector1]
        # sampled_theta_y = theta_y[selector1]

        sampled_omega   = omega[selector1]
        sampled_theta   = theta[selector1]
        sampled_phi     = phi[selector1]

        sampled_x       = x[selector1]
        sampled_y       = y[selector1]

        sampled_zeta    = np.random.normal( 0 , self.config.beam.sigmaZ  , number_photons )
        sampled_z       = +0.5 * sampled_zeta
        sampled_t       = -0.5 * sampled_zeta


        sampled_U_in  = np.asarray([pt0[selector1],px0[selector1],py0[selector1],pz0[selector1]])


        sampled_xi_peak = xi_peak[selector1]


        # photon position
        X0           = sampled_t
        X1           = sampled_x
        X2           = sampled_y
        X3           = sampled_z

        # photon momentum
        K0           = sampled_omega
        K1           = sampled_omega * np.sin(sampled_theta) * np.cos(sampled_phi)
        K2           = sampled_omega * np.sin(sampled_theta) * np.sin(sampled_phi)
        K3           = sampled_omega * np.cos(sampled_theta) 


        # electrons that have emitted
        P1           = px0[selector1] * elec_mass - K1
        P2           = py0[selector1] * elec_mass - K2
        P3           = pz0[selector1] * elec_mass - K3 
        P0           = np.sqrt( elec_mass**2 + P1**2 + P2**2 + P3**2 )


        # Calculation of Stokes Parameters of emitted Photons
        sampled_Spectrum_object = self.Compton_Spectrum( sampled_U_in , sampled_xi_peak , 
                                                    self.config.laser.omega0 , self.config.laser.sigma / self.config.control.laser.pulse_rescale_bias , 
                                                    sampled_omega, sampled_theta, sampled_phi, 
                                                    self.config.laser.poldegree, self.config.laser.polangle, 
                                                    self.config.control.detector.a0_freq_correction )
        
        S1_sp, S2_sp, S3    = sampled_Spectrum_object.StokesParameters()

        S1 = np.cos(2*sampled_phi) * S1_sp - np.sin(2*sampled_phi) * S2_sp
        S2 = np.sin(2*sampled_phi) * S1_sp + np.cos(2*sampled_phi) * S2_sp

        S0             = np.ones(S1.shape, dtype=float)




        sampled_weights    = base_photon_weight * np.ones(number_photons)
        sampled_K_photon   = np.asarray( (K0,K1,K2,K3) )             / self.config.unit.momentum_scale
        sampled_X          = np.asarray( (X0,X1,X2,X3) )
        assigned_Stokes    = np.asarray( (S0,S1,S2,S3) )
        sampled_P_electron = np.asarray( (P0,P1,P2,P3) )             / self.config.unit.momentum_scale



        self.xi_peak      = np.concatenate( (self.xi_peak,    sampled_xi_peak    ) , axis=0)
        self.W_photon     = np.concatenate( (self.W_photon,   sampled_weights    ) , axis=0)
        self.W_electron   = np.concatenate( (self.W_electron, sampled_weights    ) , axis=0)

        self.X_photon     = np.concatenate( (self.X_photon,   sampled_X          ) , axis=1)
        self.K_photon     = np.concatenate( (self.K_photon,   sampled_K_photon   ) , axis=1)
        self.S_photon     = np.concatenate( (self.S_photon,   assigned_Stokes    ) , axis=1)

        self.X_electron   = np.concatenate( (self.X_electron, sampled_X          ) , axis=1)
        self.P_electron   = np.concatenate( (self.P_electron, sampled_P_electron ) , axis=1)


        print ('   total photon number:',len(self.W_electron))




#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################
#######################################################################################################################################


"""
class AtomicSimulation(H5Writer, H5Reader, ICSAnalysis):

    def __init__(self, input_filename ):

        print ("AtomicSimulation!")


        self.filename = input_filename

        with open( self.filename + '.yml', 'r' ) as stream:
            self.input_dict = yaml.load(stream, Loader=yaml.SafeLoader)


        self.targettype = self.input_dict['control']['beam']['targettype']
        print (f'Target type is <<{self.targettype}>>')
        
        if self.targettype=='orbital':
            from .orbitals import OrbitalSampler
            self.OrbitalSamplerObject = OrbitalSampler(self.input_dict)

        self.xsection = self.input_dict['control']['xsection']


        if self.xsection=='Full':
            self.Compton_Spectrum = Compton_Spectrum_Full
        else:
            raise ValueError

        self.read_laser_parameters()
        self.read_beam_parameters()
        self.read_detector()


        self.unit_X = self.input_dict['unit']['position']
        if self.unit_X not in ('micron',):
            raise ValueError

        self.unit_P = self.input_dict['unit']['momentum']
        if self.unit_P not in ('eV','keV','MeV','GeV'):
            raise ValueError

        if self.unit_P == 'eV':
            self.momentum_scale = 1
        elif self.unit_P == 'keV':
            self.momentum_scale = 1e3
        elif self.unit_P == 'MeV':
            self.momentum_scale = 1e6
        elif self.unit_P == 'GeV':
            self.momentum_scale = 1e9



        # translate TFWHM==FWHM in fs --> sigma parameter which is dimensionless, specific for cos^2 envelope
        self.sigma = 0.25*np.pi/np.arccos(1/2**0.25)*c/hbarc * self.TFWHM * self.omega0
        # print (0.25*np.pi/np.arccos(1/2**0.25)*c/hbarc , 2.0852201339)
        # self.sigma = 2.0852201339 * self.TFWHM * self.omega0
        # the numerical factor is 0.25*pi/arccos(1/2**(1/4)) * 1.52, where the factor 1.52 comes from the transition from eV to fs


        self.number_electrons =  self.beam_charge / 1.60217653e-19


        # energy specific for cos^2 envelope
        self.total_energy    = 3/32 * self.omega0 * ( 0.5*self.a0**2 ) * elec_mass**2 / finestruct * self.w0**2 * self.sigma / hbarc**2
        self.total_energy_J  = self.total_energy * joule_per_eV


        self.pulse_rescale_bias = float( self.input_dict['control']['laser']['pulse_rescale_bias']  )
        print (f' >> pulse rescaling bias: {self.pulse_rescale_bias:02.2f}: sigma = {self.sigma:.2f}     -> {self.sigma/self.pulse_rescale_bias:.2f}    ')
        print (f'                              : TFWHM = {self.TFWHM:.2f} fs -> {self.TFWHM/self.pulse_rescale_bias:.2f} fs')



        self.a0_freq_correction = bool(self.input_dict['control']['detector']['a0_freq_correction'])


        # initialize ouput arrays
        self.X_photon       = np.empty((4,0), dtype=float)
        self.K_photon       = np.empty((4,0), dtype=float)
        self.S_photon       = np.empty((4,0), dtype=float)
        self.W_photon       = np.empty(0, dtype=float)
        self.xi_peak        = np.empty(0, dtype=float)


        self.X_electron     = np.empty((4,0), dtype=float)
        self.P_electron     = np.empty((4,0), dtype=float)
        self.W_electron     = np.empty(0, dtype=float)


    def run(self):

        self.define_batches()

        for j,n in enumerate(self.sampling_batches):
            print (j,n)
            self.MC_sampling_one_batch(j,n)

        self.write_file()




    def define_batches(self):
        self.sample_electrons_total  = int(float( self.input_dict['control']['beam']['sample_electrons'] ))

        try:
            self.sample_batch_size   = int(float(self.input_dict['control']['beam']['sample_batch_size']))
        except KeyError:
            self.sample_batch_size   = int(1e7)
            print (f'"control/sample_batch_size" not specified, set to default value {self.sample_batch_size:.1g}')

        n_samp = self.sample_electrons_total // self.sample_batch_size
        r_samp = self.sample_electrons_total % self.sample_batch_size

        if r_samp:
            self.sampling_batches = [self.sample_batch_size]*int(n_samp) + [r_samp]
        else:
            self.sampling_batches = [self.sample_batch_size]*int(n_samp)

        print (f'n_batches={n_samp}, remainder={r_samp}, number_of_particles={sum(self.sampling_batches):.4g}={self.sample_electrons_total:.4g}')

        number_electrons_int  = int( self.number_electrons )

        # sample_electrons = int(5e7) # number of elecrons used for MC sampling, determines weight of each emitted photon
        self.electron_weight  = number_electrons_int / self.sample_electrons_total
        print (f' >> number_electrons   = {self.number_electrons:.1g} -> {number_electrons_int}')
        print (f' >> sample_electrons   = {self.sample_electrons_total}')
        print (f' >> electron weight    = {self.electron_weight:.4g}')





    def generate_electron_distribution(self, number_sample_electrons):
        

        if self.targettype=='orbital':
            p_au     = self.OrbitalSamplerObject.SampleFromCDF( number_sample_electrons )
            p        = p_au * finestruct            # from atomic units to natural units, p is normalized to mc


            cos_theta = np.random.uniform(-1.,1.,number_sample_electrons)
            phi       = np.random.uniform(0,2*np.pi,number_sample_electrons)

            theta     = np.arccos(cos_theta)

            x         = np.zeros(number_sample_electrons)
            y         = np.zeros(number_sample_electrons)

            pz0          = p * cos_theta
            px0          = p * np.sin(theta) * np.cos(phi)
            py0          = p * np.sin(theta) * np.sin(phi)
            pt0          = np.sqrt( 1 + p**2 )            
            # gamma        = pt0


        elif self.targettype=='beam':

            x         = np.zeros(number_sample_electrons)
            y         = np.zeros(number_sample_electrons)
            theta_x   = np.zeros(number_sample_electrons)
            theta_y   = np.zeros(number_sample_electrons)

            gamma        = self.gamma0 * np.ones( number_sample_electrons )
            p            = np.sqrt(gamma**2-1.)

            
            pz0          = p*np.cos(theta_x)*np.cos(theta_y)
            px0          = p*np.sin(theta_x)
            py0          = p*np.sin(theta_y)
            pt0          = np.sqrt( 1 + p**2 )


        return x,y,pt0,px0,py0,pz0



    def MC_sampling_one_batch(self, jj, number_sample_electrons):



        print (f'  > batch {jj:03d}: {number_sample_electrons:.2g} macroelectrons of weight {self.electron_weight:.5g}')
        
        # define electron beam with correlated phase space
        # x,theta_x,y,theta_y,gamma,pt0,px0,py0,pz0 = self.generate_electron_distribution( number_sample_electrons )
        x,y,pt0,px0,py0,pz0 = self.generate_electron_distribution( number_sample_electrons )

        # radial position of electron at the interaction point        
        r            = np.sqrt( x**2 + y**2 )

        # peak value of laser intensity at electron radial location
        xi_peak      = self.a0 * np.exp( -r**2/self.w0**2 )


        U_in         = np.asarray([pt0,px0,py0,pz0])

        
        # random photon kinematics in detector range
        omega        = np.random.uniform(self.omega_detector[0], self.omega_detector[1], number_sample_electrons) 
        theta        = np.random.uniform(self.theta_detector[0], self.theta_detector[1], number_sample_electrons) 
        phi          = np.random.uniform(self.phi_detector[0]  , self.phi_detector[1]  , number_sample_electrons)
        

        Spectrum_object = self.Compton_Spectrum( U_in , xi_peak , self.omega0 , self.sigma / self.pulse_rescale_bias , 
                                                    omega, theta, phi, 
                                                    self.poldegree, self.polangle, self.a0_freq_correction )
        rad_int         = Spectrum_object.cross_section()


        spectrum = self.pulse_rescale_bias * rad_int * np.sin(theta)
        spec_max = np.amax(spectrum)



        s_val        = np.random.uniform(0             , spec_max,       number_sample_electrons )


        samplingvolume     = spec_max * (self.omega_detector[1]-self.omega_detector[0]) \
                                      * (self.theta_detector[1]-self.theta_detector[0]) \
                                      * (self.phi_detector[1]  -self.phi_detector[0])

        base_photon_weight = samplingvolume * self.electron_weight


        # selector1 are the particles that have been successfully sampled
        selector1          = s_val < spectrum
        number_photons     = sum(selector1)
        print ('   base photon weight :' , base_photon_weight)
        print ('   number photons     :' , number_photons)



        sampled_omega   = omega[selector1]
        sampled_theta   = theta[selector1]
        sampled_phi     = phi[selector1]

        sampled_x       = x[selector1]
        sampled_y       = y[selector1]

        sampled_zeta    = np.random.normal( 0 , self.beam_length  , number_photons )
        sampled_z       = +0.5 * sampled_zeta
        sampled_t       = -0.5 * sampled_zeta


        sampled_U_in  = np.asarray([pt0[selector1],px0[selector1],py0[selector1],pz0[selector1]])


        sampled_xi_peak = xi_peak[selector1]


        # photon position
        X0           = sampled_t
        X1           = sampled_x
        X2           = sampled_y
        X3           = sampled_z

        # photon momentum
        K0           = sampled_omega
        K1           = sampled_omega * np.sin(sampled_theta) * np.cos(sampled_phi)
        K2           = sampled_omega * np.sin(sampled_theta) * np.sin(sampled_phi)
        K3           = sampled_omega * np.cos(sampled_theta) 


        # electrons that have emitted
        P1           = px0[selector1] * elec_mass - K1
        P2           = py0[selector1] * elec_mass - K2
        P3           = pz0[selector1] * elec_mass - K3 
        P0           = np.sqrt( elec_mass**2 + P1**2 + P2**2 + P3**2 )


        # Calculation of Stokes Parameters of emitted Photons
        sampled_Spectrum_object = self.Compton_Spectrum( sampled_U_in , sampled_xi_peak , self.omega0 , self.sigma / self.pulse_rescale_bias , 
                                                    sampled_omega, sampled_theta, sampled_phi, 
                                                    self.poldegree, self.polangle, self.a0_freq_correction )
        
        S1_sp, S2_sp, S3    = sampled_Spectrum_object.StokesParameters()

        S1 = np.cos(2*sampled_phi) * S1_sp - np.sin(2*sampled_phi) * S2_sp
        S2 = np.sin(2*sampled_phi) * S1_sp + np.cos(2*sampled_phi) * S2_sp

        S0             = np.ones(S1.shape, dtype=float)




        sampled_weights    = base_photon_weight * np.ones(number_photons)
        sampled_K_photon   = np.asarray( (K0,K1,K2,K3) )             / self.momentum_scale
        sampled_X          = np.asarray( (X0,X1,X2,X3) )
        assigned_Stokes    = np.asarray( (S0,S1,S2,S3) )
        sampled_P_electron = np.asarray( (P0,P1,P2,P3) )             / self.momentum_scale



        self.xi_peak      = np.concatenate( (self.xi_peak,    sampled_xi_peak    ) , axis=0)
        self.W_photon     = np.concatenate( (self.W_photon,   sampled_weights    ) , axis=0)
        self.W_electron   = np.concatenate( (self.W_electron, sampled_weights    ) , axis=0)

        self.X_photon     = np.concatenate( (self.X_photon,   sampled_X          ) , axis=1)
        self.K_photon     = np.concatenate( (self.K_photon,   sampled_K_photon   ) , axis=1)
        self.S_photon     = np.concatenate( (self.S_photon,   assigned_Stokes    ) , axis=1)

        self.X_electron   = np.concatenate( (self.X_electron, sampled_X          ) , axis=1)
        self.P_electron   = np.concatenate( (self.P_electron, sampled_P_electron ) , axis=1)


        print ('   total photon number:',len(self.W_electron))
"""


