import numpy as np
import yaml
import h5py

from scipy.integrate import cumtrapz

from .constants import hbarc,c, elec_mass,finestruct, joule_per_eV
from .spectrum  import Compton_Spectrum_Greiner, Compton_Spectrum_Landau, Compton_Spectrum_Full
from .__init__  import __version__
from .auxiliary import beam_covariance_matrix, gaussian_sum_normalized
from .inout     import H5Writer, ParameterReader, H5Reader



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


    # rotated Stokes parameters in lab frame
    @property
    def S1_labframe(self):
        _,Stokes1,Stokes2,_ = self.S_photon
        phi = self.phi
        Stokes_rot1 = np.cos(2*phi) * Stokes1 - np.sin(2*phi) * Stokes2
        return Stokes_rot1

    @property
    def S2_labframe(self):
        _,Stokes1,Stokes2,_ = self.S_photon
        phi = self.phi
        Stokes_rot2 = np.sin(2*phi) * Stokes1 + np.cos(2*phi) * Stokes2
        return Stokes_rot2

    
    @property
    def LinearPolarizationDegree(self):
        _,Stokes1,Stokes2,_ = self.S_photon
        PolDegree = np.sqrt(Stokes1**2 + Stokes1**2 )
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
        return self.thetax * self.baseline

    @property
    def yoffset(self):
        return self.thetay * self.baseline



class ICSSimulation(H5Writer, ParameterReader, H5Reader, ICSAnalysis):

    def __init__(self, input_filename ):

        self.filename = input_filename

        with open( self.filename + '.yml', 'r' ) as stream:
            self.input_dict = yaml.load(stream, Loader=yaml.SafeLoader)


        self.xsection = self.input_dict['control']['xsection']

        if self.xsection=='Landau':
            self.Compton_Spectrum = Compton_Spectrum_Landau
        elif self.xsection=='Greiner':
            self.Compton_Spectrum = Compton_Spectrum_Greiner
        elif self.xsection=='Full':
            self.Compton_Spectrum = Compton_Spectrum_Full


        self.read_laser_parameters()
        self.read_beam_parameters()
        self.read_detector()


        # translate Tpulse==FWHM in fs --> sigma parameter which is dimensionless, specific for cos^2 envelope
        self.sigma = 0.25*np.pi/np.arccos(1/2**0.25)*c/hbarc * self.Tpulse * self.omega0
        # print (0.25*np.pi/np.arccos(1/2**0.25)*c/hbarc , 2.0852201339)
        # self.sigma = 2.0852201339 * self.Tpulse * self.omega0
        # the numerical factor is 0.25*pi/arccos(1/2**(1/4)) * 1.52, where the factor 1.52 comes from the transition from eV to fs


        self.total_energy    = 3/32 * self.omega0 * self.a0**2 * elec_mass**2 / finestruct * self.w0**2 * self.sigma / hbarc**2
        self.total_energy_J  = self.total_energy * joule_per_eV

        if self.sigma_rescale:
            self.sigma_crit =  float( self.input_dict['control']['laser']['sigma_crit']  )
            if self.sigma > self.sigma_crit:
                self.sigma_rescalefactor = self.sigma / self.sigma_crit
                print (f' >> rescale pulse duration: {self.sigma:.2f} -> {self.sigma_crit:.2f}: sigma_rescalefactor = {self.sigma_rescalefactor:.2f}')
            else:
                self.sigma_rescalefactor = 1.
        else:
            # self.sigma_crit          = 0
            self.sigma_rescalefactor = 1.


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
            # print (self.W_photon.shape)
            # print (self.K_photon.shape)
            self.MC_sampling_one_batch(j,n)

        self.write_file()




    def define_batches(self):
        self.sample_electrons_total  = int(float( self.input_dict['control']['beam']['sample_electrons'] ))

        try:
            self.sample_batch_size   = int(float(self.input_dict['control']['beam']['sample_batch_size']))
        except KeyError:
            self.sample_batch_size   = int(1e7)
            print (f'"control/sample_batch_size" not specified, set to default value {self.sample_batch_size:.1g}')


        # beam_charge  = float( self.input_dict['beam']['beam_charge'])



        n_samp = self.sample_electrons_total // self.sample_batch_size
        r_samp = self.sample_electrons_total % self.sample_batch_size

        if r_samp:
            self.sampling_batches = [self.sample_batch_size]*int(n_samp) + [r_samp]
        else:
            self.sampling_batches = [self.sample_batch_size]*int(n_samp)

        print (f'n_batches={n_samp}, remainder={r_samp}, number_of_particles={sum(self.sampling_batches):.4g}={self.sample_electrons_total:.4g}')

        number_electrons = int( self.beam_charge / 1.60217653e-19)

        # sample_electrons = int(5e7) # number of elecrons used for MC sampling, determines weight of each emitted photon
        self.electron_weight  = number_electrons / self.sample_electrons_total
        print (f' >> number_electrons   = {number_electrons}')
        print (f' >> sample_electrons   = {self.sample_electrons_total}')
        print (f' >> electron weight    = {self.electron_weight:.4g}')





    def generate_electron_beam(self, number_sample_electrons):
        
        
        # rms angles for x and y space
        rms_angle_X   = self.emittance_X / self.beam_size_X / self.gamma0
        rms_angle_Y   = self.emittance_Y / self.beam_size_Y / self.gamma0
        
        
        #Mean and Covariant matrices for bivariate gaussian
        
        mean_x = [0,0]  # x-offset and x' offset for focus 
        cov_x = beam_covariance_matrix(self.beam_size_X, rms_angle_X, self.beam_focus_z )

        mean_y = [0,0]  # y-offset and y' offset for focus
        cov_y = beam_covariance_matrix(self.beam_size_Y, rms_angle_Y, self.beam_focus_z )

        
        #Sampling (x,x') and (y,y')
        x,theta_x = np.random.multivariate_normal(mean_x, cov_x, number_sample_electrons).T
        y,theta_y = np.random.multivariate_normal(mean_y, cov_y, number_sample_electrons).T
        
        
        
        gamma        = np.random.normal( self.gamma0 , self.gamma0*self.energyspread , number_sample_electrons )
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
        xi_peak      = self.a0 * np.exp( -r**2/self.w0**2 )


        U_in         = np.asarray([pt0,px0,py0,pz0])

        
        # photon detector
        omega        = np.random.uniform(self.omega_detector[0], self.omega_detector[1], number_sample_electrons) 
        theta        = np.random.uniform(self.theta_detector[0], self.theta_detector[1], number_sample_electrons) 
        phi          = np.random.uniform(self.phi_detector[0]  , self.phi_detector[1]  , number_sample_electrons)
        

        Spectrum_object = self.Compton_Spectrum( U_in , xi_peak , self.omega0 , self.sigma / self.sigma_rescalefactor , 
                                                    omega, theta, phi, 
                                                    self.poldegree, self.polangle, self.a0_freq_correction )
        rad_int         = Spectrum_object.cross_section()


        spectrum = self.sigma_rescalefactor * rad_int * np.sin(theta)
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



        sampled_gamma   = gamma[selector1]
        sampled_theta_x = theta_x[selector1]
        sampled_theta_y = theta_y[selector1]

        sampled_omega   = omega[selector1]
        sampled_theta   = theta[selector1]
        sampled_phi     = phi[selector1]

        sampled_x       = x[selector1]
        sampled_y       = y[selector1]

        sampled_zeta    = np.random.normal( 0 , self.beam_length  , number_photons )
        sampled_z       = +0.5 * sampled_zeta
        sampled_t       = -0.5 * sampled_zeta


        sampled_xi_peak = xi_peak[selector1]


        K0           = sampled_omega
        K1           = sampled_omega * np.sin(sampled_theta) * np.cos(sampled_phi)
        K2           = sampled_omega * np.sin(sampled_theta) * np.sin(sampled_phi)
        K3           = sampled_omega * np.cos(sampled_theta) 


        # electrons that have emitted
        px_out  = px0[selector1] * elec_mass - K1
        py_out  = py0[selector1] * elec_mass - K2
        pz_out  = pz0[selector1] * elec_mass - K3 
        pt_out  = np.sqrt( elec_mass**2 + px_out**2 + py_out**2 + pz_out**2 )


        # Calculation of Stokes Parameters of emitted Photons
        sampled_U_in  = np.asarray([pt0[selector1],px0[selector1],py0[selector1],pz0[selector1]])

        sampled_Spectrum_object = self.Compton_Spectrum( sampled_U_in , sampled_xi_peak , self.omega0 , self.sigma / self.sigma_rescalefactor , 
                                                    sampled_omega, sampled_theta, sampled_phi, 
                                                    self.poldegree, self.polangle, self.a0_freq_correction )
        S1_new, S2_new, S3_new  = sampled_Spectrum_object.StokesParameters()
        S0_new                  = np.ones(S1_new.shape, dtype=float)




        sampled_weights    = base_photon_weight * np.ones(number_photons)
        sampled_K_photon   = np.asarray( (K0,K1,K2,K3) )
        sampled_X          = np.asarray( (sampled_t,sampled_x,sampled_y,sampled_z) )
        assigned_Stokes    = np.asarray( (S0_new,S1_new,S2_new,S3_new) )
        sampled_P_electron = np.asarray( (pt_out,px_out,py_out,pz_out) )



        self.xi_peak      = np.concatenate( (self.xi_peak,    sampled_xi_peak     ) , axis=0)
        self.W_photon     = np.concatenate( (self.W_photon,   sampled_weights     ) , axis=0)
        self.W_electron   = np.concatenate( (self.W_electron, sampled_weights     ) , axis=0)

        self.X_photon     = np.concatenate( (self.X_photon,   sampled_X          ) , axis=1)
        self.K_photon     = np.concatenate( (self.K_photon,   sampled_K_photon   ) , axis=1)
        self.S_photon     = np.concatenate( (self.S_photon,   assigned_Stokes    ) , axis=1)

        self.X_electron   = np.concatenate( (self.X_electron, sampled_X          ) , axis=1)
        self.P_electron   = np.concatenate( (self.P_electron, sampled_P_electron ) , axis=1)


        print ('   total photon number:',len(self.W_electron),len(self.W_photon))




