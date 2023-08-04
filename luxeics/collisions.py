import numpy as np
from numpy import *
from scipy import *
from scipy.integrate import cumtrapz
#from scipy.interpolate import RectBivariateSpline


from .constants import elec_mass,finestruct
from .spectrum import Compton_spectrum_linear,StokesParameters
from .inout import *
from .__init__ import __version__

import yaml



"""
Routines for a collision of a beam with a laser pulse
"""

"""
focused electron beam:
random.multivariate_normal(mean, cov, size=None, check_valid='warn', tol=1e-8)
to sample correlated transverse phase space ellipses, (x,px) and (y,py);
need to determine the covariance from the actual beam parameters
"""


S3 = 1

############################################
# This is where the real stuff starts ######
############################################


def Covariance_Matrix(beam_size,rms_angle,L):
    cov_at_IP_Focus = np.array([[beam_size**2, 0], [0, rms_angle**2]])
     


    #Inverse of beam transport matrix
    M = np.array([[1, -L], [0, 1]])
    #Calculating the Covariance matrix at the Interaction Point
    cov=(M@cov_at_IP_Focus@M.T)
    return(cov)


def GenerateElectronBeam(input_filename, sample_electrons):
    
    
    
    with open( input_filename, 'r' ) as stream:
        input_dict = yaml.load(stream, Loader=yaml.SafeLoader)
    
    
    
    
    
    
    # extracting beam parameters from the YML File
    gamma0       = float( input_dict['beam']['gamma'] )
    energyspread = float( input_dict['beam']['energyspread'] )
    emittance    = float( input_dict['beam']['emittance'] )
    beam_size_X  = float( input_dict['beam']['sigmaX'] )        # transverse beam size x axis in microns
    beam_size_Y  = float( input_dict['beam']['sigmaY'] )        # transverse beam size y axis in microns
    beam_length  = float( input_dict['beam']['sigmaL'] )        # longitudinal beam size in microns 
    beam_charge  = float( input_dict['beam']['beam_charge'])
    L            = float( input_dict['beam']['L'])
    
    
    
    # rms angles for x and y space
    rms_angle_X   = emittance / beam_size_X / gamma0
    rms_angle_Y   = emittance / beam_size_Y / gamma0
    
    
    #Mean and Covariant matrices for bivariate gaussian
    
    mean_x = [0,0]  # x-offset and x' offset for focus 
    cov_x = Covariance_Matrix(beam_size_X,rms_angle_X,L)

    mean_y = [0,0]  # y-offset and y' offset for focus
    cov_y = Covariance_Matrix(beam_size_Y,rms_angle_Y,L)

    
    #Sampling (x,x') and (y,y')
    x,theta_x = np.random.multivariate_normal(mean_x, cov_x, sample_electrons).T
    y,theta_y = np.random.multivariate_normal(mean_y, cov_y, sample_electrons).T
    
    
    
    gamma        = np.random.normal( gamma0 , gamma0*energyspread , sample_electrons )
    beta         = sqrt(1-1./gamma**2)

    
    pz0          = gamma*beta*cos(theta_x)*cos(theta_y)
    px0          = gamma*beta*sin(theta_x)
    py0          = gamma*beta*sin(theta_y)
    pt0          = sqrt( 1 + px0**2 + py0**2 + pz0**2 )


    return(x,theta_x,y,theta_y,gamma,pt0,px0,py0,pz0)


def gaussian_sum_normalized( mean , rms , nbeam ):
    """
    Help generate normalized Gaussian distribtions for
    """

    values        = mean + linspace(-3*rms,3*rms,nbeam)
    dvalues       = values[1] - values[0]
    beam_weights  = np.exp(-(values-mean)**2/2/rms**2) / sqrt(2*pi) / rms * dvalues

    return values, beam_weights


def integrate_spectrum_emittancemacroelectron( a0 , omega0 , sigma , sigma_rescalefactor , # plane wave laser parameters
                                               omega, theta, phi   ,                       # detector           
                                               beam_energies,weights_energy,               # beam energy spread
                                               beam_angles_x,weights_angle_x,              # beam theta_x spread
                                               beam_angles_y,weights_angle_y  ):

    # empty array for spectrum
    spectrum = zeros( (omega*theta*phi).shape )


    # iterate over the angular beam distribution and sum the spectra
    for ii,(gammab,w_gamma) in enumerate(zip(beam_energies,weights_energy)):

        print (f' >> {ii}: gamma = {gammab:.4f}, w_gamma = {w_gamma:.4g}')
        betab       = sqrt(1-1./gammab**2)


        for theta_x,w_x in zip(beam_angles_x,weights_angle_x):
            for theta_y,w_y in zip(beam_angles_y,weights_angle_y):



                # print (theta_x,theta_y,w_x,w_y),

                # set up initial electron momentum
                pz0    = gammab*betab*cos(theta_x)*cos(theta_y)
                px0    = gammab*betab*sin(theta_x)
                py0    = gammab*betab*sin(theta_y)
                pt0    = sqrt( 1 + px0**2 + py0**2 + pz0**2 )

                # print (gammab,pt0)


                U_in   = [pt0,px0,py0,pz0]

                spec_weight = w_x*w_y*w_gamma * sigma_rescalefactor

                spectrum += spec_weight * Compton_spectrum_linear( U_in , a0 , omega0 , sigma / sigma_rescalefactor , omega, theta, phi, S3 )

    return spectrum




def main_program( input_filename ):

    with open( input_filename, 'r' ) as stream:
        input_dict = yaml.load(stream, Loader=yaml.SafeLoader)


    filename   = input_dict['control']['name']
    mode       = input_dict['control']['mode']


    if mode=='full':
        print ('>>> mode == full')
        # X_photon,K_photon,W_photon,X_electron,P_electron,W_electron,xi_peak = main_program_full( input_filename )
        X_photon,K_photon,W_photon,Stokes_photon,X_electron,P_electron,W_electron,xi_peak = main_program_full( input_filename )

    # elif mode=='simple':
    #     print ('>>> mode == simple')
    #     main_program_simple( input_filename )

    else:
        NotImplementedError



    ###############
    # HDF5 output #
    ###############

    with h5py.File( filename+'.h5' , 'w' ) as ff:

        # build
        ff['build/branch']  = 'devel'
        ff['build/version'] = __version__

        # config
        recursively_save_dict_contents_to_group(ff, 'config/', input_dict )

        with open(input_filename,'r') as fff:
            inputfile =''.join( fff.readlines() )
        ff['config/input-file'] = inputfile


        # spectrum
        fs=ff.create_group('final-state')



        fs['photon/position']               = X_photon
        fs['photon/position'].attrs['unit'] = 'micron'
        fs['photon/momentum']               = K_photon
        fs['photon/momentum'].attrs['unit'] = 'eV'
        fs['photon/weight'  ]               = W_photon
        fs['photon/weight'  ].attrs['unit'] = '1'
        fs['photon/Stokes'  ]               = Stokes_photon
        fs['photon/Stokes'  ].attrs['unit'] = '1'




        fs['electron/position']               = X_electron
        fs['electron/position'].attrs['unit'] = 'micron'
        fs['electron/momentum']               = P_electron
        fs['electron/momentum'].attrs['unit'] = 'eV'
        fs['electron/weight'  ]               = W_electron
        fs['electron/weight'  ].attrs['unit'] = '1'

        fs['photon/xi_peak']                     = xi_peak
        fs['photon/xi_peak'].attrs['unit']       = '1'





def main_program_full( input_filename ):


    with open( input_filename, 'r' ) as stream:
        input_dict = yaml.load(stream, Loader=yaml.SafeLoader)

    # control parameters
    filename                = input_dict['control']['name']
    sigma_rescale           = bool(input_dict['control']['laser']['sigma_rescale'])
    sample_electrons_total  = int(float( input_dict['control']['sample_electrons'] ))

    try:
        sample_batch_size       = int(float(input_dict['control']['sample_batch_size']))
    except KeyError:
        print ('"control/sample_batch_size" not specified, set to default value 1e8')
        sample_batch_size       = int(1e7)



    n_samp = sample_electrons_total // sample_batch_size
    r_samp = sample_electrons_total % sample_batch_size

    if r_samp:
        sampling_batches = [sample_batch_size]*int(n_samp) + [r_samp]
    else:
        sampling_batches = [sample_batch_size]*int(n_samp)

    print (f'n_batches={n_samp}, remainder={r_samp}, number_of_particles={sum(sampling_batches):.4g}={sample_electrons_total:.4g}')






    # laser parameters
    omega0 = float( input_dict['laser']['omega0']  )
    a0     = float( input_dict['laser']['a0']      )
    Tpulse = float( input_dict['laser']['Tpulse']  )
    pol    = float( input_dict['laser']['pol']     )
    w0     = float( input_dict['laser']['w0']      )
    pulse  = input_dict['laser']['pulse']

    if pulse!='cos2':
        raise NotImplementedError

    # translate Tpulse==FWHM in fs --> sigma parameter which is dimensionless
    sigma = 2.0852201339 * Tpulse * omega0

    if sigma_rescale:
        sigma_crit =  float( input_dict['control']['laser']['sigma_crit']  )
        if sigma > sigma_crit:
            sigma_rescalefactor = sigma / sigma_crit
            print (f' >> rescale pulse duration: {sigma:.2f} -> {sigma_crit:.2f}: sigma_rescalefactor = {sigma_rescalefactor:.2f}')
        else:
            sigma_rescalefactor = 1.
    else:
        sigma_crit          = 0
        sigma_rescalefactor = 1.



    # extract beam parameters: Parsing of the Beam Parameters has been shifted to the function "Generate_electron_Beam"
    
#     gamma0       = float( input_dict['beam']['gamma'] )
#     energyspread = float( input_dict['beam']['energyspread'] )
#     emittance    = float( input_dict['beam']['emittance'] )
#     beam_size_X    = float( input_dict['beam']['sigmaX'] )        # transverse beam size x axis in microns
#     beam_size_Y    = float( input_dict['beam']['sigmaY'] )        # transverse beam size y axis in microns
    beam_length  = float( input_dict['beam']['sigmaL'] )        # longitudinal beam size in microns 
    beam_charge  = float( input_dict['beam']['beam_charge'])
    

    # extract detector parameters
    omega_detector = [float(w) for w in input_dict['detector']['omega']]
    theta_detector = [float(t) for t in input_dict['detector']['theta']]
    phi_detector   = [float(p) for p in input_dict['detector']['phi']]



    number_electrons = int( beam_charge / 1.60217653e-19)

    # sample_electrons = int(5e7) # number of elecrons used for MC sampling, determines weight of each emitted photon
    electron_weight  = number_electrons / sample_electrons_total
    print (f' >> number_electrons   = {number_electrons}')
    print (f' >> sample_electrons   = {sample_electrons_total}')
    print (f' >> electron weight    = {electron_weight:.4g}')


    keep_xi_peak   = []
    X_photon       = []
    K_photon       = []
    W_photon       = []

    Stokes_photon  = []

    X_electron     = []
    P_electron     = []
    W_electron     = []




    print (' >> MC sampling')

    for jj,sample_electrons in enumerate(sampling_batches):

        print (f'  > batch {jj} : {sample_electrons:.2g} macroelectrons of weight {electron_weight:.5g}')

        # print (len(W_photon))

        
        # Correlated Phase Space
        x,theta_x,y,theta_y,gamma,pt0,px0,py0,pz0 = GenerateElectronBeam(input_filename,sample_electrons)
        
        r            = np.sqrt( x**2 + y**2 )
        xi_peak      = a0 * exp( -r**2/w0**2 )
        U_in         = asarray([pt0,px0,py0,pz0])

        
        # photon detector
        omega        = np.random.uniform(omega_detector[0],     omega_detector[1],     sample_electrons) 
        theta        = np.random.uniform(theta_detector[0],     theta_detector[1],     sample_electrons) 
        phi          = np.random.uniform(phi_detector[0]*np.pi, phi_detector[1]*np.pi, sample_electrons)
        
        spec_weight  =  sigma_rescalefactor


        # print (gamma.shape,theta_x.shape,U_in.shape,)
        rad_int  = Compton_spectrum_linear( U_in , xi_peak , omega0 , sigma / sigma_rescalefactor , omega, theta, phi, S3 )

        spectrum = spec_weight * rad_int * sin(theta)

        spec_max = amax(spectrum)



        s_val        = np.random.uniform(0             , spec_max,       sample_electrons)


        samplingvolume     = spec_max * (omega_detector[1]-omega_detector[0]) * (theta_detector[1]-theta_detector[0])
        base_photon_weight = 2*pi * samplingvolume * electron_weight

        selector1          = s_val < spectrum
        number_photons     = sum(selector1)
        print ('   base photon weight :' , base_photon_weight)
        print ('   number photons     :' , number_photons)



        keep_gamma   = gamma[selector1]
        keep_theta_x = theta_x[selector1]
        keep_theta_y = theta_y[selector1]

        keep_omega   = omega[selector1]
        keep_theta   = theta[selector1]
        keep_phi     = phi[selector1]



        keep_x       = x[selector1]
        keep_y       = y[selector1]
        keep_zeta    = np.random.normal( 0 , beam_length  , number_photons )

        keep_z       = +0.5 * keep_zeta
        keep_t       = -0.5 * keep_zeta


        K0           = keep_omega
        K1           = keep_omega * sin(keep_theta) * cos(keep_phi)
        K2           = keep_omega * sin(keep_theta) * sin(keep_phi)
        K3           = keep_omega * cos(keep_theta) 



        keep_xi_peak.extend(xi_peak[selector1])

        X_photon.extend(    np.asarray( (keep_t,keep_x,keep_y,keep_z) ).T)
        K_photon.extend(    np.asarray( (K0,K1,K2,K3) ).T                )
        W_photon.extend(    base_photon_weight * ones(number_photons)    )


        # electrons that have emitted
        px_out  = px0[selector1] * elec_mass - K1
        py_out  = py0[selector1] * elec_mass - K2
        pz_out  = pz0[selector1] * elec_mass - K3 
        pt_out  = sqrt( elec_mass**2 + px_out**2 + py_out**2 + pz_out**2 )


        X_electron.extend( np.asarray( (keep_t,keep_x,keep_y,keep_z) ).T )
        P_electron.extend( np.asarray( (pt_out,px_out,py_out,pz_out) ).T )
        W_electron.extend( base_photon_weight * ones(number_photons)     )

        print ('   total photon number:',len(W_electron))

        # Stokes Parameters of emitted Photons
        U_in_keep  = asarray([pt0[selector1],px0[selector1],py0[selector1],pz0[selector1]])
        S1_new, S2_new, S3_new = StokesParameters( U_in_keep , a0 , omega0 , sigma , keep_omega, keep_theta, keep_phi , S3 )
        # print (S1_new,S2_new,S3_new)


        Stokes_photon.extend( np.asarray([S1_new,S2_new,S3_new]).T )



    # convert lists into ndarrays
    X_photon     = np.asarray(X_photon)
    K_photon     = np.asarray(K_photon)
    W_photon     = np.asarray(W_photon)
    Stokes_photon= np.asarray(Stokes_photon)
    keep_xi_peak = np.asarray(keep_xi_peak)

    X_electron   = np.asarray(X_electron)
    P_electron   = np.asarray(P_electron)
    W_electron   = np.asarray(W_electron)





    return X_photon,K_photon,W_photon,Stokes_photon,X_electron,P_electron,W_electron,keep_xi_peak








