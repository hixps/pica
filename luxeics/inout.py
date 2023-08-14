import h5py
import numpy as np
from .__init__ import __version__


"""
Routines for input/output


Output file format

build   --- 
config
final-state
    spectrum



"""

class ParameterReader():

    def __init__(self):
        pass

    def read_laser_parameters(self):

        sigma_rescale           = bool(self.input_dict['control']['laser']['sigma_rescale'])


        self.omega0 = float( self.input_dict['laser']['omega0']  )
        self.a0     = float( self.input_dict['laser']['a0']      )
        self.Tpulse = float( self.input_dict['laser']['Tpulse']  )
        self.pol    = float( self.input_dict['laser']['pol']     )
        self.w0     = float( self.input_dict['laser']['w0']      )
        self.S3     = float( self.input_dict['laser']['S3']      )
        self.pulse  = self.input_dict['laser']['pulse']

        if self.pulse!='cos2':
            raise NotImplementedError

        # translate Tpulse==FWHM in fs --> sigma parameter which is dimensionless
        self.sigma = 2.0852201339 * self.Tpulse * self.omega0

        if sigma_rescale:
            sigma_crit =  float( input_dict['control']['laser']['sigma_crit']  )
            if sigma > sigma_crit:
                self.sigma_rescalefactor = sigma / sigma_crit
                print (f' >> rescale pulse duration: {sigma:.2f} -> {sigma_crit:.2f}: sigma_rescalefactor = {sigma_rescalefactor:.2f}')
            else:
                self.sigma_rescalefactor = 1.
        else:
            sigma_crit               = 0
            self.sigma_rescalefactor = 1.


    def read_beam_parameters(self):
        # extracting beam parameters from the YML File
        self.gamma0       = float( self.input_dict['beam']['gamma'] )
        self.energyspread = float( self.input_dict['beam']['energyspread'] )
        self.emittance_X  = float( self.input_dict['beam']['emittanceX'] )
        self.emittance_Y  = float( self.input_dict['beam']['emittanceY'] )
        self.beam_size_X  = float( self.input_dict['beam']['sigmaX'] )        # transverse beam size x axis in microns
        self.beam_size_Y  = float( self.input_dict['beam']['sigmaY'] )        # transverse beam size y axis in microns
        self.beam_length  = float( self.input_dict['beam']['sigmaL'] )        # longitudinal beam size in microns 
        self.beam_charge  = float( self.input_dict['beam']['beam_charge'])
        self.beam_focus_z = float( self.input_dict['beam']['beam_focus_z'])


class H5Writer():

    def __init__(self):
        pass


    def write_build(self, ff):
        # branch
        # version
        # commit-hash?
        # features?
        # with h5py.File(filename+'.h5' , 'w' ) as ff:
        ff['build/branch']  = 'devel'
        ff['build/version'] = __version__


    def write_file(self):

        with h5py.File( self.filename+'.h5', 'w' ) as ff:

            # build
            self.write_build(ff)
            # ff['build/branch']  = 'devel'
            # ff['build/version'] = __version__

            # config
            recursively_save_dict_contents_to_group(ff, 'config/', self.input_dict )

            with open(self.filename + '.yml', 'r') as fff:
                inputfile =''.join( fff.readlines() )
            ff['config/input-file'] = inputfile


            # spectrum
            fs=ff.create_group('final-state')



            fs['photon/position']               = self.X_photon
            fs['photon/position'].attrs['unit'] = 'micron'
            fs['photon/momentum']               = self.K_photon
            fs['photon/momentum'].attrs['unit'] = 'eV'
            fs['photon/weight'  ]               = self.W_photon
            fs['photon/weight'  ].attrs['unit'] = '1'
            fs['photon/Stokes'  ]               = self.Stokes_photon
            fs['photon/Stokes'  ].attrs['unit'] = '1'

            fs['photon/xi_peak']                = self.xi_peak
            fs['photon/xi_peak'].attrs['unit']  = '1'



            fs['electron/position']               = self.X_electron
            fs['electron/position'].attrs['unit'] = 'micron'
            fs['electron/momentum']               = self.P_electron
            fs['electron/momentum'].attrs['unit'] = 'eV'
            fs['electron/weight'  ]               = self.W_electron
            fs['electron/weight'  ].attrs['unit'] = '1'


def recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....
    """
    for key, item in dic.items():
        # print(key,item)
        if isinstance(item, dict):
            recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        elif isinstance(item, (list,tuple)):
            # has some problems with arrays containing strings. thus convert into floats (hoping that works always)
            h5file[path + key] = np.asarray(item,dtype=float)
        elif isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes,float,int)):
            h5file[path + key] = item
        else:
            raise ValueError('Cannot save %s type'%type(item))



def config_writer():
    # beam
    # control
    # laser
    # input-file
    # mpi-tasks
    # output
    with h5py.File(filename+'.h5' , 'w' ) as ff:

        b=ff.create_group('config')


    return True




    return True

def write_finalstate(filename,omega,theta,spectrum):
    # final-state
        # ics_spectrum
            # grid
            # spectrum
    with h5py.File(filename+'.h5' , 'w' ) as ff:
        finalstate = ff.create_group('final-state')
        finalstate['omega']=omega
        finalstate['theta']=theta
        finalstate['spectrum']=spectrum

    return





