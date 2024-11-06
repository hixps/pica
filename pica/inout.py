import numpy as np
import h5py

from dataclasses import dataclass


from .__init__ import __version__
from .constants import hbar



"""
Routines for input/output


Output file format

build   --- 
config
final-state
    spectrum


"""


@dataclass
class Control_Beam:
    sample_electrons: float | str | int
    sample_batch_size: float | str | int = '1e7'
        
    def __post_init__(self):
        self.sample_electrons  = int(float(self.sample_electrons))
        self.sample_batch_size = int(float(self.sample_batch_size))

@dataclass
class Control_Laser:
    pulse_rescale_bias: float
        
@dataclass
class Control_Detector:
    a0_freq_correction: bool

@dataclass
class Control:
    sampling: str
    xsection: str
    beam: Control_Beam
    laser: Control_Laser
    detector: Control_Detector
        
@dataclass
class Unit:
    momentum: str
    position: str
        
    def __post_init__(self):
        if self.position not in ('micron',):
            raise ValueError

        if self.momentum == 'eV':
            self.momentum_scale = 1
        elif self.momentum == 'keV':
            self.momentum_scale = 1e3
        elif self.momentum == 'MeV':
            self.momentum_scale = 1e6
        elif self.momentum == 'GeV':
            self.momentum_scale = 1e9
        else:
            raise ValueError


@dataclass
class Beam:
    gamma: float | str
    energyspread: float | str
    emittanceX: float | str
    emittanceY: float | str
    sigmaX: float | str
    sigmaY: float | str
    sigmaL: float | str
    beam_charge: float | str
    beam_focus_z: float | str
    baseline: float | str
        
    def __post_init__(self):
        self.gamma        = float(self.gamma)
        self.energyspread = float(self.energyspread)
        self.emittanceX   = float(self.emittanceX)
        self.emittanceY   = float(self.emittanceY)
        self.sigmaX       = float(self.sigmaX)
        self.sigmaY       = float(self.sigmaY)
        self.sigmaL       = float(self.sigmaL)
        self.beam_charge  = float(self.beam_charge)
        self.beam_focus_z = float(self.beam_focus_z)
        self.baseline     = float(self.baseline)

@dataclass
class Laser:
    a0: float
    omega0: float
    TFWHM: float
    pulse: str
    w0: float
    polangle: float
    poldegree: float

    def __post_init__(self):
        # translate TFWHM==FWHM in fs --> sigma parameter which is dimensionless, specific for cos^2 envelope
        self.sigma = 0.25*np.pi/np.arccos(1/2**0.25)/hbar * self.TFWHM * self.omega0
        # print (0.25*np.pi/np.arccos(1/2**0.25)*c/hbarc , 2.0852201339)
        # self.sigma = 2.0852201339 * self.TFWHM * self.omega0
        # the numerical factor is 0.25*pi/arccos(1/2**(1/4)) * 1.52, where the factor 1.52 comes from the transition from eV to fs


@dataclass
class Detector:
    # pdim:  int
    omega: list[float | str]
    theta: list[float | str]
    phi:   list[float | str]

    def __post_init__(self):
        self.omega = [float(i) for i in self.omega]
        self.theta = [float(i) for i in self.theta]
        self.phi   = [float(i) for i in self.phi  ]


        
        
@dataclass
class PICA_Config:
    control: Control
    unit: Unit
    beam: Beam
    laser: Laser
    detector: Detector
 




# class ParameterReader():

#     def __init__(self):
#         pass

#     def read_laser_parameters(self):


#         self.omega0    = float( self.input_dict['laser']['omega0']  )
#         self.a0        = float( self.input_dict['laser']['a0']      )
#         self.TFWHM     = float( self.input_dict['laser']['TFWHM']  )
#         self.w0        = float( self.input_dict['laser']['w0']      )
#         self.poldegree = float( self.input_dict['laser']['poldegree']   )
#         self.polangle  = float( self.input_dict['laser']['polangle']    )
#         self.pulse     = self.input_dict['laser']['pulse']

#         if self.pulse!='cos2':
#             raise NotImplementedError



#     def read_beam_parameters(self):
#         # extracting beam parameters from the YML File
#         self.gamma0       = float( self.input_dict['beam']['gamma'] )
#         self.energyspread = float( self.input_dict['beam']['energyspread'] )
#         self.emittance_X  = float( self.input_dict['beam']['emittanceX'] )
#         self.emittance_Y  = float( self.input_dict['beam']['emittanceY'] )
#         self.beam_size_X  = float( self.input_dict['beam']['sigmaX'] )        # transverse beam size x axis in microns
#         self.beam_size_Y  = float( self.input_dict['beam']['sigmaY'] )        # transverse beam size y axis in microns
#         self.beam_length  = float( self.input_dict['beam']['sigmaL'] )        # longitudinal beam size in microns 
#         self.beam_charge  = float( self.input_dict['beam']['beam_charge'])
#         self.beam_focus_z = float( self.input_dict['beam']['beam_focus_z'])
#         self.baseline     = float( self.input_dict['beam']['baseline'])

#     def read_detector(self):

#         # extract detector parameters
#         self.pdim           = int(self.input_dict['detector']['pdim'])
#         self.omega_detector = [float(w) for w in self.input_dict['detector']['omega']]
#         self.theta_detector = [float(t) for t in self.input_dict['detector']['theta']]
#         self.phi_detector   = [float(p) for p in self.input_dict['detector']['phi']]


class H5Writer():

    def __init__(self):
        pass


    def write_build(self, ff):
        # branch
        # version
        # commit-hash?
        # features?
        ff['build/branch']  = 'polarization'
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



            fs['photon/position']                     = self.X_photon.T
            fs['photon/position'].attrs['unit']       = self.config.unit.position
            fs['photon/momentum']                     = self.K_photon.T
            fs['photon/momentum'].attrs['unit']       = self.config.unit.momentum
            fs['photon/weight'  ]                     = self.W_photon
            fs['photon/weight'  ].attrs['unit']       = '1'
            fs['photon/polarization'  ]               = self.S_photon.T
            fs['photon/polarization'  ].attrs['unit'] = '1'

            fs['photon/xi_peak']                      = self.xi_peak
            fs['photon/xi_peak'].attrs['unit']        = '1'



            fs['electron/position']               = self.X_electron.T
            fs['electron/position'].attrs['unit'] = self.config.unit.position
            fs['electron/momentum']               = self.P_electron.T
            fs['electron/momentum'].attrs['unit'] = self.config.unit.momentum
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
    # units
    # input-file
    with h5py.File(filename+'.h5' , 'w' ) as ff:

        b=ff.create_group('config')

    return True


class H5Reader():

    def __init__(self):
        pass


    def read(self):
        with h5py.File(self.filename+'.h5' ,'r') as ff:

            self.K_photon       = ff['final-state/photon/momentum'      ][()].T
            self.X_photon       = ff['final-state/photon/position'      ][()].T
            self.S_photon       = ff['final-state/photon/polarization'  ][()].T
            self.W_photon       = ff['final-state/photon/weight'        ][()]

            self.xi_peak        = ff['final-state/photon/xi_peak'       ][()]

            self.X_electron     = ff['final-state/electron/position'    ][()].T
            self.P_electron     = ff['final-state/electron/momentum'    ][()].T
            self.W_electron     = ff['final-state/electron/weight'      ][()]





