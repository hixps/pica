"""
Data structures for reading of config files and for reading simulation data into a simulation class



"""

import numpy as np
import h5py

from dataclasses import dataclass

from .__init__ import __version__
from .constants import hbar


@dataclass
class Control_Beam:
    sample_electrons: float 
    sample_batch_size: float = 1e7
        
    def __post_init__(self):
        self.sample_electrons  = int( self.sample_electrons)
        self.sample_batch_size = int( self.sample_batch_size)

@dataclass
class Control_Laser:
    pulse_rescale_bias: float
        
@dataclass
class Control_Detector:
    a0_freq_correction: bool

@dataclass
class Control:
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
    gamma: float
    energyspread: float 
    emittanceX: float 
    emittanceY: float
    sigmaX: float
    sigmaY: float
    sigmaZ: float
    charge: float
    focus_z: float
    baseline: float
        
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
    omega: list[float]
    theta: list[float]
    phi:   list[float]

@dataclass
class PICA_Config:
    control: Control
    unit: Unit
    beam: Beam
    laser: Laser
    detector: Detector
 


class H5Reader():

    def __init__(self):
        pass


    def read(self):
        with h5py.File( f'{self.filename}.h5', 'r') as ff:

            self.K_photon       = ff['final-state/photon/momentum'      ][()].T
            self.X_photon       = ff['final-state/photon/position'      ][()].T
            self.S_photon       = ff['final-state/photon/polarization'  ][()].T
            self.W_photon       = ff['final-state/photon/weight'        ][()]

            self.xi_peak        = ff['final-state/photon/xi_peak'       ][()]

            self.X_electron     = ff['final-state/electron/position'    ][()].T
            self.P_electron     = ff['final-state/electron/momentum'    ][()].T
            self.W_electron     = ff['final-state/electron/weight'      ][()]





