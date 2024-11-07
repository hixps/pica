"""
Routines for writing the *.h5 files


Output file format

build   --- 
config
final-state
    spectrum


"""

import numpy as np
import h5py

from .__init__ import __version__



class H5Writer():

    def __init__(self):
        pass


    def write_build(self, ff):
        # branch
        # version
        # commit-hash?
        # features?
        ff['build/branch']  = 'main'
        ff['build/version'] = __version__


    def write_file(self):

        with h5py.File( self.filename+'.h5', 'w' ) as ff:

            # build
            self.write_build(ff)
            # ff['build/branch']  = 'devel'
            # ff['build/version'] = __version__

            # config
            _recursively_save_dict_contents_to_group(ff, 'config/', self.input_dict )

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


def _recursively_save_dict_contents_to_group(h5file, path, dic):
    """
    ....
    """
    for key, item in dic.items():
        # print(key,item)
        if isinstance(item, dict):
            _recursively_save_dict_contents_to_group(h5file, path + key + '/', item)
        elif isinstance(item, (list,tuple)):
            # has some problems with arrays containing strings. thus convert into floats (hoping that works always)
            h5file[path + key] = np.asarray(item,dtype=float)
        elif isinstance(item, (np.ndarray, np.int64, np.float64, str, bytes,float,int)):
            h5file[path + key] = item
        else:
            raise ValueError('Cannot save %s type'%type(item))





