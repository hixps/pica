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

def write_build(filename):
    # branch
    # commit-hash
    # features
    # version
    with h5py.File(filename+'.h5' , 'w' ) as ff:
        ff['build/branch']  = 'devel'
        ff['build/version'] = __version__



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


if __name__ == '__main__':

    write_build('test_io.h5')
    write_config('test_io.h5')
    write_finalstate('test_io.h5')
    # print ('1')



