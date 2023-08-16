# luxeics


Simulate ICS gamma photon spectra for the LUXE experiment https://arxiv.org/abs/2102.02032.

Input files in ```yaml``` format

Output files in PTARMIGAN-hdf5 format: https://github.com/tgblackburn/ptarmigan

## Dependencies

* python >= 3.9
* numpy
* scipy
* h5py
* pyyaml

## Installation

```python
python -m build;
pip install dist/luxeics-<version number>.tar.gz 
```

or

```python
python setup.py install
```


## Useage

```python
import luxeics

SIM =  luxeics.ICSSimulation( input_filename='input_file' )
SIM.run()
```


## Changelog

* v1.2.0 
	* linear laser polarization and linear Stokes parameters of gamma rays
* v1.2.1
	* laser polarization control via polarization degree and polarization angle (only LP)
	* "Full" polarized cross section based on Jauch-Rohrlich expression with explicit polarization vectors
	* Analysis routines for rotating Stokes vectors from scattering plane to laboratory frame	
	
