# pica

Polarized ICS CAlculator


Simulate ICS gamma photon spectra for the LUXE experiment https://arxiv.org/abs/2102.02032.

Input files in ```yaml``` format

Output files in PTARMIGAN-hdf5 format: https://github.com/tgblackburn/ptarmigan

## Source

`https://github.com/hixps/pica'


## Dependencies

* `python >= 3.9`
* `numpy`
* `scipy`
* `h5py`
* `pyyaml`

## Installation

```python
python -m build;
pip install dist/pica-<version number>.tar.gz 
```

or

```python
python setup.py install
```


## Useage

```python
import pica

SIM =  pica.ICSSimulation( input_filename='input_file' )
SIM.run()
```


## Changelog

* v1.2.0 
	* linear laser polarization and linear Stokes parameters of gamma rays
* v1.2.1
	* laser polarization control via polarization angle (only LP, 100% polarization degree)
	* "Full" polarized cross section based on Jauch-Rohrlich expression with explicit polarization vectors
	* Analysis routines for rotating Stokes vectors from scattering plane to laboratory frame	
* v1.2.2
	* partially polarized laser	
* v1.2.3
	* software rename luxeics -> pica
	
