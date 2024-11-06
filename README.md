# pica

Polarized ICS CAlculator


Simulate ICS gamma photon spectra for the LUXE experiment https://arxiv.org/abs/2102.02032.

Input files in ```yaml``` format

Output files in PTARMIGAN-hdf5 format: https://github.com/tgblackburn/ptarmigan

## Source

```https://github.com/hixps/pica```


## Dependencies

* `python >= 3.10`
* `numpy`
* `scipy`
* `h5py`
* `pyyaml`

## Installation

```console
pip install <path to pica source>
	
pip install https://github.com/hixps/pica
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
* v1.2.4
	* switched Stokes parameter definition to S3=circular
	* Stokes parameter S0=1 for consistency with PTARMIGAN
	* rename 'Stokes' group in hdf5 output to 'polarization'
	* switched to 'pulse_rescale_bias' instead of 'sigma_rescaling' parameter
	* record units of output, unit section in input file
		* config/unit/position = ('micron')
		* config/unit/momentum = ('eV','keV','MeV','GeV')
	* pulse length parameter renamed to TFWHM
	* cleanup of old cross sections, only Cross_Section_Full useable
* v1.2.5
	* Stokes parameters in the x-z plane are now standard
	* Implemented methods to calculate Stokes parameters in scattering plane

* v1.2.8
	* self.number_electrons
	* cross section normalization benchmarked against PTARMIGAN
	
* v1.2.9
	* fixed the setup	
	
