# PICA

**P**olarized **I**CS **Ca**lculator

Monte Carlo code for the simulation of Inverse Compton Scattering gamma-ray spectra including the photon polarization.
The primary usage is for the LUXE experiment https://arxiv.org/abs/2102.02032.


## Source

https://github.com/hixps/pica


## Dependencies

* `python >= 3.10`
* `numpy`
* `scipy`
* `h5py`
* `pyyaml`
* `dacite`

## Installation

```console
pip install <path to pica source>
```

or

```console
pip install git+https://github.com/hixps/pica.git
```


## Useage

```python
import pica

SIM =  pica.ICSSimulation( input_filename='input_file' )
SIM.run()
```


## Input

Input files for configuring the simulation in the `yaml` format.



## Output

Simulation output files are `hdf5` in PTARMIGAN-compatible format: https://github.com/tgblackburn/ptarmigan


## Changelog


* v1.2.9
	* fixed the setup	
* v1.2.8
	* self.number_electrons
	* cross section normalization benchmarked against PTARMIGAN	
* v1.2.5
	* Stokes parameters in the x-z plane are now standard
	* Implemented methods to calculate Stokes parameters in scattering plane
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



	

	
