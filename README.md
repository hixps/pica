# luxeics


Simulate ICS gamma photon spectra for the LUXE experiment https://arxiv.org/abs/2102.02032.

Input files in ```yaml``` format

Output files in PTARMIGAN-hdf5 format: https://github.com/tgblackburn/ptarmigan

## Installation

```python
python setup.py install
```

## Useage

```python
import luxeics
luxeics.main_program( input_filename='config.yml' )
```
