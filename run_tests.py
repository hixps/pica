
import yaml
import numpy as np


import luxeics
# import matplotlib.pyplot as plt


input_filename = 'large_beam_full.yml'

with open( input_filename, 'r' ) as stream:
    input_dict = yaml.load(stream, Loader=yaml.SafeLoader)

print (input_dict)



luxeics.main_program( input_filename )



