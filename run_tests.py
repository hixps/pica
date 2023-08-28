
import yaml
import numpy as np


import pica
# import matplotlib.pyplot as plt


input_filename = 'small_beam_full.yml'

with open( input_filename, 'r' ) as stream:
    input_dict = yaml.load(stream, Loader=yaml.SafeLoader)

print (input_dict)



pica.main_program( input_filename )



