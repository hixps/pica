from .__version__ import __version__
from . import config

from .constants import *
from .collisions import *


# import re
# float_loader = yaml.SafeLoader
# float_loader.add_implicit_resolver(
#     u'tag:yaml.org,2002:float',
#     re.compile(u'''^(?:
#      [-+]?(?:[0-9][0-9_]*)\\.[0-9_]*(?:[eE][-+]?[0-9]+)?
#     |[-+]?(?:[0-9][0-9_]*)(?:[eE][-+]?[0-9]+)
#     |\\.[0-9_]+(?:[eE][-+][0-9]+)?
#     |[-+]?[0-9][0-9_]*(?::[0-5]?[0-9])+\\.[0-9_]*
#     |[-+]?\\.(?:inf|Inf|INF)
#     |\\.(?:nan|NaN|NAN))$''', re.X),
#     list(u'-+0123456789.'))

