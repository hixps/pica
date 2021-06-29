import setuptools
from distutils.core import setup


with open("README.md", "r") as fh:
    long_description = fh.read()

# setuptools.setup(
setup(
    name="qcrad",
    version="0.2",
    author="Daniel Seipt",
    author_email="d.seipt@gsi.de",
    description="ICS simulations, calculate quasi-classical radiation spectra from trajectories",
    # long_description=long_description,
    # long_description_content_type="text/markdown",
    # url          = 'https://github.com/danielseipt/rotating_cascade.git',
    packages     = ['qcrad'] ,
    package_dir  = {'qcrad' : 'qcrad'},
    install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'mpmath'
      ],
    test_suite='nose.collector',
    tests_require=['nose'],
)

# print setuptools.find_packages()