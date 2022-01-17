import setuptools
from distutils.core import setup


with open("README.md", "r") as fh:
    long_description = fh.read()

# setuptools.setup(
setup(
    name="luxeics",
    version="1.1",
    author="Daniel Seipt",
    author_email="d.seipt@hi-jena.gsi.de",
    description="ICS simulations for LUXE",
    url          = 'https://github.com/danielseipt/luxeics.git',
    packages     = ['luxeics'] ,
    package_dir  = {'luxeics' : 'luxeics'},
    install_requires=[
                'numpy',
                'scipy',
                'matplotlib'
                ],
    test_suite='nose.collector',
    tests_require=['nose'],
)

# print setuptools.find_packages()