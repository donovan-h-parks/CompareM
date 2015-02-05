from distutils.core import setup

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'aai', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='aai',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['aai'],
    scripts=['bin/aai'],
    package_data={'aai' : ['VERSION']},
    url='http://pypi.python.org/pypi/aai/',
    license='GPL3',
    description='Calculate amino acid identity (AAI) between genomes.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.8.0"],
)
