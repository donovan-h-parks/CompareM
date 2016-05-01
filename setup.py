from distutils.core import setup
from Cython.Build import cythonize

import os

def version():
    setupDir = os.path.dirname(os.path.realpath(__file__))
    versionFile = open(os.path.join(setupDir, 'comparem', 'VERSION'))
    return versionFile.read().strip()

setup(
    name='comparem',
    version=version(),
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['comparem', 'comparem.plots'],
    scripts=['bin/comparem'],
    package_data={'comparem': ['VERSION']},
    ext_modules = cythonize("./comparem/aai_calculator.pyx"),
    url='http://pypi.python.org/pypi/comparem/',
    license='GPL3',
    description='A toolbox for comparative genomics.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.9.0",
        "biolib >= 0.0.20"],
)
