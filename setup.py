from distutils.core import setup

setup(
    name='aai',
    version='0.0.1',
    author='Donovan Parks',
    author_email='donovan.parks@gmail.com',
    packages=['aai'],
    scripts=['bin/aai'],
    url='http://pypi.python.org/pypi/aai/',
    license='GPL3',
    description='Assess the quality of putative genome bins.',
    long_description=open('README.md').read(),
    install_requires=[
        "numpy >= 1.8.0",
        "scipy >= 0.9.0",
        "matplotlib >= 1.1.0"],
)
