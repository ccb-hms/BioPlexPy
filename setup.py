from setuptools import setup, find_packages

setup(
   name='BioPlexPy',
   version='0.99.0',
   description='Python-side access to PPI data from Gygi lab',
   author='Roger Vargas',
   author_email='roger_vargas@g.harvard.edu',
   packages=['BioPlexPy'],  #same as name
   install_requires=['pandas','requests','anndata','networkx','numpy','matplotlib','pypdb'], #external packages as dependencies
)