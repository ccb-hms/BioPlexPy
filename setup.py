from setuptools import setup, find_packages

setup(
   name='BioPlexPy',
   version='1.0',
   description='Python-side access to PPI data from Gygi lab',
   author='Roger Vargas',
   author_email='roger_vargas@g.harvard.edu',
   packages=['BioPlexPy'],  #same as name
   install_requires=['pandas','requests','anndata'], #external packages as dependencies
)