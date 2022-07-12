from setuptools import setup, find_packages

setup(
   name='bioplexpy',
   version='1.0.0',
   description='Python-side access to PPI data from Gygi lab',
   author='Roger Vargas, Ludwig Geistlinger, Tyrone Lee',
   author_email='roger_vargas@g.harvard.edu',
   packages=['bioplexpy'],  #same as name
   install_requires=['pandas','requests','anndata','networkx','numpy','matplotlib','biopython','scipy','pypdb'], #external packages as dependencies
)