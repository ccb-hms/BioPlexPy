# setup.py
from setuptools import setup, find_packages
from distutils.util import convert_path

main_ns = {}
ver_path = convert_path('bioplexpy/version.py')
with open(ver_path) as ver_file:
    exec(ver_file.read(), main_ns)

setup(
   name='bioplexpy',
   version=main_ns['__version__'],
   description='Python-side access to PPI data from Gygi lab',
   author='Roger Vargas, Ludwig Geistlinger, Tyrone Lee',
   author_email='roger_vargas@g.harvard.edu',
   packages=['bioplexpy'],  #same as name
   install_requires=['pandas','requests','anndata','networkx','numpy','matplotlib','biopython','scipy','pypdb'], #external packages as dependencies
)