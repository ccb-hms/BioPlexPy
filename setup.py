from setuptools import setup, find_packages

setup(
   name='BioPlexPy',
   version='0.99.0',
   description='Python-side access to PPI data from Gygi lab',
   author='Roger Vargas',
   author_email='roger_vargas@g.harvard.edu',
   packages=['BioPlexPy'],  #same as name
<<<<<<< HEAD
   install_requires=['pandas','requests','anndata','networkx','numpy','matplotlib'], #external packages as dependencies
)
=======
   install_requires=['pandas','requests','anndata','networkx'], #external packages as dependencies
)
>>>>>>> bfa64393a8a9a2653b6024435ca48a44057136aa
