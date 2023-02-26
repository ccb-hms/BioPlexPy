# BioPlexPy
-------------
[![Documentation Status](https://readthedocs.org/projects/bioplexpy/badge/?version=latest)](http://bioplexpy.readthedocs.io/?badge=latest)

Python-side access to PPI data from Gygi lab

## Citation
-------------

If you use the BioPlexPy package in published research, please cite:

Ludwig Geistlinger, Roger Vargas Jr, Tyrone Lee, Joshua Pan, Edward Huttlin, Robert Gentleman (2023) BioPlexR and BioPlexPy: integrated data products for the analysis of human protein interactions. *Bioinformatics*. doi: [10.1093/bioinformatics/btad091](https://doi.org/10.1093/bioinformatics/btad091). 

## Installation
-------------
The package can be installed through pypi or through github

`pip install bioplexpy`
or 
`pip install git+https://github.com/ccb-hms/BioPlexPy.git#egg=BioPlexPy`


## Usage
-------------
See the [BioPlex examples notebook](https://github.com/ccb-hms/BioPlexPy/blob/main/docs/BioPlex_Examples.ipynb) for basic usage and how to obtain
BioPlex datasets.

https://bioplexpy.readthedocs.io/en/latest/

## Running Tests
-------------

* `python -m doctest BioPlexPy/` to run the respective test for the file
*Warning*: By default, doctests will not print anything if tests are successful.
Add the `-v` option to print verbose output.

