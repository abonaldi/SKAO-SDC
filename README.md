# Science Data Challenge 3a: EoR foreground removal

This repository contains a collection of scripts used to generate mock SKA-Mid-observed neutral hydrogen (HI) datacubes for the second SKA Science Data Challenge ([SDC2](https://sdc2.astronomers.skatelescope.org/)). A detailed description on the simulations can be found in Section 3 of the SDC2 [paper](https://arxiv.org/abs/2303.07943).

### Prerequisites

#### Python environment and FITSIO

The Python environment has been exported to `SDC_full_dependencies.yml` which lists the Python dependencies.

Also in use is a version of [fitsio](https://pypi.org/project/fitsio/) that has been modified to allow arbitrarily large file sizes to be written directly to disk, and to allow subsets of data to be written to the file on disk (this functionality is currently unavailable in AstroPy). The modified version can be downloaded from [here](https://drive.google.com/drive/folders/15h0hE-cnqvS6xpX90qtX_Ji1wzC65V9R?usp=sharing). 

To set up the environment for running the pipeline, first create the Python environment via conda

`conda env create --file SDC_full_dependencies.yml`

Activate this environment

`conda activate SDC2`

Then install fitsio by navigating to the directory containing the modified source code, e.g.:

`cd /home/software/fitsio`

before installing via

`python setup.py install`

#### MIRIAD 

The MIRIAD package ([Sault et al., 2011](https://ui.adsabs.harvard.edu/abs/2011ascl.soft06007S/abstract)) is required for the telescope simulations. Tarballs and full installation instructions are available from [here](https://www.atnf.csiro.au/computing/software/miriad/INSTALL.html). For SKA simulations, a non-default 1024 antenna version is required. To compile this version, edit in the miriad source code the `miriad/inc/maxmid.h` and `miriad/inc/maxdimc.h` files to specify  `MAXANT = 1024`. Examples of these files can be found with a description of the steps taken to install a 1024 antenna version on Ubuntu 22 in `miriad_installation/`.

#### Source input catalogues

The pipeline uses catalogue files produced using [T-RECS](https://github.com/abonaldi/TRECS). The catalogues produced for use during SDC2 can be downloaded from [here](https://drive.google.com/drive/folders/15h0hE-cnqvS6xpX90qtX_Ji1wzC65V9R?usp=sharing). For use with these scripts, place the `SDC2_catalogues` directory inside the directory that contains this repository.

### Basic usage

`python run_SDC3_pipeline.py`

### Detailed usage

#### How to run the SDC2 simulation pipeline

An end-to-end pipeline, `run_SDC3_pipeline.py`, can be used to simulate the data products produced for SDC3. 


### Changes 

These scripts does not reproduce the exact same output as the data products used in SDC3. This is due to a change in the code to run in parallel that has modified the was that random seed is handles with respect to the version used for SDC3a.






