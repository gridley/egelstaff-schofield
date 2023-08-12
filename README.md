# Direct Egelstaff-Schofield Law Sampling Demo

This code is meant to accompany the paper submitted to Nuclear
Science and Engineering, "Direct Sampling of the Egelstaff-Schofield
Scattering Law and the Phonon Sampling Method for Liquids". In summary,
we show how to sample from the Egelstaff-Schofield law using a few exact
formulas, which allows the use of Trainer's phonon sampling method (PSM)
with liquids. The file sample_egelstaff.C includes C++ code that implements
our sampling techniques. The Jupyter notebook shows how to generate the
requisite input data using a modified version of NJOY and OpenMC, and
uses the C++ code to test our new method against conventional techniques
for sampling the thermal scattering law.

In order to run this code, you must have a functional copy of OpenMC installed,
and a copy of NJOY. These can be found respectively here:

https://github.com/openmc-dev/openmc

https://github.com/njoy/NJOY2016

Patches must be applied to each piece of software to carry out the calculations
presented in our paper. The git-apply command can be used to apply the patches
contained in this repo to both NJOY and OpenMC. This can be done by running the
following:

```
cd <openmc location>/openmc
git checkout 382bcb2e8ea63b5ffd9ead818c651e5aa9d4c7ce
git apply <location of this repo>/0001-read-debye-waller-lambda-from-modded-njoy.patch

cd <njoy location>/njoy
git checkout 715681a37164e797daa86e037590cf3ee0f955ea
git apply <location of this repo>/0001-leapr-patch.patch
```
you must then rebuild NJOY, and possibly reinstall the OpenMC python library.

The Jupyter notebook here shows how we convert from a vibrational density of states
to the scattering law, subtracting out the translational component. This allows
us to test our exact sampling technique for the Egelstaff-Schofield law. The Python
code in the Jupyter notebook calls our C++ implementations of the sampling routines.
To prepare the C++ code to be called by Python, the following must be ran:

```
g++ -fPIC -c sample_egelstaff.C
g++ sample_egelstaff.o -shared -o sample_egelstaff.so
```
