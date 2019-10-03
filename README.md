# CAM: The Community Atmosphere Model

CAM code is stored in this repository on branches other than master.  The details are explained below.
Please see the [wiki](https://github.com/ESCOMP/CAM/wiki) for complete documentation on CAM.

## How to checkout and use CAM:

The instructions below assume you have cloned this repository and are in the repository directory. For example:
```
git clone https://github.com/ESCOMP/CAM
cd CAM
```

### To run CAM compatible with the CESM2.1 release:
```
git checkout cam_cesm2_1_rel_33
./manage_externals/checkout_externals
```
NOTE -- to view this branch in github, go to the "Branch:master" pulldown menu and select cam_cesm2_1_rel

### To use unsupported CAM **development** code:

## NOTE: This is **unsupported** development code and is subject to the [CESM developer's agreement](http://www.cgd.ucar.edu/cseg/development-code.html).
```
git checkout cam6_2_000
./manage_externals/checkout_externals
```
### CAM Documentation - https://ncar.github.io/CAM/doc/build/html/index.html

### CAM6 namelist settings - http://www.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html

