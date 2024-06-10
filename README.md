# CAM: The Community Atmosphere Model

CAM code is stored in this repository on branches other than main.  The details are explained below.

Please see the [wiki](https://github.com/ESCOMP/CAM/wiki) for complete documentation on CAM, getting started with git and how to contribute to CAM's development.

This repository has four officially-supported branches:
* main - contains this readme and the Code of Conduct information
* cam_cesm2_2_rel - contains the current CESM2.2 CAM code
* cam_cesm2_1_rel - contains the current CESM2.1 CAM code
* cam_development - contains the current CAM development code (see the important note below before using this branch)

There are other branches in this repo that have a `zz_` prefix which are used for special projects.  Please note that use of these branches is not-supported and may very likely produce un-scientific results in certain configurations.

## How to checkout and use CAM:

The instructions below assume you have cloned this repository and are in the repository directory. For example:
```
git clone https://github.com/ESCOMP/CAM
cd CAM
```
### To run CAM compatible with the CESM2.2 release:
```
git checkout cam_cesm2_2_rel_02
./manage_externals/checkout_externals
```
### To run CAM compatible with the CESM2.1 release:
```
git checkout cam_cesm2_1_rel_41
./manage_externals/checkout_externals
```
### To view the release branches in Github, go to the "Branch:main" pulldown menu and select the appropriate branch.

### To use unsupported CAM **development** code:

## NOTE: This is **unsupported** development code and is subject to the [CESM developer's agreement](https://www.cgd.ucar.edu/sections/cseg/policies).
```
git checkout cam6_3_162
./bin/git-fleximod update
```
### CAM Documentation - https://ncar.github.io/CAM/doc/build/html/index.html

### CAM6 namelist settings - https://docs.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html

