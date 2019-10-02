# CAM: The Community Atmosphere Model

CAM Documentation - https://ncar.github.io/CAM/doc/build/html/index.html

CAM6 namelist settings - http://www.cesm.ucar.edu/models/cesm2/settings/current/cam_nml.html

**CAM Code coming soon**

Please see the [wiki](https://github.com/ESCOMP/CAM/wiki) for information.

## How to checkout and use CAM:

The instructions below assume you have cloned this repository and are in the repository directory. For example:
```
git clone https://github.com/ESCOMP/CAM
cd CAM
```

### To run CAM compatible with the CESM2.1 release:
```
git checkout cam_cesm2_1_rel_32
./manage_externals/checkout_externals
```

### To use unsupported CAM **development** code:

## NOTE: This is **unsupported** development code and is subject to the [CESM developer's agreement](http://www.cgd.ucar.edu/cseg/development-code.html).
```
git checkout cam_development
./manage_externals/checkout_externals
```
