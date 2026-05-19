# Instructions on how to use COSP-RTTOV in CAM
### Jonah Shaw, 2026/02/18 (jonah.shaw@colorado.edu)

The document describes how to use COSP-RTTOV (Shaw et al., 2025) in CAM. COSP-RTTOV is an extension of the COSPv2.0 satellite simulator (Swales et al., 2018) that produces satellite-like spectral radiation fields (e.g. radiances and brightness temperatures) using the CESM state vector and the RTTOV v13.2 radiative transfer model (Saunders et al., 2018). COSP and COSP-RTTOV are open-source software available from the Cloud Feedback Model Intercomparison Project (CFMIP) via github (https://github.com/CFMIP/COSPv2.0).

Shaw, J. K., Swales, D. J., DeSouza-Machado, S., Turner, D. D., Kay, J. E., and Schneider, D. P.: COSP-RTTOV-1.0: flexible radiation diagnostics to enable new science applications in model evaluation, climate change detection, and satellite mission design, Geosci. Model Dev., 18, 4935–4950, https://doi.org/10.5194/gmd-18-4935-2025, 2025.

Swales, D. J., Pincus, R., and Bodas-Salcedo, A.: The Cloud Feedback Model Intercomparison Project Observational Simulator Package: Version 2, Geosci. Model Dev., 11, 77–81, https://doi.org/10.5194/gmd-11-77-2018, 2018.

Saunders, R., Hocking, J., Turner, E., Rayer, P., Rundle, D., Brunel, P., Vidot, J., Roquet, P., Matricardi, M., Geer, A., Bormann, N., and Lupu, C., 2018: An update on the RTTOV fast radiative transfer model (currently at version 12), Geosci. Model Dev., 11, 2717-2737, https://doi.org/10.5194/gmd-11-2717-2018.

## Installation of RTTOV
RTTOV is a free-to-use software, but the license requires each user to register to access the source code. Users should follow official instructions to access and download RTTOV v13.2 here: https://nwp-saf.eumetsat.int/site/software/rttov/download/#Software

Once you have downloaded RTTOV, follow RTTOV's instructions to compile the libraries. Building RTTOV can take ~30 minutes, so a pre-existing build is linked to CESM for each case using COSP-RTTOV. The file "compile_RTTOV.sh" can be used as a template to compile RTTOV and download example radiative transfer coefficient files. You may need to change your compile and update paths to existing NETCDF and HDF5 libraries. NETCDF and HDF5 are used to efficiently read RTTOV coefficient files. If these libraries are not available, RTTOV can also read coefficients from .txt and binary formats, but these files are large and may take up large amounts of memeory. See RTTOV v13.2 documentation for details.

## Linking RTTOV into COSP and CAM
Once RTTOV has been built, update Makefile.rttov to point to your installation and include appropriate flags. Specifically:
1. Set to $RTTOVDIR to point to your installation.
2. If needed, update FFLAGS and LDFLAGS_ARCH to match settings used to build RTTOV.
3. If needed, update HDF5_PREFIX, FFLAGS_HDF5, LDFLAGS_HDF5, NETCDF_PREFIX, FFLAGS_NETCDF, and LDFLAGS_NETCDF to match the settings used to build RTTOV.
4. If you don't want to force a compilation every time, comment out the line "_run_cmd(rttov_command, rttov_src_dir)" in the CAM cime_config/buildlib python script.

When building CAM/CESM, make the following changes to link RTTOV.
1. Set RTTOV as an environment variable (e.g. export RTTOV=TRUE, setenv RTTOV=TRUE)
2. Add COSP and RTTOV to CAM_CONFIG_OPTS before running ./case.setup by running './xmlchange --append CAM_CONFIG_OPTS="-cosp -rttov"'
3. Run './xmlchange FORCE_BUILD_SMP=TRUE' before running ./case.setup.

## Producing satellite-like outputs

### CAM namelist settings
COSP-RTTOV allows you to specify an arbitrary number of satellite instruments to simulate. Trigger COSP-RTTOV by editing the cam namelist file (user_nl_cam) within your case. For example:

cosp_lrttov_sim = .true.,
cosp_rttov_Ninstruments = 1,
cosp_rttov_instrument_namelists = '/glade/u/home/jonahshaw/Scripts/git_repos/COSPv2.0/driver/run/instrument_nls/cosp2_rttov_inst_CESM.txt',

These lines turn on RTTOV, specify that 1 instrument will be simulated, and provide an absolute path to a namelist file specifying how outputs should be produced. In the case of multiple instruments, "cosp_rttov_Ninstruments" should be increased and "cosp_rttov_instrument_namelists" should contain multiple strings with character limit 256.

### Instrument coefficients.
RTTOV performs fast radiative transfer calculations using instrument-specific coefficients derived from line-by-line calculations. To simulate an instrument, coefficients must be available and you must download them to be read by COSP-RTTOV. Clear-sky calculations use Optical Depth coefficients. All-sky calculations use cloud coefficients. COSP does not handle the radiative effects of aerosols. Coefficients are available at https://nwp-saf.eumetsat.int/site/software/rttov/download/coefficients/coefficient-download/ and example code to download them is at the end of "compile_RTTOV.sh".

### Instrument namelists.
Each instrument requires a namelist (example "cosp2_rttov_inst_CESM.txt") to specify:
1. The instrument itself (e.g. AIRS, IASI) and which trace gases are used in radiative transfer calculations.
Set "rttov_srcDir" to point to your RTTOV build or the parent directory of where your RTTOV coefficients are stored.
Set "rttov_coefDir" as the directory RTTOV coefficients are stored. "rttov_srcDir" and "rttov_coefDir" are used together when load coefficients.
Set "OD_coef_filepath" to the optical depth (clear-sky) RTTOV coefficients file. This file will be loaded using the combined path: $rttov_srcDir/$rttov_coefDir/$OD_coef_filepath
If you are performing all-sky calculations, set "cld_coef_filepath" to the RTTOV cloud coefficients file. This file will be loaded using the combined path: $rttov_srcDir/$rttov_coefDir/$OD_coef_filepath
Set the trace gas logicals (SO2_data, N2O_data, CO_data, CO2_data, CH4_data, ozone_data) to match your OD coefficient files. The coefficient file name indicates which trace gases are active.

2. Which fields are simulated (radiances, brightness temperatures, clear-sky vs. all-sky).
Set Lrttov_bt, Lrttov_rad, Lrttov_refl, and Lrttov_cld to determine which outputs will be produced. Note that SW and NIR reflectances have not been validated.
If "Lrttov_cld" is true, both all-sky and clear-sky fields will be computed. Otherwise only clear-sky fields will be computed.

3. Which instrument channels are simulated.
Set "nchannels_rec" to determine how many fields will be simulated.
If "Lchannel_filepath" is .false., then the first "nchannels_rec" will be simulated. If .true., then the entry for "channel_filepath" will be used to simulated requested channels.
Set "channel_filepath" as the absolute path to a .csv file following the format of "rttov_channel_input_AIRSL1C_subset.csv". The first column contains the indices of the channels that should be simulated. Channels are 1-indexed, and only the first "nchannels_rec" will be simulated if the .csv exceeds this length (in the example files, only the first 5 channels are simulated). Note that if you are using a coefficient file containing a subset of channels that channels will be selected from this subset.

4. How variables from CAM will be passed into RTTOV.
Some settings specify how the CAM state will be interpreted by RTTOV. Changing these fields is not recommended as it will produce physically-inconsistent output fields.
See: Lrttov_gridbox_cldmmr, rttov_gas_units, rttov_extendatmos
The radiative properties of clouds can use multiple parametrizations internal to RTTOV by changing [rttov_clw_scheme, rttov_ice_scheme, rttov_icede_param]. The values for these fields set in "cosp2_rttov_inst_CESM.txt" gave the best agreement with CAM's RRTMG-LW during validation of COSP-RTTOV using single-column simulations (see Shaw (2025)). Modifications are not advised. See RTTOV documentation for specifics.

5. Satellite viewing zenith angle.
"rttov_ZenAng" determines a constant SZA in units of degrees.

6. Simplified orbit swathing is set by adjusting: ["rttov_Nlocaltime", "rttov_localtime", "rttov_localtime_width"]. For example, to simulate a swathing pattern centered at 1:30am and 1:30pm local times with 1500km swath width, these fields would be:
  rttov_Nlocaltime=2,
  rttov_localtime=1.5,13.5,
  rttov_localtime_width=1500,1500,
An arbitrary number of swaths can be simulated. See Shaw (2025) for details.

*Remaining settings are mostly not relevant to running COSP-RTTOV in CAM. PC-RTTOV works but the memory output may practically be too large for CAM simulations. Thousands of channels can be efficiently produced from single-column experiments though.*
