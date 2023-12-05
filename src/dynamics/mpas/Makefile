CPPFLAGS := -D_MPI -DMPAS_NATIVE_TIMERS -DMPAS_GIT_VERSION=unknown -DMPAS_NAMELIST_SUFFIX=atmosphere -DMPAS_CAM_DYCORE -DMPAS_PIO_SUPPORT -DMPAS_NO_ESMF_INIT
ifdef PIODEF
  CPPFLAGS += $(PIODEF)
endif
ifeq ($(strip $(COMP_INTERFACE)),nuopc)
  CPPFLAGS += -DMPAS_EXTERNAL_ESMF_LIB
endif

REGISTRY_FILE := $(MPAS_SRC_ROOT)/dycore/src/core_atmosphere/Registry.xml

VPATH := $(MPAS_SRC_ROOT)/dycore/src/core_atmosphere/diagnostics:$(MPAS_SRC_ROOT)/dycore/src/core_atmosphere/dynamics:
VPATH += $(MPAS_SRC_ROOT)/dycore/src/core_atmosphere:$(MPAS_SRC_ROOT)/dycore/src/tools/registry:$(MPAS_SRC_ROOT)/dycore/src/tools/input_gen:
VPATH += $(MPAS_SRC_ROOT)/dycore/src/operators:$(MPAS_SRC_ROOT)/dycore/src/framework:$(MPAS_SRC_ROOT)/dycore/src/external/ezxml:
VPATH += $(MPAS_SRC_ROOT)/driver

INTERFACE_OBJS = \
	mpas_atm_core.o \
	mpas_atm_core_interface.o \
	mpas_atm_dimensions.o \
	mpas_atm_threading.o \
	mpas_atm_halos.o \
	cam_mpas_subdriver.o

DYN_OBJS = \
	mpas_atm_time_integration.o \
	mpas_atm_boundaries.o \
	mpas_atm_iau.o

DIAG_OBJS = \
	mpas_atm_diagnostics_manager.o \
	mpas_atm_diagnostics_utils.o

DIAGNOSTICS = \
	mpas_atm_diagnostic_template.o \
	mpas_isobaric_diagnostics.o \
	mpas_cloud_diagnostics.o \
	mpas_convective_diagnostics.o \
	mpas_pv_diagnostics.o \
	mpas_soundings.o

REG_OBJS = \
	parse.o \
	dictionary.o \
	gen_inc.o \
	fortprintf.o \
	utility.o

STREAMS_GEN_OBJS = \
	streams_gen.o \
	test_functions.o

OPERATOR_OBJS = \
	mpas_vector_operations.o \
	mpas_matrix_operations.o \
	mpas_tensor_operations.o \
	mpas_rbf_interpolation.o \
	mpas_vector_reconstruction.o \
	mpas_spline_interpolation.o \
	mpas_tracer_advection_helpers.o \
	mpas_tracer_advection_mono.o \
	mpas_tracer_advection_std.o \
	mpas_geometry_utils.o

FRAME_OBJS = \
	mpas_kind_types.o \
	mpas_framework.o \
	mpas_timer.o \
	mpas_timekeeping.o \
	mpas_constants.o \
	mpas_attlist.o \
	mpas_hash.o \
	mpas_sort.o \
	mpas_block_decomp.o \
	mpas_block_creator.o \
	mpas_dmpar.o \
	mpas_abort.o \
	mpas_decomp.o \
	mpas_threading.o \
	mpas_io.o \
	mpas_io_streams.o \
	mpas_bootstrapping.o \
	mpas_io_units.o \
	mpas_stream_manager.o \
	mpas_stream_list.o \
	mpas_forcing.o \
	mpas_c_interfacing.o \
	random_id.o \
	pool_hash.o \
	mpas_derived_types.o \
	mpas_domain_routines.o \
	mpas_field_routines.o \
	mpas_pool_routines.o \
	xml_stream_parser.o \
	regex_matching.o \
	mpas_log.o \
	mpas_halo.o \
	mpas_string_utils.o


UTIL_OBJS = \
	ezxml.o


all:
	( $(MAKE) xml )
	( $(MAKE) framework )
	( $(MAKE) operators )
	( $(MAKE) registry )
	( $(MAKE) streams_gen )
	( $(MAKE) incs )
	( $(MAKE) dycore )
	( ar -ru libmpas.a $(INTERFACE_OBJS) $(DYN_OBJS) $(DIAG_OBJS) $(DIAGNOSTICS) $(OPERATOR_OBJS) $(FRAME_OBJS) $(UTIL_OBJS) )


#
# EZXML
#
xml: ezxml.o


#
# Framework dependencies
#
framework: $(FRAME_OBJS)

mpas_framework.o: mpas_dmpar.o \
                  mpas_io.o \
                  mpas_derived_types.o \
                  mpas_domain_routines.o \
                  mpas_field_routines.o \
                  mpas_pool_routines.o \
                  mpas_timer.o \
                  mpas_sort.o \
                  mpas_io_units.o \
                  mpas_block_decomp.o \
                  mpas_stream_manager.o \
                  mpas_c_interfacing.o \
                  mpas_halo.o

mpas_abort.o: mpas_kind_types.o mpas_io_units.o mpas_threading.o

mpas_constants.o: mpas_kind_types.o

mpas_log.o: mpas_derived_types.o mpas_io_units.o mpas_abort.o mpas_threading.o mpas_c_interfacing.o

mpas_attlist.o: mpas_kind_types.o mpas_io_units.o mpas_derived_types.o

mpas_derived_types.o: mpas_kind_types.o mpas_constants.o $(TYPE_DEPS)

mpas_domain_routines.o: mpas_derived_types.o mpas_pool_routines.o mpas_dmpar.o

mpas_field_routines.o: mpas_derived_types.o duplicate_field_array.inc duplicate_field_scalar.inc mpas_threading.o mpas_attlist.o

mpas_pool_routines.o: mpas_derived_types.o mpas_field_routines.o mpas_threading.o mpas_log.o

mpas_decomp.o: mpas_derived_types.o mpas_stream_manager.o mpas_log.o

mpas_hash.o : mpas_derived_types.o

mpas_dmpar.o: mpas_sort.o mpas_kind_types.o mpas_derived_types.o mpas_hash.o mpas_threading.o mpas_pool_routines.o mpas_log.o

mpas_sort.o: mpas_kind_types.o mpas_log.o

mpas_timekeeping.o: mpas_string_utils.o mpas_kind_types.o mpas_derived_types.o mpas_dmpar.o mpas_threading.o mpas_log.o

mpas_timer.o: mpas_kind_types.o mpas_dmpar.o mpas_threading.o mpas_log.o

mpas_block_decomp.o: mpas_derived_types.o mpas_hash.o mpas_io_units.o mpas_dmpar.o

mpas_block_creator.o: mpas_dmpar.o mpas_hash.o mpas_sort.o mpas_io_units.o mpas_block_decomp.o mpas_stream_manager.o mpas_decomp.o mpas_abort.o $(DEPS)

mpas_io.o: mpas_dmpar.o mpas_attlist.o mpas_log.o

mpas_io_streams.o: mpas_attlist.o mpas_derived_types.o mpas_timekeeping.o mpas_io.o mpas_pool_routines.o add_field_indices.inc mpas_log.o $(DEPS)

mpas_bootstrapping.o: mpas_derived_types.o mpas_dmpar.o mpas_block_decomp.o mpas_block_creator.o mpas_sort.o mpas_timekeeping.o mpas_io_streams.o mpas_stream_manager.o random_id.o mpas_log.o $(DEPS)

mpas_io_units.o: mpas_kind_types.o

mpas_threading.o: mpas_kind_types.o

mpas_stream_list.o: mpas_derived_types.o mpas_kind_types.o mpas_io_streams.o mpas_timekeeping.o regex_matching.o mpas_log.o

mpas_stream_manager.o: mpas_io_streams.o mpas_timekeeping.o mpas_derived_types.o mpas_kind_types.o mpas_c_interfacing.o mpas_stream_list.o mpas_dmpar.o mpas_io.o mpas_threading.o mpas_log.o

mpas_forcing.o: mpas_derived_types.o mpas_timekeeping.o mpas_stream_manager.o mpas_log.o mpas_io_units.o

mpas_halo.o: mpas_derived_types.o mpas_pool_routines.o mpas_log.o

#
# Operator dependencies
#
operators: $(OPERATOR_OBJS)

mpas_tensor_operations.o: mpas_vector_operations.o mpas_matrix_operations.o

mpas_rbf_interpolation.o: mpas_vector_operations.o

mpas_vector_reconstruction.o: mpas_rbf_interpolation.o

mpas_tracer_advection_helpers.o: mpas_geometry_utils.o

mpas_tracer_advection_mono.o: mpas_tracer_advection_helpers.o

mpas_tracer_advection_std.o: mpas_tracer_advection_helpers.o

mpas_geometry_utils.o: mpas_vector_operations.o mpas_matrix_operations.o


#
# Registry
#
registry: $(REG_OBJS) ezxml.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $(REG_OBJS) ezxml.o


#
# Default stream file generator
#
streams_gen: $(STREAMS_GEN_OBJS) ezxml.o
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ $(STREAMS_GEN_OBJS) ezxml.o


#
# Include files manufactured by the registry code generator
#
incs: $(REGISTRY_FILE)
	( cpp -P -traditional $(CPPFLAGS) -I$(MPAS_SRC_ROOT)/dycore/src/core_atmosphere/diagnostics $(REGISTRY_FILE) > Registry_processed.xml )
	( ./registry Registry_processed.xml )

#
# Dycore
#
dycore: $(DYN_OBJS) $(DIAGNOSTICS) $(DIAG_OBJS) $(INTERFACE_OBJS)

mpas_atm_time_integration.o: mpas_atm_iau.o mpas_atm_dimensions.o mpas_atm_boundaries.o


#
# Diagnostics
#
mpas_atm_diagnostics_manager.o: mpas_atm_diagnostics_utils.o $(DIAGNOSTICS)

mpas_cloud_diagnostics.o: mpas_atm_diagnostics_utils.o
mpas_convective_diagnostics.o: mpas_atm_diagnostics_utils.o
mpas_isobaric_diagnostics.o: mpas_atm_diagnostics_utils.o
mpas_pv_diagnostics.o: mpas_atm_diagnostics_utils.o


#
# Core interface objects
#
mpas_atm_core_interface.o: mpas_atm_core.o incs

mpas_atm_core.o: mpas_atm_threading.o mpas_atm_time_integration.o mpas_atm_diagnostics_manager.o mpas_atm_halos.o

cam_mpas_subdriver.o: mpas_atm_core_interface.o mpas_derived_types.o mpas_framework.o mpas_domain_routines.o mpas_pool_routines.o


clean:
	rm registry Registry_processed.xml *.inc *.o *.mod lib*.a

# Cancel the built-in implicit rule for Modula-2 files (.mod) to avoid having 'make'
# try to create .o files from Fortran .mod files
%.o : %.mod

%.o : %.c
	$(CC) -c $(CFLAGS) -DUNDERSCORE $(CPPFLAGS) $(CPPINCLUDES) -I$(MPAS_SRC_ROOT)/dycore/src/external/ezxml $<

%.o : %.F
	$(FC) -c $(CPPFLAGS) $(FFLAGS) $(CPPINCLUDES) $(FCINCLUDES) $<

%.o : %.F90
	$(FC) -c $(CPPFLAGS) $(FFLAGS) $(CPPINCLUDES) $(FCINCLUDES) $<

