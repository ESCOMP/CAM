"""
build_cam_namelists.py

Using MPAS-Atmosphere's Registry.xml file as input, generates content to be included
in CAM's namelist definition files:
* namelist_definitions.mpas should be included in cam/bld/namelist_files/namelist_definition.xml
* namelist_defaults.mpas should be included in cam/bld/namelist_files/namelist_defaults_cam.xml
* build-namelist.mpas should be included in cam/bld/build-namelist
* namelist_reads.mpas should be included in cam/src/dynamics/mpas/dyn_comp.F90
"""

import xml.etree.ElementTree as ET
import re


def cam_namespace(config):
	return re.sub('config_', 'mpas_', config)


def char_to_string(nmltype):
	return re.sub('character', 'char*512', nmltype)


def build_namelists():
	tree = ET.parse('Registry.xml')
	root = tree.getroot()


	#
	# Section to be included in CAM namelist_definitions.xml file
	#
	f = open('namelist_definitions.mpas', 'w')
	for child in root:
		if child.tag == 'nml_record':
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					f.write('<entry id="{}" type="{}" category="mpas" group="{}" valid_values="">\n'.format(cam_namespace(nmlopt.attrib['name']), char_to_string(nmlopt.attrib['type']), child.attrib['name']))
					f.write('{}\n'.format(nmlopt.attrib['description']))
					f.write('Default: {}\n'.format(nmlopt.attrib['default_value']))
					f.write('</entry>\n\n')
	f.close()


	#
	# Section to be included in CAM namelist_defaults.xml file
	#
	f = open('namelist_defaults.mpas', 'w')
	for child in root:
		if child.tag == 'nml_record':
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					f.write('<{}>{}</{}>\n'.format(cam_namespace(nmlopt.attrib['name']), nmlopt.attrib['default_value'], cam_namespace(nmlopt.attrib['name'])))
	f.close()


	#
	# Section to be included in CAM bulid-namelist file
	#
	f = open('build-namelist.mpas', 'w')
	for child in root:
		if child.tag == 'nml_record':
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					f.write('    add_default($nl, \'{}\');\n'.format(cam_namespace(nmlopt.attrib['name'])));
	f.close()


	#
	# Auto-generated Fortran for reading namelists
	#
	f = open('namelist_reads.mpas', 'w')
	f.write('!-----------------------------------------------------------------------\n')
	f.write('!  routine cam_mpas_namelist_read\n')
	f.write('!\n')
	f.write('!> \\brief Reads MPAS-A dycore namelists and adds the namelists to the MPAS configPool\n')
	f.write('!> \details\n')
	f.write('!>  Given the name of a file containing namelists and an MPAS pool, reads the dycore\n')
	f.write('!>  namelists from that file and adds the namelist options to the pool.\n')
	f.write('!\n')
	f.write('!>  Only the CAM masterproc actually opens and reads from the specified file. Upon return,\n')
	f.write('!>  if no errors were encountered, all MPI ranks have valid namelists in their configPool.\n')
	f.write('!\n')
	f.write('!>  A value of zero is returned if no errors were encountered, and a non-zero value is returned\n')
	f.write('!>  if any errors were encountered in reading the namelist file.\n')
	f.write('!\n')
	f.write('!   WARNING: This routine was auto-generated based on the MPAS-Atmosphere Registry.xml file\n')
	f.write('!\n')
	f.write('!            Rather than editing this function directly, edit core_atmosphere/Registry.xml\n')
	f.write('!            and re-run the build_cam_namelist.py script to regenerate this function.\n')
	f.write('!\n')
	f.write('!-----------------------------------------------------------------------\n')
	f.write('function cam_mpas_namelist_read(namelistFilename, configPool) result(ierr)\n')
	f.write('\n')
	f.write('   use units, only : getunit, freeunit\n')
	f.write('   use spmd_utils, only : mpicom, masterproc, masterprocid, &\n')
	f.write('                          mpi_integer, mpi_real8,  mpi_logical, mpi_character, mpi_success\n')
	f.write('   use shr_kind_mod, only : shr_kind_r8\n')
	f.write('   use cam_logfile, only : iulog\n')
	f.write('   use namelist_utils, only : find_group_name\n')
	f.write('\n')
	f.write('   use mpas_derived_types, only : mpas_pool_type\n')
	f.write('   use mpas_kind_types, only : StrKIND\n')
	f.write('   use mpas_pool_routines, only : mpas_pool_add_config\n')
	f.write('\n')
	f.write('   implicit none\n')
	f.write('\n')
	f.write('   character(len=*), intent(in) :: namelistFilename\n')
	f.write('   type (mpas_pool_type), intent(inout) :: configPool\n')
	f.write('\n')
	f.write('   integer :: ierr   ! Return value\n')
	f.write('\n')
	f.write('   integer :: unitNumber\n')
	f.write('\n')
	f.write('   integer :: mpi_ierr\n')
	
	f.write('\n')
	for child in root:
		if child.tag == 'nml_record':
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					if nmlopt.attrib['type'] == 'integer':
						f.write('   integer                 :: {} = {}\n'.format(cam_namespace(nmlopt.attrib['name']), nmlopt.attrib['default_value']))
					elif nmlopt.attrib['type'] == 'real':
						f.write('   real (kind=shr_kind_r8) :: {} = {}\n'.format(cam_namespace(nmlopt.attrib['name']), nmlopt.attrib['default_value']))
					elif nmlopt.attrib['type'] == 'character':
						f.write('   character (len=StrKIND) :: {} = \'{}\'\n'.format(cam_namespace(nmlopt.attrib['name']), nmlopt.attrib['default_value']))
					elif nmlopt.attrib['type'] == 'logical':
						f.write('   logical                 :: {} = .{}.\n'.format(cam_namespace(nmlopt.attrib['name']), nmlopt.attrib['default_value']))

	f.write('\n')
	for child in root:
		if child.tag == 'nml_record' and 'used_by' in child.attrib and re.search('cam', child.attrib['used_by']):
			first = True
			f.write('   namelist /{}/ &\n'.format(child.attrib['name']))
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					if not first:
						f.write(', &\n')
					f.write('           {}'.format(cam_namespace(nmlopt.attrib['name'])))
					first = False
			f.write('\n\n')
	
	f.write('   if (masterproc) then\n')
	f.write('      write(iulog,*) \'Reading MPAS-A dycore namelist from \', trim(namelistFilename)\n')
	f.write('      unitNumber = getunit()\n')
	f.write('      open(unit=unitNumber, file=trim(namelistFilename), status=\'old\', form=\'formatted\')\n')
	f.write('   end if\n')
	
	for child in root:
		if child.tag == 'nml_record' and 'used_by' in child.attrib and re.search('cam', child.attrib['used_by']):
			f.write('\n')
			f.write('   !\n')
			f.write('   ! Read namelist group &{}\n'.format(child.attrib['name']))
			f.write('   !\n')
			f.write('   if (masterproc) then\n')
			f.write('      rewind(unitNumber)\n')
			f.write('      call find_group_name(unitNumber, \'{}\', status=ierr)\n'.format(child.attrib['name']))
			f.write('      if (ierr == 0) then\n')
			f.write('         read(unitNumber, {}, iostat=ierr)\n'.format(child.attrib['name']))
			f.write('      else\n')
			f.write('         close(unit=unitNumber)\n')
			f.write('         call freeunit(unitNumber)\n')
			f.write('      end if\n')
			f.write('   end if\n')
			f.write('   call mpi_bcast(ierr, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)\n');
			f.write('   if (ierr /= 0) then\n')
			f.write('      if (masterproc) then\n')
			f.write('         write(iulog,*) \'Failed to read namelist group &{}\'\n'.format(child.attrib['name']))
			f.write('      end if\n')
			f.write('      return\n')
			f.write('   end if\n')
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					if nmlopt.attrib['type'] == 'integer':
						f.write('   call mpi_bcast({}, 1, mpi_integer, masterprocid, mpicom, mpi_ierr)\n'.format(cam_namespace(nmlopt.attrib['name'])));
					elif nmlopt.attrib['type'] == 'real':
						f.write('   call mpi_bcast({}, 1, mpi_real8, masterprocid, mpicom, mpi_ierr)\n'.format(cam_namespace(nmlopt.attrib['name'])));
					elif nmlopt.attrib['type'] == 'character':
						f.write('   call mpi_bcast({}, StrKIND, mpi_character, masterprocid, mpicom, mpi_ierr)\n'.format(cam_namespace(nmlopt.attrib['name'])));
					elif nmlopt.attrib['type'] == 'logical':
						f.write('   call mpi_bcast({}, 1, mpi_logical, masterprocid, mpicom, mpi_ierr)\n'.format(cam_namespace(nmlopt.attrib['name'])));
					f.write('   if (mpi_ierr /= mpi_success) then\n');
					f.write('      if (masterproc) then\n');
					f.write('         write(iulog,*) \'MPI_Bcast failed for namelist option {}\'\n'.format(cam_namespace(nmlopt.attrib['name'])));
					f.write('      end if\n');
					f.write('      return\n');
					f.write('   end if\n');
	
			f.write('\n')
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					f.write('   call mpas_pool_add_config(configPool, \'{}\', {})\n'.format(nmlopt.attrib['name'], cam_namespace(nmlopt.attrib['name'])))

	f.write('\n')
	f.write('   if (masterproc) then\n')
	f.write('      close(unit=unitNumber)\n')
	f.write('      call freeunit(unitNumber)\n')
	f.write('   end if\n')
	f.write('\n')
	f.write('   if (masterproc) then\n')
	f.write('      write(iulog,*) \'MPAS-A dycore configuration:\'\n');
	for child in root:
		if child.tag == 'nml_record' and 'used_by' in child.attrib and re.search('cam', child.attrib['used_by']):
			for nmlopt in child:
				if 'used_by' in nmlopt.attrib and re.search('cam', nmlopt.attrib['used_by']):
					if nmlopt.attrib['type'] == 'character':
						f.write('      write(iulog,*) \'   {} = \', trim({})\n'.format(cam_namespace(nmlopt.attrib['name']), cam_namespace(nmlopt.attrib['name'])))
					else:
						f.write('      write(iulog,*) \'   {} = \', {}\n'.format(cam_namespace(nmlopt.attrib['name']), cam_namespace(nmlopt.attrib['name'])))
	f.write('   end if\n')
	f.write('\n')
	f.write('end function cam_mpas_namelist_read\n')
	f.close()

if __name__ == "__main__":
    build_namelists()
