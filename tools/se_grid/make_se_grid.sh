#! /bin/bash

##
## General syntax help function
## Usage: help <exit status>
##
help () {
  local sname="`basename ${0}`"
  local hname="Usage: ${sname}"
  local hprefix="`echo ${hname} | tr '[!-~]' ' '`"
  echo "${hname} [ --ne <ne> ]"
  echo "${hprefix} [ --pg <#> ] physgrid FVM size (default no physgrid)"
  echo "${hprefix} [ --filename <new SE grid name> ]"
  echo "${hprefix} [ --np <#> ] Set number of GLL points (default 4)"
  echo "${hprefix} [ --fine-ne <#> ] Effective NE for refined mesh files"
  echo "${hprefix} [ --mesh-file <mesh filename> ]"
  echo "${hprefix} [ --mesh-name <mesh name> ]"
  echo "${hprefix} [ --root <cam_root_dir> ]"
  echo "${hprefix} [ --project <project account number> ]"
  if [ $1 -eq 0 ]; then
    echo ""
    echo "${sname} is used to create SCRIP and domain files for new SE grids."
    echo -e "\nTo create a normal (fixed size) GLL 1 degree grid:"
    echo "${sname} --ne 30"
    echo -e "\nTo create a 1 degree grid with a 3x3 FVM physics grid:"
    echo "${sname} --ne 30 --pg 3"
    echo -e "\nTo create a refined mesh grid with refinement from 1 to 1/4 deg:"
    echo "${sname} --ne 30 --fine-ne 120 --mesh-name <RAREA> --mesh-file <fname>.g"
    echo "  where <fname>.g is the mesh file read by the SE dycore and"
    echo "  <RAREA> is the name of the refined area (e.g., CONUS)"
    echo -e "\nNote that before running this script, you must enter the new"
    echo "grid in two places:"
    echo "In <root>/components/cam/bld/config_files/horiz_grid.xml and"
    echo "<root>/cime/config/cesm/config_grids.xml"
  fi
  exit $1
}

##
## Error output function (should be handed a string)
##
perr() {
  echo -e "\nERROR: ${@}\n"
  help 1
}

##
## Calculate the number of points in a CAM-SE grid
## Inputs are <ne> <pg> <np>
##
npoints() {
  local np
  local npt

  if [ $2 -gt 0 ]; then
    npt=$(( $1 * $1 * $2 * $2 * 6 ))
  else
    np=$(( $3 - 1 ))
    npt=$(( $1 * $1 * $np * $np * 6 + 2 ))
  fi
  echo "${npt}"
}

##
## Calculate the approximate resolution of a CAM-SE grid
## Inputs are <ne> <np>
##
resolution() {
  local np
  local nres
  local ndeg

  np=$(( $2 - 1 ))
  nres=$(( $1 * $np * 4 ))
  if [ $nres -gt 360 ]; then
    # Fraction of a degree, assume 1/integer
    nres=$(( $nres / 360 ))
    echo "1/${nres}"
  else
    # Resolution of less than a degree, round up
    ndeg=$(( 360 + ( $nres / 2 ) ))
    nres=$(( $ndeg / $nres ))
    echo "${nres}"
  fi
}

##
## Calculate the number of points in a CAM-SE refined mesh grid
## Inputs are <meshfile> <np>
##
nptRefined() {
  local np
  local nelem
  local npt=0
  # Check for ncdump
  if [ -n "`type ncdump 2> /dev/null`" ]; then
    np=$(( $2 - 1 ))
    nelem=`ncdump -h ${1} | grep num_elem | head -1 | cut -f3 -d' '`
    npt=$(( $nelem * $np * $np + 2 ))
  fi
  echo "${npt}"
}

##
## Return 0 status iff $1 is found in $2
##
find_in_file() {
  local fstr="`grep ${1} ${2} 2> /dev/null`"
  if [ -n "${fstr}" ]; then
    return 0
  else
    return 1
  fi
}

camroot=""
domaingrid=""
filename=""
finene=0
gridalias=""
gridname=""
meshfile=""
meshname=""
ne=0
np=4
pg=0
project=""
scriptroot=""

## Process our input arguments
while [ $# -gt 0 ]; do
  case $1 in
    --h | -h | --help | -help)
      help 0
      ;;
    --ne | -ne)
      if [ $# -lt 2 ]; then
        perr "${1} requires NE, the number of elements across a cube edge"
      fi
      ne=${2}
      shift
      ;;
    --np | -np)
      if [ $# -lt 2 ]; then
        perr "${1} requires NP, the number of gll points across an element edge"
      fi
      np=${2}
      shift
      ;;
    --pg | -pg)
      if [ $# -lt 2 ]; then
        perr "${1} requires PG, the number of FVM cols across an element edge"
      fi
      pg=${2}
      shift
      ;;
    --fine-ne | -fine-ne)
      if [ $# -lt 2 ]; then
        perr "${1} requires PG, the number of FVM cols across an element edge"
      fi
      finene=${2}
      shift
      ;;
    --filename | -filename)
      if [ $# -lt 2 ]; then
        perr "${1} requires the name of the output SCRIP file"
      fi
      filename="${2}"
      shift
      ;;
    --mesh-file | -mesh-file)
      if [ $# -lt 2 ]; then
        perr "${1} requires the name of the refined grid mesh file"
      fi
      meshfile="${2}"
      shift
      ;;
    --mesh-name | -mesh-name)
      if [ $# -lt 2 ]; then
        perr "${1} requires a name for this refined grid"
      fi
      meshname="${2}"
      shift
      ;;
    --project | -project)
      if [ $# -lt 2 ]; then
        perr "${1} requires a project number"
      fi
      project="--project ${2}"
      shift
      ;;
    --root | -root)
      if [ $# -lt 2 ]; then
        perr "${1} requires a cam root directory"
      fi
#  We need to make sure that the install directory is a full path
#  First, we see if it looks like a full path (not sure Windows will like this)
      case $1 in
        /*)
          camroot="$2"
        ;;
        *)
          camroot="`(cd $2; pwd -P)`"
          retval=$?
          if [ $retval -ne 0 ]; then
            perr "CAM root must exist"
          fi
      esac
      if [ ! -d "${camroot}" ]; then
        perr "The specified CAM directory, \"${2}\", does not exist."
        exit 1
      fi
      if [ ! -d "${camroot}/cime/scripts" ]; then
        perr "The specified CAM directory, \"${2}\", does not appear to be a CAM root."
        exit 1
      fi
      shift
      ;;
    *)
      perr "Unrecognized option, \"${1}\""
      ;;
  esac
  shift
done

scriptdir="$( cd $( dirname $0 ); pwd -P )"
trialroot="$( dirname $( dirname $( dirname ${scriptdir} ) ) )"
if [ -z "${camroot}" ]; then
  # Are we already in a CAM root or CIME scripts directory?
  if [ -f "./create_newcase" ]; then
    camroot="$( dirname $( dirname $( pwd -P ) ) )"
    if [ ! -d "${camroot}/components/cam" -o ! -d "${camroot}/cime" ]; then
      perr "We seem to be in a CIME directory but not in a CAM sandbox"
    fi
  elif [ -d "./components/cam" -a -d "./cime" ]; then
    camroot="`pwd -P`"
  elif [ -d "${trialroot}/components/cam" -a -d "${trialroot}/cime" ]; then
    camroot="${trialroot}"
  else
    perr "We don't seem to be in a CAM sandbox, use --root option"
  fi
fi

scriptroot="${camroot}/cime/scripts"

# Some basic checks
if [ $np -le 0 ]; then
  perr "NP must be a positive integer"
fi
if [ $ne -le 0 ]; then
  perr "Need an --ne <#> entry to construct a new grid filename"
fi

# Construct a gridname
if [ $pg -gt 0 ]; then
  if [ -n "${meshfile}" ]; then
    perr "FVM physics grid not supported for refined mesh grids"
  fi
  if [ -n "${meshname}" ]; then
    perr "FVM physics grid not supported for refined mesh grids"
  fi
  if [ $finene -gt 0 ]; then
    perr "FVM physics grid not supported for refined mesh grids"
  fi
  gridname="ne${ne}np${np}.pg${pg}"
  domaingrid="ne${ne}pg${pg}"
  ncol=`npoints $ne $pg $np`
  resolution="`resolution ${ne} ${np}`"
  if [ -n "`echo ${resolution} | grep '/'`" ]; then
    mask="tx0.1v2"
    malias="mt12"
  elif [ $resolution -gt 2 ]; then
    mask="gx3v7"
    malias="mg37"
  else
    mask="gx1v7"
    malias="mg17"
  fi
  gridalias="${domaingrid}_${domaingrid}_${malias}"
  desc="${gridname} is a Spectral Elem ${resolution}-deg grid with a ${pg}x${pg} FVM physics grid:"
elif [ -n "${meshname}" -o -n "${meshfile}" -o $finene -gt 0 ]; then
  # Refined mesh, make sure we have all the info
  if [ $finene -le 0 ]; then
    perr "Refined mesh grids require a --fine-ne entry"
  fi
  if [ -z "${meshfile}" ]; then
    perr "Refined mesh grids require a --mesh-file entry"
  fi
  if [ -z "${meshname}" ]; then
    perr "Refined mesh grids require a --mesh-name entry"
  fi
  gridname="ne0np${np}${meshname}.ne${ne}x$(( $finene / $ne ))"
  domaingrid="ne0${meshname}ne${ne}x$(( $finene / $ne ))"
  resolution="`resolution ${finene} ${np}`"
  if [ -n "`echo ${resolution} | grep '/'`" ]; then
    mask="tx0.1v2"
    malias="mt12"
  elif [ $resolution -gt 2 ]; then
    mask="gx3v7"
    malias="mg37"
  else
    mask="gx1v7"
    malias="mg17"
  fi
  gridalias="${domaingrid}_${domaingrid}_${malias}"
  ncol=`nptRefined ${meshfile} ${np}`
  desc="${gridname} is a Spectral Elem refined mesh grid:"
else
  # Standard GLL grid
  gridname="ne${ne}np${np}"
  domaingrid="ne${ne}"
  ncol=`npoints $ne $pg $np`
  resolution="`resolution ${ne} ${np}`"
  if [ -n "`echo ${resolution} | grep '/'`" ]; then
    mask="tx0.1v2"
    malias="mt12"
  elif [ $resolution -gt 2 ]; then
    mask="gx3v7"
    malias="mg37"
  else
    mask="gx1v7"
    malias="mg17"
  fi
  gridalias="${domaingrid}_${domaingrid}_${malias}"
  ncol=`npoints $ne $pg $np`
  desc="${gridname} is a Spectral Elem ${resolution}-deg grid:"
fi

datestr="`date '+%y%m%d'`"
if [ -z "${filename}" ]; then
  # We were not passed a filename so we have to construt one
  filename="${gridname}_SCRIP_${datestr}.nc"
fi

echo "gridname = ${gridname}"
echo "domaingrid = ${domaingrid}"
echo "gridalias = ${gridalias}"
echo "ncol = ${ncol}"
echo "filename = ${filename}"

## Check for an entry for this grid in horiz_grid.xml
hgridfile="${camroot}/components/cam/bld/config_files/horiz_grid.xml"

if [ ! -f "${hgridfile}" ]; then
  echo $hgridfile
  perr "Cannot find horiz_grid.xml -- THIS SHOULD NOT HAPPEN!"
fi

if ! find_in_file "hgrid=\\\"${gridname}\\\"" "${hgridfile}"; then
  echo "Did not find ${gridname} in `basename ${hgridfile}`"
  echo -e "Add the following grid entry to ${hgridfile}:\n"
  echo "<horiz_grid dyn=\"se\" hgrid=\"${gridname}\" ncol=\"${ncol}\" csne=\"${ne}\" csnp=\"${np}\" npg=\"${pg}/>"
  echo -e "\nThen repeat command"
  exit 1
fi

cgridfile="${camroot}/cime/config/cesm/config_grids.xml"
find_in_file "alias=\\\"${gridalias}\\\"" ${cgridfile}
resa=$?
find_in_file "domain name=\\\"${gridname}\\\"" ${cgridfile}
resd=$?

if [ $resa -ne 0 -o $resd -ne 0 ]; then
  if [ $resa -ne 0 ]; then
    echo -e "\nDid not find a model_grid entry for ${gridname}"
    echo -e "Add the following model_grid entry to ${cgridfile}:\n"
    echo "<model_grid alias=\"${gridalias}\" not_compset=\"_POP|_CLM\">"
    echo "  <grid name=\"atm\">${gridname}</grid>"
    echo "  <grid name=\"lnd\">${gridname}</grid>"
    echo "  <grid name=\"ocnice\">${gridname}</grid>"
    echo "  <mask>${mask}</mask>"
    echo "</model_grid>"
  fi
  if [ $resd -ne 0 ]; then
    echo -e "\nDid not find a domain entry for ${gridname}"
    echo -e "\nAdd the following domain entry to ${cgridfile}:\n"
    echo "<domain name=\"${gridname}\">"
    echo "  <nx>${ncol}</nx> <ny>1</ny>"
    echo "  <file grid=\"atm|lnd\" mask=\"${mask}\">\$DIN_LOC_ROOT/share/domains/domain.lnd.${gridname}_${mask}.${datestr}.nc</file>"
    echo "  <file grid=\"ocnice\"  mask=\"${mask}\">\$DIN_LOC_ROOT/share/domains/domain.ocn.${gridname}_${mask}_${datestr}.nc</file>"
    echo "  <desc>${desc}</desc>"
    echo "</domain>"
  fi
  echo -e "\nThen repeat command"
  exit 1
fi

# If we get this far, it is time to create the files
cset="FADIAB"
if [ -d "/glade/scratch" ]; then
  scrdir="/glade/scratch/${USER}"
elif [ -d "/scratch/cluster" ]; then
  scrdir="/scratch/cluster/${USER}"
else
  scrdir="`pwd -P`"
fi
if [ ! -d "${scrdir}" ]; then
  mkdir -p ${scrdir}
fi
cdir="${scrdir}/${cset}_${datestr}"
rm -rf ${cdir}
cd ${scriptroot}
res=$?
if [ $res -ne 0 ]; then
  perr "Unable to cd to script root, ${scriptroot}"
fi
./create_newcase --case ${cdir} --compset ${cset} --run-unsupported --res ${gridalias} ${project} --walltime 00:20:00
res=$?
if [ $res -ne 0 ]; then
  perr "create_newcase FAILED"
fi
cd ${cdir}
res=$?
if [ $res -ne 0 ]; then
  perr "Unable to change to case directory, ${cdir}"
fi
./xmlchange DOUT_S=FALSE,DEBUG=FALSE,ATM_GRID=${gridname}
res=$?
if [ $res -ne 0 ]; then
  perr "xmlchange 1 FAILED"
fi
./xmlchange --append CAM_CONFIG_OPTS=-analytic_ic
res=$?
if [ $res -ne 0 ]; then
  perr "xmlchange 2 FAILED"
fi
./xmlchange STOP_OPTION=nsteps,STOP_N=1
res=$?
if [ $res -ne 0 ]; then
  perr "xmlchange 3 FAILED"
fi
./case.setup
res=$?
if [ $res -ne 0 ]; then
  perr "case.setup FAILED"
fi
echo "se_write_grid_file = 'SCRIP'" >> user_nl_cam
res=$?
if [ $res -ne 0 ]; then
  perr "Unable to add se_write_grid_file to user_nl_cam"
fi
echo "se_grid_filename = '${domaingrid}_scrip_${datestr}.nc'" >> user_nl_cam
res=$?
if [ $res -ne 0 ]; then
  perr "Unable to add se_grid_filename to user_nl_cam"
fi
echo "ncdata = '\$DIN_LOC_ROOT/atm/cam/inic/homme/cami-mam3_0000-01-01_ne30np4_L30_c130424.nc'" >> user_nl_cam
res=$?
if [ $res -ne 0 ]; then
  perr "Unable to add ncdata to user_nl_cam"
fi
if [ $pg -gt 0 ]; then
  echo "se_fv_nphys = ${pg}" >> user_nl_cam
  res=$?
  if [ $res -ne 0 ]; then
    perr "Unable to add se_fv_nphys to user_nl_cam"
  fi
fi

if [ -n "${meshname}" -o -n "${meshfile}" -o $finene -gt 0 ]; then
  echo "se_refined_mesh = .true. ">> user_nl_cam
  res=$?
  if [ $res -ne 0 ]; then
    perr "Unable to add se_ne to user_nl_cam"
  fi
  echo "se_fine_ne=${finene} ">> user_nl_cam
  res=$?
  if [ $res -ne 0 ]; then
    perr "Unable to add se_fine_ne to user_nl_cam"
  fi
  echo "se_mesh_file = '${meshfile}' ">> user_nl_cam
  res=$?
  if [ $res -ne 0 ]; then
    perr "Unable to add se_mesh_file to user_nl_cam"
  fi
fi

if [ -n "`type execca 2> /dev/null`" ]; then
  execca ./case.build
else
  ./case.build
fi
res=$?
if [ $res -ne 0 ]; then
  perr "case.build FAILED"
fi
tname="temp_${datestr}_`date '+%H%M%S'`"
./case.submit 2>&1 | tee ${tname}
res=$?
if [ $res -ne 0 ]; then
  perr "case.submit FAILED"
fi

# Wait until the file is created.
jobno="`grep Submitted ${tname}`"
if [[ "${jobno}" =~ ^[^0-9]*([0-9]*)[^0-9].*$ ]]; then
  jobno="${BASH_REMATCH[1]}"
else
  jobno=""
fi
rm ${tname}
if [ -x "`which bjobs 2> /dev/null`" ]; then
  qs="`which bjobs`"
elif [ -x "`which qstat 2> /dev/null`" ];then
  qs="`which qstat`"
else
  qs=""
fi
if [ -n "${qs}" -a -n "${jobno}" ]; then
  while [ -n "`${qs} | grep ${jobno}`" ]; do
    echo "Waiting for job ${jobno} to finish"
    sleep 20
  done
fi

# Print out instructions for finishing up
echo "You probably need to create the domain files in the <domain> entry"
echo "you added to ${cgridfile}."
echo "To do this, you need to create an ocean to atmosphere conservative"
echo "mapping file using the script:"
echo "${camroot}/cime/tools/mapping/gen_mapping_files/gen_cesm_maps.sh"

# Last but not least, create the domain files
gendomdir="${camroot}/cime/tools/mapping/gen_domain_files"
gendomname="gen_domain"
gendom="${gendomdir}/${gendomname}"
echo "You also need to create the domain files so that the correct ocean"
echo "mask can be accessed."
if [ ! -x "${gendom}" ]; then
  echo "Before creating domain files, you must generate the ${gendomname} program."
  echo "Follow the directions in ${gendomdir}/INSTALL"
fi
echo "To create domain files, follow the directions in ${gendomdir}/README"
