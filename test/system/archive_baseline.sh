#!/bin/sh -f

echo

if [ $# -ne 1 ]; then
    echo "Invoke archive_baseline.sh -help for usage."
    exit 1
fi

if [ $1 == "-help" ] || [ $1 == "--help" ]; then
cat << EOF1
NAME

	archive_baseline.sh - archive pretag baselines to set locations on
                              izumi and derecho.


SYNOPSIS

	archive_baseline.sh TAGNAME
	  [-help]


ENVIROMENT VARIABLES

	CESM_TESTDIR - Directory that contains the CESM finished results you wish to archive.
	CAM_FC      - Compiler used,  used on derecho (INTEL, NVHPC) and izumi (GNU,NAG), where the compiler
                      name is appended to the archive directory.


BASELINE ARCHIVED LOCATION

	izumi:    /fs/cgd/csm/models/atm/cam/pretag_bl/TAGNAME_gnu
	          /fs/cgd/csm/models/atm/cam/pretag_bl/TAGNAME_nag
        derecho:  /glade/campaign/cesm/community/amwg/cam_baselines/TAGNAME_intel
                  /glade/campaign/cesm/community/amwg/cam_baselines/TAGNAME_nvhpc



HOW TO USE ARCHIVE BASELINES

	on izumi:
          env CESM_TESTDIR=/scratch/cluster/YourName/aux_cam_gnu_yyyymmddsssss CAM_FC=GNU ./archive_baseline.sh cam6_4_XXX
          env CESM_TESTDIR=/scratch/cluster/YourName/aux_cam_nag_yyyymmddsssss CAM_FC=NAG ./archive_baseline.sh cam6_3_XXX

        on derecho:
          env CESM_TESTDIR=/glade/derecho/scratch/YourName/aux_cam_intel_yyyymmddsssss CAM_FC=INTEL ./archive_baseline.sh cam6_4_XXX
          env CESM_TESTDIR=/glade/derecho/scratch/YourName/aux_cam_nvhpc_yyyymmddsssss CAM_FC=NVHPC ./archive_baseline.sh cam6_4_XXX


WARNING

	System changes can cause answer changes. So you may need to create new baselines
        if you are getting unexpected baseline failures.

EOF1
exit
fi

hostname=`hostname`
case $hostname in

  iz*)
    echo "server: izumi"
    if [ -z "$CAM_FC" ]; then
       echo "Must specify CAM_FC"
    fi
    test_file_list="tests_pretag_izumi_${CAM_FC,,}"
    cam_tag=$1_${CAM_FC,,}
    baselinedir="/fs/cgd/csm/models/atm/cam/pretag_bl/$cam_tag"
  ;;

  de*)
    echo "server: derecho"
    if [ -z "$CAM_FC" ]; then
      echo "Must specify CAM_FC"
    fi
    test_file_list="tests_pretag_derecho_${CAM_FC,,}"
    cam_tag=$1
    baselinedir="/glade/campaign/cesm/community/amwg/cam_baselines/$cam_tag"
  ;;

  * ) echo "ERROR: machine $hostname not currently supported"; exit 1 ;;
esac

if [ -d ${baselinedir} ]; then
   echo "ERROR: Baseline $baselinedir already exists."
   exit 1
fi

if [ -n "$CESM_TESTDIR" ]; then

    echo " "
    mkdir -p $baselinedir
    root_baselinedir=`dirname $baselinedir`
    echo "CESM Archiving to $root_baselinedir/$cam_tag"
    if [ -d $CESM_TESTDIR/baselines ]; then
      echo "Using cp to archive baselines."
      cp -r $CESM_TESTDIR/baselines/. $root_baselinedir/$cam_tag
    else
      echo "Using bless_test_results to archive baselines."
      ../../cime/CIME/Tools/bless_test_results -p -t '' -c '' -r $CESM_TESTDIR --baseline-root $root_baselinedir -b $cam_tag -f -s
    fi

    echo " "
fi

case $hostname in

    de* | izumi)
	if [ -z "$CESM_TESTDIR" ]; then
	    echo '***********************************************************************************'
	    echo 'INFO: The aux_cam and test_cam tests were NOT archived'
	    echo "INFO: Must set CESM_TESTDIR (test-root in the create_test) to archive aux_cam tests"
	    echo '***********************************************************************************'
	fi
	;;

esac
