#!/bin/sh 
# Test for missing r8
#


# Check physics
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/cam
rc=$?
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/camrt
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/rrtmg -s aer_src
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/simple
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/waccm
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/waccmx
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/cam
rc=$?
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/camrt
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/rrtmg -s aer_src
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/simple
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/waccm
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/waccmx
rc=`expr $? + $rc`

fi

#Check Ionosphere
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/ionosphere
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/ionosphere
rc=`expr $? + $rc`

fi

#Check Chemistry
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/chemistry
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/chemistry
rc=`expr $? + $rc`

fi

#Check Dynamics
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/dynamics/se
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/dynamics/fv
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/dynamics/eul
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/dynamics/se
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/dynamics/fv
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/dynamics/eul
rc=`expr $? + $rc`

fi

#Check other
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/advection
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/control
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/utils
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/advection
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/control
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/utils
rc=`expr $? + $rc`

fi

#Check coupler
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/cpl
rc=`expr $? + $rc`



else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/cpl
rc=`expr $? + $rc`

fi

echo $rc

if [ $rc = 255 ]; then
   rc=1
fi



echo $rc
exit $rc
