#!/bin/sh
# Test for missing r8, plus a check that every constant that DOES carry
# the _r8 kind suffix is actually a floating point literal (i.e. it
# contains a decimal point or an exponent).
#
# This is TR8.sh with an extra pass added: ADDREALKIND_EXE (below) finds
# real constants that are missing the "_r8" kind suffix, but it does not
# notice the opposite mistake -- an integer literal that was kind-tagged
# with "_r8" without ever being given a decimal point/exponent (e.g.
# "3_r8" instead of "3.0_r8" or "3._r8"). That produces an INTEGER
# literal of kind r8, not the intended REAL literal, so it is checked
# for separately here with check_r8_has_dot().
#

# check_r8_has_dot <directory> [comma,separated,skip,list]
#
# Greps the Fortran source tree rooted at <directory> for numeric
# literals suffixed with "_r8" whose digits are not preceded by a "."
# and are not preceded by an exponent letter (e/E/d/D). Such literals
# are integer constants, not real constants, despite the "_r8" suffix.
# The optional second argument is a comma-separated list of
# sub-directories/files to skip, using the same skip semantics as
# ADDREALKIND_EXE's "-s" option.
check_r8_has_dot() {
  dir="$1"
  skiplist="$2"

  [ -d "$dir" ] || return 0

  filelist=$(cd "$dir" && find . -type f \( -name '*.F90' -o -name '*.f90' -o -name '*.F' -o -name '*.f' \) 2>/dev/null)

  if [ -n "$skiplist" ]; then
    old_ifs="$IFS"
    IFS=','
    for s in $skiplist; do
      filelist=$(printf '%s\n' "$filelist" | grep -vE "(^|/)${s}(/|\$)")
    done
    IFS="$old_ifs"
  fi

  [ -z "$filelist" ] && return 0

  # A literal tagged "_r8" is only a legitimate real constant if its
  # mantissa contains a decimal point OR it has an exponent part (e.g.
  # "3.0_r8", "3._r8", ".5_r8", "1e10_r8", "1.0E+10_r8", "5d0_r8"). If
  # neither is present (e.g. "3_r8", "100_r8", "-5_r8") it is really an
  # INTEGER literal of kind r8. A single-character grep lookbehind can't
  # tell an exponent sign ("E+0_r8", valid) apart from an operator sign
  # ("a + 3_r8", not valid), so this is parsed with perl instead.
  # The lookbehind also must not treat every preceding "." as "mid
  # number" -- old-style relational operators butt a "." directly
  # against the literal (e.g. "res.gt.1.e-14_r8"), so only a dot that is
  # itself preceded by a digit (i.e. an actual decimal point of the same
  # number, as in "...3.5_r8") disqualifies a fresh match start.
  bad=$(cd "$dir" && printf '%s\n' "$filelist" | tr '\n' '\0' | xargs -0 perl -ne '
    my $bad = 0;
    while (/(?<![0-9A-Za-z_])(?<!\d\.)(\d+\.?\d*|\.\d+)([eEdD][+-]?\d+)?_[rR]8\b/g) {
      if (!defined($2) && $1 !~ /\./) { $bad = 1; last }
    }
    print "$ARGV:$.:$_" if $bad;
    close ARGV if eof;
  ')
  scan_rc=$?

  # check exit status for the scan actually completing
  # e.g., missing perl, an xargs without "-0", and any future typo in the perl above.
  if [ $scan_rc -ne 0 ]; then
    echo "TR8_enhanced: ERROR: the _r8 scan of \"${dir}\" failed to run (rc = ${scan_rc})"
    return 1
  fi

  if [ -n "$bad" ]; then
    echo "TR8_enhanced:  found _r8-tagged literal(s) with no decimal point/exponent under \"${dir}\":"
    echo "$bad"
    return 1
  fi
  return 0
}

# Check physics
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/cam
rc=$?
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/cam
rc=`expr $? + $rc`
#ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/atmos_phys/schemes -s mmm,musica
#rc=`expr $? + $rc`
#check_r8_has_dot $CAM_ROOT/components/cam/src/atmos_phys/schemes
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/cam7
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/cam7
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/camrt
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/camrt
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/carma
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/carma
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/rrtmg -s aer_src
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/rrtmg aer_src
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/rrtmgp -s data,ext
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/rrtmgp data,ext
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/simple
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/simple
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/waccm
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/waccm
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/physics/waccmx
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/physics/waccmx
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/cam
rc=$?
check_r8_has_dot $CAM_ROOT/src/physics/cam
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/cam7
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/cam7
#ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/atmos_phys/schemes -s mmm,musica
#rc=`expr $? + $rc`
#check_r8_has_dot $CAM_ROOT/src/atmos_phys/schemes
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/camrt
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/camrt
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/carma
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/carma
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/rrtmg -s aer_src
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/rrtmg aer_src
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/rrtmgp -s data,ext
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/rrtmgp data,ext
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/simple
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/simple
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/waccm
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/waccm
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/physics/waccmx
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/physics/waccmx
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/infrastructure
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/infrastructure
rc=`expr $? + $rc`

fi

#Check Ionosphere
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/ionosphere
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/ionosphere
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/ionosphere
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/ionosphere
rc=`expr $? + $rc`

fi

#Check Chemistry
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/chemistry -s geoschem/geoschem_src,cloud_j,hetp
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/chemistry geoschem/geoschem_src,cloud_j,hetp
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/chemistry -s geoschem/geoschem_src,cloud_j,hetp
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/chemistry geoschem/geoschem_src,cloud_j,hetp
rc=`expr $? + $rc`

fi

#Check Dynamics
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/dynamics/fv3 -s atmos_cubed_sphere,microphys,src_override
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/dynamics/fv3 atmos_cubed_sphere,microphys,src_override
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/dynamics/se
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/dynamics/se
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/dynamics/fv
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/dynamics/fv
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/dynamics/mpas -s dycore
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/dynamics/mpas dycore
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/dynamics/fv3 -s atmos_cubed_sphere,microphys,src_override
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/dynamics/fv3 atmos_cubed_sphere,microphys,src_override
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/dynamics/se
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/dynamics/se
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/dynamics/fv
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/dynamics/fv
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/dynamics/mpas -s dycore
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/dynamics/mpas dycore
rc=`expr $? + $rc`

fi

#Check other
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/control
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/control
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/utils
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/utils
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/control
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/control
rc=`expr $? + $rc`
ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/utils
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/utils
rc=`expr $? + $rc`

fi

#Check coupler
if [ -d "${CAM_ROOT}/components/cam" ]; then

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/components/cam/src/cpl
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/components/cam/src/cpl
rc=`expr $? + $rc`

else

ruby $ADDREALKIND_EXE -r r8 -l 1 -d $CAM_ROOT/src/cpl
rc=`expr $? + $rc`
check_r8_has_dot $CAM_ROOT/src/cpl
rc=`expr $? + $rc`

fi

echo $rc

if [ $rc -ne  0 ]; then
   rc=1
fi



echo $rc
exit $rc
