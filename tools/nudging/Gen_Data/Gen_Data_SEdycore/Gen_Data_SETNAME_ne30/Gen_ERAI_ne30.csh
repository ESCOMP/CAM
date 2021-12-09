#! /bin/csh -f

#---------------------------------------
# Script to Generate YOTC Nudging data
#---------------------------------------

#===================================================================
# CONFIGRATION SECTION:
#===================================================================

# Set a REFERENCE (Starting) Date and the numbe rof days to process
#-------------------------------------------------------------------
set RUNNUM   = 01 
set REF_DATE = '19790101'
set NUM_DAYS = 365

# Set INPUT/OUTPUT/TMP directories, 
#-----------------------------------
set NAMELIST = './Config/Config_makeIC-'$RUNNUM'.nl'
set MYLOGDIR  = './LOG/LOG_002.'$RUNNUM'/'
set MYTMPDIR  = './TMP/TMP_002.'$RUNNUM'/'
set MYOUTDIR  = '/PATH-FOR-OUTPUT-DATA/'
set INPUTDIR  = '/PATH-TO-REANALYSES-INPUT/'

# Set ESMF options
#------------------------
set ESMF_interp = 'conserve'
set ESMF_pole   = 'none'
set ESMF_clean  = 'False'
set TMP_clean   = 'True'

# Set Processing options
#------------------------
set CASE                   = 'ERAI_ne30np4_L30'
set DYCORE                 = 'se'
set PRECISION              = 'float'
set VORT_DIV_TO_UV         = 'True'
set SST_MASK               = 'False'
set ICE_MASK               = 'False'
set OUTPUT_PHIS            = 'True'
set REGRID_ALL             = 'True'
set ADJUST_STATE_FROM_TOPO = 'True'
set MASS_FIX               = 'True'

# Set files containig OUTPUT Grid structure and topography
#---------------------------------------------------------
set fname_grid_info        = '/glade/p/cesm/cseg/inputdata/atm/cam/inic/se/cami_0000-01-01_ne30np4_L30_c120315.nc'
set fname_phis_output      = '/glade/p/cesm/cseg/inputdata/atm/cam/topo/USGS-gtopo30_ne30np4_16xdel2.nc'
#set fname_grid_info        = '/glade/p/cesm/cseg/inputdata/atm/cam/inic/se/cami-mam3_0000-01-ne120np4_L30_c110928.nc'
#set fname_phis_output      = '/glade/p/cesm/cseg/inputdata/atm/cam/topo/USGS-gtopo30_ne120np4_16xdel2-PFC-consistentSGH.nc'
set ftype_phis_output      = 'SE_TOPOGRAPHY'

# Set INPUT filename format, type of file, and number of time records per file
#------------------------------------------------------------------------------
set fname   = ( none none none none none none )
set ftype   = ( none none none none none none )
set ftime   = ( none none none none none none )
set hoursec = ( 00000 21600 43200 64800 )
set hourstr = (   00    06    12    18  )

set fname[1] = "YYYYMM/ei.oper.an.ml.regn128sc.YYYYMMDDHH"
set ftype[1] = "Era_Interim_627.0_sc"
set ftime[1] = "1X"

# Set the OUTPUT fileds and the index of the file containing them
#------------------------------------------------------------------
set fields       = ( U V T PS Q )
set field_findex = ( 1 1 1 1  1 )

# Set the index of the file containing the INPUT topography
#------------------------------------------------------------------
set phis_findex  = 1


#===================================================================
# PROCESSING SECTION:
#   Loop over range of times, 
#   generate inpute namelist for each time, 
#   then generate the desired Nudging data
#===================================================================

# Construct the fields string
#-----------------------------
set fieldstr  = $fields[1]
set jj        = $field_findex[1]
@ jj = $jj - 1 
set findexstr = $jj
set ii = 1
set fnum = $#fields
while( $ii <  $fnum )
  @ ii = $ii + 1
  set jj = $field_findex[$ii]
  @ jj = $jj - 1
  set fieldstr  = `echo $fieldstr','$fields[$ii]`
  set findexstr = `echo $findexstr','$jj`
end

set jj = $phis_findex
@ jj = $jj - 1 
set phisstr = $jj

# Build NCL fortran library
#---------------------------
#/usr/local/bin/WRAPIT         MAKEIC.stub MAKEIC.f90
/contrib/ncl-6.1.0/bin/WRAPIT MAKEIC.stub MAKEIC.f90

#
# Loop over range of times to process
#  For the first time true, force 
#  generation of ESMF weight datasets
#=====================================
set ESMF_clean  = 'True'
set NDAYS = 0
while ( $NDAYS < $NUM_DAYS )

  # Set current time values
  #------------------------
  set Yearstr = "`date --date=$REF_DATE+${NDAYS}day +%Y`"
  set Monstr  = "`date --date=$REF_DATE+${NDAYS}day +%m`"
  set Daystr  = "`date --date=$REF_DATE+${NDAYS}day +%d`"
  echo $Yearstr ' / ' $Monstr ' / ' $Daystr

  # Loop over the hourly values (4X daily)
  #========================================
  foreach hnum ( 1 2 3 4 )
    
    # Set Values dependend upon $hnum, clean 
    # up TMP files at the end of each DAY
    #----------------------------------------------------------------------
    set datestr = $Yearstr$Monstr$Daystr$hoursec[$hnum]
    set LOGFILE = $MYLOGDIR'/LogNCL.'$Yearstr$Monstr$Daystr$hoursec[$hnum]
    if( $hnum == 4 ) then
      set TMP_clean   = 'True'
    else
      set TMP_clean   = 'False'
    endif

    # Create Namelist file for current time
    #----------------------------------------
    echo '                                       '             > $NAMELIST
    echo '! Generated Namelist for maekIC_se.ncl '            >> $NAMELIST
    echo '!-------------------                   '            >> $NAMELIST
    echo '&makeic_nl'                                         >> $NAMELIST
    echo ' MYTMPDIR="'$MYTMPDIR'"'                            >> $NAMELIST
    echo ' MYOUTDIR="'$MYOUTDIR'"'                            >> $NAMELIST
    echo ' INPUTDIR="'$INPUTDIR'"'                            >> $NAMELIST
    echo ' INPUTDIR="'$INPUTDIR'"'                            >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' ESMF_interp="'$ESMF_interp'"'                      >> $NAMELIST
    echo ' ESMF_pole  ="'$ESMF_pole'"'                        >> $NAMELIST
    echo ' ESMF_clean ="'$ESMF_clean'"'                       >> $NAMELIST
    echo ' TMP_clean  ="'$TMP_clean'"'                        >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' REF_DATE              ="'$REF_DATE'"'              >> $NAMELIST
    echo ' CASE                  ="'$CASE'"'                  >> $NAMELIST
    echo ' DYCORE                ="'$DYCORE'"'                >> $NAMELIST
    echo ' PRECISION             ="'$PRECISION'"'             >> $NAMELIST
    echo ' VORT_DIV_TO_UV        ="'$VORT_DIV_TO_UV'"'        >> $NAMELIST
    echo ' SST_MASK              ="'$SST_MASK'"'              >> $NAMELIST
    echo ' ICE_MASK              ="'$ICE_MASK'"'              >> $NAMELIST
    echo ' OUTPUT_PHIS           ="'$OUTPUT_PHIS'"'           >> $NAMELIST
    echo ' REGRID_ALL            ="'$REGRID_ALL'"'            >> $NAMELIST
    echo ' ADJUST_STATE_FROM_TOPO="'$ADJUST_STATE_FROM_TOPO'"'>> $NAMELIST
    echo ' MASS_FIX              ="'$MASS_FIX'"'              >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' fname_phis_output     ="'$fname_phis_output'"'     >> $NAMELIST
    echo ' ftype_phis_output     ="'$ftype_phis_output'"'     >> $NAMELIST
    echo ' fname_grid_info       ="'$fname_grid_info'"'       >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    echo ' fields       ="'$fieldstr'"'                       >> $NAMELIST
    echo ' source_files ="'$findexstr'"'                      >> $NAMELIST
    echo ' fname_phis_in='$phisstr                            >> $NAMELIST
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      set file = $fname[$ii]
      set file = `echo $file | sed 's/YYYY/'$Yearstr'/g'`
      set file = `echo $file | sed 's/MM/'$Monstr'/g'`
      set file = `echo $file | sed 's/DD/'$Daystr'/g'`
      if( $ftime[$ii] == "4X" ) then
        set file = `echo $file | sed 's/HH/'$hourstr[1]'/'`
        set file = `echo $file | sed 's/HH/'$hourstr[4]'/'`
      else
        set file = `echo $file | sed 's/HH/'$hourstr[$hnum]'/g'`
      endif
      echo " fname$jj="'"'$file'"'                      >> $NAMELIST
    end
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      echo " ftype$jj="'"'$ftype[$ii]'"'                      >> $NAMELIST
    end
    echo '                                    '               >> $NAMELIST
    set ii = 0
    set fnum = $#fname
    while( $ii <  $fnum )
      set jj = $ii
      @ ii = $ii + 1
      if( $fname[$ii] == "none" ) then
        echo " fdate$jj="'"-1"'                               >> $NAMELIST
      else
        echo " fdate$jj="'"'$datestr'"'                       >> $NAMELIST
      endif
    end
    echo "/"                                                  >> $NAMELIST

    # Execute NCL program to process the data
    #----------------------------------------
    cat $NAMELIST                            >& $LOGFILE.cfg
    ncl makeIC_se_002.ncl NameNumber=$RUNNUM >& $LOGFILE

  # End Loop over the hourly values (4X daily)
  #===========================================
   set ESMF_clean  = 'False'
  end # foreach hnum ( 1 2 3 4 )

# End loop over time
#=======================
@ NDAYS = $NDAYS + 1
end # while ( $NDAYS < $NUM_DAYS )

