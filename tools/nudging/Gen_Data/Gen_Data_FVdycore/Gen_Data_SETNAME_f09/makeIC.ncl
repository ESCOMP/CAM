;
;-------------------------------------------------------------------
; Written in Jun 2012 by Jerry Olson
;-----------------------------------------------------------------------
; Create Initial Conditions file for CAM from other Centers' analyses
;-----------------------------------------------------------------------

external MAKEIC "./MAKEIC.so"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./makeIC_extract_analyses_info.ncl"
load "./makeIC_Create_field_Master_List.ncl"
load "./makeIC_procedures.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

begin

;  NameNumber = 1
  print("--------")
  if(isinteger(NameNumber)) then
   print("Namelist file number to use:")
   print(NameNumber)
  else
   print("NameNumber: The Namelist file number MUST be given on the command line")
   exit
  end if

  ;--------------
  ; Set Constants
  ;--------------
  dummy                  = 0.0
  gravit                 = 9.80616D
  ps0                    = 100000.D
  tmass0                 = 98222.0D/gravit

  ;----------------------------------------------
  ; Initialize field meta-data with user settings
  ;----------------------------------------------
  ;field                  = Initialize_data_processing_settings
  field                  = Initialize_data_processing_namelist(NameNumber)

  ;--------------------------------------------------
  ; Initialize local variables with "field" meta-data
  ;--------------------------------------------------
  ifield                 = 0

  nfiles                 = field@nfiles
  dycore                 = field@dycore
  precision              = field@precision
  plono                  = field@plono
  plato                  = field@plato
  plevo                  = field@plevo
  vort_div_to_UV         = field@vort_div_to_UV
  sst_mask               = field@sst_mask
  ice_mask               = field@ice_mask
  fname                  = field@fname
  ftype                  = field@ftype
  fdate                  = field@fdate
  ref_date               = field@ref_date
  fname_lev_info         = field@fname_lev_info
  mytmpdir               = field@mytmpdir
  mytmpdiro              = field@mytmpdiro
  adjust_state_from_topo = field@adjust_state_from_topo

  mass_fix               = False
  if(dycore .eq. "eul") then
    mass_fix             = True
  end if

;------------------------
; Define disk directories
;------------------------

  print ("/bin/mkdir -p "+mytmpdir)
  system("/bin/mkdir -p "+mytmpdir)
  print ("/bin/mkdir -p "+mytmpdiro)
  system("/bin/mkdir -p "+mytmpdiro)

;-----------------------
; Build output file name
;-----------------------

  fileout = build_filename(field, fdate(0))

;---------------------------------------------------------------
; define nc file for output (including variables and attributes)
;---------------------------------------------------------------

  print ("/bin/rm    -f "+mytmpdiro+"/"+fileout)
  system("/bin/rm    -f "+mytmpdiro+"/"+fileout)

  print (" ")
  print ("Initialize "+mytmpdiro+"/"+fileout)
  print ("   (file will be '"+precision+"' precision)")
  print (" ")

  cdf = addfile(mytmpdiro+"/"+fileout,"c")

;---------------------------------------------
; Read ALL input files from HPSS or local disk
;---------------------------------------------

  print ("")
  print ("Stage all input files to work directory")
  print ("")

  lfname             = fname
  file_prefix_str    = ""

  do ifile           = 0,nfiles-1

;--------------------------------------------------------------
; Check if file is just a duplicate that's already been read in
;--------------------------------------------------------------

    duplicate = False
    if (ifile .gt. 0) then
      do ifile_tmp = 0,ifile-1
        if(fname(ifile) .eq. fname(ifile_tmp) ) then
          duplicate = True
          break
        end if
      end do
    end if

    if(duplicate) then
      lfname(ifile)   = lfname(ifile_tmp)
    else
      file_prefix_str = ifile
      lfname(ifile)   = get_file(file_prefix_str,mytmpdir,fname(ifile))
    end if

  end do

  lfname_lev_info    = get_file("_lev_info"   ,mytmpdir,fname_lev_info)

  field@lfname = lfname

;---------------------
; Open ALL input files
;---------------------

  ifile = 0
  if(ifile .le. nfiles-1) then
    print ("Opening file "+lfname(ifile))
    data0 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")

;----------------------------------
; Make place-holders out of data1-7
;----------------------------------

    data1 = data0
    data2 = data0
    data3 = data0
    data4 = data0
    data5 = data0
    data6 = data0
    data7 = data0

  end if
  ifile = 1
  if(ifile .le. nfiles-1) then
    delete(data1)
    print ("Opening file "+lfname(ifile))
    data1 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")
  end if
  ifile = 2
  if(ifile .le. nfiles-1) then
    delete(data2)
    print ("Opening file "+lfname(ifile))
    data2 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")
  end if
  ifile = 3
  if(ifile .le. nfiles-1) then
    delete(data3)
    print ("Opening file "+lfname(ifile))
    data3 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")
  end if
  ifile = 4
  if(ifile .le. nfiles-1) then
    delete(data4)
    print ("Opening file "+lfname(ifile))
    data4 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")
  end if
  ifile = 5
  if(ifile .le. nfiles-1) then
    delete(data5)
    print ("Opening file "+lfname(ifile))
    data5 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")
  end if
  ifile = 6
  if(ifile .le. nfiles-1) then
    delete(data6)
    print ("Opening file "+lfname(ifile))
    data6 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")
  end if
  ifile = 7
  if(ifile .le. nfiles-1) then
    delete(data7)
    print ("Opening file "+lfname(ifile))
    data7 = addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")
  end if

  print ("Opening file "+lfname_lev_info)
  data_lev_info = addfile(mytmpdir+"/"+lfname_lev_info, "r")

;-----------------------------------------------------------------------------
; Initialize list of input field names that maps to list of output field names
;-----------------------------------------------------------------------------

  field_modified        = initialize_input_field_names(field, data0, data1, data2, data3, data4, data5, data6, data7)
  delete(field)
  field                 = field_modified

  nfields               = field@nfields
  source_file           = field@source_file
  dimensions            = field@dimensions
  horz_interp           = field@horz_interp
  vert_interp           = field@vert_interp
  spec_trunc            = field@spec_trunc
  loutput               = field@loutput
  nfields_out           = num(loutput)
  ps_adjust             = any(field@post_process .eq. "ps_adjust"            )
  q_adjust              = any(field@post_process .eq. "q_adjust"             )
  cwat_adjust           = any(field@post_process .eq. "cwat_adjust"          )
  cloud_fraction_adjust = any(field@post_process .eq. "cloud_fraction_adjust")
  ts_adjust             = any(field@post_process .eq. "ts_adjust"            )

;--------------
; Sanity checks
;--------------

  exit_script = False

  if(dycore .ne. "eul" .and. dycore .ne. "fv" .and. dycore .ne. "homme" ) then
     print("Error:  'dycore' must be:  'eul', 'fv', or 'homme'")
     print("        'dycore' is currently:  "+dycore)
     exit_script = True
  end if
     
  if(precision .ne. "float" .and. precision .ne. "double") then
     print("Error:  'precision' must be:  'float' or 'double' ")
     print("        'precision' is currently:  "+precision)
     exit_script = True
  end if
     
  is_SST_cpl           = False
  is_ice_cov           = False
  is_LANDFRAC_analysis = False
  is_ICEFRAC_analysis  = False

  any_vert_interp  = False
  vert_quad_interp = False

  do i = 0,nfields-1
      
    if(field(i) .eq. "SST_cpl") then
       is_SST_cpl       = True
    end if
    if(field(i) .eq. "ice_cov") then
       is_ice_cov       = True
    end if
    if(field(i) .eq. "LANDFRAC_analysis") then
       is_LANDFRAC_analysis = True
    end if
    if(field(i) .eq. "ICEFRAC_analysis") then
       is_ICEFRAC_analysis = True
    end if

    if(vert_interp(i) .ne. "no_interp" .and. dimensions(i) .eq."3D") then
      any_vert_interp = True
      if(vert_interp(i) .eq."quad_opt1") then
        vert_quad_interp = True
      end if
      if(vert_interp(i) .eq."quad_logp_opt1") then
        vert_quad_interp = True
      end if
    end if

    do ii = i+1,nfields-1
      if(field(i) .eq. field(ii)) then
        i1  = i + 1
        ii1 = ii + 1
        print("Error:  Specifying output field more than once")
        print("        Field "+i1 +" = "+field(i ))
        print("        Field "+ii1+" = "+field(ii))
        exit_script = True
      end if
    end do

    if(dimensions(i) .ne."3D" .and. dimensions(i) .ne."2D") then
      print("Error:  attribute 'dimensions' must be set to '2D' or '3D'")
      print("     For field "+field(i)+", 'dimensions' is set to "+ \
             dimensions(i))
      exit_script = True
    end if
    if(horz_interp(i) .ne."no_interp" .and. \
       horz_interp(i) .ne."cubic"     .and. \
       horz_interp(i) .ne."spectral"  .and. \
       horz_interp(i) .ne."binning"   .and. \
       horz_interp(i) .ne."cubic_sp"     ) then
         print("Error:  attribute 'horz_interp' must be set to "+ \
               "'no_interp', 'cubic', 'binning', 'spectral', or 'cubic_sp'")
         print("     For field "+field(i)+", 'horz_interp' is set to "+ \
                horz_interp(i))
       exit_script = True
    end if
    if(vert_interp(i) .ne."no_interp" .and. \
       vert_interp(i) .ne."linear"    .and. \
       vert_interp(i) .ne."linear_opt1"    .and. \
       vert_interp(i) .ne."linear_opt2"    .and. \
       vert_interp(i) .ne."quad_opt1"    .and. \
       vert_interp(i) .ne."quad_logp_opt1"    .and. \
       vert_interp(i) .ne."linear_logp"    .and. \
       vert_interp(i) .ne."linear_logp_opt1"    .and. \
       vert_interp(i) .ne."linear_logp_opt2"     ) then
         print("Error:  attribute 'vert_interp' must be set to: ")
         print("'no_interp', 'linear', 'quad_opt1', 'quad_logp_opt1', 'linear_logp', ")
         print("'linear_opt1', 'linear_logp_opt1', ")
         print("'linear_opt2', or 'linear_logp_opt2'")
         print("     For field "+field(i)+", 'vert_interp' is set to "+ \
                vert_interp(i))
       exit_script = True
    end if

  end do

;-------------------------
; Exit if any errors found
;-------------------------

  if(exit_script) then
    exit
  end if

;-----------------------
; Initialize output file
;-----------------------

  initialize_IC(cdf, field, nfields, nfields_out, fname)

printVarSummary(cdf)

;-----------------------------------------------
; annotate output file and do consistency checks
;-----------------------------------------------

  if(vort_div_to_UV) then
    cdf@vort_div_to_UV = "True; U/V derived from vort/div from input file"
  end if
  cdf@ps_adjust = "True; Adjust Ps per ECMWF algorithm"
  if(mass_fix) then
    cdf@mass_fix = "True; Fix mass of atmosphere to be consistent with model mass fixer"
  end if
  if(q_adjust) then
    cdf@q_adjust = "True; Eliminate super-saturation"
  end if
  if(cwat_adjust) then
    cdf@cwat_adjust = "True; Remove negatives from cloud water species"
  end if
  if(cloud_fraction_adjust) then
    cdf@cloud_fraction_adjust = "True; Bracket between 0. and 1."
  end if
  if(ts_adjust) then
    cdf@ts_adjust = "True; Adjust TS for topography per ECMWF algorithm"
  end if
  if(sst_mask .and. is_SST_cpl) then
    cdf@sst_mask = "True; Use landfrac and icefrac input masks during SST_cpl interpolation"
    if(.not. is_LANDFRAC_analysis .or. .not. is_ICEFRAC_analysis) then
      print ("Error:     if 'sst_mask' set to True, then LANDFRAC_analysis and ICEFRAC_analysis must also be specified in the field list")
      print ("           (Else, set 'sst_mask' to False)")
      exit_script = True
    end if
  end if
  if(ice_mask .and. is_ice_cov) then
    cdf@ice_mask = "True; Use landfrac during 'ice_cov' interpolation"
    if(.not. is_LANDFRAC_analysis) then
      print ("Error:     if 'ice_mask' set to True, then LANDFRAC_analysis must also be specified in the field list")
      print ("           (Else, set 'ice_mask' to False)")
      exit_script = True
    end if
  end if
  cdf@interp_log_ps = "True; Interpolate ln(Ps) rather than Ps"

;-------------------------
; Exit if any errors found
;-------------------------

  if(exit_script) then
    exit
  end if

;-------------------------------------------
; Fill output dataset with level information
;-------------------------------------------

  if(plevo .gt. 0) then
     initialize_IC_levels(cdf  , field, fname_lev_info, data_lev_info )
  end if

;-----------------------------------------------
; Determine date info for Input and Output files
;-----------------------------------------------

  initialize_IC_dates(cdf  , nfiles, mytmpdir, lfname, ftype, fdate, ref_date, \
                      data0, data1   , data2 , data3, data4, data5, data6, data7 )
  tmp         = doubletointeger(fdate)
  delete(fdate)
  delete(field@fdate)
  fdate       = tmp
  field@fdate = fdate
  delete(tmp)


  count = 0


;--------------------------------------------------------------
; Input fields that will be used to process output fields later
;--------------------------------------------------------------


;-----
; PHIS
;-----

  find_field_index(field, "PHIS_input", nfields, ifield, ifile)
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)
  print (" ")
  print ("Input PHIS_input from file "+mytmpdir+"/"+lfname(ifile)+" for later processing")
  print (" ")

  xfield = read_reformat(data, field, ifield)

;------------
; Interpolate
;------------

  phis_in = horz_process(xfield , cdf, field, ifield, plono, plato, cdf->lon, cdf->lat)
  tmp     = spectral_tr (phis_in, cdf, field, ifield)
  delete(phis_in)
  phis_in = tmp
  delete(tmp)

  output2CAM(phis_in, cdf, field, ifield, count, nfields_out )

  delete(xfield)

;--------------------------------------------------
; Landfrac and icefrac mask (for SST determination)
;--------------------------------------------------

  if(sst_mask .and. is_SST_cpl) then

    find_field_index(field, "LANDFRAC_analysis", nfields, ifield, ifile)
    if( isvar("data") ) then
      delete (data)
    end if
    data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)
    print (" ")
    print ("Input LANDFRAC_analysis from file "+mytmpdir+"/"+lfname(ifile)+" for later processing")
    print (" ")

    xfield = read_reformat(data, field, ifield)

    dims        = dimsizes( xfield(lat|: ,lon|: ) )
    landfrac_in = new( (/dims(0), dims(1) /), double)
    landfrac_in = xfield

    delete(xfield)

    find_field_index(field, "ICEFRAC_analysis", nfields, ifield, ifile)
    if( isvar("data") ) then
      delete (data)
    end if
    data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)
    print (" ")
    print ("Input ICEFRAC_analysis from file "+mytmpdir+"/"+lfname(ifile)+" for later processing")
    print (" ")

    xfield = read_reformat(data, field, ifield)

    dims       = dimsizes( xfield(lat|: ,lon|: ) )
    icefrac_in = new( (/dims(0), dims(1) /), double)
    icefrac_in = xfield

    delete(xfield)
    delete(dims)
  end if

;------------------------------------------
; Landfrac (mask for ice_cov determination)
;------------------------------------------

  if(ice_mask .and. is_ice_cov .and. .not. isvar("landfrac_in") ) then

    find_field_index(field, "LANDFRAC_analysis", nfields, ifield, ifile)
    if( isvar("data") ) then
      delete (data)
    end if
    data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)
    print (" ")
    print ("Input LANDFRAC_analysis from file "+mytmpdir+"/"+lfname(ifile)+" for later processing")
    print (" ")

    xfield = read_reformat(data, field, ifield)

    dims        = dimsizes( xfield(lat|: ,lon|: ) )
    landfrac_in = new( (/dims(0), dims(1) /), double)
    landfrac_in = xfield

    delete(dims)
    delete(xfield)

  end if


;-------------------
; Process PHIS first
;-------------------


  find_field_index(field, "PHIS", nfields, ifield, ifile)
  if( isvar("data") ) then
    delete (data)
  end if
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

  print (" ")
  print ("Begin processing PHIS from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

  xfield = read_reformat(data, field, ifield)

;-------------------------------------------------------------
; Pass through horizontal processor.  No interpolation is done
;-------------------------------------------------------------
printVarSummary(xfield)
printVarSummary(plono)
printVarSummary(plato)
printVarSummary(cdf->lat)
printVarSummary(cdf->lat)
  phis_out = horz_process(xfield  , cdf, field, ifield, plono, plato, cdf->lon, cdf->lat)
  tmp      = spectral_tr (phis_out, cdf, field, ifield)
  delete(phis_out)
  phis_out = tmp
  delete(tmp)

;-------
; Output
;-------

  output2CAM(phis_out, cdf, field, ifield, count, nfields_out )
  delete(xfield)

;-----------
; Process PS
;-----------


  find_field_index(field, "PS", nfields, ifield, ifile)
  ifield_Ps = ifield
  if( isvar("data") ) then
    delete (data)
  end if
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

  print (" ")
  print ("Begin processing PS from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

  xfield = read_reformat(data, field, ifield_Ps)

;------------
; Interpolate
;------------

  ps_xx = xfield
  print ("     Converting field from Ps to log(Ps) prior to interpolation")
  xfield = log(xfield)

  ps_yy = horz_process(xfield, cdf, field, ifield_Ps, plono, plato, cdf->lon, cdf->lat)
  ps_yy_no_spect = ps_yy
  if(.not. ps_adjust .and. horz_interp(ifield_Ps) .ne. "spectral") then
    tmp   = spectral_tr (ps_yy , cdf, field, ifield_Ps)
    delete(ps_yy)
    ps_yy = tmp
    delete(tmp)
  end if

  print ("     Converting field from log(Ps) to Ps post interpolation")
  ps_yy = exp(ps_yy)
  ps_yy_no_spect = exp(ps_yy_no_spect)
  ps_yy_adjust = ps_yy

  delete(xfield)


;------------------------------
; Process T, Q, and Ps together
;------------------------------


  T_horz_interp_only = 0.
  T_input            = 0.
  T_output           = 0.


  delete(T_horz_interp_only)
  delete(T_input)

  find_field_index(field, "T", nfields, ifield, ifile)
  ifield_T = ifield
  if( isvar("data") ) then
    delete (data)
  end if
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

;-------
; Read T
;-------

  print (" ")
  print ("Begin processing T from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

  T_input = read_reformat(data, field, ifield_T)

;---------------------------
; Horizontally Interpolate T
;---------------------------

  yy  = horz_process(T_input , cdf, field, ifield_T, plono, plato, cdf->lon, cdf->lat)
  if( vert_interp(ifield_T) .eq. "no_interp" .and. horz_interp(ifield_T) .ne. "spectral") then
    tmp = spectral_tr (yy, cdf, field, ifield_T)
    T_horz_interp_only = tmp
    delete(tmp)
  else
    T_horz_interp_only = yy
  end if

  delete(yy)

;----------------------
; Open input file for Q
;----------------------

  find_field_index(field, "Q", nfields, ifield, ifile)
  ifield_Q = ifield
  if( isvar("data") ) then
    delete (data)
  end if
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

  print (" ")
  print ("Begin processing Q from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

;-----------------
; Read Q from file
;-----------------

  xfield   = read_reformat(data, field, ifield_Q)

;---------------------------
; Horizontally Interpolate Q
;---------------------------

  Q_horz_interp_only = process_Q_horz(cdf                , data               , xfield   , field   , \
                                      ifile              , ifield_Q           , ftype    , ps_xx   , \
                                      ps_yy_no_spect     , T_input            , T_horz_interp_only )
  if(q_adjust) then
    print ("     Adjusting Q to eliminate RH < 0. or > 1.")
    dims = dimsizes( Q_horz_interp_only )
    press_in_yy  = build_press(data, ftype(ifile), ps_yy_no_spect, "FULL")
    MAKEIC::q2rh(dims(1) ,dims(0) ,dims(2) ,Q_horz_interp_only ,T_horz_interp_only ,press_in_yy )
    MAKEIC::rh2q(dims(1) ,dims(0) ,dims(2) ,Q_horz_interp_only ,T_horz_interp_only ,press_in_yy )
    delete(dims)
    delete(press_in_yy)
  end if
  delete(xfield)

;----------
; Adjust Ps
;----------

  print ("     Adjusting PS")
  press_in_yy       = build_press(data, ftype(ifile), ps_yy, "FULL")
  press_in_yy_half  = build_press(data, ftype(ifile), ps_yy, "HALF")

  dims = dimsizes( T_horz_interp_only )
  MAKEIC::psadj(dims(0) ,dims(0)+1, dims(1) ,dims(2) ,T_horz_interp_only, \
                press_in_yy ,press_in_yy_half ,phis_in ,phis_out ,ps_yy , \
                ps_yy_adjust)

  if(spec_trunc(ifield_Ps) .gt. 0) then
    print ("     Converting field from Ps to log(Ps) prior to truncation")
    ps_yy_adjust = log(ps_yy_adjust)
    tmp          = spectral_tr(ps_yy_adjust, cdf, field, ifield_Ps)
    delete(ps_yy_adjust)
    ps_yy_adjust = tmp
    delete(tmp)
    print ("     Converting field from log(Ps) to Ps post truncation")
    ps_yy_adjust = exp(ps_yy_adjust)
  end if

  delete(dims)
  delete(press_in_yy)
  delete(press_in_yy_half)


;-----------------------------------------------
; Atmospheric mass correction via CAM mass fixer
;  (iterate between mass fixer and Q/T interps)
;-----------------------------------------------

  if(mass_fix) then

    hyai = cdf->hyai
    hybi = cdf->hybi
    gw   = cdf->gw

    do iter = 1,2

      print (" ")
      print ("     ******* Atmospheric Mass Fixer, iteration "+iter+" ******* ")
      print (" ")

;-------------------------------
; Vertically interpolate T and Q
;-------------------------------

      if( vert_interp(ifield_T) .ne. "no_interp") then
        print ("     Vertical interpolation of T")
      end if
      delete(T_output)
      T_output = process_T_vert(cdf                , data         , field , ifile   , \
                                ifield_T           , ftype        , ps_yy_no_spect  , phis_in , \
                                T_horz_interp_only , ps_yy_adjust )

      if( vert_interp(ifield_Q) .ne. "no_interp") then
        print ("     Vertical interpolation of Q")
      end if
      Q_output = process_Q_vert(cdf                , data     , field          , ifile            , \
                                ifield_Q           , ftype    , ps_yy_no_spect ,                    \
                                T_horz_interp_only , T_output , ps_yy_adjust   , Q_horz_interp_only )

      dims = dimsizes( Q_output )
      if(q_adjust) then
        print ("     Adjusting Q to eliminate RH < 0. or > 1.")
        press_out_yy = build_press(cdf , "CAM"       , ps_yy_adjust, "FULL")
        MAKEIC::q2rh(dims(1) ,dims(0) ,dims(2) ,Q_output ,T_output ,press_out_yy )
        MAKEIC::rh2q(dims(1) ,dims(0) ,dims(2) ,Q_output ,T_output ,press_out_yy )
        delete(press_out_yy)
      end if

      MAKEIC::mass_fixer(dims(0) , dims(0)+1, dims(1), dims(2), Q_output ,hyai ,hybi ,gw ,gravit ,ps0 , tmass0, ps_yy_adjust )
    end do

    delete(dims)
    delete(hyai)
    delete(hybi)
    delete(gw)

  end if

;-------------------------------
; Vertically interpolate T and Q
;-------------------------------

  if( vert_interp(ifield_T) .ne. "no_interp") then
    print ("     Vertical interpolation of T")
  end if
  delete(T_output)
  T_output = process_T_vert(cdf                , data         , field , ifile   , \
                            ifield_T           , ftype        , ps_yy_no_spect  , phis_in , \
                            T_horz_interp_only , ps_yy_adjust )
  if( vert_interp(ifield_Q) .ne. "no_interp") then
    print ("     Vertical interpolation of Q")
  end if
  Q_output = process_Q_vert(cdf                , data     , field          , ifile            , \
                            ifield_Q           , ftype    , ps_yy_no_spect ,                    \
                            T_horz_interp_only , T_output , ps_yy_adjust   , Q_horz_interp_only )
  if(q_adjust) then
    print ("     Adjusting Q to eliminate RH < 0. or > 1.")
    dims = dimsizes( Q_output )
    press_out_yy = build_press(cdf , "CAM"       , ps_yy_adjust, "FULL")
    MAKEIC::q2rh(dims(1) ,dims(0) ,dims(2) ,Q_output ,T_output ,press_out_yy )
    MAKEIC::rh2q(dims(1) ,dims(0) ,dims(2) ,Q_output ,T_output ,press_out_yy )
    delete(press_out_yy)
    delete(dims)
  end if

;----------------
; Output T, Q, Ps
;----------------

  output2CAM (ps_yy_adjust , cdf, field, ifield_Ps, count, nfields_out )
  output2CAM (T_output     , cdf, field, ifield_T , count, nfields_out )
  output2CAM (Q_output     , cdf, field, ifield_Q , count, nfields_out )


;-----------------
; Process U and V
;-----------------


;-----------------------------------
; Search for "U" or "US"/"V" or "VS"
;-----------------------------------

  find_field_index(field, "U", nfields, ifield, ifile)
  ifieldu = ifield
  ifileu  = ifile

  find_field_index(field, "V", nfields, ifield, ifile)
  ifieldv = ifield
  ifilev  = ifile

  print ("Begin processing U from file "+mytmpdir+"/"+lfname(ifileu))
  print ("Begin processing V from file "+mytmpdir+"/"+lfname(ifilev))
  print (" ")

  ifile = ifileu

  if(ifileu .ne. ifilev) then
    print ("   Error:  U/V inputs must come from the same file")
    print ("           U info is coming from file "+lfname(ifileu))
    print ("           V info is coming from file "+lfname(ifilev))
    exit
  end if

;----------------
; Open input file
;----------------

  if( isvar("data") ) then
    delete (data)
  end if
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

  xfieldu = read_reformat(data, field, ifieldu)

  xfieldv = read_reformat(data, field, ifieldv)

;-------------------------------------------------------------------------------
; Compute 3D array of INPUT/OUTPUT pressures (for use in vertical interpolation)
;-------------------------------------------------------------------------------

  press_in_yy  = build_press(data, ftype(ifile), ps_yy_no_spect, "FULL")
  press_out_yy = build_press(cdf , "CAM"       , ps_yy_adjust  , "FULL")

;-----------------------------------------------------------
; If inputs are vorticity and divergence, convert to U and V
;-----------------------------------------------------------

  if(vort_div_to_UV) then

;------------------
; Define U/V arrays
;------------------

    u = xfieldu
    v = xfieldv

;---------------
; Convert to U/V
;---------------

    print ("     Converting vorticity and divergence to U/V")
    print (" ")

    vrdv2uvg(xfieldu, xfieldv, u, v)

    delete(xfieldu)
    delete(xfieldv)

    xfieldu = u
    xfieldv = v

    delete(u)
    delete(v)

  end if

;----------------
; Interpolate U/V
;----------------

  u_lat_out = plato
  u_lev_out = dimsizes(xfieldu&lev)
  u_lon_out = plono
  v_lat_out = plato
  v_lev_out = dimsizes(xfieldv&lev)
  v_lon_out = plono
  u_flat    = cdf->lat
  u_flon    = cdf->lon
  v_flat    = cdf->lat
  v_flon    = cdf->lon

  yyu  = new( (/u_lev_out, u_lat_out, u_lon_out /), double)
  yyv  = new( (/v_lev_out, v_lat_out, v_lon_out /), double)
  horz_process_v(xfieldu, xfieldv, cdf, field, ifieldu, ifieldv, u_lon_out, u_lat_out, v_lon_out, v_lat_out, \
                 u_flon, u_flat, v_flon, v_flat, yyu, yyv)

  u_lev_out = plevo
  v_lev_out = plevo

  zzu  = vert_process(yyu    , cdf, field, ifieldu, press_in_yy     , press_out_yy)
  zzv  = vert_process(yyv    , cdf, field, ifieldv, press_in_yy     , press_out_yy)
  tmpu = zzu
  tmpv = zzv
  delete(zzu)
  delete(zzv)
  zzu  = new( (/u_lev_out, u_lat_out, u_lon_out /), double)
  zzv  = new( (/v_lev_out, v_lat_out, v_lon_out /), double)
  spectral_v_tr(tmpu, tmpv, cdf, field, ifieldu, ifieldv, zzu, zzv)
  delete(tmpu)
  delete(tmpv)

  delete(u_lat_out)
  delete(u_lev_out)
  delete(u_lon_out)
  delete(v_lat_out)
  delete(v_lev_out)
  delete(v_lon_out)
  delete(u_flat)
  delete(u_flon)
  delete(v_flat)
  delete(v_flon)

;-------
; Output
;-------
  UU = zzu
  VV = zzv

  output2CAM(zzu, cdf, field, ifieldu, count, nfields_out )
  output2CAM(zzv, cdf, field, ifieldv, count, nfields_out )

  delete(xfieldu)
  delete(xfieldv)
  delete(yyu)
  delete(yyv)
  delete(press_in_yy)
  delete(press_out_yy)

  print ("Done processing U from file ")
  print ("Done processing V from file ")

;-----------------------------------------------------
; Process US/VS: just interpolate the U/V values
;----------------------------------------------------
  find_field_index(field, "US", nfields, ifield, ifile)
  ifieldu = ifield
  ifileu  = ifile
  find_field_index(field, "VS", nfields, ifield, ifile)
  ifieldv = ifield
  ifilev  = ifile
  ifile = ifileu

  if( isvar("data") ) then
    delete (data)
  end if
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

  xfieldu = read_reformat(data, field, ifieldu)
  xfieldv = read_reformat(data, field, ifieldv)

  us_lon  = cdf->lon
  us_slat = cdf->slat
  vs_slon = cdf->slon
  vs_lat  = cdf->lat

;printVarSummary(UU)
;printVarSummary(VV)
  Ve = edge_reformat(VV)
  Ue = edge_reformat(UU)
;  Ve = VV
;  Ue = UU
;printVarSummary(Ue)
;printVarSummary(Ve)

  UUgridName       = mytmpdir+"/UUgrid.nc"
  USgridName       = mytmpdir+"/USgrid.nc"
  UU2USwgtFileName = mytmpdir+"/UU2USwgt.nc"

  VVgridName       = mytmpdir+"/VVgrid.nc"
  VSgridName       = mytmpdir+"/VSgrid.nc"
  VV2VSwgtFileName = mytmpdir+"/VV2VSwgt.nc"

  if(isfilepresent(VVgridName)) then
    print("Found VVgridName description "+VVgridName)
  else
    print("Creating VVgridName description "+VVgridName)
    VVopt                = True
    VVopt@ForceOverwrite = True
    rectilinear_to_SCRIP(VVgridName, Ve&lat,Ve&lon,VVopt)
  end if

  if(isfilepresent(VSgridName)) then
    print("Found VSgridName description "+VSgridName)
  else
    print("Creating VSgridName description "+VSgridName)
    VSopt                = True
    VSopt@ForceOverwrite = True
    rectilinear_to_SCRIP(VSgridName, vs_lat,vs_slon,VSopt)
  end if

  if(isfilepresent(VV2VSwgtFileName)) then
    print("Found VV2VSwgtFileName description "+VV2VSwgtFileName)
  else
    print("Creating VV2VSwgtFileName description "+VV2VSwgtFileName)
    VV2VSopt                      = True
    ;VV2VSopt@InterpMethod         = "conserve"
    VV2VSopt@InterpMethod         = "bilinear"
    VV2VSopt@IgnoreUnmappedPoints = False
    VV2VSopt@ForceOverwrite       = True
    ESMF_regrid_gen_weights(VVgridName,VSgridName,VV2VSwgtFileName,VV2VSopt)
  end if

  if(isfilepresent(UUgridName)) then
    print("Found UUgridName description "+UUgridName)
  else
    print("Creating UUgridName description "+UUgridName)
    UUopt                = True
    UUopt@ForceOverwrite = True
    rectilinear_to_SCRIP(UUgridName, Ue&lat, Ue&lon,UUopt)
  end if

  if(isfilepresent(USgridName)) then
    print("Found USgridName description "+USgridName)
  else
    print("Creating USgridName description "+USgridName)
    USopt                = True
    USopt@ForceOverwrite = True
    rectilinear_to_SCRIP(USgridName,us_slat, us_lon,USopt)
  end if

  if(isfilepresent(UU2USwgtFileName)) then
    print("Found UU2USwgtFileName description "+UU2USwgtFileName)
  else
    print("Creating UU2USwgtFileName description "+UU2USwgtFileName)
    UU2USopt                      = True
    ;UU2USopt@InterpMethod         = "conserve"
    UU2USopt@InterpMethod         = "bilinear"
    UU2USopt@IgnoreUnmappedPoints = False
    UU2USopt@ForceOverwrite       = True
    ESMF_regrid_gen_weights(UUgridName,USgridName,UU2USwgtFileName,UU2USopt)
  end if

  Opt = True
  US = ESMF_regrid_with_weights(Ue,UU2USwgtFileName,Opt)
  Opt = True
  VS = ESMF_regrid_with_weights(Ve,VV2VSwgtFileName,Opt)

printVarSummary(US)
printVarSummary(VS)
 
  output2CAM(US, cdf, field, ifieldu, count, nfields_out )
  output2CAM(VS, cdf, field, ifieldv, count, nfields_out )

  delete(xfieldu)
  delete(xfieldv)

;----------------------------------
; Begin processing all other fields
;----------------------------------



  do ifile = 0,nfiles-1

;----------------
; Open input file
;----------------

    if( isvar("data") ) then
      delete (data)
    end if
    data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

;------------------------------------------
; Begin processing all fields for this file
;------------------------------------------

    print (" ")
    print ("Begin processing fields from file "+mytmpdir+"/"+lfname(ifile))
    print (" ")

    do ifield = 0,nfields-1

      if( loutput(ifield) .and. source_file(ifield) .eq. ifile .and. .not. field@processed(ifield) ) then
        xfield = read_reformat(data, field, ifield)

        if(sst_mask .and. field(ifield) .eq. "SST_cpl") then
          print ("     Masking Ts using icefrac and landfrac for interpolation of "+ \
                 " SSTs (and converting to deg C)")
          dims = dimsizes(xfield(lat|: ,lon|: ) )
          tmp  = new( (/dims(0), dims(1) /), double)
          tmp  = xfield
          delete(xfield)
          xfield = tmp
          MAKEIC::mask_sst(dims(0) , dims(1),           xfield     (lat|: ,lon|: ), \
                           landfrac_in (lat|: ,lon|: ), icefrac_in (lat|: ,lon|: )  )
          delete(dims)
          delete(tmp )
        end if
        if(ice_mask .and. field(ifield) .eq. "ice_cov") then
          print ("     Masking icefrac using landfrac before interpolation of "+ \
                 " ice_cov")
          dims = dimsizes(xfield(lat|: ,lon|: ) )
          tmp  = new( (/dims(0), dims(1) /), double)
          tmp  = xfield
          delete(xfield)
          xfield = tmp
          MAKEIC::mask_ice(dims(0) , dims(1),           xfield     (lat|: ,lon|: ), \
                                                       landfrac_in (lat|: ,lon|: )  )
          delete(dims)
          delete(tmp )
        end if

;-------------------------------------------------------------------------------
; Compute 3D array of INPUT/OUTPUT pressures (for use in vertical interpolation)
;-------------------------------------------------------------------------------

        if(dimensions (ifield) .eq. "3D" .and. \
           vert_interp(ifield) .ne. "no_interp") then
             press_in_yy      = build_press(data, ftype(ifile), ps_yy_no_spect, "FULL")
             press_in_yy_half = build_press(data, ftype(ifile), ps_yy_no_spect, "HALF")
             press_out_yy     = build_press(cdf , "CAM"       , ps_yy_adjust, "FULL")
        else
             press_in_yy      = dummy
             press_in_yy_half = dummy
             press_out_yy     = dummy
        end if

;------------
; Interpolate
;------------

        yy  = horz_process(xfield, cdf, field, ifield, plono, plato, cdf->lon, cdf->lat)
        if( vert_interp(ifield) .eq. "no_interp" .or. dimensions(ifield) .eq. "2D") then
          if( horz_interp(ifield) .ne. "spectral") then
            tmp = spectral_tr (yy, cdf, field, ifield)
            delete(yy)
            yy  = tmp
            delete(tmp)
          end if
        end if

        if(field@post_process(ifield) .eq. "ts_adjust") then
          print ("     Adjusting "+field(ifield)+" for change in PHIS")
          dims = dimsizes( yy(lat|: ,lon|: ) )
          MAKEIC::tsadj(dims(0) ,dims(1), phis_in (lat|: ,lon|: ), \
                                          phis_out(lat|: ,lon|: ), yy(lat|: ,lon|: ) )
          delete(dims)
        end if
        if(field@post_process(ifield) .eq. "cwat_adjust") then
          print ("     Adjusting "+field(ifield)+" to eliminate negative values")
          dims = dimsizes( yy )
          MAKEIC::myminmax(dims(0) ,dims(1) ,dims(2) ,yy, 0., 1.e36)
          delete(dims)
        end if
        if(field@post_process(ifield) .eq. "cloud_fraction_adjust") then
          print ("     Adjusting "+field(ifield)+" between 0. and 1.")
          dims = dimsizes( yy )
          MAKEIC::myminmax(dims(0) ,dims(1) ,dims(2) ,yy, 0., 1.   )
          delete(dims)
        end if
             
        if(vert_interp(ifield) .eq. "quad_opt1"     .and. dimensions  (ifield) .eq. "3D" .or. \
          vert_interp(ifield)  .eq. "quad_logp_opt1".and. dimensions  (ifield) .eq. "3D") then
          zz      = new( (/plevo , plato  , plono /), double)
          zz!0    = "lev"
          zz!1    = "lat"
          zz!2    = "lon"
          dims    = dimsizes( yy )
          if(vert_interp(ifield) .eq. "quad_opt1") then
            print("     QUADRATIC (in pressure) vertical interpolation for this field (OPTION 1)")
            loglin = 1
          end if
          if(vert_interp(ifield) .eq. "quad_logp_opt1") then
            print("     QUADRATIC (in log pressure) vertical interpolation for this field (OPTION 1)")
            loglin = 0
          end if
          MAKEIC::vert_quad_opt1(dims(0) ,dims(0)+1, plevo,  dims(1) ,dims(2) , \
                                 yy(lat|: ,lon|: ,lev|: ), \
                                 press_in_yy      (lat|: ,lon|: ,lev|: ), \
                                 press_in_yy_half (lat|: ,lon|: ,lev|: ), \
                                 press_out_yy     (lat|: ,lon|: ,lev|: ), \
                                 phis_in(lat|: ,lon|: ), ps_yy_no_spect  (lat|: ,lon|: ), \
                                 zz(lat|: ,lon|: ,lev|: ), loglin )
          delete(dims)
          zz&lat = cdf->lat
          zz&lev = cdf->lev
          zz&lon = cdf->lon
        else
          zz = vert_process(yy    , cdf, field, ifield, press_in_yy, press_out_yy)
        end if

        if(field@post_process(ifield) .eq. "cwat_adjust") then
          print ("     Adjusting "+field(ifield)+" to eliminate negative values")
          dims = dimsizes( zz )
          MAKEIC::myminmax(dims(0) ,dims(1) ,dims(2) ,zz, 0., 1.e36)
          delete(dims)
        end if
        if(field@post_process(ifield) .eq. "cloud_fraction_adjust") then
          print ("     Adjusting "+field(ifield)+" between 0. and 1.")
          dims = dimsizes( zz )
          MAKEIC::myminmax(dims(0) ,dims(1) ,dims(2) ,zz, 0., 1.   )
          delete(dims)
        end if

        if( vert_interp(ifield) .ne. "no_interp" .and. dimensions(ifield) .eq. "3D") then
          tmp = spectral_tr (zz, cdf, field, ifield)
          delete(zz)
          zz  = tmp
          delete(tmp)
        end if


;-------
; Output
;-------

        output2CAM(zz, cdf, field, ifield, count, nfields_out )

        delete(yy)
        delete(zz)
        delete(xfield)
        delete(press_in_yy)
        delete(press_in_yy_half)
        delete(press_out_yy)

      end if

    end do

  end do

;---------------
; Fat lady sings
;---------------
  print(" ")
  print(" ")
  print("SUCCESSFULLY COMPLETED PROCESSING")
  print(" ")
  print(" ")
  if(field@tmp_clean .eq. "True") then
    print (" ")
    print ("Cleaning up TEMP files:")
    print ("-----------------------")
    print ("/bin/rm -f "+mytmpdir+"/file*")
    system("/bin/rm -f "+mytmpdir+"/file*")
  end if

exit

end

