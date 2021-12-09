;
;-------------------------------------------------------------------
; Original Version Written in Jun 2012 by Jerry Olson for generating
; data on a uniform grid.
;
; Adapted for nonuniform SE grid. - Patrick Callaghan Sept 2012
;
;-----------------------------------------------------------------------
; Create Initial Conditions file for CAM from other Centers' analyses
;-----------------------------------------------------------------------

external MAKEIC "./MAKEIC.so"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "./makeIC_extract_analyses_info.ncl"
load "./makeIC_Create_field_Master_List.ncl"
load "./makeIC_se_procedures.ncl"
load "./ESMF_regridding.ncl"

begin

;  NameNumber=1
  print("--------")
  if(isinteger(NameNumber)) then
   print("Namelist file number to use:")
   print(NameNumber)
  else
   print("NameNumber: The Namelist file number MUST be given on the command line")
   exit
  end if

  ;----------
  ; Set Constants
  ;----------
  dummy =0.0
  gravit=9.80616D
  ps0   =100000.D
  tmass0=98222.0D/gravit

  ;----------------------------------------------
  ; Initialize field meta-data with user settings
  ;----------------------------------------------
  ; field=Initialize_data_processing_settings
  field=Initialize_diag_settings(NameNumber)

  ;--------------------------------------------------
  ; Initialize local variables with "field" meta-data
  ;--------------------------------------------------
  ifield                =0
  ifield_phis           =0
  nfiles                =field@nfiles
  dycore                =field@dycore
  precision             =field@precision
  vort_div_to_UV        =field@vort_div_to_UV
  sst_mask              =field@sst_mask
  ice_mask              =field@ice_mask
  fname                 =field@fname
  ftype                 =field@ftype
  fdate                 =field@fdate
  ref_date              =field@ref_date
  fname_grid_info       =field@fname_grid_info
  mytmpdir              =field@mytmpdir
  myoutdir              =field@myoutdir 
  adjust_state_from_topo=field@adjust_state_from_topo
  mass_fix              =field@mass_fix

  ;------------------------
  ; Define disk directories
  ;------------------------
  print ("/bin/mkdir -p "+mytmpdir)
  system("/bin/mkdir -p "+mytmpdir)
  print ("/bin/mkdir -p "+myoutdir)
  system("/bin/mkdir -p "+myoutdir)

  ;-----------------------
  ; Build output file name
  ;-----------------------
  fileout=build_filename(field,fdate(0))

  ;---------------------------------------------------------------
  ; define nc file for output (including variables and attributes)
  ;---------------------------------------------------------------
  print ("/bin/rm -f "+myoutdir+"/"+fileout)
  system("/bin/rm -f "+myoutdir+"/"+fileout)

  print (" ")
  print ("Initialize "+myoutdir+"/"+fileout)
  print ("   (file will be '"+precision+"' precision)")
  print (" ")

  cdf = addfile(myoutdir+"/"+fileout,"c")

  ;--------------------------------------------------------
  ; Loop over ALL input files, read from HPSS or local disk
  ;--------------------------------------------------------
  print ("")
  print ("Stage all input files to work directory")
  print ("")

  lfname         =fname
  file_prefix_str=""
  do ifile=0,nfiles-1

    ;--------------------------------------------------------------
    ; Check if file is just a duplicate that's already been read in
    ;--------------------------------------------------------------
    duplicate=False
    if(ifile .gt. 0) then
      do ifile_tmp=0,ifile-1
        if(fname(ifile) .eq. fname(ifile_tmp) ) then
          duplicate = True
          break
        end if
      end do ; ifile_tmp=0,ifile-1
    end if

   if(duplicate) then
     lfname(ifile)  =lfname(ifile_tmp)
   else
     file_prefix_str=ifile
     lfname(ifile)  =get_file(file_prefix_str,mytmpdir,fname(ifile))
   end if

  end do ; ifile=0,nfiles-1

  lfname_grid_info=get_file("_grid_info"   ,mytmpdir,fname_grid_info)
  field@lfname   =lfname

  ;---------------------
  ; Open ALL input files
  ;---------------------
  ifile = 0
  if(ifile .le. nfiles-1) then
    print ("Opening file "+lfname(ifile))
    data0=addfile(mytmpdir+"/"+lfname(ifile)+file_suffix(ftype(ifile)), "r")

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

  print ("Opening file "+lfname_grid_info)
  data_grid_info = addfile(mytmpdir+"/"+lfname_grid_info, "r")

  ;-----------------------------------------------------------------------------
  ; Initialize list of input field names that maps to list of output field names
  ;-----------------------------------------------------------------------------
  field_modified       =initialize_input_field_names(field,data0,data1,data2,data3, \
                                                           data4,data5,data6,data7)
  delete(field)
  field                =field_modified

  nfields              =field@nfields
  source_file          =field@source_file
  dimensions           =field@dimensions
  horz_interp          =field@horz_interp
  vert_interp          =field@vert_interp
  spec_trunc           =field@spec_trunc
  loutput              =field@loutput
  nfields_out          =num(loutput)
  ps_adjust            =any(field@post_process .eq. "ps_adjust"            )
  q_adjust             =any(field@post_process .eq. "q_adjust"             )
  cwat_adjust          =any(field@post_process .eq. "cwat_adjust"          )
  cloud_fraction_adjust=any(field@post_process .eq. "cloud_fraction_adjust")
  ts_adjust            =any(field@post_process .eq. "ts_adjust"            )

  ;-------------------------------------------
  ; Sanity checks: Look for erroreous inputs
  ;-------------------------------------------
   exit_script = False

   if(dycore .ne. "se" ) then
     print("Error:  'dycore' must be:  'se'")
     print("        'dycore' is currently:  "+dycore)
     exit_script=True
   end if

   if(precision .ne. "float" .and. precision .ne. "double") then
     print("Error:  'precision' must be:  'float' or 'double' ")
     print("        'precision' is currently:  "+precision)
     exit_script=True
   end if

   is_SST_cpl          =False
   is_ice_cov          =False
   is_LANDFRAC_analysis=False
   is_ICEFRAC_analysis =False
   any_vert_interp     =False
   vert_quad_interp    =False

   ;----------------------------------------
   ; Loop over fields: check field values
   ;----------------------------------------
   do ii=0,(nfields-1)

     if(field(ii) .eq. "SST_cpl") then
       is_SST_cpl=True
     end if
     if(field(ii) .eq. "ice_cov") then
       is_ice_cov=True
     end if
     if(field(ii) .eq. "LANDFRAC_analysis") then
       is_LANDFRAC_analysis=True
     end if
     if(field(ii) .eq. "ICEFRAC_analysis") then
       is_ICEFRAC_analysis=True
     end if

     if(vert_interp(ii) .ne. "no_interp" .and. dimensions(ii) .eq."3D") then
       any_vert_interp=True
       if(vert_interp(ii) .eq."quad_opt1") then
         vert_quad_interp=True
       end if
       if(vert_interp(ii) .eq."quad_logp_opt1") then
         vert_quad_interp=True
       end if
     end if

     do jj=(ii+1),(nfields-1)
       if(field(ii) .eq. field(jj)) then
         ii1=ii+1
         jj1=jj+1
         print("Error:  Specifying output field more than once")
         print("        Field "+ii1+" = "+field(ii))
         print("        Field "+jj1+" = "+field(jj))
         exit_script=True
       end if
     end do ; jj=ii+1,nfields-1

     if(dimensions(ii) .ne."3D" .and. dimensions(ii) .ne."2D") then
       print("Error:  attribute 'dimensions' must be set to '2D' or '3D'")
       print("     For field "+field(i)+", 'dimensions' is set to "+ dimensions(ii))
       exit_script=True
     end if
     if((horz_interp(ii).ne."no_interp") .and. \
        (horz_interp(ii).ne."bilinear" ) .and. \
        (horz_interp(ii).ne."patch"    ) .and. \
        (horz_interp(ii).ne."conserve" )       ) then
       print("Error:  attribute 'horz_interp' of ESMF must be set to "+ \
             "'no_interp', 'bilinear', 'patch', or 'conserve'")
       print("     For field "+field(ii)+", 'horz_interp' is set to "+horz_interp(ii))
       exit_script=True
     end if
     if((vert_interp(ii).ne."no_interp"       ) .and. \
        (vert_interp(ii).ne."linear"          ) .and. \
        (vert_interp(ii).ne."linear_opt1"     ) .and. \
        (vert_interp(ii).ne."linear_opt2"     ) .and. \
        (vert_interp(ii).ne."quad_opt1"       ) .and. \
        (vert_interp(ii).ne."quad_logp_opt1"  ) .and. \
        (vert_interp(ii).ne."linear_logp"     ) .and. \
        (vert_interp(ii).ne."linear_logp_opt1") .and. \
        (vert_interp(ii).ne."linear_logp_opt2")       ) then
       print("Error:  attribute 'vert_interp' must be set to: ")
       print("'no_interp', 'linear', 'quad_opt1', 'quad_logp_opt1', 'linear_logp', ")
       print("'linear_opt1', 'linear_logp_opt1', ")
       print("'linear_opt2', or 'linear_logp_opt2'")
       print("     For field "+field(ii)+", 'vert_interp' is set to " \
                        +vert_interp(ii))
       exit_script=True
     end if

   end do ; ii=0,(nfields-1)

  ;-------------------------
  ; Exit if any errors found
  ;-------------------------
  if(exit_script) then
    exit
  end if

  ;-----------------------
  ; Initialize output file
  ;-----------------------
  initialize_IC(cdf,field,nfields,nfields_out,fname,fname_grid_info,data_grid_info)
  ncol  = field@ncol
  plono = field@plono
  plato = field@plato
  plevo = field@plevo

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
      print ("Error:     if 'sst_mask' set to True, then LANDFRAC_analysis ")
      print (" and ICEFRAC_analysis must also be specified in the field list")
      print ("           (Else, set 'sst_mask' to False)")
      exit_script = True
    end if
  end if
  if(ice_mask .and. is_ice_cov) then
    cdf@ice_mask = "True; Use landfrac during 'ice_cov' interpolation"
    if(.not. is_LANDFRAC_analysis) then
      print ("Error:     if 'ice_mask' set to True, then LANDFRAC_analysis ")
      print (" must also be specified in the field list")
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

  ;-----------------------------------------------
  ; Determine date info for Input and Output files
  ;-----------------------------------------------
  initialize_IC_dates(cdf,nfiles,mytmpdir,lfname,ftype,fdate,ref_date, \
                      data0,data1,data2,data3,data4,data5,data6,data7)
  tmp=doubletointeger(fdate)
  delete(fdate)
  delete(field@fdate)
  fdate=tmp
  field@fdate=fdate
  delete(tmp)

  ;------------------------------------------------------
  ; Input PHIS: Select the file containing the data,
  ; Read in and Interpolate
  ;------------------------------------------------------
  find_field_index(field, "PHIS_input", nfields, ifield, ifile)
  data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)
  print (" ")
  print ("Input PHIS_input from file "+mytmpdir+"/"+lfname(ifile)+" for later processing")
  print (" ")
  xfield = read_reformat(data, field, ifield)

  ;------------------------------------------------------
  ; Generate ESMF mapping files. ASSUME that all input data will 
  ; be on the same horizonatal grid.
  ;------------------------------------------------------
  print (" ")
  print (" Preparing ESMF regridding datasets:")
  print ("------------------------------------")
  plati=xfield&lat
  ploni=xfield&lon

  srcLab="_"+sprinti("%0.3i",dimsizes(plati(lat|:)))+"x" \
            +sprinti("%0.3i",dimsizes(ploni(lon|:)))
  dstLab="_"+sprinti("%0.8i",ncol)
  srcGridName=mytmpdir+"/SrcGrid"+srcLab+".nc"
;***********************************************************
; WORK AROUND: the 'unstructured_to_ESMF' routine is
;              not good for non-graphic interpolation.
;              It's okay for plotting values, but not
;              for the mapping we want to do. For now a 
;              SCRIP file for ne30 will be sliced in so 
;              we can get some work done.
;----------------------------------------------------------
;  dstGridName=mytmpdir+"/DstGrid"+dstLab+".nc"
;  dstGridName="/glade/p/cesmdata/cseg/mapping/grids/ne120np4_pentagons_100310.nc"
  dstGridName="/glade/p/cesmdata/cseg/mapping/grids/ne30np4_091226_pentagons.nc"

;***********************************************************
  wgtGridName=mytmpdir+"/WgtGrid"+srcLab+dstLab+"_" \
                      +field@esmf_interp+"_"+field@esmf_pole+".nc"

  ; If flagged, remove any pre-existing ESMF weight files
  ;-------------------------------------------------------
  if(field@esmf_clean .eq. "True") then
   print (" Cleaning ESMF files:")
   print ("/bin/rm -f "+srcGridName)
   system("/bin/rm -f "+srcGridName)
;***********************************************************
; WORK AROUND: the 'unstructured_to_ESMF' routine
;  print ("/bin/rm -f "+dstGridName)
;  system("/bin/rm -f "+dstGridName)
;***********************************************************
   print ("/bin/rm -f "+wgtGridName)
   system("/bin/rm -f "+wgtGridName)
  end if

  ESMFopt=True
  ESMFopt@ForceOverwrite=True
  ESMFopt@SrcGridLat    =plati
  ESMFopt@SrcGridLon    =ploni
  ESMFopt@DstGridLat    =plato
  ESMFopt@DstGridLon    =plono
;  ESMFopt@InterpMethod  ="patch" ; "conserve" ; "bilinear"
;  ESMFopt@Pole          ="teeth"    ; "all" ; "teeth"
  ESMFopt@InterpMethod  =field@esmf_interp
  ESMFopt@Pole          =field@esmf_pole
  ESMFopt@SrcESMF       =False
;***********************************************************
; WORK AROUND: the 'unstructured_to_ESMF' routine
;  ESMFopt@DstESMF       =True
  ESMFopt@DstESMF       =False
;***********************************************************
  ESMFopt@Debug         =False ; True

  ; Use existing source grid map, otherwise create one
  ;---------------------------------------------------
  if(isfilepresent(srcGridName)) then
    print(" Using source grid:"+srcGridName)
  else
    print(" Creating source grid:"+srcGridName)
    rectilinear_to_SCRIP(srcGridName,plati,ploni,ESMFopt)
  end if

  ; Use existing destination grid map, otherwise create one
  ;---------------------------------------------------------
;***********************************************************
; WORK AROUND: the 'unstructured_to_ESMF' routine
;  if(isfilepresent(dstGridName)) then
    print(" Using destination grid:"+dstGridName)
;  else
;    print(" Creating destination grid:"+dstGridName)
;    unstructured_to_ESMF(dstGridName,plato,plono,ESMFopt)
;  end if
;***********************************************************

  ; Use existing ESMF weights, otherwise create one
  ;---------------------------------------------------------
  if(isfilepresent(wgtGridName)) then
    print(" Using ESMF weights:"+wgtGridName)
  else
    print(" Creating ESMF weights:"+wgtGridName)
    ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtGridName,ESMFopt)
  end if
  print (" ")

  ;------------------------------------
  ; Interpolate PHIS to output grid
  ;------------------------------------
;DIAG
 print("PHIS grid")
 print(dimsizes(xfield))
;DIAG
  phis_in = ESMF_regrid_with_weights(xfield,wgtGridName,ESMFopt)
  phis_in!0="ncol"
  delete(xfield)

  ;------------------------------------------------
  ; Initialize Field counter for output2CAM routine 
  ; and write out PHIS values.
  ;------------------------------------------------
  count = 0
  output2CAM(phis_in , cdf, field, ifield     , count, nfields_out )

  ;---------------------------------------------------------------
  ; Read in "PHIS" output topography, save values to output file 
  ;---------------------------------------------------------------
  find_field_index(field, "PHIS", field@nfields, ifield_phis, ifile)
  if( isvar("data") ) then
    delete (data)
  end if
  data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)
  print (" ")
  print ("Begin processing PHIS from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

  phis_out = data->$field@source_field(ifield_phis)$
  field@processed(ifield_phis) = True
  output2CAM(phis_out, cdf, field, ifield_phis, count, nfields_out )

  ;--------------------------------------------------
  ; Landfrac and icefrac mask (for SST determination)
  ;--------------------------------------------------
  if(sst_mask .and. is_SST_cpl) then

    ; Select the file containing the data
    ;-----------------------------------------
    find_field_index(field, "LANDFRAC_analysis", nfields, ifield, ifile)
    if( isvar("data") ) then
      delete (data)
    end if
    data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)
    print (" ")
    print ("Input LANDFRAC_analysis from file "+mytmpdir+"/"+lfname(ifile))
    print (" for later processing")
    print (" ")

    xfield = read_reformat(data, field, ifield)

    dims       =dimsizes(xfield(lat|: ,lon|: ))
    landfrac_in=new( (/dims(0), dims(1) /), double)
    landfrac_in=xfield
    delete(xfield)

    find_field_index(field, "ICEFRAC_analysis", nfields, ifield, ifile)
    if( isvar("data") ) then
      delete (data)
    end if
    data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)
    print (" ")
    print ("Input ICEFRAC_analysis from file "+mytmpdir+"/"+lfname(ifile))
    print (" for later processing")
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
    data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)
    print (" ")
    print ("Input LANDFRAC_analysis from file "+mytmpdir+"/"+lfname(ifile))
    print (" for later processing")
    print (" ")

    xfield=read_reformat(data,field,ifield)

    dims       =dimsizes( xfield(lat|: ,lon|: ) )
    landfrac_in=new( (/dims(0), dims(1) /), double)
    landfrac_in=xfield

    delete(dims)
    delete(xfield)
  end if

  ;--------------------------
  ; Process PS: Interpolate
  ;--------------------------
  find_field_index(field, "PS", nfields, ifield, ifile)
  ifield_Ps = ifield
  if( isvar("data") ) then
    delete (data)
  end if
  data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)
  print (" ")
  print ("Begin processing PS from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

  xfield = read_reformat(data, field, ifield_Ps)

  ps_xx = xfield
  print ("     Converting field from Ps to log(Ps) prior to interpolation")
  xfield = log(xfield)

  ; Interpolate PS to output grid    # DIAG-FIX #??? are 'no_spect/adjust' versions needed??? 
  ;------------------------------------
;DIAG
 print("PS grid")
 print(dimsizes(xfield))
;DIAG
  ps_yy = ESMF_regrid_with_weights(xfield,wgtGridName,ESMFopt)
  ps_yy!0="ncol"
  ps_yy_no_spect = ps_yy

  print ("     Converting field from log(Ps) to Ps post interpolation")
  ps_yy         =exp(ps_yy)
  ps_yy_no_spect=exp(ps_yy_no_spect)
  ps_yy_adjust  =ps_yy
  delete(xfield)

;------------------------------
; Process T, Q, and Ps together
;------------------------------
  T_horz_interp_only = 0.
  T_input            = 0.
  T_output           = 0.

  ;------------------------------------------
  ; Read T and Horizontally Interpolate T
  ;------------------------------------------
  delete(T_horz_interp_only)
  delete(T_input)

  find_field_index(field, "T", nfields, ifield, ifile)
  ifield_T = ifield
  if( isvar("data") ) then
    delete (data)
  end if
  data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)

  print (" ")
  print ("Begin processing T from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

  ; Interpolate TEMP to output grid    # DIAG-FIX #??? are 'interp_only' version needed??? 
  ;------------------------------------
  T_input = read_reformat(data, field, ifield_T)
;DIAG
 print("TEMP grid")
 print(dimsizes(T_input))
;DIAG
  yy = ESMF_regrid_with_weights(T_input,wgtGridName,ESMFopt)
  yy!0=T_input!0
  yy!1="ncol"
  T_horz_interp_only = yy
  delete(yy)

;*******************************************************************************
;DIAG 
;; if(False) then
;;   printVarSummary(ploni)
;;   printVarSummary(plono)
;;   printVarSummary(ps_xx)
;;   printVarSummary(ps_yy)
;;   printVarSummary(T_input)
;;   printVarSummary(T_horz_interp_only)
;;   wks1=gsn_open_wks("X11","Input Data")
;;   wks2=gsn_open_wks("X11","Interpolated Data")
;;   gsn_define_colormap(wks1,"BlAqGrYeOrReVi200")   ; choose colormap
;;   gsn_define_colormap(wks2,"BlAqGrYeOrReVi200")   ; choose colormap
;; 
;;   lat1d = plato          ; Pull off lat/lon variable and
;;   lon1d = plono          ; convert both to 1D.
;; 
;;   res2                  = True         ; plot modifications desired
;;   res2@gsnMaximize      = True         ; Maximize size of plot in frame
;;   res2@gsnSpreadColors  = True         ; Use full colormap, but start
;;   res2@sfXArray         = lon1d        ; Required to tell NCL where to
;;   res2@sfYArray         = lat1d        ; overlay data on globe.
;;   res2@mpCenterLonF     = 180.
;;   res2@cnFillOn         = True         ; Turn on contour fill
;;   res2@cnFillMode       = "RasterFill"
;;   res2@cnLinesOn        = False        ; Turn off contour lines
;;   res2@cnLineLabelsOn   = False        ; Turn off contour line labels
;;   res2@lbLabelAutoStride= True         ; Clean up labelbar labels.
;;   res2@tiMainString     = "Interpolated Data"
;; 
;;   res1                  = True         ; plot modifications desired
;;   res1@gsnMaximize      = True         ; Maximize size of plot in frame
;;   res1@gsnSpreadColors  = True         ; Use full colormap, but start
;; ;  res1@sfXArray         = ploni        ; Required to tell NCL where to
;; ;  res1@sfYArray         = plati        ; overlay data on globe.
;;   res2@mpCenterLonF     = 180.
;;   res1@cnFillOn         = True         ; Turn on contour fill
;;   res1@cnFillMode       = "RasterFill"
;;   res1@cnLinesOn        = False        ; Turn off contour lines
;;   res1@cnLineLabelsOn   = False        ; Turn off contour line labels
;;   res1@lbLabelAutoStride= True         ; Clean up labelbar labels.
;;   res1@tiMainString     = "Input Data"
;; 
;;   Xvar1 = T_input(0,:,:)
;;   Xvar2 = T_horz_interp_only(0,:)
;;   map1 = gsn_csm_contour_map_ce(wks1,Xvar1,res1)
;;   map1 = gsn_csm_contour_map_ce(wks1,Xvar1,res1)
;;   map2 = gsn_csm_contour_map_ce(wks2,Xvar2,res2)
;;   Xvar1 = T_input(90,:,:)
;;   Xvar2 = T_horz_interp_only(90,:)
;;   map1 = gsn_csm_contour_map_ce(wks1,Xvar1,res1)
;;   map2 = gsn_csm_contour_map_ce(wks2,Xvar2,res2)
;;   Xvar1 = T_input(45,:,:)
;;   Xvar2 = T_horz_interp_only(45,:)
;;   map1 = gsn_csm_contour_map_ce(wks1,Xvar1,res1)
;;   map2 = gsn_csm_contour_map_ce(wks2,Xvar2,res2)
;;   map1 = gsn_csm_contour_map_ce(wks1,ps_xx,res1)
;;   map2 = gsn_csm_contour_map_ce(wks2,ps_yy,res2)
;; end if
;DIAG 
;*******************************************************************************

  ;----------------------------------------
  ; Read Q and Horizontally Interpolate Q
  ;----------------------------------------
  find_field_index(field, "Q", nfields, ifield, ifile)
  ifield_Q = ifield
  if( isvar("data") ) then
    delete (data)
  end if
  data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)

  print (" ")
  print ("Begin processing Q from file "+mytmpdir+"/"+lfname(ifile))
  print (" ")

  xfield   = read_reformat(data, field, ifield_Q)

;DIAG
 print("Q grid")
 print(dimsizes(xfield))
;DIAG
  Q_horz_interp_only = ESMF_regrid_with_weights(xfield,wgtGridName,ESMFopt)
  Q_horz_interp_only!0=xfield!0
  Q_horz_interp_only!1="ncol"

  if(q_adjust) then
    ; Force surface pressure into 2D array so build_press() works
    ;--------------------------------------------------------------
    ps_dim = dimsizes(ps_yy_no_spect)
    ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_no_spect,0)
    ps_conform!0="ncol"
    ps_conform!1="null"
    press_conform  = build_press(data, ftype(ifile), ps_conform, "FULL")
    press_in_yy  = press_conform(:,:,0)

    ; Make adjustments to Q
    ;----------------------------
    print ("     Adjusting Q to eliminate RH < 0. or > 1.")
    dims = dimsizes( Q_horz_interp_only )
    MAKEIC::q2rh_se(dims(0) ,dims(1) ,Q_horz_interp_only , \
                    T_horz_interp_only ,press_in_yy )
    MAKEIC::rh2q_se(dims(0) ,dims(1) ,Q_horz_interp_only , \
                    T_horz_interp_only ,press_in_yy )
    delete(dims)
    delete(ps_dim)
    delete(ps_conform)
    delete(press_conform)
    delete(press_in_yy)
  end if
  delete(xfield)

  ;----------
  ; Adjust Ps
  ;----------
  print ("     Adjusting PS")
  ; Force surface pressure into 2D array so build_press() works
  ;--------------------------------------------------------------
  ps_dim = dimsizes(ps_yy)
  ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy,0)
  ps_conform!0="ncol"
  ps_conform!1="null"
  press_conform  = build_press(data,ftype(ifile),ps_conform, "FULL")
  press_in_yy  = press_conform(:,:,0)
  delete(press_conform)
  press_conform  = build_press(data,ftype(ifile),ps_conform, "HALF")
  press_in_yy_half=press_conform(:,:,0)
  delete(press_conform)
  delete(ps_conform)
  delete(ps_dim)

  dims = dimsizes( T_horz_interp_only )
  MAKEIC::psadj_se(dims(0) ,dims(0)+1, dims(1) ,T_horz_interp_only, \
                press_in_yy ,press_in_yy_half ,phis_in ,phis_out ,ps_yy , \
                ps_yy_adjust)
  delete(dims)
  delete(press_in_yy)
  delete(press_in_yy_half)

  ;-----------------------------------------------
  ; Atmospheric mass correction via CAM mass fixer
  ;  (iterate between mass fixer and Q/T interps)
  ;-----------------------------------------------
  if(mass_fix) then
    hyai=cdf->hyai
    hybi=cdf->hybi
    area=cdf->area

    do iter=1,2
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
      T_output=process_T_vert(cdf               ,data        ,field , ifile ,         \
                              ifield_T          ,ftype       ,ps_yy_no_spect,phis_in, \
                              T_horz_interp_only,ps_yy_adjust)

      if( vert_interp(ifield_Q) .ne. "no_interp") then
        print ("     Vertical interpolation of Q")
      end if
      Q_output=process_Q_vert(cdf               ,data    ,field         ,ifile,            \
                              ifield_Q          ,ftype   ,ps_yy_no_spect,                  \
                              T_horz_interp_only,T_output,ps_yy_adjust  ,Q_horz_interp_only)

      dims=dimsizes( Q_output )
      if(q_adjust) then
        print ("     Adjusting Q to eliminate RH < 0. or > 1.")
        ps_dim = dimsizes(ps_yy_adjust)
        ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_adjust,0)
        ps_conform!0="ncol"
        ps_conform!1="null"
        press_conform  = build_press(cdf , "CAM" ,ps_conform, "FULL")
        press_out_yy  = press_conform(:,:,0)
        delete(press_conform)
        delete(ps_conform)
        delete(ps_dim)
        MAKEIC::q2rh_se(dims(0) ,dims(1) ,Q_output , T_output ,press_out_yy )
        MAKEIC::rh2q_se(dims(0) ,dims(1) ,Q_output , T_output ,press_out_yy )
        delete(press_out_yy)
      end if

      MAKEIC::mass_fixer_se(dims(0),(dims(0)+1),dims(1),Q_output, \
                            hyai,hybi,area,gravit,ps0,tmass0,ps_yy_adjust)
    end do ; iter=1,2

    delete(dims)
    delete(hyai)
    delete(hybi)
    delete(area)
  end if ; (mass_fix) then

  ;-------------------------------
  ; Vertically interpolate T and Q
  ;-------------------------------
  if( vert_interp(ifield_T) .ne. "no_interp") then
    print ("     Vertical interpolation of T")
  end if
  delete(T_output)
  T_output=process_T_vert(cdf               ,data        ,field , ifile ,          \
                          ifield_T          ,ftype       ,ps_yy_no_spect,phis_in , \
                          T_horz_interp_only,ps_yy_adjust)
  if( vert_interp(ifield_Q) .ne. "no_interp") then
    print ("     Vertical interpolation of Q")
  end if
  Q_output=process_Q_vert(cdf               ,data    ,field         ,ifile,            \
                          ifield_Q          ,ftype   ,ps_yy_no_spect,                  \
                          T_horz_interp_only,T_output,ps_yy_adjust  ,Q_horz_interp_only)
  if(q_adjust) then
    print ("     Adjusting Q to eliminate RH < 0. or > 1.")
    dims = dimsizes( Q_output )
    ps_dim = dimsizes(ps_yy_adjust)
    ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_adjust,0)
    ps_conform!0="ncol"
    ps_conform!1="null"
    press_conform  = build_press(cdf , "CAM" ,ps_conform, "FULL")
    press_out_yy  = press_conform(:,:,0)
    delete(press_conform)
    delete(ps_conform)
    delete(ps_dim)
    MAKEIC::q2rh_se(dims(0) ,dims(1) ,Q_output , T_output ,press_out_yy )
    MAKEIC::rh2q_se(dims(0) ,dims(1) ,Q_output , T_output ,press_out_yy )
    delete(press_out_yy)
    delete(dims)
  end if

  ;----------------
  ; Output T, Q, Ps
  ;----------------
  output2CAM (ps_yy_adjust, cdf, field, ifield_Ps, count, nfields_out )
  output2CAM (T_output    , cdf, field, ifield_T , count, nfields_out )
  output2CAM (Q_output    , cdf, field, ifield_Q , count, nfields_out )

;------------------------------
; Now, Process U and V
;------------------------------

  ;-----------------------------------
  ; Search for "U" / "V" 
  ;-----------------------------------
  find_field_index(field,"U" ,nfields,ifield,ifile)
  ifieldu=ifield
  ifileu =ifile

  find_field_index(field, "V", nfields, ifield, ifile)
  ifieldv=ifield
  ifilev =ifile

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

  ;-------------------------------
  ; Open input file, read in U/V
  ;-------------------------------
  if( isvar("data") ) then
    delete (data)
  end if
  data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)

  xfieldu=read_reformat(data,field,ifieldu)
  xfieldv=read_reformat(data,field,ifieldv)

  ;-------------------------------------------------------------------------------
  ; Compute 3D array of INPUT/OUTPUT pressures (for use in vertical interpolation)
  ;  -- Force surface pressure into 2D array so build_press() works
  ;--------------------------------------------------------------
  ps_dim = dimsizes(ps_yy_no_spect)
  ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_no_spect,0)
  ps_conform!0="ncol"
  ps_conform!1="null"
  press_conform  = build_press(data, ftype(ifile), ps_conform, "FULL")
  press_in_yy  = press_conform(:,:,0)
  delete(press_conform)
  delete(ps_conform)
  delete(ps_dim)

  ps_dim = dimsizes(ps_yy_adjust)
  ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_adjust,0)
  ps_conform!0="ncol"
  ps_conform!1="null"
  press_conform  = build_press(cdf , "CAM" ,ps_conform, "FULL")
  press_out_yy  = press_conform(:,:,0)
  delete(press_conform)
  delete(ps_conform)
  delete(ps_dim)

  ;-----------------------------------------------------------
  ; If inputs are vorticity and divergence, convert to U and V
  ;-----------------------------------------------------------
  if(vort_div_to_UV) then
    ; Define U/V arrays, Convert to U/V
    ;-----------------------------------
    u=xfieldu
    v=xfieldv
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
;DIAG
 print("U grid")
 print(dimsizes(xfieldu))
;DIAG
  yyu = ESMF_regrid_with_weights(xfieldu,wgtGridName,ESMFopt)
  yyu!0=xfieldu!0
  yyu!1="ncol"

;DIAG
 print("V grid")
 print(dimsizes(xfieldv))
;DIAG
  yyv = ESMF_regrid_with_weights(xfieldv,wgtGridName,ESMFopt)
  yyv!0=xfieldv!0
  yyv!1="ncol"

  zzu = vert_process(yyu,cdf,field,ifieldu,press_in_yy,press_out_yy)
  zzv = vert_process(yyv,cdf,field,ifieldv,press_in_yy,press_out_yy)

  ;-------
  ; Output
  ;-------
  output2CAM(zzu, cdf, field, ifieldu, count, nfields_out )
  output2CAM(zzv, cdf, field, ifieldv, count, nfields_out )
  delete(xfieldu)
  delete(xfieldv)
  delete(yyu)
  delete(yyv)
  delete(press_in_yy)
  delete(press_out_yy)


;----------------------------------
; Now, Process all other fields
;----------------------------------
  do ifile=0,(nfiles-1)

    ;----------------
    ; Open input file
    ;----------------
    if( isvar("data") ) then
      delete (data)
    end if
    data=choose_file(data0,data1,data2,data3,data4,data5,data6,data7,ifile)

    ;------------------------------------------
    ; Begin processing all fields for this file
    ;------------------------------------------
    print (" ")
    print ("Begin processing fields from file "+mytmpdir+"/"+lfname(ifile))
    print (" ")
    do ifield=0,(nfields-1)

      if((    loutput(ifield)          ).and. \
         (source_file(ifield).eq.ifile ).and. \
         (.not. field@processed(ifield))      ) then
        xfield = read_reformat(data, field, ifield)

        if(sst_mask .and. field(ifield) .eq. "SST_cpl") then
          print ("     Masking Ts using icefrac and landfrac for interpolation of "+ \
                 " SSTs (and converting to deg C)")
          dims=dimsizes(xfield(lat|: ,lon|: ) )
          tmp =new( (/dims(0), dims(1) /), double)
          tmp =xfield
          delete(xfield)
          xfield=tmp
          MAKEIC::mask_sst(dims(0),dims(1),               xfield(lat|: ,lon|: ), \
                           landfrac_in(lat|: ,lon|: ),icefrac_in(lat|: ,lon|: )  )
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
          MAKEIC::mask_ice(dims(0),dims(1),xfield(lat|: ,lon|: ),landfrac_in(lat|: ,lon|: ))
          delete(dims)
          delete(tmp )
        end if

        ;-------------------------------------------------------------------------------
        ; Compute 3D array of INPUT/OUTPUT pressures (for use in vertical interpolation)
        ;-------------------------------------------------------------------------------
        if((dimensions (ifield) .eq. "3D"       ) .and. \
           (vert_interp(ifield) .ne. "no_interp")       ) then
          ; Force surface pressure into 2D array so build_press() works
          ;--------------------------------------------------------------
          ps_dim = dimsizes(ps_yy_no_spect)
          ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_no_spect,0)
          ps_conform!0="ncol"
          ps_conform!1="null"
          press_conform  = build_press(data,ftype(ifile),ps_conform, "FULL")
          press_in_yy  = press_conform(:,:,0)
          delete(press_conform)
          press_conform  = build_press(data,ftype(ifile),ps_conform, "HALF")
          press_in_yy_half=press_conform(:,:,0)
          delete(press_conform)
          delete(ps_conform)
          delete(ps_dim)
          ps_dim = dimsizes(ps_yy_adjust)
          ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_adjust,0)
          ps_conform!0="ncol"
          ps_conform!1="null"
          press_conform  = build_press(cdf , "CAM" ,ps_conform, "FULL")
          press_out_yy  = press_conform(:,:,0)
          delete(press_conform)
          delete(ps_conform)
          delete(ps_dim)
        else
          press_in_yy     =dummy
          press_in_yy_half=dummy
          press_out_yy    =dummy
        end if

        ;------------
        ; Interpolate
        ;------------
;DIAG
 print("XXX grid")
 print(dimsizes(xfield))
;DIAG
        yy = ESMF_regrid_with_weights(xfield,wgtGridName,ESMFopt)
        yy!0=xfield!0
        yy!1="ncol"

        if(field@post_process(ifield) .eq. "ts_adjust") then
          print ("     Adjusting "+field(ifield)+" for change in PHIS")
          dims = dimsizes( yy(ncol|: ))
          MAKEIC::tsadj_se(dims(0),phis_in (ncol|: ), phis_out(ncol|: ), yy(ncol|: ))
          delete(dims)
        end if
        if(field@post_process(ifield) .eq. "cwat_adjust") then
          print ("     Adjusting "+field(ifield)+" to eliminate negative values")
          dims = dimsizes( yy )
          MAKEIC::myminmax_se(dims(0) ,dims(1) ,yy, 0., 1.e36)
          delete(dims)
        end if
        if(field@post_process(ifield) .eq. "cloud_fraction_adjust") then
          print ("     Adjusting "+field(ifield)+" between 0. and 1.")
          dims = dimsizes( yy )
          MAKEIC::myminmax_se(dims(0) ,dims(1) ,yy, 0., 1.   )
          delete(dims)
        end if

        if((vert_interp(ifield).eq."quad_opt1"     ) .and. \
           ( dimensions(ifield).eq."3D"            ) .or.  \
           (vert_interp(ifield).eq."quad_logp_opt1") .and. \
           ( dimensions(ifield).eq."3D"            )      ) then
          zz  =new( (/plevo , ncol  /), double)
          zz!0="lev"
          zz!1="ncol"
          dims=dimsizes( yy )
          if(vert_interp(ifield) .eq. "quad_opt1") then
            print("     QUADRATIC (in pressure) vertical interpolation for ")
            print("               this field (OPTION 1)")
            loglin = 1
          end if
          if(vert_interp(ifield) .eq. "quad_logp_opt1") then
            print("     QUADRATIC (in log pressure) vertical interpolation for ")
            print("               this field (OPTION 1)")
            loglin = 0
          end if
          MAKEIC::vert_quad_opt1_se(dims(0),(dims(0)+1),plevo,dims(1),  \
                                               yy(ncol|: ,lev|: ),      \
                                      press_in_yy(ncol|: ,lev|: ),      \
                                 press_in_yy_half(ncol|: ,lev|: ),      \
                                     press_out_yy(ncol|: ,lev|: ),      \
                                          phis_in(ncol|: ),             \
                                   ps_yy_no_spect(ncol|: ),             \
                                               zz(ncol|: ,lev|: ),loglin)
          delete(dims)
          zz&lev=cdf->lev
        else
          zz=vert_process(yy,cdf,field,ifield,press_in_yy,press_out_yy)
        end if

        if(field@post_process(ifield) .eq. "cwat_adjust") then
          print ("     Adjusting "+field(ifield)+" to eliminate negative values")
          dims = dimsizes( zz )
          MAKEIC::myminmax_se(dims(0) ,dims(1) ,zz, 0., 1.e36)
          delete(dims)
        end if
        if(field@post_process(ifield) .eq. "cloud_fraction_adjust") then
          print ("     Adjusting "+field(ifield)+" between 0. and 1.")
          dims = dimsizes( zz )
          MAKEIC::myminmax_se(dims(0) ,dims(1) ,zz, 0., 1.   )
          delete(dims)
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

    end do ; ifield=0,(nfields-1)

  end do ; ifile=0,(nfiles-1)

  ;----------------------------
  ; Misc. Cleanup
  ;----------------------------
  delete(ESMFopt)

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
