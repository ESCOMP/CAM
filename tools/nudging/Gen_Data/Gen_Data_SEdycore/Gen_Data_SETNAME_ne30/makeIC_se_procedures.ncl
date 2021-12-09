undef("comma_delimit")
undef("stringtological")
undef("date2YMDH_str")
undef("choose_file")
undef("date_to_time_365")
undef("date_to_time")
undef("vert_process")
undef("build_press")
undef("Initialize_diag_settings")
undef("build_filename")
undef("get_file")
undef("initialize_input_field_names")
undef("initialize_IC")
undef("initialize_IC_levels")
undef("initialize_IC_dates")
undef("find_field_index")
undef("flip_longitudes")
undef("flip_latitudes")
undef("flip_levels")
undef("read_reformat")
undef("horz_process")
undef("spectral_tr")
undef("spectral_v_tr")
undef("output2CAM")
undef("process_Q_horz")
undef("process_T_vert")
undef("process_Q_vert")
undef("horz_process_v")
undef("ESMF_pole_fix")


;=====================================================================================
function comma_delimit(strng:string)
;
;************************************************************************;
;                                                                        ;
; comma_delimit                                                          ;
;   Add comma at beginning and end of string that is intended to be a    ;
;   list of strings                                                      ;
;                                                                        ;
;   strng:  string;  string which is intended to be a list of strings    ;
;                                                                        ;
;************************************************************************;

local string_out

begin
  string_out = ","+strng+","
  return (string_out)
end
;=====================================================================================


;=====================================================================================
function stringtological(strng:string)
;
;************************************************************************;
;                                                                        ;
; stringtological                                                        ;
;   Convert string containing "False" or "True" to logical               ;
;   strng:  string;  string containing "False" or "True"                 ;
;                                                                        ;
;************************************************************************;

local logical_out

begin
  ;--------------------------------------------------
  ; Check that string contains only "False" or "True"
  ;--------------------------------------------------
  if(strng .ne. "False" .and. strng .ne. "True") then
    print("Error in 'stringtological'")
    print("String needs to be set to 'False' or 'True'")
    print("Currently, string is set to:  "+strng)
    exit
  end if

  logical_out = where(strng .eq. "True",True,False)

  return (logical_out)
end
;=====================================================================================


;=====================================================================================
procedure date2YMDH_str(fdate:double  , year_str:string, month_str:string, \
                        day_str:string, hour_str:string                    )
;
;************************************************************************;
;                                                                        ;
; date2YMDH_str                                                          ;
;                                                                        ;
;            fdate     : double ; date yyyymmddhhhhh                     ;
;            year_str  : string ; date yyyy                              ;
;            month_str : string ; date mm                                ;
;            day_str   : string ; date dd                                ;
;            hour_str  : string ; date hhhhh                             ;
;                                                                        ;
; Split "fdate" into year, month, day, hour strings                      ;
;                                                                        ;
;************************************************************************;

local year, month, day, hour, year_char, month_char, day_char, hour_char, \
      tmp1, tmp2, tmp3, tmp4, tmp5, dyear, dmonth, dday

begin
  dyear  = fdate
  dmonth = fdate
  dday   = fdate

  ;----------------------------------------------------------
  ; Split date record into year, month, day, and hour strings
  ;----------------------------------------------------------
  year   = doubletointeger(fdate/1000000000.D)
  dyear  = year
  month  = doubletointeger(fdate - dyear*1000000000.D)
  month  = month/10000000
  dmonth = month
  day    = doubletointeger(fdate - dyear*1000000000.D - dmonth*10000000.D)
  day    = day/100000
  dday   = day
  hour   = doubletointeger(fdate - dyear*1000000000.D - dmonth*10000000.D - dday*100000.D)

  ;---------------
  ; Create strings
  ;---------------
  year_char     = new( (/4/) , character)
  tmp1          = year/1000
  year_char(0)  = integertochar(tmp1+48)
  tmp2          = (year - tmp1*1000)/100
  year_char(1)  = integertochar(tmp2+48)
  tmp3          = (year - tmp1*1000 - tmp2*100)/10
  year_char(2)  = integertochar(tmp3+48)
  tmp4          =  year - tmp1*1000 - tmp2*100 - tmp3*10
  year_char(3)  = integertochar(tmp4+48)

  month_char    = new( (/2/) , character)
  tmp1          = month/10
  month_char(0) = integertochar(tmp1+48)
  tmp2          = month - tmp1*10
  month_char(1) = integertochar(tmp2+48)
  day_char      = new( (/2/) , character)
  tmp1          = day/10
  day_char(0)   = integertochar(tmp1+48)
  tmp2          = day - tmp1*10
  day_char(1)   = integertochar(tmp2+48)

  hour_char     = new( (/5/) , character)
  tmp1          = hour/10000
  hour_char(0)  = integertochar(tmp1+48)
  tmp2          = (hour - tmp1*10000)/1000
  hour_char(1)  = integertochar(tmp2+48)
  tmp3          = (hour - tmp1*10000 - tmp2*1000)/100
  hour_char(2)  = integertochar(tmp3+48)
  tmp4          = (hour - tmp1*10000 - tmp2*1000 - tmp3*100)/10
  hour_char(3)  = integertochar(tmp4+48)
  tmp5          =  hour - tmp1*10000 - tmp2*1000 - tmp3*100 - tmp4*10
  hour_char(4)  = integertochar(tmp5+48)

  year_str      = chartostring(year_char )
  month_str     = chartostring(month_char)
  day_str       = chartostring(day_char  )
  hour_str      = chartostring(hour_char )

  return
end
;=====================================================================================


;=====================================================================================
function choose_file(data0:file, data1:file, data2:file,  \
                     data3:file, data4:file, data5:file,  \
                     data6:file, data7:file, ifile:integer)
;
;************************************************************************;
;                                                                        ;
; choose_file                                                            ;
;                                                                        ;
;     data0-data7 : file    ; Set of input files                         ;
;     ifile       : integer ; input file number                          ;
;                                                                        ;
; Choose file from list                                                  ;
;                                                                        ;
;************************************************************************;

local data, processed

begin
  processed = False

  if( ifile .eq. 0) then
    data      = data0
    processed = True
  end if

  if( ifile .eq. 1) then
    data      = data1
    processed = True
  end if

  if( ifile .eq. 2) then
    data      = data2
    processed = True
  end if

  if( ifile .eq. 3) then
    data      = data3
    processed = True
  end if

  if( ifile .eq. 4) then
    data      = data4
    processed = True
  end if

  if( ifile .eq. 5) then
    data      = data5
    processed = True
  end if

  if( ifile .eq. 6) then
    data      = data6
    processed = True
  end if

  if( ifile .eq. 7) then
    data      = data7
    processed = True
  end if

  if(.not. processed) then
    print("Error:  No file chosen from list")
    exit
  end if

  return(data)
end
;=====================================================================================


;=====================================================================================
function date_to_time_365(date:double)
;
;************************************************************************;
;                                                                        ;
; date_to_time                                                           ;
;                                                                        ;
;            date: double;  date (YYYYMMDDHHHHH)                         ;
;                                                                        ;
; Construct "absolute" time (in days) based on input date                ;
;                                                                        ;
;                                                                        ;
;************************************************************************;

local  m, dayspm, year, month, day, fraction, tmp, time

begin
  fraction = new(  1, double )
  dayspm   = new( 12, double )
  dayspm   = (/ 31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31. /)

  ;----------------------------------
  ; Compute "absolute" time (in days)
  ;----------------------------------
  time = new( 1, double)
  tmp  = new( 1, double)

  if(date .gt. 2100000000000. ) then
     print("Error:  Date must be in form of YYYYMMDDHHHHH")
     print("Date is: "+date)
     exit
  end if

  year     = doubletointeger(date/1000000000)
  tmp      = year*100
  month    = doubletointeger(date/10000000 - tmp)
  tmp      = tmp *100    + month*100
  day      = doubletointeger(date/100000   - tmp)
  tmp      = tmp *100000 + day  *100000
  fraction = (date - tmp)/86400.

  if(year .gt. 0) then
    if(month .lt. 1 .or. month .gt. 12) then
      print("Error:  range of months must be between 1 and 12")
      print("month = "+month)
      exit
    end if
  end if

  if(month .gt. 0) then
    if(day .lt. 1 .or. day .gt. dayspm(month-1)) then
      print("Error:  day must be between 1 and "+dayspm(month-1))
      print("        for month "+month)
      print("day = "+day)
      exit
    end if
  end if

  time = year*365
  if(month .gt. 1) then
     do m = 0,month-2
        time = time + dayspm(m)
     end do
  end if
  time = time + day
  time = time + fraction

  return(time)
end
;=====================================================================================


;=====================================================================================
function date_to_time(date:double)
;
;************************************************************************;
;                                                                        ;
; date_to_time                                                           ;
;                                                                        ;
;            date: double;  date (YYYYMMDDHHHHH)                         ;
;                                                                        ;
; Construct "absolute" time (in days) based on input date                ;
;                                                                        ;
;                                                                        ;
;************************************************************************;

local  m, dayspm, year, month, day, fraction, tmp, time

begin
  fraction = new(  1, double )
  dayspm   = new( 12, double )
  dayspm   = (/ 31.,28.,31.,30.,31.,30.,31.,31.,30.,31.,30.,31. /)

  ;----------------------------------
  ; Compute "absolute" time (in days)
  ;----------------------------------
  time = new( 1, double)
  tmp  = new( 1, double)

  if(date .gt. 2100000000000. ) then
     print("Error:  Date must be in form of YYYYMMDDHHHHH")
     print("Date is: "+date)
     exit
  end if

  year     = doubletointeger(date/1000000000)
  tmp      = year*100
  month    = doubletointeger(date/10000000 - tmp)
  tmp      = tmp *100    + month*100
  day      = doubletointeger(date/100000   - tmp)
  tmp      = tmp *100000 + day  *100000
  fraction = (date - tmp)/86400.

  if(year .eq. (4*doubletointeger(int2dble(year)/4.))) then
   dayspm(1) = 29.
  end if

  if(year .gt. 0) then
    if(month .lt. 1 .or. month .gt. 12) then
      print("Error:  range of months must be between 1 and 12")
      print("month = "+month)
      exit
    end if
  end if

  if(month .gt. 0) then
    if(day .lt. 1 .or. day .gt. dayspm(month-1)) then
      print("Error:  day must be between 1 and "+dayspm(month-1))
      print("        for month "+month)
      print("day = "+day)
      exit
    end if
  end if

  time = year*365.25
  if(month .gt. 1) then
     do m = 0,month-2
        time = time + dayspm(m)
     end do
  end if
  time = time + day
  time = time + fraction

  return(time)
end
;=====================================================================================


;=====================================================================================
function vert_process(xfield:numeric  , cdf:file       , \
                      field:string    , ifield:numeric , \
                      pres_in:numeric , pres_out:numeric )
;
;***********************************************************************;
;                                                                       ;
; vert_process                                                          ;
;                                                                       ;
;            xfield  : numeric; INPUT field                             ;
;            cdf     : file   ; NCL Netcdf OUTPUT file                  ;
;            field   : string ; Names of of OUTPUT field variables      ;
;            ifield  : integer; OUTPUT field index                      ;
;            pres_in : numeric; pressure field of input field to be     ;
;                               vertically interpolated                 ;
;            pres_out: numeric; pressure field to which input field must;
;                               be vertically interpolated              ;
;                                                                       ;
; Vertically interpolate/copy input field to output field               ;
;                                                                       ;
;***********************************************************************;

local yy, zz, l_interpv, vert_interp, dims, dim_lev, loglin

begin
  ;------------------
  ; Check if 3D field
  ;------------------
  if(field@dimensions(ifield) .ne. "3D") then
    return(xfield)
  end if

  l_interpv   = False
  vert_interp = field@vert_interp(ifield)
  yy          = xfield

  ;-----------------------------------
  ; Determine dimensions of input file
  ;-----------------------------------
  dims = dimsizes(yy)
  dim_lev = dims(0)

  ;----------------------------------------------------
  ; Process input data according to interpolation flags
  ;----------------------------------------------------
  zz      = new( (/field@plevo, dims(1) /), double)
  zz!0    = "lev"
  zz!1    = "ncol"

  ;--------------------------
  ; No vertical interpolation
  ;--------------------------
  if(vert_interp .eq. "no_interp") then
    zz    = yy
    l_interpv = True
  end if

  ;------------------------------------------------
  ; Linear ( in P or ln(P) ) vertical interpolation
  ;------------------------------------------------
  if(vert_interp .eq. "linear" .or. vert_interp .eq. "linear_logp") then
    if(vert_interp .eq. "linear") then
      print("     LINEAR (in pressure) vertical interpolation for this field")
      loglin = -1
    end if
    if(vert_interp .eq. "linear_logp") then
      print("     LINEAR (in log pressure) vertical interpolation for this field")
      loglin = -2
    end if
    zz(ncol|: ,lev|field@plevo-1:0 ) = \
                                  int2p(pres_in (ncol|: ,lev|dim_lev-1:0     ) , \
                                        yy      (ncol|: ,lev|dim_lev-1:0     ) , \
                                        pres_out(ncol|: ,lev|field@plevo-1:0 ) ,loglin)
    print("       Warning:  int2p can possibly return missing values")
    l_interpv = True
  end if

  if(vert_interp .eq. "linear_opt1" .or. vert_interp .eq. "linear_logp_opt1") then
    if(vert_interp .eq. "linear_opt1") then
      print("     LINEAR (in pressure) vertical interpolation for this field (OPTION 1)")
      loglin = 1
    end if
    if(vert_interp .eq. "linear_logp_opt1") then
      print("     LINEAR (in log pressure) vertical interpolation for this field (OPTION 1)")
      loglin = 0
    end if
    MAKEIC::vert_int_opt1_se(dims(1),dim_lev,field@plevo, \
                             pres_in (ncol|: ,lev|: ), pres_out (ncol|: ,lev|: ), \
                             yy      (ncol|: ,lev|: ), zz       (ncol|: ,lev|: ), \
                             loglin )
    l_interpv = True
  end if

  if(vert_interp .eq. "linear_opt2" .or. vert_interp .eq. "linear_logp_opt2") then
    if(vert_interp .eq. "linear_opt2") then
      print("     LINEAR (in pressure) vertical interpolation for this field (OPTION 2)")
      loglin = 1
    end if
    if(vert_interp .eq. "linear_logp_opt2") then
      print("     LINEAR (in log pressure) vertical interpolation for this field (OPTION 2)")
      loglin = 0
    end if
    MAKEIC::vert_int_opt2_se(dims(1),dim_lev,field@plevo, \
                             pres_in (ncol|: ,lev|: ), pres_out (ncol|: ,lev|: ), \
                             yy      (ncol|: ,lev|: ), zz       (ncol|: ,lev|: ), \
                             loglin )
    l_interpv = True
  end if

  ;----------------------
  ; Add coordinate values
  ;----------------------
  zz&lev = field@lev

  ;------------------
  ; Modify attributes
  ;------------------
  if(isatt(zz,"_FillValue")) then
    delete(zz@_FillValue)
  end if

  if(.not. l_interpv) then
    print("       Error:  field not processed; '"+vert_interp+ \
                          "' is not a valid interpolation flag")
    exit
  end if

  delete(yy)

  return(zz)
end
;=====================================================================================


;=====================================================================================
function build_press(data:file, ftype:string, ps:numeric, level_type:string)
;
;************************************************************************;
;                                                                        ;
; build_press                                                            ;
;                                                                        ;
;            data    : file   ; Analysis dataset                         ;
;            ftype   : string ; file type                                ;
;            ps      : numeric; Surface pressure                         ;
;            level_type: string; half/full levels                        ;
;                                                                        ;
; Build 3-D pressure field based upon Ps and A's and B's                 ;
;                                                                        ;
;************************************************************************;

local press, hya, hyb, p0, dims, press_i, analysis_info

begin
  p0 = 100000.

  if(level_type .ne. "FULL" .and. level_type .ne. "HALF") then
    print("Error (function build_press):  Must specify 'FULL' or 'HALF' in the last argument of this function")
    exit
  end if

  ;------------------------------
  ; Read in A's and B's from file
  ;------------------------------
  analysis_info = extract_analyses_info(data, ftype)
  if(level_type .eq. "FULL") then
    hya = analysis_info@hyam
    hyb = analysis_info@hybm
  end if
  if(level_type .eq. "HALF") then
    hya = analysis_info@hyai
    hyb = analysis_info@hybi
  end if

  ;---------------------------------
  ; Compute pressures for all levels
  ;---------------------------------
  press = pres_hybrid_ccm(ps, p0, hya, hyb)

  ;------------------------------------------------------------
  ; Special case:  for JRA_25, full-level pressures are derived
  ; from interface pressures
  ;------------------------------------------------------------
  if(ftype .eq. "JRA_25" .and. level_type .eq. "FULL") then
    dims    = dimsizes(press)
    delete(press)
    delete(hya)
    delete(hyb)
    press   = new( (/ dims(0)  , dims(1), dims(2) /), double)
    press_i = new( (/ dims(0)+1, dims(1), dims(2) /), double)
    hya     = analysis_info@hyai
    hyb     = analysis_info@hybi
    press_i = pres_hybrid_ccm(ps, p0, hya, hyb)
    MAKEIC::jra_25_press_full_levels(dims(0)+1, dims(0), dims(1), dims(2), press_i, press  )
    delete (dims)
    delete (press_i)
  end if

  press!0 = "lev"
  press!1 = ps!0
  press!2 = ps!1

  return(press)
end
;=====================================================================================


;=====================================================================================
function Initialize_diag_settings(NameNumber:integer)
;
;************************************************************************;
;                                                                        ;
; Initialize_diag_settings                                               ;
;   Diagnostic Settings for processing all fields                        ;
;   "field" is a master variable carrying meta-data needed for file      ;
;   processing                                                           ;
;                                                                        ;
;************************************************************************;

local MYTMPDIR,MYOUTDIR,INPUTDIR, \
      fields,source_files,REF_DATE, \
      fname,fname0,fname1,fname2,fname3,fname4,fname5, \
      ftype,ftype0,ftype1,ftype2,ftype3,ftype4,ftype5, \
      fdate,fdate0,fdate1,fdate2,fdate3,fdate4,fdate5, \
      fname_phis_input,ftype_phis_input,fname_phis_output,ftype_phis_output,fname_grid_info, \
      CASE,DYCORE,PRECISION,PTRM,PTRN,PTRK,MASS_FIX, \
      VORT_DIV_TO_UV,SST_MASK,ICE_MASK,OUTPUT_PHIS,REGRID_ALL,ADJUST_STATE_FROM_TOPO, \
      ii,field,dim_fname,end_of_array,nfiles,field_dim,tmp, \
      tmp_source_file, tmp_fname, tmp_ftype, tmp_fdate

begin
  ;--------------------------------------
  ; Set values rather than from getenv()
  ;--------------------------------------
  DMYSTR=oneDtostring(new(180,character))
  MYTMPDIRx              =DMYSTR
  MYOUTDIRx              =DMYSTR
  INPUTDIRx              =DMYSTR
  ESMF_interpx           =DMYSTR
  ESMF_polex             =DMYSTR
  ESMF_cleanx            =DMYSTR
  TMP_cleanx             =DMYSTR
  REF_DATEx              =DMYSTR
  CASEx                  =DMYSTR
  DYCOREx                =DMYSTR
  PRECISIONx             =DMYSTR
  VORT_DIV_TO_UVx        =DMYSTR
  SST_MASKx              =DMYSTR
  ICE_MASKx              =DMYSTR
  OUTPUT_PHISx           =DMYSTR
  REGRID_ALLx            =DMYSTR
  ADJUST_STATE_FROM_TOPOx=DMYSTR
  MASS_FIXx              =DMYSTR
  fname_phis_outputx     =DMYSTR
  ftype_phis_outputx     =DMYSTR
  fname_grid_infox       =DMYSTR
  fieldsx                =DMYSTR
  source_filesx          =DMYSTR
  fname_phis_inputx      =DMYSTR
  ftype_phis_inputx      =DMYSTR
  fname0x                =DMYSTR
  fname1x                =DMYSTR
  fname2x                =DMYSTR
  fname3x                =DMYSTR
  fname4x                =DMYSTR
  fname5x                =DMYSTR
  ftype0x                =DMYSTR
  ftype1x                =DMYSTR
  ftype2x                =DMYSTR
  ftype3x                =DMYSTR
  ftype4x                =DMYSTR
  ftype5x                =DMYSTR
  fdate0x                =DMYSTR
  fdate1x                =DMYSTR
  fdate2x                =DMYSTR
  fdate3x                =DMYSTR
  fdate4x                =DMYSTR
  fdate5x                =DMYSTR
  MAKEIC::read_namelist(NameNumber,MYTMPDIRx,MYOUTDIRx,INPUTDIRx,              \
                        ESMF_interpx, ESMF_polex, ESMF_cleanx, TMP_cleanx,     \
                        REF_DATEx, CASEx, DYCOREx, PRECISIONx, VORT_DIV_TO_UVx,\
                        SST_MASKx, ICE_MASKx, OUTPUT_PHISx, REGRID_ALLx,       \
                        ADJUST_STATE_FROM_TOPOx, MASS_FIXx, fname_phis_outputx,\
                        ftype_phis_outputx, fname_grid_infox, fieldsx,         \
                        source_filesx, fname_phis_inputx, ftype_phis_inputx,   \
                        fname0x, fname1x, fname2x, fname3x, fname4x, fname5x,  \
                        ftype0x, ftype1x, ftype2x, ftype3x, ftype4x, ftype5x,  \
                        fdate0x, fdate1x, fdate2x, fdate3x, fdate4x, fdate5x)
  MYTMPDIR              =str_squeeze(MYTMPDIRx)
  MYOUTDIR              =str_squeeze(MYOUTDIRx)
  INPUTDIR              =str_squeeze(INPUTDIRx)
  ESMF_interp           =str_squeeze(ESMF_interpx)
  ESMF_pole             =str_squeeze(ESMF_polex)
  ESMF_clean            =str_squeeze(ESMF_cleanx)
  TMP_clean             =str_squeeze(TMP_cleanx)
  REF_DATE              =str_squeeze(REF_DATEx )
  CASE                  =str_squeeze(CASEx)
  DYCORE                =str_squeeze(DYCOREx)
  PRECISION             =str_squeeze(PRECISIONx)
  VORT_DIV_TO_UV        =str_squeeze(VORT_DIV_TO_UVx)
  SST_MASK              =str_squeeze(SST_MASKx)
  ICE_MASK              =str_squeeze(ICE_MASKx)
  OUTPUT_PHIS           =str_squeeze(OUTPUT_PHISx)
  REGRID_ALL            =str_squeeze(REGRID_ALLx)
  ADJUST_STATE_FROM_TOPO=str_squeeze(ADJUST_STATE_FROM_TOPOx)
  MASS_FIX              =str_squeeze(MASS_FIXx)
  fname_phis_output     =str_squeeze(fname_phis_outputx)
  ftype_phis_output     =str_squeeze(ftype_phis_outputx)
  fname_grid_info       =str_squeeze(fname_grid_infox)
  fields                =str_squeeze(fieldsx)
  source_files          =str_squeeze(source_filesx)
  fname_phis_input      =str_squeeze(fname_phis_inputx)
  ftype_phis_input      =str_squeeze(ftype_phis_inputx)
  fname0                =str_squeeze(fname0x)
  fname1                =str_squeeze(fname1x)
  fname2                =str_squeeze(fname2x)
  fname3                =str_squeeze(fname3x)
  fname4                =str_squeeze(fname4x)
  fname5                =str_squeeze(fname5x)
  ftype0                =str_squeeze(ftype0x)
  ftype1                =str_squeeze(ftype1x)
  ftype2                =str_squeeze(ftype2x)
  ftype3                =str_squeeze(ftype3x)
  ftype4                =str_squeeze(ftype4x)
  ftype5                =str_squeeze(ftype5x)
  fdate0                =str_squeeze(fdate0x)
  fdate1                =str_squeeze(fdate1x)
  fdate2                =str_squeeze(fdate2x)
  fdate3                =str_squeeze(fdate3x)
  fdate4                =str_squeeze(fdate4x)
  fdate5                =str_squeeze(fdate5x)


;  print(MYTMPDIR)
;  print(MYOUTDIR)
;  print(INPUTDIR)
;  print(CASE)
;  print(PRECISION)
;  print(fname0)
;  print(fname1)
;  print(fname2)
;  print(fname3)
;  print(fname4)
;  print(fname5)
;
;
;;  MYTMPDIR              ="./TMP/"
;;  MYOUTDIR              ="./YOTC-Output/"
;;  INPUTDIR              ="../YOTC-Input/"
;;  REF_DATE              ="20080501"
;;  CASE                  ="YOTC2CAM_ne30np4_L30"
;;  DYCORE                ="se"
;;  PRECISION             ="float"
;;  VORT_DIV_TO_UV        ="False"
;;  SST_MASK              ="False"
;;  ICE_MASK              ="False"
;;  OUTPUT_PHIS           ="True"
;;  REGRID_ALL            ="False"
;;  ADJUST_STATE_FROM_TOPO="True"
;;  MASS_FIX              ="True"

;;  fname_phis_output     ="../Grid-Reference/USGS-gtopo30_ne30np4_16xdel2-PFC-consistentSGH.nc"
;;  ftype_phis_output     ="SE_TOPOGRAPHY"
;;  fname_grid_info       ="../Grid-Reference/cami_0000-01-01_ne30np4_L30_c120315.nc"

;;  fname0="2008050100_000_018_130128_132128_799.yt.oper.an.ml"
;;  fname1="2008050100_000_018_129128_152128_799.yt.oper.an.ml"
;;  fname2="2008050100_018_018_133128_248128_400.yt.oper.an.ml"
;;  fname3="none"
;;  fname4="none"
;;  fname5="none"

;;  ftype0="YOTC_sh"
;;  ftype1="YOTC_PS_Z"
;;  ftype2="ECMWF_gg"
;;  ftype3="none"
;;  ftype4="none"
;;  ftype5="none"

;;  fdate0="2008050164800"
;;  fdate1="2008050164800"
;;  fdate2="2008050164800"
;;  fdate3="-1"
;;  fdate4="-1"
;;  fdate5="-1"

;;  fields      ="U,V,T,PS,Q"
;;  source_files="0,0,0,1,2"
;;  fname_phis_input=fname1
;;  ftype_phis_input=ftype1

  fname=fname0+","+fname1+","+fname2+","+ fname3+","+fname4+","+fname5
  ftype=ftype0+","+ftype1+","+ftype2+","+ftype3+","+ftype4+","+ftype5
  ftype=ftype0+","+ftype1+","+ftype2+","+ftype3+","+ftype4+","+ftype5
  fdate=fdate0+","+fdate1+","+fdate2+","+fdate3+","+fdate4+","+fdate5

;printVarSummary(fname)
;print(fname)

  PTRM     ="42"
  PTRN     ="42"
  PTRK     ="42"

  ;--------------------------------------------------
  ; Now Process input values into 'field' variable
  ;--------------------------------------------------
  field                  =cshstringtolist(comma_delimit(fields))
  field@fname            =cshstringtolist(comma_delimit(fname))
  field@ftype            =cshstringtolist(comma_delimit(ftype))
  field@fdate            =stringtodouble(cshstringtolist(comma_delimit(fdate)))
  field@source_file      =stringtointeger(cshstringtolist(comma_delimit(source_files)))
  field@fname_phis_input =INPUTDIR+fname_phis_input
  field@ftype_phis_input =ftype_phis_input
  field@fname_phis_output=fname_phis_output
  field@ftype_phis_output=ftype_phis_output
  field@fname_grid_info   =fname_grid_info

  field@ref_date         =stringtointeger(REF_DATE)
  field@mytmpdir         =MYTMPDIR
  field@myoutdir         =MYOUTDIR
  field@inputdir         =INPUTDIR
  field@esmf_interp      =ESMF_interp
  field@esmf_pole        =ESMF_pole
  field@esmf_clean       =ESMF_clean
  field@tmp_clean        =TMP_clean

  field@dycore                =DYCORE
  field@case                  =CASE
  field@precision             =PRECISION
  field@ptrmo                 =stringtointeger(PTRM)
  field@ptrno                 =stringtointeger(PTRN)
  field@ptrko                 =stringtointeger(PTRK)
  field@vort_div_to_UV        =stringtological(VORT_DIV_TO_UV)
  field@sst_mask              =stringtological(SST_MASK)
  field@ice_mask              =stringtological(ICE_MASK)
  field@output_phis           =stringtological(OUTPUT_PHIS)
  field@regrid_all_cam_fields =stringtological(REGRID_ALL)
  field@adjust_state_from_topo=stringtological(ADJUST_STATE_FROM_TOPO)
  field@mass_fix              =stringtological(MASS_FIX)

  field@spectral_trunc_setting= field@ptrmo

  ;--------------------------------------------------------------
  ; Only SE dycore is allowed for this version
  ;--------------------------------------------------------------
  if( field@dycore .ne. "se" ) then
    print("Error in 'Initialize_data_processing_settings'")
    print("  Only DYCORE = 'se' is allowed for this version ")
    exit
  end if

  ;--------------------------------------------------------------
  ; Make sure at least U, V, T, Q, and PS are in the "field" list
  ;--------------------------------------------------------------
  if( (.not. any(field .eq. "U" )) .or. (.not. any(field .eq. "V")) .or. \
      (.not. any(field .eq. "T" )) .or. (.not. any(field .eq. "Q")) .or. \
      (.not. any(field .eq. "PS"))                                       ) then
    print("Error in 'Initialize_data_processing_settings'")
    print("  List of 'FIELDS' does not contain the minimum required fields:  U, V, T, Q, and PS")
    print("  Current list contains:  "+field)
    exit
  end if

  ;--------------------------
  ; Find end of list of files
  ;--------------------------
  dim_fname    = dimsizes(field@fname)
  end_of_array = False
  ii           = -1
  nfiles       = dim_fname
  do while((.not. end_of_array).and.(ii .le. dim_fname-2))
    ii = ii + 1
    if(field@fname(ii) .eq. "none") then
      nfiles = ii
      end_of_array = True
    end if
    if(field@fname(ii) .ne. "none") then
     field@fname(ii)=field@inputdir+field@fname(ii)
     print(field@fname(ii))
    end if
  end do

  if(nfiles .lt. 1) then
    print("Error in 'Initialize_data_processing_settings'")
    print("Asking for no input files:  "+nfiles)
    print("Script requires at least 1 input file")
    exit
  end if

  if(nfiles .gt. 6) then
    print("Error in 'Initialize_data_processing_settings'")
    print("Asking for too many input files:  "+nfiles)
    print("Max is 6")
    exit
  end if

  ;------------------------------------------------------------------
  ; Check that dimension of lists of each file attribute is at 
  ; least as large as the number of input files requested
  ;------------------------------------------------------------------
  if(dimsizes(field@ftype) .lt. nfiles) then
    print("Error in 'Initialize_data_processing_settings'")
    print("Dimension of 'ftype' needs to be at least as large as 'fname'")
    print("Dimension of 'ftype' = "+dimsizes(field@ftype))
    print("Dimension of 'fname'  = "+nfiles)
    exit
  end if

  if(dimsizes(field@fdate) .lt. nfiles) then
    print("Error in 'Initialize_data_processing_settings'")
    print("Dimension of 'fdate' needs to be at least as large as 'fname'")
    print("Dimension of 'fdate' = "+dimsizes(field@fdate))
    print("Dimension of 'fname'  = "+nfiles)
    exit
  end if

  if(dimsizes(field@source_file) .ne. dimsizes(field)) then
    print("Error in 'Initialize_data_processing_settings'")
    print("Dimension of 'source_files' needs to be the same as 'fields'")
    print("Dimension of 'source_files' = "+dimsizes(field@source_file))
    print("Dimension of 'fields'       = "+dimsizes(field      ))
    exit
  end if

  ;--------------------------------------------------------------------------
  ; Add space for SE topography file which provides the horizontal 
  ; output grid as well as "PHIS_output" if needed
  ;--------------------------------------------------------------------------
  if(ismissing(field@fname_phis_output) .or. ismissing(field@ftype_phis_output)) then
    print("")
    print("Error in 'Initialize_data_processing_settings'")
    print("")
    print("  Output Topography file is required:")
    print("        'fname_phis_output' = "+field@fname_phis_output)
    print("        'ftype_phis_output' = "+field@ftype_phis_output)
    exit
  end if

  field_dim                              = dimsizes(field)
  tmp                                    = field
  tmp_source_file                        = field@source_file
  tmp_fname                              = field@fname(0:nfiles-1)
  tmp_ftype                              = field@ftype(0:nfiles-1)
  tmp_fdate                              = field@fdate(0:nfiles-1)
  delete(field)
  field                                  = new( (/field_dim+1/), string )
  field(0:field_dim-1)                   = tmp
  field(field_dim:field_dim)             = (/ "PHIS" /)

  delete(field@source_file)
  delete(field@fname)
  delete(field@ftype)
  delete(field@fdate)

  field@source_file                      = new( (/field_dim+1/), integer)
  field@fname                            = new( (/nfiles   +1/), string )
  field@ftype                            = new( (/nfiles   +1/), string )
  field@fdate                            = new( (/nfiles   +1/), double )

  field@source_file(0:field_dim-1)       = tmp_source_file
  field@fname      (0:nfiles-1   )       = tmp_fname
  field@ftype      (0:nfiles-1   )       = tmp_ftype
  field@fdate      (0:nfiles-1   )       = tmp_fdate

  field@source_file(field_dim:field_dim) = (/ nfiles                  /)
  field@fname      (nfiles:nfiles      ) = (/ field@fname_phis_output /)
  field@ftype      (nfiles:nfiles      ) = (/ field@ftype_phis_output /)
  field@fdate      (nfiles:nfiles      ) = (/ 0                       /)

  field@nfiles                           = nfiles+1
  nfiles=nfiles+1

  ;--------------------------------------------------------------------------
  ; Add space for "PHIS_input" to field and file lists if required
  ;--------------------------------------------------------------------------
  if(field@adjust_state_from_topo) then
    if(ismissing(field@fname_phis_input ) .or. ismissing(field@ftype_phis_input)) then
      print("")
      print("Error in 'Initialize_data_processing_settings'")
      print("")
      print("        If 'adjust_state_from_topo' is True, then ")
      print("           'fname_phis_input' and 'ftype_phis_input' must be set")
      print("")
      print("        'fname_phis_input' = "+field@fname_phis_input)
      print("        'ftype_phis_input' = "+field@ftype_phis_input)
      exit
    end if

    delete(field_dim)
    delete(tmp)
    delete(tmp_source_file)
    delete(tmp_fname)
    delete(tmp_ftype)
    delete(tmp_fdate)

    field_dim                              = dimsizes(field)
    tmp                                    = field
    tmp_source_file                        = field@source_file
    tmp_fname                              = field@fname(0:nfiles-1)
    tmp_ftype                              = field@ftype(0:nfiles-1)
    tmp_fdate                              = field@fdate(0:nfiles-1)
    delete(field)
    field                                  = new( (/field_dim+1/), string )
    field(0:field_dim-1)                   = tmp
    field(field_dim:field_dim)             = (/ "PHIS_input" /)

    delete(field@source_file)
    delete(field@fname)
    delete(field@ftype)
    delete(field@fdate)

    field@source_file                      = new( (/field_dim+1/), integer)
    field@fname                            = new( (/nfiles   +1/), string )
    field@ftype                            = new( (/nfiles   +1/), string )
    field@fdate                            = new( (/nfiles   +1/), double )

    field@source_file(0:field_dim-1)       = tmp_source_file
    field@fname      (0:nfiles-1   )       = tmp_fname
    field@ftype      (0:nfiles-1   )       = tmp_ftype
    field@fdate      (0:nfiles-1   )       = tmp_fdate

    field@source_file(field_dim:field_dim) = (/ nfiles                 /)
    field@fname      (nfiles:nfiles      ) = (/ field@fname_phis_input /)
    field@ftype      (nfiles:nfiles      ) = (/ field@ftype_phis_input /)
    field@fdate      (nfiles:nfiles      ) = (/ 0                      /)

    field@nfiles                             = nfiles+1
  end if

  return (field)
end
;=====================================================================================


;=====================================================================================
function build_filename(field:string, fdate:double )
;
;************************************************************************;
;                                                                        ;
; build_filename                                                         ;
;                                                                        ;
;            field   : string ; name of field to be read from input file ;
;            fdate   : double ; date yyyymmddhhhhh                       ;
;                                                                        ;
; Create file name from date and case information                        ;
;                                                                        ;
;************************************************************************;

local name, year_str, month_str, day_str, hour_str

begin
  ;----------------------------------------------------------
  ; Split date record into year, month, day, and hour strings
  ;----------------------------------------------------------
  year_str  = new( 1 , string )
  month_str = new( 1 , string )
  day_str   = new( 1 , string )
  hour_str  = new( 1 , string )
  date2YMDH_str(fdate, year_str, month_str, day_str, hour_str)
  name = field@case+".cam2.i."+year_str+"-"+month_str+"-"+day_str+"-"+hour_str+".nc"
  return(name)
end
;=====================================================================================


;=====================================================================================
function get_file( file_prefix_str:string, mytmpdir:string, fname:string )
;
;************************************************************************;
;                                                                        ;
; get_file                                                               ;
;   Copy file from local disk or HPSS                                    ;
;                                                                        ;
;   file_prefix_str  : string ; File prefix to add to local filename     ;
;   mytmpdir         : string ; Work directory of input files            ;
;   fname            : string ; Full pathname of file                    ;
;                                                                        ;
;************************************************************************;

local string_tmp, check_remote_file, check_local_file, ichar_ind, ichars, lfname

begin
  ;-------------------------------------
  ; Extract file name from full pathname
  ;-------------------------------------
  ichar_ind = str_index_of_substr(fname, "/", -1)
  ichars    = stringtochar(fname)
  lfname    = "file"+file_prefix_str+"_"+chartostring(ichars(ichar_ind+1:dimsizes(ichars)-1))

  ;-------------------------------------------------
  ; Check if file is already in input work directory
  ;-------------------------------------------------
  check_local_file = systemfunc("ls "+mytmpdir+"/"+lfname)

  ;--------------------------------
  ; If not, then access remote file
  ;--------------------------------
  if(ismissing(check_local_file) ) then

    ;----------------------------------------------------------------------------
    ; Check if remote file is on disk.  If not, then try accessing file from HPSS
    ;----------------------------------------------------------------------------
    check_remote_file = systemfunc("ls "+fname)

    if(ismissing(check_remote_file) ) then
      string_tmp = "hsi -q get "+          mytmpdir+"/"+lfname+" : "+fname
    else
      string_tmp = "/bin/cp -i "+fname+" "+mytmpdir+"/"+lfname
    end if

    print ("------------------------------------------------------------")
    print (" "+string_tmp)
    print ("------------------------------------------------------------")
    system(string_tmp)

    ;----------------------------------------------
    ; Check if file made it to input work directory
    ;----------------------------------------------
    delete(check_local_file)
    check_local_file = systemfunc("ls "+mytmpdir+"/"+lfname)
    if(ismissing(check_local_file) ) then
      print (" ")
      print (" ")
      print(" Attempt to copy file "+fname+" was unsuccessful")
      print (" ")
      print (" ")
      exit
    end if

  else
    print(" Copy not done; file "+mytmpdir+"/"+lfname+" already exists")
  end if

  print (" ")

  return (lfname)
end
;=====================================================================================


;=====================================================================================
function initialize_input_field_names(field_old:string, data0:file, data1:file, \
                                                        data2:file, data3:file, \
                                                        data4:file, data5:file, \
                                                        data6:file, data7:file )
;
;**************************************************************************;
;                                                                          ;
; initialize_input_field_names                                             ;
;                                                                          ;
;            field_old   : string ; Array of field variables to go in      ;
;                                   output file                            ;
;            data0-data7 : file   ; Set of input files                     ;
;                                                                          ;
; Modify list of output fields to include all fields of CAM file if CAM    ;
; file is "file0"                                                          ;
;                                                                          ;
; Initialize a list of input fields that maps to the list of ouput fields  ;
;                                                                          ;
;**************************************************************************;

local ii, mm, vv, ifield, ifile, lfname, ftype, data, Found, exit_script, input_fields, \
      field_dim, tmp, Master_List, var_atts, CAM_fields, CAM_select_fields,          \
      field_exclude_list, ind_extra_fields, dim_extra_fields, field, source_file

begin
  ;-------------------------------------------
  ; Copy "field_old" to local variable "field"
  ;-------------------------------------------
  field = field_old

  ;-----------------------------------------------------------------------------------
  ; List of fields NOT to include on list of *extra* CAM fields to be regridded if the
  ; input CAM file is "file0" (state variables are already on "field")
  ;-----------------------------------------------------------------------------------
  field_exclude_list = (/ "T","U","V","US","VS","Q","PS","PHIS","LANDM_COSLAT","LANDM",  \
                          "ICEFRAC","LANDFRAC","DELP","SGH","SGH30","lat","lon","area"  /)

  ;-----------------------------------------------------------
  ; If file 0 is a CAM file, then regrid ALL the 2-D and 3-D 
  ; fields in addition to the state fields
  ;-----------------------------------------------------------
  if(field@ftype(0) .eq. "CAM" .and. field@regrid_all_cam_fields) then
    CAM_fields         = getfilevarnames(data0)
    CAM_select_fields  = new( dimsizes(CAM_fields), logical )
    CAM_select_fields  = False

    do ii = 0,(dimsizes(CAM_fields)-1)
      if( .not. any(field_exclude_list .eq. CAM_fields(ii)) .and. \
          .not. any(field              .eq. CAM_fields(ii)) ) then
        if(  ( any(getfilevardims(data0, CAM_fields(ii)) .eq. "lon" )   .and. \
               any(getfilevardims(data0, CAM_fields(ii)) .eq. "lat" )   .and. \
         .not. any(getfilevardims(data0, CAM_fields(ii)) .eq. "ilev") ) .or.  \
             ( any(getfilevardims(data0, CAM_fields(ii)) .eq. "ncol")   .and. \
         .not. any(getfilevardims(data0, CAM_fields(ii)) .eq. "ilev") ) ) then
          CAM_select_fields(ii) = True
        end if
      end if
    end do ; ii = 0,(dimsizes(CAM_fields)-1)
    ind_extra_fields = ind(CAM_select_fields)
    dim_extra_fields = dimsizes(ind_extra_fields)

    ;------------------------------------------------------------
    ; Add extra input CAM fields onto already existing field list
    ;------------------------------------------------------------
    field_dim                                                 = dimsizes(field)
    tmp                                                       = field
    source_file                                               = field@source_file
    delete(field)
    field                                                     = new( field_dim+dim_extra_fields, string )
    field(0:field_dim-1)                                      = tmp
    field(field_dim:field_dim+dim_extra_fields-1)             = CAM_fields(ind_extra_fields)

    delete(field@source_file)
    field@source_file                                         = new( field_dim+dim_extra_fields, integer)
    field@source_file(0:field_dim-1)                          = source_file
    field@source_file(field_dim:field_dim+dim_extra_fields-1) = 0
    delete(tmp)

  end if

  ;--------------------------------------------------------------------------------------------------
  ; Declarations:  Set up processing meta-data for each variable to be processed in the "field" array
  ;--------------------------------------------------------------------------------------------------
  field_dim         =dimsizes(field)
  field@source_field=new( (/field_dim/), string )   ; input (analysis) field name
  field@dimensions  =new( (/field_dim/), string )   ; "2D" or "3D" designation extracted
  field@horz_interp =new( (/field_dim/), string )   ; flag to indicate which type of horizontal 
                                                    ; interpolation to use
  field@vert_interp =new( (/field_dim/), string )   ; flag to indicate which type of vertical
                                                    ; interpolation to use
  field@spec_trunc  =new( (/field_dim/), integer)   ; Spectral truncation to use (if Eulerian dycore)
  field@bin_factor  =new( (/field_dim/), double )   ; bin factor to use if horizontal 
                                                    ; interpolation is "binning" (typical value:  1.)
  field@processed   =new( (/field_dim/), logical)   ; Set to False until field is processed
  field@loutput     =new( (/field_dim/), logical)   ; Flag to output field (if False, field is 
                                                    ; used only for internal processing only)
  field@long_name   =new( (/field_dim/), string )   ; long name of field
  field@units       =new( (/field_dim/), string )   ; units of field
  field@min         =new( (/field_dim/), float  )   ; minimum possible value for this field 
                                                    ; (checked against final interpolated field)
  field@max         =new( (/field_dim/), float  )   ; maximum possible value for this field (
                                                    ; checked against final interpolated field)
  field@post_process=new( (/field_dim/), string )   ; post-processing flags

  do ii= 0,(field_dim-1)
    field@processed(ii)=False
  end do

  ;----------------------------------------------------
  ; Map input field names to list of output field names
  ;----------------------------------------------------
  exit_script = False
  do ifield = 0,(field_dim-1)
    Found        = False
    ifile        = field@source_file(ifield)
    lfname       = field@lfname(ifile)
    ftype        = field@ftype(ifile)
    data         = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)
    input_fields = get_input_fields(ftype)

    if(ftype .eq. "CAM") then
      tmp = field(ifield)

      ;--------------
      ; Special cases
      ;--------------
      if((tmp.eq."U").or.(tmp.eq."US").or.                    \
         (tmp.eq."V").or.(tmp.eq."VS").or.(tmp.eq."PHIS_input")) then

        if(tmp .eq. "U" .or. tmp .eq. "US") then
          if (isfilevar(data, "U")) then
            field@source_field(ifield) = "U"
            Found                      = True
          end if
          if (isfilevar(data, "US")) then
            field@source_field(ifield) = "US"
            Found                      = True
          end if
          if(.not. Found) then
            print("Error in 'initialize_input_field_names'")
            print("Could not find the fields U or US in file "+lfname)
            exit_script = True
          end if
        end if

        if(tmp .eq. "V" .or. tmp .eq. "VS") then
          if (isfilevar(data, "V")) then
            field@source_field(ifield) = "V"
            Found                      = True
          end if
          if (isfilevar(data, "VS")) then
            field@source_field(ifield) = "VS"
            Found                      = True
          end if
          if(.not. Found) then
            print("Error in 'initialize_input_field_names'")
            print("Could not find the fields V or VS in file "+lfname)
            exit_script = True
          end if
        end if

        if(tmp .eq. "PHIS_input") then
          if (isfilevar(data, "PHIS")) then
            field@source_field(ifield) = "PHIS"
          else
            print("Error in 'initialize_input_field_names'")
            print("Could not find the field PHIS in file "+lfname)
            exit_script = True
          end if
        end if
      else
        ;---------------------------------------------
        ; Input names are same as output names for CAM
        ;---------------------------------------------
        if (isfilevar(data, tmp)) then
          field@source_field(ifield) = tmp
        else
          print("Error in 'initialize_input_field_names'")
          print("Could not find the field "+tmp+" in file "+lfname)
          exit_script                = True
        end if
      end if
    else
      ;------------------------------------------------------------------
      ; For foreign analyses, find input names to map to output CAM names
      ;------------------------------------------------------------------
      tmp = input_fields@$field(ifield)$
      if (isfilevar(data, tmp)) then
        field@source_field(ifield) = tmp
      else
        print("Error in 'initialize_input_field_names'")
        print("Could not find the field "+tmp+" in file "+lfname)
        exit_script                = True
      end if
    end if

    delete(tmp)
    delete(Found)
    delete(ifile)
    delete(lfname)
    delete(ftype)
    delete (data)
    delete (input_fields)
  end do ; ifield = 0,(field_dim-1)

  ;-------------------------
  ; Exit if any errors found
  ;-------------------------
  if(exit_script) then
    exit
  end if

  field@nfields = dimsizes(field)

  ;----------------------------------------
  ; Fill "field" meta-data from Master list
  ;----------------------------------------
  Master_List = Create_field_Master_List(field@spectral_trunc_setting, field)
  var_atts = getvaratts(Master_List)

  do ii=0,(field@nfields-1)

    mm=ind(field(ii).eq.Master_List)
    ;-------------------------------------------------
    ; If no meta-data from Master List for this field:
    ;-------------------------------------------------
    if(ismissing(mm)) then

      ;---------------------------------------------------------------------------
      ; If file0 is a CAM file, create a meta-data template and fill in the blanks
      ;---------------------------------------------------------------------------
      if(field@ftype(0) .eq. "CAM") then
        print("Warning in 'Initialize_data_processing_settings'")
        print("Requested output field variable '"+field(ii)+"' not defined in 'Create_field_Master_List'")
        print("Using generic regridding on this field (binning in the horizontal, linear in ln(P) in vertical)")

        mm = ind("generic" .eq. Master_List)

        do vv=0,(dimsizes(var_atts)-1)
          field@$var_atts(vv)$(ii) = Master_List@$var_atts(vv)$(mm)
        end do

        ;-------------------------------------------------------------------------------
        ; Fill in blanks from the "generic" meta-data with meta-data from the input file
        ;-------------------------------------------------------------------------------
        if(isfilevaratt(data0, field(ii), "units"    )) then
          field@units    (ii)  = data0->$field(ii)$@units
        end if
        if(isfilevaratt(data0, field(ii), "long_name")) then
          field@long_name(ii)  = data0->$field(ii)$@long_name
        end if
        field@dimensions(ii)   = "2D"
        if( any(getfilevardims(data0, field(ii)) .eq. "lev" ) ) then
          field@dimensions(ii) = "3D"
        end if

        ;------------
        ; Else, abort
        ;------------
      else
        print("Error in 'Initialize_data_processing_settings'")
        print("Requested output field variable '"+field(ii)+"' not defined in 'Create_field_Master_List'")
        exit
      end if

      ;----------------------------------------------
      ; Otherwise, fill in meta-data from Master List
      ;----------------------------------------------
    else
      do vv=0,(dimsizes(var_atts)-1)
        field@$var_atts(vv)$(ii) = Master_List@$var_atts(vv)$(mm)
      end do
    end if
  end do ; ii=0,(field@nfields-1)

  return (field)
end
;=====================================================================================


;=====================================================================================
procedure initialize_IC(cdf:file           , field:string,   nfields:integer,      \
                        nfields_out:integer, fname:string, fname_grid_info:string, \
                                                            data_grid_info:file    )
;
;**************************************************************************;
;                                                                          ;
; initialize_IC                                                            ;
;                                                                          ;
;            cdf             : file   ; NCL Netcdf OUTPUT file             ;
;            field           : string ; Array of field variables to go in  ;
;                                       output file                        ;
;            nfields         : integer; number of fields                   ;
;            nfields_out     : integer; number of fields output to HT      ;
;            fname           : string ; Array of source (INPUT) file names ;
;            fname_grid_info : string ; Name of Grid Input file            ;
;            data_grid_info  : file   ; Grid Info File                     ;
;                                                                          ;
; Initializes a CAM IC file                                                ;
;                                                                          ;
;**************************************************************************;

local ii, jj, count, ifile, name, nfiles,   \
      tmps, tmpn, time_dim, lev_dim, ncol_dim, pi
begin
  pi     = 4.D*atan(1.D)
  nfiles = field@nfiles

  ;------------------------------------------------------------------
  ; Read in Grid Dimensions and Lat/Lon values for each gridpoint
  ;------------------------------------------------------------------
  print("OUTPUT grid based on file "+fname_grid_info)
  file_dim_names = getvardims     (data_grid_info)
  file_dim_sizes = getfiledimsizes(data_grid_info)
  field@ncol     = file_dim_sizes( ind(file_dim_names .eq. "ncol") )
  field@plevo    = file_dim_sizes( ind(file_dim_names .eq. "lev") )
  field@plono    = data_grid_info->lon
  field@plato    = data_grid_info->lat
  field@area     = data_grid_info->area
  delete(file_dim_names)
  delete(file_dim_sizes)

  ;------------------------------------------------------
  ; Get hybrid level information from the "grid_info" file
  ;------------------------------------------------------
  if(field@plevo .gt. 0) then
    field@hyam = (/ data_grid_info->hyam /)
    field@hybm = (/ data_grid_info->hybm /)
    field@hyai = (/ data_grid_info->hyai /)
    field@hybi = (/ data_grid_info->hybi /)
    field@lev  = (/ field@hyam + field@hybm /)*1000.
    field@ilev = (/ field@hyai + field@hybi /)*1000.
  end if

  ;-------------------------
  ; Define global dimensions
  ;-------------------------
  if(field@plevo .gt. 0) then
    if(field@dycore .eq. "se") then
      filedimdef(cdf,(/"ncol"    ,"lev"      ,"ilev"       ,"time","scalar"/),\
                     (/field@ncol,field@plevo,field@plevo+1, 1    , 1      /),\
                     (/ False    , False     , False       , True , False  /) )
    end if
  else
    if(field@dycore .eq. "se") then
      filedimdef(cdf,(/"ncol"    ,"time","scalar"/),\
                     (/field@ncol, 1    , 1      /),\
                     (/False     , True , False  /) )
    end if
  end if

  ;-------------------------
  ; Define global attributes
  ;-------------------------
  cdf@Conventions = "NCAR-CSM"
  cdf@Description = "Re-gridded CAM-SE Initial Conditions from the following source(s):"
  do ifile=0,(nfiles-1)
    ii = ifile+1
    name = "source"+ii
    cdf@$name$  = fname(ifile)
  end do
  cdf@source_grid_info = field@fname_grid_info

  if(field@plevo .gt. 0) then
    cdf@hybrid_sigma_pressure = "p(i,j,k) = A(k)*PO + B(k)*PS(i,j)"
  end if
  cdf@logname     =getenv("LOGNAME")
  cdf@host        =getenv("HOST")
  cdf@hosttype    =getenv("HOSTTYPE")
  cdf@Command     ="NCL script"
  cdf@Date_Written=systemfunc("date")
  cdf@case        =field@case

  ;--------------------------------------------
  ; Define universal parameters and coordinates
  ;--------------------------------------------
  filevardef (cdf,"time","double",(/"time"/))
   cdf->time@long_name    = "time"
   cdf->time@units        = "days since base date (nbdate)"
   cdf->time@calendar     = "noleap"

  filevardef (cdf,"ndbase","integer",(/"scalar"/))
   cdf->ndbase@long_name  = "base day"

  filevardef (cdf,"nsbase","integer",(/"scalar"/))
   cdf->nsbase@long_name  = "seconds of base day"

  filevardef (cdf,"nbdate","integer",(/"scalar"/))
   cdf->nbdate@long_name  = "base date (YYYYMMDD)"

  filevardef (cdf,"nbsec","integer",(/"scalar"/))
   cdf->nbsec@long_name   = "seconds of base date"

  filevardef (cdf,"ndcur","integer",(/"time"/))
   cdf->ndcur@long_name   = "current day (from base day)"

  filevardef (cdf,"nscur","integer",(/"time"/))
   cdf->nscur@long_name   = "current seconds of current day"

  filevardef (cdf,"date","integer",(/"time"/))
   cdf->date@long_name    = "current date (YYYYMMDD)"

  filevardef (cdf,"datesec","integer",(/"time"/))
   cdf->datesec@long_name = "current seconds of current date"

  if(field@plevo .gt. 0) then
    filevardef (cdf,"P0","double",(/"scalar"/))
     cdf->P0@long_name      = "reference pressure"
     cdf->P0@units          = "Pa"

    filevardef (cdf,"lev","double",(/"lev"/))
     cdf->lev@long_name     = "hybrid level at midpoints (1000*(A+B))"
     cdf->lev@units         = "level"
     cdf->lev@positive      = "down"
     cdf->lev@standard_name = "atmosphere_hybrid_sigma_pressure_coordinate" ;
     cdf->lev@formula_terms = "a: hyam b: hybm p0: P0 ps: PS" ;

    filevardef (cdf,"ilev","double",(/"ilev"/))
     cdf->ilev@long_name    = "hybrid level at interfaces (1000*(A+B))"
     cdf->ilev@units        = "level"
     cdf->ilev@positive     = "down"
     cdf->ilev@standard_name = "atmosphere_hybrid_sigma_pressure_coordinate" ;
     cdf->ilev@formula_terms = "a: hyai b: hybi p0: P0 ps: PS" ;

    filevardef (cdf,"hyai","double",(/"ilev"/))
     cdf->hyai@long_name    = "hybrid A coefficient at layer interfaces"

    filevardef (cdf,"hybi","double",(/"ilev"/))
     cdf->hybi@long_name    = "hybrid B coefficient at layer interfaces"

    filevardef (cdf,"hyam","double",(/"lev"/))
     cdf->hyam@long_name    = "hybrid A coefficient at layer midpoints"

    filevardef (cdf,"hybm","double",(/"lev"/))
     cdf->hybm@long_name    = "hybrid B coefficient at layer midpoints"
  end if

  ;-----------------------
  ; Pre-set some variables
  ;-----------------------
  cdf->time    = 0.
  cdf->ndbase  = 0
  cdf->nsbase  = 0
  cdf->nbsec   = 0
  cdf->ndcur   = 0
  cdf->nscur   = 0
  cdf->datesec = 0

  if(field@plevo .gt. 0) then
    cdf->P0      = 100000.D
    cdf->lev     = -999.
    cdf->ilev    = -999.
    cdf->hyai    = -999.
    cdf->hybi    = -999.
    cdf->hyam    = -999.
    cdf->hybm    = -999.
  end if

;--------------------------------------------------
; Define dycore-specific parameters and coordinates
;--------------------------------------------------

  ;---------
  ; CAM-SE
  ;---------
  if(field@dycore .eq. "se") then

    filevardef (cdf,"lat","double",(/"ncol"/))
     cdf->lat@long_name     = "latitude"
     cdf->lat@units         = "degrees_north"

    filevardef (cdf,"lon","double",(/"ncol"/))
     cdf->lon@long_name     = "longitude"
     cdf->lon@units         = "degrees_east"

    filevardef (cdf,"area","double",(/"ncol"/))
     cdf->area@long_name     = "gridpoint area"
     cdf->area@units         = "normalized"

    filevardef (cdf,"ntrm","integer",(/"scalar"/))
     cdf->ntrm@long_name    = "spectral truncation parameter M"

    filevardef (cdf,"ntrn","integer",(/"scalar"/))
     cdf->ntrn@long_name    = "spectral truncation parameter N"

    filevardef (cdf,"ntrk","integer",(/"scalar"/))
     cdf->ntrk@long_name    = "spectral truncation parameter K"

    ;-----------------------
    ; Pre-set some variables
    ;-----------------------

    cdf->lon =field@plono
    cdf->lat =field@plato
    cdf->area=field@area

    cdf->ntrm=field@ptrmo
    cdf->ntrn=field@ptrno
    cdf->ntrk=field@ptrko

    field@lat=cdf->lat
    field@lon=cdf->lon

    ;---------------------------------------------------
    ; Define field variables, dimensions, and attributes
    ;---------------------------------------------------
    time_dim = "time"
    lev_dim  = "lev"
    count    = 0
    do ii=0,(nfields-1)
      if(field@loutput(ii)) then
        count = count + 1

        print ("   Declare space for output field:  "+field(ii)+ \
               " (field "+count+" of "+nfields_out+")")

        ncol_dim  = "ncol"
        if(field@dimensions(ii) .eq. "2D") then
          filevardef (cdf,field(ii),field@precision,(/time_dim,ncol_dim/))
        end if
        if(field@dimensions(ii) .eq. "3D") then
          filevardef (cdf,field(ii),field@precision,(/time_dim,lev_dim,ncol_dim /))
        end if

        cdf->$field(ii)$@long_name           =       field@long_name   (ii)
        cdf->$field(ii)$@units               =       field@units       (ii)
        cdf->$field(ii)$@source_file         = fname(field@source_file (ii))
        cdf->$field(ii)$@source_field_name   =       field@source_field(ii)
        cdf->$field(ii)$@horz_interp_flag    =       field@horz_interp (ii)
        if(field@dimensions(ii) .eq. "3D") then
          cdf->$field(ii)$@vert_interp_flag  =       field@vert_interp (ii)
        end if
      end if
    end do ; ii=0,(nfields-1)
  end if

  ;-------------------------------------------
  ; Fill output dataset with level information
  ;-------------------------------------------
  if(field@plevo .gt. 0) then
    cdf->hyam@source_file = fname_grid_info
    cdf->hybm@source_file = fname_grid_info
    cdf->hyai@source_file = fname_grid_info
    cdf->hybi@source_file = fname_grid_info
    cdf->lev@source_file  = fname_grid_info
    cdf->ilev@source_file = fname_grid_info
    cdf->hyam  = (/ field@hyam /)
    cdf->hybm  = (/ field@hybm /)
    cdf->lev   = (/ field@lev  /)
    cdf->hyai  = (/ field@hyai /)
    cdf->hybi  = (/ field@hybi /)
    cdf->ilev  = (/ field@ilev /)
    delete(cdf->hyam@_FillValue)
    delete(cdf->hybm@_FillValue)
    delete(cdf->hyai@_FillValue)
    delete(cdf->hybi@_FillValue)
    delete(cdf->lev@_FillValue )
    delete(cdf->ilev@_FillValue)
  end if

  return
end
;=====================================================================================


; ;=====================================================================================
; procedure initialize_IC_levels(cdf:file, field:string,fname_grid_info:string, \
;                                                       data_grid_info :file    )
; ;
; ;************************************************************************;
; ;                                                                        ;
; ; initialize_IC_levels                                                   ;
; ;                                                                        ;
; ;            cdf     : file   ; NCL Netcdf OUTPUT file                   ;
; ;            field   : string ; Array of field variables to go in file   ;
; ;            fname_grid_info: string ; INPUT file name with output level  ;
; ;                                     info                               ;
; ;            data_grid_info : file   ; Input file                         ;
; ;                                                                        ;
; ; Fill output dataset with level information:                            ;
; ;                                                                        ;
; ;   - Fill with level info defined in file_grid_info                      ;
; ;                                                                        ;
; ;************************************************************************;
; 
; begin
;   print("")
;   if(field@plevo .gt. 0) then
;     ;------------------------------------------------------
;     ; Get hybrid level information from the "grid_info" file
;     ;------------------------------------------------------
;     print("OUTPUT levels based on file "+fname_grid_info)
; 
;     field@hyam = (/ data_grid_info->hyam /)
;     field@hybm = (/ data_grid_info->hybm /)
;     field@hyai = (/ data_grid_info->hyai /)
;     field@hybi = (/ data_grid_info->hybi /)
; 
;     cdf->hyam@source_file = fname_grid_info
;     cdf->hybm@source_file = fname_grid_info
;     cdf->hyai@source_file = fname_grid_info
;     cdf->hybi@source_file = fname_grid_info
;     cdf->lev@source_file  = fname_grid_info
;     cdf->ilev@source_file = fname_grid_info
; 
;     field@lev  = (/ field@hyam + field@hybm /)*1000.
;     field@ilev = (/ field@hyai + field@hybi /)*1000.
; 
;     cdf->hyam  = (/ field@hyam /)
;     cdf->hybm  = (/ field@hybm /)
;     cdf->lev   = (/ field@lev  /)
;     cdf->hyai  = (/ field@hyai /)
;     cdf->hybi  = (/ field@hybi /)
;     cdf->ilev  = (/ field@ilev /)
;     delete(cdf->hyam@_FillValue)
;     delete(cdf->hybm@_FillValue)
;     delete(cdf->hyai@_FillValue)
;     delete(cdf->hybi@_FillValue)
;     delete(cdf->lev@_FillValue )
;     delete(cdf->ilev@_FillValue)
;   end if
; 
; end
; ;=====================================================================================


;=====================================================================================
procedure initialize_IC_dates(cdf:file     , nfiles:integer, mytmpdir:string, lfname:string  , \
                              ftype:string , fdate:double   ,                                  \
                              ref_date:integer              , data0:file     ,                 \
                              data1:file   , data2:file     , data3:file     ,                 \
                              data4:file   , data5:file     , data6:file     , data7:file      )
;
;************************************************************************;
;                                                                        ;
; initialize_IC_dates                                                    ;
;                                                                        ;
;            cdf     : file   ; NCL Netcdf OUTPUT file                   ;
;            nfiles  : integer; number of input files                    ;
;            mytmpdir: string ; local disk pathname                      ;
;            lfname  : string ; Array of local input file names          ;
;            ftype   : string ; Array of input file types                ;
;            fdate   : double ; Array of dates and/or time slice pointers;
;            ref_date: integer; Base year/month/day                      ;
;            data0-data7 :file; Set of input files                       ;
;                                                                        ;
;   - fdate is an input/output array of dates for each input file.       ;
;   - On input, fdate may be filled with a combination of either dates   ;
;     or time slice pointers                                             ;
;        a:  if filled with dates, replace each date with the            ;
;            appropriate time slice pointer for each file                ;
;        b:  The date of the FIRST input file will be the date stamped   ;
;            on the OUTPUT file                                          ;
;                                                                        ;
;************************************************************************;

local ii, ifile, data, datei, houri, date, hour, dims, tmp, date_keep,   \
      hour_keep, found, hourf, time_keep, time_ref, ddate, date_int,     \
      ichar_ind1, ichar_ind2, ichars, year_str, month_str, day_str,      \
      hour_str, analysis_info

begin
  print("")
  do ifile=0,(nfiles-1)

    print("Checking date for "+mytmpdir+"/"+lfname(ifile) )
    data = choose_file(data0, data1, data2, data3, data4, data5, data6, data7, ifile)

    if(fdate(ifile) .lt. 0) then
      print("Error:  'fdate' must be a positive-definite integer")
      exit
    end if

    analysis_info = extract_analyses_info(data, ftype(ifile))
    datei = analysis_info@datei
    houri = analysis_info@houri
    hourf = analysis_info@hourf
    dims  = dimsizes(datei)

    found = False
    if(dims .gt. 1) then
      do ii =0,(dims-1)
        date = (/ datei(ii) /)
        hour = (/ houri(ii) /)

        if(fdate(ifile) .lt. 10000000) then
          if(fdate(ifile) .gt. dims(0)-1) then
            print("ERROR:  'fdate' pointer is out of bounds of date array")
            exit
          end if
          if(fdate(ifile) .eq. ii .and. houri(ii) .eq. hourf(ii)) then
            date_keep = date
            hour_keep = hour
            found     = True
          end if
        else
          tmp = fdate(ifile)
          tmp = date
          tmp = tmp*100000 + hour
          if(fdate(ifile) .eq. tmp .and. houri(ii) .eq. hourf(ii)) then
            date_keep    = date
            hour_keep    = hour
            fdate(ifile) = ii
            found        = True
          end if
          delete(tmp)
        end if
      end do ; ii =0,(dims-1)

      if( .not. found) then
        print("ERROR:  desired date not found on file")
        exit
      end if
    else
      date = (/ datei /)
      hour = (/ houri /)

      if(fdate(ifile) .lt. 10000000) then
        if(fdate(ifile) .gt. 0) then
          print("ERROR:  'fdate' pointer is out of bounds of date array")
          exit
        end if
        if(houri .eq. hourf) then
          date_keep = date
          hour_keep = hour
          found     = True
        end if
      else
        tmp = fdate(ifile)
        tmp = date
        tmp = tmp*100000 + hour
        if(fdate(ifile) .eq. tmp .and. houri .eq. hourf) then
          date_keep    = date
          hour_keep    = hour
          fdate(ifile) = 0
          found        = True
        end if
        delete(tmp)
      end if

      if( .not. found) then
        print("ERROR:  desired date not found on file")
        exit
      end if
    end if

    delete(datei)
    delete(houri)
    delete(hourf)
    delete(dims)

    ;------------------------------------------------
    ; Compute time in days relative to reference date
    ;------------------------------------------------
    print("date = "+date_keep+"; "+hour_keep+"Z; time slice pointer = "+fdate(ifile))
    print("")
    if(ifile .eq. 0) then

      if(ref_date .gt. 31000000) then
        print("Error:  Reference date must be in form of YYYYMMDD")
        print("Reference date = "+ref_date)
        exit
      end if

      ddate     = new( 1, double)
      ddate     = date_keep
      ddate     = ddate*100000
      ddate     = ddate + hour_keep
      time_keep = date_to_time(ddate)
      ddate     = ref_date
      ddate     = ddate*100000
      time_ref  = date_to_time(ddate)

      if(time_ref .gt. time_keep) then
        print("ERROR:  reference time gt current time")
        print("reference date = "+ref_date )
        print("current   date = "+date_keep+"; current hour = "+hour_keep)
        exit
      end if

      cdf->ndbase     = 0
      cdf->nsbase     = 0
      cdf->nbdate     = ref_date
      cdf->nbsec      = 0
      cdf->ndcur      = doubletointeger(time_keep - time_ref)
      cdf->nscur      = hour_keep
      cdf->date       = date_keep
      cdf->datesec    = hour_keep
      cdf->time       = time_keep - time_ref

      ;-----------------------------------------------
      ; If making SST dataset, modify "time" meta-data
      ;-----------------------------------------------
      tmp       = fdate(0)
      tmp       = ref_date*100000.D
      year_str  = new( 1 , string )
      month_str = new( 1 , string )
      day_str   = new( 1 , string )
      hour_str  = new( 1 , string )
      date2YMDH_str(tmp, year_str, month_str, day_str, hour_str)
      if (isfilevar(cdf, "SST_cpl")) then
        cdf->time@units    = "days since "+year_str+"-"+month_str+"-"+day_str+" 00:00:00"
        cdf->time@calendar = "365_day"
      end if
      delete(tmp)

      print("")
      print("----------------------------------------------------------------------- ")
      print("OUTPUT file time info:")
      print("")
      print("    time    = "+cdf->time)
      print("    ndbase  = "+cdf->ndbase)
      print("    nsbase  = "+cdf->nsbase)
      print("    nbdate  = "+cdf->nbdate)
      print("    nbsec   = "+cdf->nbsec)
      print("    ndcur   = "+cdf->ndcur)
      print("    nscur   = "+cdf->nscur)
      print("    date    = "+cdf->date)
      print("    datesec = "+cdf->datesec)
      print("")
      print("----------------------------------------------------------------------- ")
      print("")
    end if

  end do ; ifile=0,(nfiles-1)

  return
end
;=====================================================================================


;=====================================================================================
procedure find_field_index(field:string, field_of_interest:string, nfields:integer, \
                           field_index:integer, file_index:integer )
;
;************************************************************************;
;                                                                        ;
; find_file_index:  find field and file index for "field_of_interest"    ;
;                                                                        ;
;   field            : string ; Array of field variables                 ;
;   field_of_interest: string ; field name searching for                 ;
;   nfields          : integer; total number of fields                   ;
;   field_index      : integer; index in "field" that contains           ;
;                               "field_of_interest"                      ;
;   file_index       : integer; file index in list of files that         ;
;                               contains "field_of_interest"             ;
;                                                                        ;
;************************************************************************;

local  ii, file_index, field_index

begin
  file_index  = -1
  field_index = -1
  do ii=0,(nfields-1)
    if(field(ii) .eq. field_of_interest) then
      file_index   = field@source_file(ii)
      field_index  = ii
    end if
  end do

  if(field_index .lt. 0) then
    print("Error:  Looking for '"+field_of_interest+"' but did not find in 'field' list")
    exit
  end if
  return
end
;=====================================================================================


;=====================================================================================
function flip_longitudes(xx:numeric)
;
;*****************************************************************************;
;                                                                             ;
; flip_longitudes:  flip longitudes of array so that Greenwich is at index 0  ;
;                   (assumes index 0 of incoming array is at dateline (-180)) ;
;                                                                             ;
;                                                                             ;
;*****************************************************************************;

local  kk, dims, dimdims, yy, dimx, dim2, dateline, error

begin
  dateline = -180.
  print("     flip longitudes of array so that Greenwich is at index 0")

  ;-------------------------------------------
  ; Check that incoming data start at dateline
  ;-------------------------------------------
  error    = abs(dateline - xx&lon(0))
  if(error .gt. 1.e-4) then
    print("Error  function 'flip_longitudes':  index 0 of longitude dimension of ")
    print("                               incoming data array is not at dateline.")
    print("index 0 begins at: "+xx&lon(0))
    exit
  end if

  dims    = dimsizes(xx)
  dimdims = dimsizes(dims)
  dimx    = dims(dimdims-1)
  dim2    = dimx/2

  if(dimdims .gt. 3) then
    print("Error  function 'flip_longitudes':  not expecting array of greater ")
    print("                                       than 3 dimensions for input.")
    print("Dimension of input array = "+dimdims)
    exit
  end if

  yy = xx
  if(dimdims .eq. 2) then
    yy(  :,   0:dim2-1) = xx(  :,dim2:dimx-1)
    yy(  :,dim2:dimx-1) = xx(  :,   0:dim2-1)
  end if
  if(dimdims .eq. 3) then
    yy(:,:,   0:dim2-1) = xx(:,:,dim2:dimx-1)
    yy(:,:,dim2:dimx-1) = xx(:,:,   0:dim2-1)
  end if
  do kk=0,(dimx-1)
    if(yy&lon(kk) .lt. 0.) then
      yy&lon(kk) = yy&lon(kk)-dateline*2.
    end if
  end do
  delete(dims)

  return(yy)
end
;=====================================================================================


;=====================================================================================
function flip_latitudes(xx:numeric)
;
;************************************************************************;
;                                                                        ;
; flip_latitudes:  flip latitudes of arrays to run S --> N               ;
;                                                                        ;
;************************************************************************;

local  kk, jj, dimdims, yy, dimy

begin
  print("     flip latitudes to run S --> N")

  ;--------------------------------------------
  ; Check that incoming data are indeed N --> S
  ;--------------------------------------------
  if(xx&lat(dimsizes(xx&lat)-1) .gt. xx&lat(0)) then
    print("Error  function 'flip_latitudes':  data already oriented S --> N.")
    exit
  end if

  dimdims  = dimsizes(dimsizes(xx))
  dimy     = dimsizes(xx&lat)

  if(dimdims .gt. 3) then
    print("Error  function 'flip_latitudes':  not expecting array of greater ")
    print("                                      than 3 dimensions for input.")
    print("Dimension of input array = "+dimdims)
    exit
  end if

  yy = xx
  do kk = 0,dimy-1
    jj = dimy - 1 - kk
    if(dimdims .eq. 2) then
      yy(kk, :) = xx(jj, :)
    end if
    if(dimdims .eq. 3) then
      yy(: ,kk, :) = xx(: ,jj, :)
    end if
  end do

  return(yy)
end
;=====================================================================================


;=====================================================================================
function flip_levels(xx:numeric)
;
;************************************************************************;
;                                                                        ;
; flip_levels:  flip levels of arrays to run Top --> Bot                 ;
;                                                                        ;
;************************************************************************;

local  jj, kk, yy, dimz

begin
  print("     flip levels to run Top --> Bot")

  ;------------------------------------------------
  ; Check that incoming data are indeed Bot --> Top
  ;------------------------------------------------
  if(xx&lev(dimsizes(xx&lev)-1) .gt. xx&lev(0)) then
    print("Error  function 'flip_levels':  data already oriented Top --> Bot")
    exit
  end if

  dimz     = dimsizes(xx&lev)

  yy = xx
  do jj=0,(dimz-1)
    kk = dimz - 1 - jj
    yy(jj, :, :) = xx(kk, :, :)
  end do

  return(yy)
end
;=====================================================================================


;=====================================================================================
function read_reformat(data:file, field:string, ifield:integer)
;
;************************************************************************;
;                                                                        ;
; read_reformat                                                          ;
;                                                                        ;
;            data    : file   ; INPUT file                               ;
;            field   : string ; Master list of field names and attributes;
;            ifield  : integer; field index                              ;
;                                                                        ;
; Read data and reformat to a CAM format                                 ;
;                                                                        ;
;************************************************************************;

local xx, xfield, dims, dimdims, processed, tmp, lat_tmp, lon_tmp, pmid,         \
      ps_in, geos_airmw, geos_h2omw, geos_vireps, sphu, rocp, tv, te, kk,ii,jj,  \
      dnames, source_file_index, source_field, fdate, time_dimname,              \
      real_imaginary_dimname, lev_dimname, lat_dimname, lon_dimname, dim_names,  \
      lat_dimname_alt1, lon_dimname_alt1, analysis_info, lev_dimname_alt1,       \
      dim_offset, dim_index, pi, coslat, ftype, tmpu, tmpv, tmpu_extents,        \
      tmpv_extents, extent_dim, lon_deg, lat_deg, fgrid

begin
  source_file_index      = field@source_file (ifield)
  source_field           = field@source_field(ifield)
  fdate                  = field@fdate       (source_file_index)
  ftype                  = field@ftype(source_file_index)

  analysis_info          = extract_analyses_info(data, ftype)

  time_dimname           = analysis_info@time_dimname
  real_imaginary_dimname = analysis_info@real_imaginary_dimname
  lev_dimname            = analysis_info@lev_dimname
  lev_dimname_alt1       = analysis_info@lev_dimname_alt1
  lat_dimname            = analysis_info@lat_dimname
  lat_dimname_alt1       = analysis_info@lat_dimname_alt1
  lon_dimname            = analysis_info@lon_dimname
  lon_dimname_alt1       = analysis_info@lon_dimname_alt1

  ;--------------------------
  ; Read field from data file
  ;--------------------------
  processed = False
  print ("   Reading input field:  "+source_field )
  xfield    = data->$source_field$

  ;-----------------------------------------------------------
  ; Build dimension set in preparation to reformat input field
  ;-----------------------------------------------------------
  dim_names = new( 5, string )
  dnames    = getvardims(xfield)

  kk = -1

  if( any(dnames .eq. time_dimname          ) ) then
    kk = kk + 1
    dim_names(kk) = time_dimname
  end if
  if( any(dnames .eq. real_imaginary_dimname) ) then
    kk = kk + 1
    dim_names(kk) = real_imaginary_dimname
  end if
  if( any(dnames .eq. lev_dimname           ) ) then
    kk = kk + 1
    dim_names(kk) = lev_dimname
  end if
  if( any(dnames .eq. lev_dimname_alt1      ) ) then
    kk = kk + 1
    dim_names(kk) = lev_dimname_alt1
  end if
  if( any(dnames .eq. lat_dimname           ) ) then
    kk = kk + 1
    dim_names(kk) = lat_dimname
  end if
  if( any(dnames .eq. lat_dimname_alt1      ) ) then
    kk = kk + 1
    dim_names(kk) = lat_dimname_alt1
  end if
  if( any(dnames .eq. lon_dimname           ) ) then
    kk = kk + 1
    dim_names(kk) = lon_dimname
  end if
  if( any(dnames .eq. lon_dimname_alt1      ) ) then
    kk = kk + 1
    dim_names(kk) = lon_dimname_alt1
  end if

  if(dimsizes(dnames) .ne. kk+1)
    print("Error:  did not find all of the dimension names in 'dnames'")
    print("Dimensions in current field = "+dnames)
    exit
  end if
  delete(dnames)

  ;------------------------------------------------------
  ; reorder input field into standard CAM dimension order
  ;------------------------------------------------------
  if( dim_names(0) .eq. time_dimname) then
    if (kk .eq. 2) then
      xx = xfield ($dim_names(0)$ | fdate, $dim_names(1)$ | : , $dim_names(2)$ | : )
      processed = True
    end if
    if (kk .eq. 3) then
      xx = xfield ($dim_names(0)$ | fdate, $dim_names(1)$ | : , \
                   $dim_names(2)$ | :    , $dim_names(3)$ | : )
      processed = True
    end if
    if (kk .eq. 4) then
      xx = xfield ($dim_names(0)$ | fdate, $dim_names(1)$ | : , \
                   $dim_names(2)$ | :    , $dim_names(3)$ | : , $dim_names(4)$ | : )
      processed = True
    end if
  else
    if (kk .eq. 1) then
      xx = xfield ($dim_names(0)$ | : , $dim_names(1)$ | : )
      processed = True
    end if
    if (kk .eq. 2) then
      xx = xfield ($dim_names(0)$ | : , $dim_names(1)$ | : , $dim_names(2)$ | : )
      processed = True
    end if
    if (kk .eq. 3) then
      xx = xfield ($dim_names(0)$ | : , $dim_names(1)$ | : , \
                   $dim_names(2)$ | : , $dim_names(3)$ | : )
      processed = True
    end if
  end if

  if(.not. processed) then
    print("Error:  field '"+source_field+"' not processed.")
    exit
  end if

  dnames = getvardims(xx)

  ;-------------------------------------------------------------------
  ; If field is expressed in spherical harmonics, convert back to grid
  ;-------------------------------------------------------------------
  if( any(dnames .eq. real_imaginary_dimname) ) then
    print ("     Converting '"+source_field+"' from spectral space to gridpoint space")
    fgrid    = shsgC(xx, dimsizes(analysis_info@lon))
    delete(xx)
    xx       = fgrid
    delete(fgrid)
    dim_offset = 1
  else
    dim_offset = 0
  end if

  ;----------------------------------------
  ; Convert coordinate arrays to CAM format
  ;----------------------------------------
  do ii = dim_offset,dimsizes(dnames)-1
    dim_index = ii - dim_offset
    if( dnames(ii) .eq. lev_dimname     ) then
      xx!dim_index = "lev"
      if( iscoord(xx,"lev") ) then
        delete(xx&lev)
      end if
      xx&lev       = new( dimsizes(analysis_info@lev     ), double )
      xx&lev       = analysis_info@lev
    end if
    if( dnames(ii) .eq. lev_dimname_alt1) then
      xx!dim_index = "lev"
      if( iscoord(xx,"lev") ) then
        delete(xx&lev)
      end if
      xx&lev       = new( dimsizes(analysis_info@lev_alt1), double )
      xx&lev       = analysis_info@lev_alt1
    end if
    if( dnames(ii) .eq. lat_dimname     ) then
      xx!dim_index = "lat"
      if( iscoord(xx,"lat") ) then
        delete(xx&lat)
      end if
      xx&lat       = new( dimsizes(analysis_info@lat     ), double )
      xx&lat       = analysis_info@lat
    end if
    if( dnames(ii) .eq. lat_dimname_alt1) then
      xx!dim_index = "lat"
      if( iscoord(xx,"lat") ) then
        delete(xx&lat)
      end if
      xx&lat       = new( dimsizes(analysis_info@lat_alt1), double )
      xx&lat       = analysis_info@lat_alt1
    end if
    if( dnames(ii) .eq. lon_dimname     ) then
      xx!dim_index = "lon"
      if( iscoord(xx,"lon") ) then
        delete(xx&lon)
      end if
      xx&lon       = new( dimsizes(analysis_info@lon     ), double )
      xx&lon       = analysis_info@lon
    end if
    if( dnames(ii) .eq. lon_dimname_alt1) then
      xx!dim_index = "lon"
      if( iscoord(xx,"lon") ) then
        delete(xx&lon)
      end if
      xx&lon       = new( dimsizes(analysis_info@lon_alt1), double )
      xx&lon       = analysis_info@lon_alt1
    end if
  end do

  dims    = dimsizes(xx)
  dimdims = dimsizes(dims)

  ;------------------------
  ; Remove cosine weighting
  ;------------------------
  if(analysis_info@remove_cos_wgt .and.                       \
     any(analysis_info@cos_wgted_var_names .eq. source_field) ) then
    print ("     Removing cosine weighting from "+source_field)
    pi     = 4.D*atan(1.D)
    coslat = cos( xx&lat*pi/180.D )
    do jj=0,(dimsizes(coslat)-1)
      if(isfloat(xx)) then
        xx(:,jj,:) = doubletofloat( xx(:,jj,:)/coslat(jj) )
      else
        xx(:,jj,:) =                xx(:,jj,:)/coslat(jj)
      end if
    end do
  end if

  ;-----------------------------------------------
  ; Flip longitudes so that dateline is at index 0
  ;-----------------------------------------------
  if(xx&lon(0) .lt. -179.) then
    tmp = xx
    tmp = flip_longitudes(xx)
    delete(xx)
    xx  = tmp
    delete(tmp)
  end if

  ;--------------------------------------
  ; Flip latitudes to orient data S --> N
  ;--------------------------------------
  if(xx&lat(0) .ge. xx&lat(dimsizes(xx&lat)-1)) then
    tmp = xx
    tmp = flip_latitudes(xx)
    delete(xx)
    xx  = tmp
    delete(tmp)
  end if

  ;---------------------------------------
  ; Flip levels to orient data Top --> Bot
  ;---------------------------------------
  if(dimdims .eq. 3) then
    if(xx&lev(0) .ge. xx&lev(dimsizes(xx&lev)-1)) then
      tmp = xx
      tmp = flip_levels(xx)
      delete(xx)
      xx  = tmp
      delete(tmp)
    end if
  end if

  ;----------------------------------
  ; If in log(Ps) form, convert to Ps
  ;----------------------------------
  if( analysis_info@logPs_to_Ps .and.                      \
      any(analysis_info@logPs_var_names .eq. source_field) ) then
    print ("     Converting "+source_field+" from log")
    xx = xx
    xx = exp(xx)
  end if

  ;------------------------------------------------------------
  ; If surface topography is in meters, convert to geopotential
  ;------------------------------------------------------------
  if( analysis_info@Z_to_PHIS .and.                    \
      any(analysis_info@Z_var_names .eq. source_field) ) then
    print ("     Converting "+source_field+" to geopotential")
    xx = xx*9.80616
  end if

  ;-------------------
  ; Multiply by factor
  ;-------------------
  if( analysis_info@factor .ne. 0. .and.                    \
      any(analysis_info@factor_var_names .eq. source_field) ) then
    print ("     Multiplying "+source_field+" by "+analysis_info@factor)
    xx = xx*analysis_info@factor
  end if

  ;-----------------------
  ; If theta, convert to T
  ;-----------------------
  if( analysis_info@convert_theta_to_T .and.               \
      any(analysis_info@theta_var_names .eq. source_field) ) then
    print ("     Converting 'theta' to 'T'")

    geos_airmw  = 28.97D
    geos_h2omw  = 18.01D
    geos_vireps = geos_airmw/geos_h2omw-1.D
    rocp        = 2.D/7.D

    tmp         = data->ps
    ps_in       = tmp(time|fdate, {lat| -90.:90.} ,{lon| 0.:360.} )
    pmid        = build_press(data, ftype, ps_in, "FULL")
    delete(tmp)
    tmp         = data->sphu
    sphu        = tmp   (time|fdate, {lev|0.:1000.} ,{lat| -90.:90.} ,{lon| 0.:360.} )
    delete(tmp)
    tmp         = xx
    tv          = tmp * ( pmid^rocp )
    te          = tv  / ( 1.D + geos_vireps * sphu )

    xx          = doubletofloat(te)

    delete(te)
    delete(tv)
    delete(tmp)
    delete(ps_in)
    delete(pmid)
    delete(sphu)
  end if

  ;-------------------------------------------------------------------------------
  ; Fix bug in NASA_MERRA data where u-field is corrupted at S. pole (replace u/v at
  ; S. pole with wave 1 u/v)
  ;-------------------------------------------------------------------------------
  if( analysis_info@PREVOCA_u_v_bugfix .and.                     \
      any(analysis_info@PREVOCA_u_v_var_names .eq. source_field) ) then
    tmpu                  = new( (/dims(0), dims(1)-2 , dims(2) /), double)
    tmpv                  = tmpu
    tmpu(:,0:dims(1)-3,:) = data->uwnd(time|fdate, lev|:, lat|1:dims(1)-2, lon|:)
    tmpv(:,0:dims(1)-3,:) = data->vwnd(time|fdate, lev|:, lat|1:dims(1)-2, lon|:)

    extent_dim            = 6
    tmpu_extents          = new( (/dims(0), extent_dim, dims(2) /), double)
    tmpv_extents          = tmpu_extents

    lon_deg = new( (/ dims(2)   /), double)
    lat_deg = new( (/ dims(1)-2 /), double)
    lon_deg = tmpu&lon
    lat_deg = tmpu&lat
    MAKEIC::extyv(dims(2)  ,dims(0)  ,dims(1)-2 ,extent_dim  , lon_deg,    \
                  lat_deg  ,tmpu     ,tmpv      ,tmpu_extents, tmpv_extents)

    if(source_field .eq. "uwnd") then
      xx(:,0,:) = doubletofloat((/ tmpu_extents(:,extent_dim/2-1,:) /))
    end if
    if(source_field .eq. "vwnd") then
      xx(:,0,:) = doubletofloat((/ tmpv_extents(:,extent_dim/2-1,:) /))
    end if

    delete(tmpu_extents)
    delete(tmpv_extents)
    delete(tmpu)
    delete(tmpv)
  end if

  ;----------------------------------------------
  ; Convert coordinate arrays to double precision
  ;----------------------------------------------
  if(isfloat(xx&lat)) then
    lat_tmp = new( (/dimsizes(xx&lat) /) , double)
    lat_tmp = xx&lat
    delete(xx&lat)
    xx&lat  = new( (/dimsizes(lat_tmp)/) , double)
    xx&lat  = lat_tmp
    delete(lat_tmp)
  end if

  if(isfloat(xx&lon)) then
    lon_tmp = new( (/dimsizes(xx&lon) /) , double)
    lon_tmp = xx&lon
    delete(xx&lon)
    xx&lon  = new( (/dimsizes(lon_tmp)/) , double)
    xx&lon  = lon_tmp
    delete(lon_tmp)
  end if

  ;-------------
  ; Check format
  ;-------------
  if(dimdims .eq. 2) then
    if(xx!0 .ne. "lat") then
      print ("   Error:  1st dimension of 2-D field is not named 'lat'")
      exit
    end if
    if(xx!1 .ne. "lon") then
      print ("   Error:  2nd dimension of 2-D field is not named 'lon'")
      exit
    end if
    if(xx&lat(0) .ge. xx&lat(dims(0)-1)) then
      print ("   Error:  latitude dimension oriented backwards")
      print ("           Should be S --> N")
      print(xx&lat)
      exit
    end if
    if(xx&lon(0) .ge. xx&lon(dims(1)-1)) then
      print ("   Error:  longitude dimension oriented incorrectly")
      print ("           Should be eastward beginning at Greenwich meridian")
      print(xx&lon)
      exit
    end if
    if( any(xx&lon(0) .ge. xx&lon(1:dims(1)-1)) ) then
      print ("   Error:  longitude array does not start at Greenwich meridian")
      print(xx&lon)
      exit
    end if
  else
    if(xx!0 .ne. "lev") then
      print ("   Error:  1st dimension of 3-D field is not named 'lev'")
      exit
    end if
    if(xx!1 .ne. "lat") then
      print ("   Error:  2nd dimension of 3-D field is not named 'lat'")
      exit
    end if
    if(xx!2 .ne. "lon") then
      print ("   Error:  3rd dimension of 3-D field is not named 'lon'")
      exit
    end if
    if(xx&lev(0) .ge. xx&lev(dims(0)-1)) then
      print ("   Error:  level dimension oriented upside down")
      print ("           Should be Top --> Bot")
      print(xx&lev)
      exit
    end if
    if(xx&lat(0) .ge. xx&lat(dims(1)-1)) then
      print ("   Error:  latitude dimension oriented backwards")
      print ("           Should be S --> N")
      print(xx&lat)
      exit
    end if
    if(xx&lon(0) .ge. xx&lon(dims(2)-1)) then
      print ("   Error:  longitude dimension oriented incorrectly")
      print ("           Should be eastward beginning at Greenwich meridian")
      print(xx&lon)
      exit
    end if
    if( any(xx&lon(0) .ge. xx&lon(1:dims(2)-1)) ) then
      print ("   Error:  longitude array does not start at Greenwich meridian")
      print(xx&lon)
      exit
    end if
  end if

  delete(dims)
  delete(xfield)

  field@processed(ifield) = True

  return(xx)
end
;=====================================================================================


;=====================================================================================
function horz_process(xfield:numeric, cdf:file, field:string, ifield:numeric,  \
                       plono:numeric, plato:numeric, flon:numeric, flat:numeric)
;
;***********************************************************************;
;                                                                       ;
; horz_process                                                          ;
;                                                                       ;
;            xfield  : numeric; INPUT field                             ;
;            cdf     : file   ; NCL Netcdf OUTPUT file                  ;
;            field   : string ; Names of of OUTPUT field variables      ;
;            ifield  : integer; OUTPUT field index                      ;
;            plono   : numeric; OUTPUT longitude dimension              ;
;            plato   : numeric; OUTPUT latitude  dimension              ;
;            flon    : numeric; OUTPUT longitudes (deg)                 ;
;            flat    : numeric; OUTPUT latitudes  (deg)                 ;
;                                                                       ;
; Horizontally interpolate/copy input field to output field             ;
;                                                                       ;
;***********************************************************************;

local kk, dims, xx, x1, yy, l_interph, horz_interp, dim_lat,       \
      dim_lev, dim_lon, sum1, sum2, extent_dim, x1_extents, limdr, \
      opt, tmp, bin_factor

begin
  l_interph = False
  horz_interp = field@horz_interp(ifield)
  xx          = xfield

  ;----------------------------------------------------------
  ; Copy "from_hires" attribute (if it exists) to output file
  ;----------------------------------------------------------
  if(isatt(xx,"from_hires")) then
    cdf->$field(ifield)$@from_hires = xx@from_hires
  end if

  ;------------------------------------
  ; Special handling for certain fields
  ;------------------------------------
  if(field(ifield) .eq. "U"  .or. field(ifield) .eq. "V"  .or. \
     field(ifield) .eq. "US" .or. field(ifield) .eq. "VS")  then
    print("     Error:  Procedure 'horz_process' does not handle the pole point")
    print("             properly for vector quantities.  Call 'horz_process_v'")
    exit
  end if

  dims = dimsizes(dimsizes(xx))

  ;-----------------------------------
  ; Determine dimensions of input file
  ;-----------------------------------
  dim_lat = dimsizes(xx&lat)
  dim_lon = dimsizes(xx&lon)
  if(dims .eq. 2) then
    dim_lev = 1
  else
    dim_lev = dimsizes(xx&lev)
  end if

  ;-----------------------
  ; Initialize work arrays
  ;-----------------------
  x1        = new( (/dim_lev, dim_lat, dim_lon /), double)
  x1!0      = "lev"
  x1!1      = "lat"
  x1!2      = "lon"
  yy        = new( (/dim_lev, plato  , plono   /), double)
  yy!0      = "lev"
  yy!1      = "lat"
  yy!2      = "lon"

  if(dims .eq. 2) then
    x1(lev|0,lat|:,lon|:) = xx(lat|:,lon|:)
    x1&lev                = 0.
  else
    x1                    = xx
  end if

;--------------------
; Begin interpolation
;--------------------

  ;----------------------------
  ; No horizontal interpolation
  ;----------------------------
  if(horz_interp .eq. "no_interp") then
    yy        = x1
    l_interph = True
  end if

  ;--------
  ; Binning
  ;--------
  if(horz_interp .eq. "binning") then
    bin_factor = field@bin_factor(ifield)
    if(bin_factor .lt. 0.) then
      bin_factor = (plato*plono*1.)/(dim_lat*dim_lon*1.)
      bin_factor = max( (/ bin_factor,1. /))
    end if
    print("     Binning grid boxes.  Bin factor = "+bin_factor)
    if( field@loutput(ifield) .and. .not. isatt(cdf->$field(ifield)$,"bin_factor") ) then
      cdf->$field(ifield)$@bin_factor = bin_factor
    end if
    opt            = True
    opt@bin_factor = bin_factor
    yy             = area_conserve_remap (xx&lon, xx&lat, x1, flon, flat, opt)

    l_interph      = True
  end if

  ;-------------------------------------------------
  ; If grid has pole points, fill with zonal average
  ;-------------------------------------------------
  if(horz_interp .ne. "no_interp") then
    if (fabs( flat(0) + 90.) .lt. 1.e-7 ) then
      print("     Filling South Pole points with zonal average")
      do kk = 0,dim_lev-1
        sum1 = dim_avg(yy(lat|0 ,lev|kk ,lon|: ))
        yy(lat|0 ,lev|kk ,lon|: ) = sum1
      end do
    end if
    if (fabs( flat(plato-1) - 90.) .lt. 1.e-7 ) then
      print("     Filling North Pole points with zonal average")
      do kk = 0,dim_lev-1
        sum2 = dim_avg(yy(lat|plato-1 ,lev|kk ,lon|: ))
        yy(lat|plato-1 ,lev|kk ,lon|: ) = sum2
      end do
    end if
  end if

  ;------------------------------------------------
  ; if 2-D, strip "lev" dimension from output array
  ;------------------------------------------------
  if(dims .eq. 2) then
    tmp = new( (/plato  , plono /), double)
    tmp = yy(lev|0,lat|:,lon|:)
    delete(yy)
    yy  = tmp
    delete(tmp)
  end if

  ;----------------------
  ; Add coordinate values
  ;----------------------
  yy&lat = flat
  yy&lon = flon
  if(dims .eq. 3) then
    yy&lev = xx&lev
  end if

  ;------------------
  ; Modify attributes
  ;------------------
  if(isatt(yy,"_FillValue")) then
    delete(yy@_FillValue)
  end if

  if(.not. l_interph) then
    print("       Error:  field not processed; '"+horz_interp+ \
                          "' is not a valid interpolation flag")
    exit
  end if

  delete(xx)
  delete(x1)

  return(yy)
end
;=====================================================================================


;=====================================================================================
function spectral_tr(xfield:numeric, cdf:file, field:string, ifield:numeric)
;
;***********************************************************************;
;                                                                       ;
; spectral_tr                                                           ;
;                                                                       ;
;            xfield  : numeric; INPUT field                             ;
;            cdf     : file   ; NCL Netcdf OUTPUT file                  ;
;            field   : string ; Names of of OUTPUT field variables      ;
;            ifield  : integer; OUTPUT field index                      ;
;                                                                       ;
; Spectrally truncate/interpolate input field to output field           ;
;                                                                       ;
;***********************************************************************;

local yy, dims, dim_lev

begin
  ;-----------------------
  ; No spectral truncation
  ;-----------------------
  if(field@spec_trunc(ifield) .lt. 0) then
    yy = xfield
    return(yy)
  end if

  dims = dimsizes(dimsizes(xfield))
  if(dims .eq. 3) then
    dim_lev = dimsizes(xfield&lev)
    yy      = new( (/dim_lev, field@plato  ,field@plono /), double)
    yy!0    = "lev"
    yy!1    = "lat"
    yy!2    = "lon"
  else
    yy      = new( (/         field@plato  ,field@plono /), double)
    yy!0    = "lat"
    yy!1    = "lon"
  end if

  ;----------------------------------------------------------
  ; Copy "from_hires" attribute (if it exists) to output file
  ;----------------------------------------------------------
  if(isatt(xfield,"from_hires")) then
    cdf->$field(ifield)$@from_hires = xfield@from_hires
  end if

  ;-------------------
  ; Process input data
  ;-------------------
  print("     SPECTRAL horizontal TRUNCATION for "+field(ifield)+ \
        " (truncation at wave"+" "+field@spec_trunc(ifield) +")")

  yy      = g2gsh(xfield, (/field@plato,field@plono/), field@spec_trunc(ifield))

  cdf->$field(ifield)$@triangular_truncation = field@spec_trunc(ifield)
  yy&lat  = field@lat
  yy&lon  = field@lon
  if(dims .eq. 3) then
    yy&lev = xfield&lev
  end if

  ;------------------
  ; Modify attributes
  ;------------------
  if(isatt(yy,"_FillValue")) then
    delete(yy@_FillValue)
  end if

  return(yy)
end
;=====================================================================================


;=====================================================================================
procedure spectral_v_tr(xfieldu:numeric, xfieldv:numeric, cdf:file, field:string, \
                        ifieldu:numeric, ifieldv:numeric, yy_u:numeric, yy_v:numeric)
;
;***********************************************************************;
;                                                                       ;
; vector spectral truncation                                            ;
;                                                                       ;
;            xfieldu : numeric; INPUT U field                           ;
;            xfieldv : numeric; INPUT V field                           ;
;            cdf     : file   ; NCL Netcdf OUTPUT file                  ;
;            field   : string ; Names of of OUTPUT field variables      ;
;            ifieldu : integer; OUTPUT U-field index                    ;
;            ifieldv : integer; OUTPUT V-field index                    ;
;            yy_u    : numeric; OUTPUT U field                          ;
;            yy_v    : numeric; OUTPUT V field                          ;
;                                                                       ;
; Spectrally truncate u/v fields as vector pair                         ;
;                                                                       ;
;***********************************************************************;

local xx_u, xx_v, dim_lev, dim_lat, dim_lon

begin
  ;-----------------------
  ; No spectral truncation
  ;-----------------------
  if(field@spec_trunc(ifieldu) .lt. 0) then
    yy_u = xfieldu
    yy_v = xfieldv
    return
  end if
  xx_u   = xfieldu
  xx_v   = xfieldv

  ;-----------------------------------
  ; Determine dimensions of input file
  ;-----------------------------------
  dim_lev = dimsizes(xx_u&lev)
  dim_lat = dimsizes(xx_u&lat)
  dim_lon = dimsizes(xx_u&lon)

  ;---------------------------
  ; Vector SPECTRAL truncation
  ;---------------------------
  print("     Vector SPECTRAL horizontal truncation for U/V (truncation at wave"+ \
        " "+field@spec_trunc(ifieldu) +")")

  yy_u!0    = "lev"
  yy_u!1    = "lat"
  yy_u!2    = "lon"

  yy_v!0    = "lev"
  yy_v!1    = "lat"
  yy_v!2    = "lon"

  g2gshv(xx_u, xx_v, yy_u, yy_v, field@spec_trunc(ifieldu))

  cdf->$field(ifieldu)$@triangular_truncation = field@spec_trunc(ifieldu)
  cdf->$field(ifieldv)$@triangular_truncation = field@spec_trunc(ifieldv)

  ;----------------------
  ; Add coordinate values
  ;----------------------
  yy_u&lev = xx_u&lev
  yy_u&lat = field@lat
  yy_u&lon = field@lon

  yy_v&lev = xx_v&lev
  yy_v&lat = field@lat
  yy_v&lon = field@lon

  ;------------------
  ; Modify attributes
  ;------------------
  if(isatt(yy_u,"_FillValue")) then
    delete(yy_u@_FillValue)
  end if
  if(isatt(yy_v,"_FillValue")) then
    delete(yy_v@_FillValue)
  end if

  delete(xx_u)
  delete(xx_v)

  return
end
;=====================================================================================


;=====================================================================================
procedure output2CAM(xfield:numeric, cdf:file, field:string, ifield:numeric, \
                     count:integer, nfields_out:integer )
;
;***********************************************************************;
;                                                                       ;
; output2CAM                                                            ;
;                                                                       ;
;            xfield  : numeric; OUTPUT field                            ;
;            cdf     : file   ; NCL Netcdf OUTPUT file                  ;
;            field   : string ; Names of of OUTPUT field variables      ;
;            ifield  : integer; OUTPUT field index                      ;
;            count   : integer; Field counter                           ;
;            nfields_out : integer; number of fields output to HT       ;
;                                                                       ;
; Dump field to output file                                             ;
;                                                                       ;
;***********************************************************************;

local dims

begin
  if(.not. field@loutput(ifield)) then
    return
  end if
  count = count + 1

  ;---------------------------------------------------------
  ; Check for missing data and abort if any data are missing
  ;---------------------------------------------------------
  xfield@_FillValue = -9999
  if(any(ismissing(xfield))) then
    print ("     Error:  Resultant interpolated field contains missing data")
    exit
  end if
  xfield@_FillValue = -999
  if(any(ismissing(xfield))) then
    print ("     Error:  Resultant interpolated field contains missing data")
    exit
  end if

  delete(xfield@_FillValue)

  ;-------------------------------
  ; Check if field exceeds max/min
  ;-------------------------------
  if(min(xfield) .lt. field@min(ifield) ) then
    print ("     Error:  Resultant interpolated field contains value less than expected minimum")
    print ("        field minimum = "+min(xfield) )
    print ("        expected minimum for this field = "+field@min(ifield) )
    print ("")
    print ("If the field, "+field(ifield)+", is not in error, modify the min value in 'makeIC_Create_field_Master_List.ncl'")
    print ("")
    print ("")
    exit
  end if
  if(max(xfield) .gt. field@max(ifield) ) then
    print ("     Error:  Resultant interpolated field contains value greater than expected maximum")
    print ("        field maximum = "+max(xfield) )
    print ("        expected maximum for this field = "+field@max(ifield) )
    print ("")
    print ("If the field, "+field(ifield)+", is not in error, modify the max value in 'makeIC_Create_field_Master_List.ncl'")
    print ("")
    print ("")
    exit
  end if

  dims = dimsizes(dimsizes(xfield))

  if(dims .ne. 2 .and. dims .ne. 3 .and. dims .ne. 1) then
    print("Error:  field to be output is neither 2-D nor 3-D")
    exit
  end if

  if(dims .eq. 3) then
     if(field@precision .eq. "double") then
        cdf->$field(ifield)$(0,:,:,:) = (/ xfield /)
     end if
     if(field@precision .eq. "float") then
        cdf->$field(ifield)$(0,:,:,:) = doubletofloat(xfield)
     end if
  end if
  if(dims .eq. 2) then
     if(field@precision .eq. "double") then
        cdf->$field(ifield)$(0,:,:  ) = (/ xfield /)
     end if
     if(field@precision .eq. "float") then
        cdf->$field(ifield)$(0,:,:  ) = doubletofloat(xfield)
     end if
  end if
  if(dims .eq. 1) then
     if(field@precision .eq. "double") then
        cdf->$field(ifield)$(0,:  ) = (/ xfield /)
     end if
     if(field@precision .eq. "float") then
        cdf->$field(ifield)$(0,:  ) = doubletofloat(xfield)
     end if
  end if

  print ("     Output field "+field(ifield)+" completed (field "+count+" of "+nfields_out+")")
  print ("")

  return
end
;=====================================================================================


;=====================================================================================
function process_Q_horz(cdf:file        , data:file     , xfield:numeric , \
                        field:string    , ifile:integer , ifield:integer , \
                        ftype:string    , ps_xx:numeric , ps_yy:numeric  , \
                        T_input:numeric , T_horz_interp_only:numeric )
;
;*******************************************************************************;
;                                                                               ;
; process_Q_horz                                                                ;
;                                                                               ;
;     cdf                : file    ; Output file                                ;
;     data               : file    ; Input file                                 ;
;     xfield             : numeric ; input Q                                    ;
;     field              : string  ; name of field to be read from input file   ;
;     ifile              : integer ; input file number                          ;
;     ifield             : integer ; field counter                              ;
;     ftype              : string  ; type of file                               ;
;     ps_xx              : numeric ; input surface pressure                     ;
;     ps_yy              : numeric ; unadjusted interpolated Ps                 ;
;     T_input            : numeric,; input T                                    ;
;     T_horz_interp_only : numeric ; T that is only Horizontally                ;
;                                    interpolated                               ;
;                                                                               ;
; Input Q and interpolate horizontally                                          ;
;                                                                               ;
;                                                                               ;
;*******************************************************************************;

local dummy, press_in_xx, press_in_yy, dims, tmp, tmp1, tmp2, tmp3, yy, xfield1

begin
  dummy   = 0.0
  xfield1 = xfield

  ;-------------------------------------------------------------------------------
  ; Compute 3D array of INPUT/OUTPUT pressures (for use in vertical interpolation)
  ;-------------------------------------------------------------------------------
  if(field@vert_interp(ifield) .ne. "no_interp") then
    press_in_xx  = build_press(data, ftype(ifile), ps_xx, "FULL")
    press_in_yy  = build_press(data, ftype(ifile), ps_yy, "FULL")
  else
    press_in_xx  = dummy
    press_in_yy  = dummy
  end if

  ;------------
  ; Interpolate
  ;------------
  yy  = horz_process(xfield1, cdf, field, ifield, field@plono, field@plato, field@lon, field@lat)
  tmp = yy
  tmp = spectral_tr (yy    , cdf, field, ifield)
  yy  = tmp
  delete(tmp)
  delete(xfield1)
  delete(press_in_xx)
  delete(press_in_yy)

  return(yy)
end
;=====================================================================================


;=====================================================================================
function process_T_vert(cdf:file       , data:file      , field:string  , \
                        ifile:integer  , ifield:integer , ftype:string  , \
                        ps_yy:numeric  , phis_in:numeric,                 \
                        T_horz_interp_only:numeric , ps_yy_adjust:numeric )
;
;*******************************************************************************;
;                                                                               ;
; process_T_vert                                                                ;
;                                                                               ;
;     cdf                : file    ; Output file                                ;
;     data               : file    ; Input file                                 ;
;     field              : string  ; name of field to be read from input file   ;
;     ifile              : integer ; input file number                          ;
;     ifield             : integer ; field counter                              ;
;     ftype              : string  ; type of file                               ;
;     ps_yy              : numeric ; unadjusted interpolated Ps                 ;
;                                    0 < RH < 100                               ;
;     phis_in            : numeric ; INPUT PHIS interpolated to output grid     ;
;     T_horz_interp_only : numeric ; T that is only Horizontally                ;
;                                    interpolated                               ;
;     ps_yy_adjust       : numeric ; Ps adjusted for differences between input/ ;
;                                    output PHIS                                ;
; Vertically interpolate T                                                      ;
;                                                                               ;
;                                                                               ;
;*******************************************************************************;

local dummy, press_in_yy, press_in_yy_half, press_out_yy, zz, dims, tmp, loglin

begin
  dummy   = 0.0

  ;-------------------------------------------------------------------------------
  ; Compute 3D array of INPUT/OUTPUT pressures (for use in vertical interpolation)
  ;-------------------------------------------------------------------------------
  if((field@dimensions (ifield).eq."3D") .and.  \
     (field@vert_interp(ifield).ne."no_interp") ) then

    ; Force surface pressure into 2D array so build_press() calls work
    ;-------------------------------------------------------------------
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
    ps_dim = dimsizes(ps_yy_adjust)
    ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_adjust,0)
    ps_conform!0="ncol"
    ps_conform!1="null"
    press_conform = build_press(cdf , "CAM" ,ps_conform, "FULL")
    press_out_yy  = press_conform(:,:,0)
    delete(press_conform)
    delete(ps_conform)
    delete(ps_dim)

;******************************************
; DIAG
;   printVarSummary(ps_yy)
;   print(ps_yy(0:10))
;   printVarSummary(ps_yy_adjust)
;   print(ps_yy_adjust(0:10))
;   print(press_in_yy(:,0))
;   print(press_in_yy_half(:,0))
;   print(press_out_yy(:,0))
; DIAG
;******************************************

  else
    press_in_yy      = dummy
    press_in_yy_half = dummy
    press_out_yy     = dummy
  end if

  dims = dimsizes( T_horz_interp_only )

  if(field@vert_interp(ifield) .eq. "quad_opt1" .or. \
    field@vert_interp(ifield) .eq. "quad_logp_opt1") then
    zz      = new( (/field@plevo , field@ncol /), double)
    zz!0    = "lev"
    zz!1    = "ncol"
    if(field@vert_interp(ifield) .eq. "quad_opt1") then
      print("     QUADRATIC (in pressure) vertical interpolation for this field (OPTION 1)")
      loglin = 1
    end if
    if(field@vert_interp(ifield) .eq. "quad_logp_opt1") then
      print("     QUADRATIC (in log pressure) vertical interpolation for this field (OPTION 1)")
      loglin = 0
    end if
    MAKEIC::vert_quad_opt1_se(dims(0) ,dims(0)+1, field@plevo,  dims(1) ,    \
                              T_horz_interp_only(ncol|: ,lev|: ),            \
                              press_in_yy       (ncol|: ,lev|: ),            \
                              press_in_yy_half  (ncol|: ,lev|: ),            \
                              press_out_yy      (ncol|: ,lev|: ),            \
                              phis_in           (ncol|: ), ps_yy  (ncol|: ), \
                              zz                (ncol|: ,lev|: ), loglin )
    zz&lev = field@lev
  else
    zz = vert_process(T_horz_interp_only , cdf, field, ifield, press_in_yy, press_out_yy)
  end if

  delete(press_in_yy)
  delete(press_in_yy_half)
  delete(press_out_yy)
  delete(dims)

  return(zz)
end
;=====================================================================================


;=====================================================================================
function process_Q_vert(cdf:file      , data:file      , field:string ,    \
                        ifile:integer , ifield:integer , ftype:string ,    \
                        ps_yy:numeric ,                                    \
                        T_horz_interp_only:numeric , T_output:numeric ,    \
                        ps_yy_adjust:numeric , Q_horz_interp_only:numeric  )
;
;*******************************************************************************;
;                                                                               ;
; process_Q_vert                                                                ;
;                                                                               ;
;     cdf                : file    ; Output file                                ;
;     data               : file    ; Input file                                 ;
;     field              : string  ; name of field to be read from input file   ;
;     ifile              : integer ; input file number                          ;
;     ifield             : integer ; field counter                              ;
;     ftype              : string  ; type of file                               ;
;     ps_yy              : numeric ; unadjusted interpolated Ps                 ;
;     T_horz_interp_only : numeric ; T that is only Horizontally                ;
;                                    interpolated                               ;
;     T_output           : numeric ; horizontally and vertically interpolated T ;
;     ps_yy_adjust       : numeric ; Ps adjusted for differences between input/ ;
;                                    output PHIS                                ;
;     Q_horz_interp_only : numeric ; Horizontally interpolated Q only           ;
;                                                                               ;
; Vertically interpolate Q                                                      ;
;                                                                               ;
;                                                                               ;
;*******************************************************************************;

local dummy, press_in_yy, press_out_yy, dims, tmp, tmp1, tmp2, tmp3, zz, xfield1

begin
  dummy   = 0.0
  xfield1 = Q_horz_interp_only

  ;-------------------------------------------------------------------------------
  ; Compute 3D array of INPUT/OUTPUT pressures (for use in vertical interpolation)
  ;-------------------------------------------------------------------------------
  if(field@vert_interp(ifield) .ne. "no_interp") then

    ; Force surface pressure into 2D array so build_press() calls work
    ;-------------------------------------------------------------------
    ps_dim = dimsizes(ps_yy)
    ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy,0)
    ps_conform!0="ncol"
    ps_conform!1="null"
    press_conform  = build_press(data,ftype(ifile),ps_conform, "FULL")
    press_in_yy  = press_conform(:,:,0)
    delete(press_conform)
    delete(ps_dim)
    ps_dim = dimsizes(ps_yy_adjust)
    ps_conform = conform_dims( (/ps_dim(0),1/),ps_yy_adjust,0)
    ps_conform!0="ncol"
    ps_conform!1="null"
    press_conform = build_press(cdf , "CAM" ,ps_conform, "FULL")
    press_out_yy  = press_conform(:,:,0)
    delete(press_conform)
    delete(ps_conform)
    delete(ps_dim)
  else
    press_in_yy  = dummy
    press_out_yy = dummy
  end if

  ;------------
  ; Interpolate
  ;------------
  zz  = vert_process(xfield1 , cdf, field, ifield, press_in_yy, press_out_yy)
  delete(xfield1)
  delete(press_in_yy)
  delete(press_out_yy)

  return(zz)
end
;=====================================================================================


;=====================================================================================
procedure horz_process_v(xfieldu:numeric, xfieldv:numeric, cdf:file, field:string, \
                         ifieldu:numeric, ifieldv:numeric,                         \
                         u_plono:numeric, u_plato:numeric,                         \
                         v_plono:numeric, v_plato:numeric,                         \
                         u_flon:numeric , u_flat:numeric ,                         \
                         v_flon:numeric , v_flat:numeric , u:numeric , v:numeric)
;
;***********************************************************************;
;                                                                       ;
; horz_process_v                                                        ;
;                                                                       ;
;            xfieldu : numeric; INPUT u-field                           ;
;            xfieldv : numeric; INPUT v-field                           ;
;            cdf     : file   ; NCL Netcdf OUTPUT file                  ;
;            field   : string ; Names of of OUTPUT field variables      ;
;            ifieldu : integer; OUTPUT u-field index                    ;
;            ifieldv : integer; OUTPUT v-field index                    ;
;            u_plono : numeric; OUTPUT longitude dimension  (U)         ;
;            u_plato : numeric; OUTPUT latitude  dimension  (U)         ;
;            v_plono : numeric; OUTPUT longitude dimension  (V)         ;
;            v_plato : numeric; OUTPUT latitude  dimension  (V)         ;
;            u_flon  : numeric; OUTPUT longitudes (deg)     (U)         ;
;            u_flat  : numeric; OUTPUT latitudes  (deg)     (U)         ;
;            v_flon  : numeric; OUTPUT longitudes (deg)     (V)         ;
;            v_flat  : numeric; OUTPUT latitudes  (deg)     (V)         ;
;            u       : numeric; OUTPUT u-field                          ;
;            v       : numeric; OUTPUT v-field                          ;
;                                                                       ;
; Horizontally interpolate/copy input vector fields to output fields    ;
;                                                                       ;
;***********************************************************************;

local xxu, xxv, x1u, x1v, l_interph, horz_interp, dim_lev,      \
      dim_lat_u, dim_lon_u, dim_lat_v, dim_lon_v, x1u_extents,  \
      x1v_extents, extent_dim, limdr, opt, bin_factor

begin
  l_interph = False
  horz_interp = field@horz_interp(ifieldu)
  xxu         = xfieldu
  xxv         = xfieldv

  ;-----------------------------------
  ; Determine dimensions of input file
  ;-----------------------------------
  dim_lat_u = dimsizes(xxu&lat)
  dim_lon_u = dimsizes(xxu&lon)
  dim_lat_v = dimsizes(xxv&lat)
  dim_lon_v = dimsizes(xxv&lon)
  dim_lev   = dimsizes(xxu&lev)

  ;----------------------------------------------------------
  ; Only "binning" can handle U and V of different dimensions
  ;----------------------------------------------------------
  if(horz_interp .ne. "no_interp" .and. horz_interp .ne. "binning") then
    if(dim_lat_u .ne. dim_lat_v  .or. dim_lon_u .ne. dim_lon_v .or. \
       any(xxu&lat .ne. xxv&lat) .or. any(xxu&lon .ne. xxv&lon) )then
       print("       Error:  U/V not processed; '"+horz_interp+ \
             "' interpolation cannot interpolate U and V of different dimensions")
       exit
    end if
  end if

  ;-----------------------
  ; Initialize work arrays
  ;-----------------------
  x1u      = new( (/dim_lev, dim_lat_u, dim_lon_u /), double)
  x1u!0    = "lev"
  x1u!1    = "lat"
  x1u!2    = "lon"
  x1v      = new( (/dim_lev, dim_lat_v, dim_lon_v /), double)
  x1v!0    = "lev"
  x1v!1    = "lat"
  x1v!2    = "lon"

  u!0      = "lev"
  v!0      = "lev"
  u!1      = "lat"
  u!2      = "lon"
  v!1      = "lat"
  v!2      = "lon"

  x1u      = xxu
  x1v      = xxv

;--------------------
; Begin interpolation
;--------------------

  ;----------------------------
  ; No horizontal interpolation
  ;----------------------------
  if(horz_interp .eq. "no_interp") then
    u         = x1u
    v         = x1v
    l_interph = True
  end if

  ;---------
  ; Spectral
  ;---------
  if(horz_interp .eq. "spectral") then
    print("     Spectral Interpolation")
    spectral_v_tr(x1u, x1v, cdf, field, ifieldu, ifieldv, u, v)
    l_interph = True
  end if

  ;--------
  ; Binning
  ;--------
  if(horz_interp .eq. "binning") then
    bin_factor = field@bin_factor(ifieldu)
    if(bin_factor .lt. 0.) then
      bin_factor = (v_plato*v_plono*1.)/(dim_lat_v*dim_lon_v*1.)
      bin_factor = max( (/ bin_factor,1. /))
    end if
    print("     Binning grid boxes.  Bin factor = "+bin_factor)
    if(field@loutput(ifieldu)) then
      cdf->$field(ifieldu)$@bin_factor = bin_factor
    end if
    if(field@loutput(ifieldv)) then
      cdf->$field(ifieldv)$@bin_factor = bin_factor
    end if
    opt            = True
    opt@bin_factor = bin_factor
    u              = area_conserve_remap (xxu&lon, xxu&lat, x1u, u_flon, u_flat, opt)
    v              = area_conserve_remap (xxv&lon, xxv&lat, x1v, v_flon, v_flat, opt)
    l_interph      = True
  end if

  ;----------------------
  ; Add coordinate values
  ;----------------------
  u&lat = u_flat
  u&lon = u_flon
  v&lat = v_flat
  v&lon = v_flon
  u&lev = xxu&lev
  v&lev = xxv&lev

  ;------------------
  ; Modify attributes
  ;------------------
  if(isatt(u,"_FillValue")) then
    delete(u@_FillValue)
  end if
  if(isatt(v,"_FillValue")) then
    delete(v@_FillValue)
  end if

  if(.not. l_interph) then
    print("       Error:  field not processed; '"+horz_interp+ \
                          "' is not a valid interpolation flag")
    exit
  end if

  delete(xxu)
  delete(xxv)
  delete(x1u)
  delete(x1v)

  return
end
;=====================================================================================


;=====================================================================================
function ESMF_pole_fix(xfield:numeric, plato:numeric, plati:numeric , ncol:numeric)
;
;**********************************************************;
;                                                          ;
; ESMF_pole_fix                                            ;
;                                                          ;
;            xfield : numeric; INPUT field                 ;
;            plato  : numeric; INPUT latitudes of xfield   ;
;            plati  : numeric; INPUT latitudes of source   ;
;            ncol   : integer; INPUT number of ploto()     ;
;**********************************************************;

local yfield, LatDlt, LatMin, LatMax, nn, AvgMin, AvgMax, xdim, rank, NumMin, NumMax


begin
  yfield = xfield
  xdim = dimsizes(xfield)
  rank = dimsizes(xdim)

  ; check if there is nothing to do
  ;---------------------------------
  if(.not.any(ismissing(yfield))) then
    return(yfield)
  end if

  ; Set a range of latitudes around the poles
  ;-------------------------------------------
  LatDlt = abs(plati(0)-plati(1))
  LatMax = max(plati)-(2*LatDlt)
  LatMin = min(plati)+(2*LatDlt)
  print(" LatMin="+LatMin+" LatMax="+LatMax+" rank="+rank)
  if(rank.gt.2) then
    print(" Error: ESMF_pole_fix: Only set up for rank 1 or 2")
    exit
  end if

  ; 2D (Lat,Lon) field
  ;----------------------
  if(rank.eq.1) then
    AvgMax = yfield(0)
    AvgMax = 0.
    NumMax = AvgMax
    AvgMin = AvgMax
    NumMin = AvgMax
    do nn=0,(ncol-1)
      if(plato(nn).gt.LatMax) then
        if(.not.ismissing(yfield(nn))) then
          NumMax=NumMax+1
          AvgMax=AvgMax+yfield(nn)
        end if
      end if
      if(plato(nn).lt.LatMin) then
        if(.not.ismissing(yfield(nn))) then
          NumMin=NumMin+1
          AvgMin=AvgMin+yfield(nn)
        end if
      end if
    end do
    if((NumMax.eq.0).or.(NumMin.eq.0)) then
      print(" Error: NumMax="+NumMax+" NumMin="+NumMin)
      exit
    end if
    AvgMax=AvgMax/NumMax
    AvgMin=AvgMin/NumMin

    do nn=0,(ncol-1)
      if(plato(nn).gt.LatMax) then
        if(ismissing(yfield(nn))) then
          yfield(nn)=AvgMax
        end if
      end if
      if(plato(nn).lt.LatMin) then
        if(ismissing(yfield(nn))) then
          yfield(nn)=AvgMin
        end if
      end if
    end do
  end if
  
  ; 3D (Lat,Lon,Lev) field
  ;------------------------
  if(rank.eq.2) then
    AvgMax = yfield(:,0)
    AvgMax(:) = 0.
    NumMax = AvgMax
    AvgMin = AvgMax
    NumMin = AvgMax
    do nn=0,(ncol-1)
      if(plato(nn).gt.LatMax) then
        if(.not.any(ismissing(yfield(:,nn)))) then
         NumMax(:)=NumMax(:)+1
         AvgMax(:)=AvgMax(:)+yfield(:,nn)
        end if
      end if
      if(plato(nn).lt.LatMin) then
        if(.not.any(ismissing(yfield(:,nn)))) then
         NumMin(:)=NumMin(:)+1
         AvgMin(:)=AvgMin(:)+yfield(:,nn)
        end if
      end if
    end do
    AvgMax=AvgMax/NumMax
    AvgMin=AvgMin/NumMin
    do nn=0,(ncol-1)
      if(plato(nn).gt.LatMax) then
        if(any(ismissing(yfield(:,nn)))) then
          yfield(:,nn)=AvgMax(:)
        end if
      end if
      if(plato(nn).lt.LatMin) then
        if(any(ismissing(yfield(:,nn)))) then
          yfield(:,nn)=AvgMin(:)
        end if
      end if
    end do
  end if

  return(yfield)
end
;=====================================================================================
