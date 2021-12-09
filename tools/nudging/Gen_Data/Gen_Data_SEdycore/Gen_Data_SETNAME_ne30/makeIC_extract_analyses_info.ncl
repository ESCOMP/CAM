undef("get_input_fields")
undef("file_suffix")
undef("extract_analyses_info")

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------

function get_input_fields(ftype:string)

;
;****************************************************************************;
;                                                                            ;
; get_input_fields                                                           ;
;  Initialize list of input fields that live on each "type" of input file    ;
;                                                                            ;
;            ftype   : string ; Analysis file type                           ;
;                                                                            ;
;****************************************************************************;

local input_fields

begin

  input_fields              = ""

  if(ftype .eq. "SE_TOPOGRAPHY") then
    input_fields@PHIS         = "PHIS"
    input_fields@LANDFRAC     = "LANDFRAC"
    input_fields@SGH          = "SGH"
    input_fields@SGH30        = "SGH30"
    input_fields@LANDM_COSLAT = "LANDM_COSLAT"
  end if
  if(ftype .eq. "YOTC_PS_Z") then
    input_fields@PS         = "LNSP_GDS50_HYBL"
    input_fields@PHIS_input = "Z_GDS50_HYBL"
  end if
  if(ftype .eq. "YOTC_sfc") then
    input_fields@SST_cpl    = "SSTK_GDS4_SFC"
    input_fields@TCW        = "TCW_GDS4_SFC"
    input_fields@TMQ        = "TCWV_GDS4_SFC"
    input_fields@SLP        = "MSL_GDS4_SFC"
    input_fields@CLDTOT     = "TCC_GDS4_SFC"
    input_fields@CLDLOW     = "LCC_GDS4_SFC"
    input_fields@CLDMED     = "MCC_GDS4_SFC"
    input_fields@CLDHGH     = "HCC_GDS4_SFC"
  end if

  if(ftype .eq. "YOTC_sh") then
    input_fields@U          = "U_GDS50_HYBL"
    input_fields@V          = "V_GDS50_HYBL"
    input_fields@US         = "U_GDS50_HYBL"
    input_fields@VS         = "V_GDS50_HYBL"
    input_fields@T          = "T_GDS50_HYBL"
    input_fields@OMEGA      = "W_GDS50_HYBL"
  end if
  if(ftype .eq. "ECMWF_gg") then
    input_fields@Q          = "Q_GDS4_HYBL"
    input_fields@RELHUM     = "RELHUM"
    input_fields@CLOUD      = "CC_GDS4_HYBL"
    input_fields@CLDLIQ     = "CLWC_GDS4_HYBL"
    input_fields@CLDICE     = "CIWC_GDS4_HYBL"
  end if
  if(ftype .eq. "ECMWF_sh") then
    input_fields@VORT       = "VO_GDS50_HYBL"
    input_fields@DIV        = "D_GDS50_HYBL"
    input_fields@U          = "VO_GDS50_HYBL"
    input_fields@V          = "D_GDS50_HYBL"
    input_fields@US         = "VO_GDS50_HYBL"
    input_fields@VS         = "D_GDS50_HYBL"
    input_fields@T          = "T_GDS50_HYBL"
    input_fields@PS         = "LNSP_GDS50_HYBL"
    input_fields@PHIS_input = "Z_GDS50_HYBL"
  end if
  if(ftype .eq. "NASA_MERRA") then
    input_fields@U          = "U"
    input_fields@V          = "V"
    input_fields@US         = "U"
    input_fields@VS         = "V"
    input_fields@T          = "T"
    input_fields@Q          = "QV"
    input_fields@PS         = "PS"
    input_fields@PHIS_input = "PHIS"
  end if
  if(ftype .eq. "NASA_MERRA_PREVOCA") then
    input_fields@U          = "uwnd"
    input_fields@V          = "vwnd"
    input_fields@US         = "uwnd"
    input_fields@VS         = "vwnd"
    input_fields@T          = "theta"
    input_fields@Q          = "sphu"
    input_fields@PS         = "ps"
    input_fields@PHIS_input = "phis"
  end if
  if(ftype .eq. "JRA_25") then
    input_fields@U          = "U_GRD_GDS4_HYBL"
    input_fields@V          = "V_GRD_GDS4_HYBL"
    input_fields@US         = "U_GRD_GDS4_HYBL"
    input_fields@VS         = "V_GRD_GDS4_HYBL"
    input_fields@T          = "TMP_GDS4_HYBL"
    input_fields@Q          = "SPF_H_GDS4_HYBL"
    input_fields@PS         = "PRES_GDS4_SFC"
    input_fields@PHIS_input = "DIST_GDS4_SFC"
  end if
  if(ftype .eq. "Era_Interim_627.0_sc") then
    input_fields@VORT       = "VO_GDS4_HYBL"
    input_fields@DIV        = "D_GDS4_HYBL"
    input_fields@U          = "VO_GDS4_HYBL"
    input_fields@V          = "D_GDS4_HYBL"
    input_fields@US         = "VO_GDS4_HYBL"
    input_fields@VS         = "D_GDS4_HYBL"
    input_fields@T          = "T_GDS4_HYBL"
    input_fields@Q          = "Q_GDS4_HYBL"
    input_fields@PS         = "LNSP_GDS4_HYBL"
    input_fields@PHIS_input = "Z_GDS4_HYBL"
    input_fields@CLOUD      = "CC_GDS4_HYBL"
    input_fields@CLDICE     = "CIWC_GDS4_HYBL"
    input_fields@CLDLIQ     = "CLWC_GDS4_HYBL"
    input_fields@OMEGA      = "W_GDS4_HYBL"
    input_fields@O3         = "O3_GDS4_HYBL"
  end if
  if(ftype .eq. "ERA40_ds117.2") then
    input_fields@VORT       = "VO_GDS50_HYBL"
    input_fields@DIV        = "D_GDS50_HYBL"
    input_fields@U          = "VO_GDS50_HYBL"
    input_fields@V          = "D_GDS50_HYBL"
    input_fields@US         = "VO_GDS50_HYBL"
    input_fields@VS         = "D_GDS50_HYBL"
    input_fields@T          = "T_GDS50_HYBL"
    input_fields@Q          = "Q_GDS4_HYBL"
    input_fields@PS         = "LNSP_GDS50_HYBL"
    input_fields@PHIS_input = "Z_GDS50_HYBL"
  end if

  if(ftype .eq. "YOTC_sfc_fcst") then
  end if

  return(input_fields)

end

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------

function file_suffix(ftype:string)

;
;****************************************************************************;
;                                                                            ;
; file_suffix                                                                ;
;  Initialize suffix string depending upon type of file being read           ;
;  (this is to add a suffix to filename when using "addfile" if suffix does  ;
;   not already exist on the filename (e.g., ".nc",".grb",".hdf", etc) )     ;
;                                                                            ;
;            ftype   : string ; Analysis file type                           ;
;                                                                            ;
;****************************************************************************;

local suffix

begin

  suffix = ""
  if(ftype .eq. "CAM") then
    suffix = ""
  end if
  if(ftype .eq. "SE_TOPOGRAPHY") then
    suffix = ""
  end if
  if(ftype .eq. "YOTC_PS_Z") then
    suffix = ".grb"
  end if
  if(ftype .eq. "YOTC_sfc") then
    suffix = ".grb"
  end if
  if(ftype .eq. "YOTC_sfc_fcst") then
    suffix = ".grb"
  end if
  if(ftype .eq. "YOTC_sh") then
    suffix = ".grb"
  end if
  if(ftype .eq. "ECMWF_gg") then
    suffix = ".grb"
  end if
  if(ftype .eq. "ECMWF_sh") then
    suffix = ""
  end if
  if(ftype .eq. "NASA_MERRA") then
    suffix = ""
  end if
  if(ftype .eq. "NASA_MERRA_PREVOCA") then
    suffix = ""
  end if
  if(ftype .eq. "JRA_25") then
    suffix = ".grb"
  end if
  if(ftype .eq. "Era_Interim_627.0_sc") then
    suffix = ".grb"
  end if
  if(ftype .eq. "ERA40_ds117.2") then
    suffix = ".grb"
  end if

  return(suffix)

end

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------

function extract_analyses_info(data:file, ftype:string)

;
;****************************************************************************;
;                                                                            ;
; extract_analyses_info                                                      ;
;  extract x, y, z dimension and other pertinent info from analysis file     ;
;  and add to "att_array" meta-data                                          ;
;                                                                            ;
;            data    : file   ; Analysis file                                ; 
;            ftype   : string ; Analysis file type                           ;
;                                                                            ;
;****************************************************************************;

local k, native_dimname, processed, att_array, file_dim_sizes, file_dim_names, \
      dim_lat, dim_lon, del_lon, gau_info, date_int, datei, houri, hourf,      \
      dims, nlev, date_read, datef

begin

  processed                        = False
  att_array                        = "Coordinates"

;---------------------------
; Initialize attribute space
;---------------------------

  att_array@time_dimname           = ""       ; Input time dimension name
  att_array@real_imaginary_dimname = ""       ; Input R/I  dimension name
  att_array@lev_dimname            = ""       ; Input lev  dimension name
  att_array@lev_dimname_alt1       = ""       ; Input lev  dimension name (alternate level     (.e.g., "ilev" in CAM files) )
  att_array@lat_dimname            = ""       ; Input lat  dimension name
  att_array@lat_dimname_alt1       = ""       ; Input lat  dimension name (alternate latitude  (.e.g., "slat" in CAM files) )
  att_array@lon_dimname            = ""       ; Input lon  dimension name
  att_array@lon_dimname_alt1       = ""       ; Input lon  dimension name (alternate longitude (.e.g., "slon" in CAM files) )

  att_array@logPs_to_Ps            = False    ; If input is in ln(Ps) form, convert to Ps
  att_array@logPs_var_names        = ""       ; Name of ln(Ps) variable(s) to be converted
  att_array@convert_theta_to_T     = False    ; Convert theta to T
  att_array@theta_var_names        = ""       ; Name of theta variable(s) to be converted
  att_array@Z_to_PHIS              = False    ; Convert Z to PHIS
  att_array@Z_var_names            = ""       ; Name of Z variable(s) to be converted
  att_array@remove_cos_wgt         = False    ; Remove cosine weighting from variable
  att_array@cos_wgted_var_names    = ""       ; Name of cosine weighted variable(s) to be converted
  att_array@PREVOCA_u_v_bugfix     = False    ; Fix bug in NASA_MERRA_PREVOCA u/v data
  att_array@PREVOCA_u_v_var_names  = ""       ; u/v name(s) in NASA_MERRA_PREVOCA datasets to be fixed
  att_array@factor                 = 0.       ; multiplying factor
  att_array@factor_var_names       = ""       ; Name of variable(s) to be modified

;-----------------------------------
; Extract info from an analysis file
;-----------------------------------

  if(ftype .eq. "CAM") then
 
    ; Standard CAM Eulerian or FV file

    file_dim_names                   = getvardims     (data)

    native_dimname                   = "time"
    att_array@time_dimname           = native_dimname
    att_array@time                   = data->$native_dimname$

    native_dimname                   = "lev"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$

    native_dimname                   = "ilev"
    att_array@lev_dimname_alt1       = native_dimname
    att_array@lev_alt1               = data->$native_dimname$

    native_dimname                   = "lat"
    att_array@lat_dimname            = native_dimname
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "slat"
    if( any(file_dim_names .eq. native_dimname) ) then
      att_array@lat_dimname_alt1     = native_dimname
      att_array@lat_alt1             = data->$native_dimname$
    end if

    native_dimname                   = "lon"
    att_array@lon_dimname            = native_dimname
    att_array@lon                    = data->$native_dimname$

    native_dimname                   = "slon"
    if( any(file_dim_names .eq. native_dimname) ) then
      att_array@lon_dimname_alt1     = native_dimname
      att_array@lon_alt1             = data->$native_dimname$
    end if

    att_array@hyam                   = data->hyam
    att_array@hybm                   = data->hybm
    att_array@hyai                   = data->hyai
    att_array@hybi                   = data->hybi

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    att_array@datei                  = data->date
    att_array@houri                  = data->datesec*1
    att_array@hourf                  = data->datesec*1.

    processed                        = True

  end if

  if(ftype .eq. "SE_TOPOGRAPHY") then

    ; CAM-SE Topography file

    native_dimname                   = "ncol"
    att_array@ncol_dimname           = native_dimname
    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)
    att_array@ncol                   = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )

    att_array@lon                    = data->lon
    att_array@lat                    = data->lat

    processed                        = True

  end if

  if(ftype .eq. "YOTC_PS_Z") then

    ; YOTC website   :  http://dss.ucar.edu/datasets/ds629.0/
    ; example dataset:  mss:/DSS/DS629.0/yt.oper.an.ml/2008/2008050100_000_018_129128_152128_799.yt.oper.an.ml

    native_dimname                   = "initial_time0_hours"
    att_array@time_dimname           = native_dimname
    att_array@time                   = data->$native_dimname$

    native_dimname                   = "real_imaginary"
    att_array@real_imaginary_dimname = native_dimname

    native_dimname                   = "g50_lat_2"
    att_array@lat_dimname            = native_dimname
    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)
    dim_lat                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    gau_info                         = gaus(dim_lat/2)
    att_array@lat                    = new( (/dim_lat  /) , double)
    att_array@lat                    = gau_info(:,0)

    native_dimname                   = "g50_lon_3"
    att_array@lon_dimname            = native_dimname
    dim_lon                          = dim_lat*2
    del_lon                          = (360.D/dim_lon)
    att_array@lon                    = new( (/dim_lon/) , double)
    att_array@lon                    = ispan(0,dim_lon-1,1)*del_lon

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    date_int                         = chartoint( stringtochar( data->initial_time0 ) ) - 48
    dims                             = dimsizes(date_int)
    if(dims(0) .ne. 4) then
      print("ERROR:  Expected number of time slices is 4")
      print("        Got "+dims(0)+" instead")
      print("date_int = "+date_int)
      exit
    end if

    datei                            = new( 4, integer)
    houri                            = new( 4, integer)

    do i                             = 0,3
      datei(i)                       =            date_int(i, 6)*10000000
      datei(i)                       = datei(i) + date_int(i, 7)*1000000
      datei(i)                       = datei(i) + date_int(i, 8)*100000
      datei(i)                       = datei(i) + date_int(i, 9)*10000
      datei(i)                       = datei(i) + date_int(i, 0)*1000
      datei(i)                       = datei(i) + date_int(i, 1)*100
      datei(i)                       = datei(i) + date_int(i, 3)*10
      datei(i)                       = datei(i) + date_int(i, 4)
      
      houri(i)                       =            date_int(i,12)*3600*10
      houri(i)                       = houri(i) + date_int(i,13)*3600
      houri(i)                       = houri(i) + date_int(i,15)*60*10
      houri(i)                       = houri(i) + date_int(i,16)*60
    end do
    hourf                            = houri * 1.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    att_array@logPs_to_Ps            = True
    att_array@logPs_var_names        = "LNSP_GDS50_HYBL"

    processed                        = True

  end if

  if(ftype .eq. "YOTC_sfc") then

    ; YOTC website   :  http://dss.ucar.edu/datasets/ds629.0/
    ; example dataset:  mss:/DSS/DS629.0/yt.oper.an.sfc/2009/2009011500_000_018_015128_238128_400.yt.oper.an.sfc

    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)

    if( any(file_dim_names .eq. "initial_time0_hours") ) then
      native_dimname                 = "initial_time0_hours"
      att_array@time_dimname         = native_dimname
      att_array@time                 = data->$native_dimname$
    end if

    if( any(file_dim_names .eq. "g4_lat_0") ) then
      native_dimname                 = "g4_lat_0"
    else
      native_dimname                 = "g4_lat_1"
    end if
    att_array@lat_dimname            = native_dimname
    dim_lat                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    att_array@lat                    = new( (/dim_lat/) , double)
    att_array@lat                    = data->$native_dimname$

    if( any(file_dim_names .eq. "g4_lat_0") ) then
      native_dimname                 = "g4_lon_1"
    else
      native_dimname                 = "g4_lon_2"
    end if
    att_array@lon_dimname            = native_dimname
    dim_lon                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    att_array@lon                    = new( (/dim_lon/) , double)
    att_array@lon                    = data->$native_dimname$

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    if( any(file_dim_names .eq. "initial_time0_hours") ) then
      date_int                         = chartoint( stringtochar( data->initial_time0 ) ) - 48
      dims                             = dimsizes(date_int)
      if(dims(0) .ne. 4) then
        print("ERROR:  Expected number of time slices is 4")
        print("        Got "+dims(0)+" instead")
        print("date_int = "+date_int)
        exit
      end if

      datei                            = new( 4, integer)
      houri                            = new( 4, integer)

      do i                             = 0,3
        datei(i)                       =            date_int(i, 6)*10000000
        datei(i)                       = datei(i) + date_int(i, 7)*1000000
        datei(i)                       = datei(i) + date_int(i, 8)*100000
        datei(i)                       = datei(i) + date_int(i, 9)*10000
        datei(i)                       = datei(i) + date_int(i, 0)*1000
        datei(i)                       = datei(i) + date_int(i, 1)*100
        datei(i)                       = datei(i) + date_int(i, 3)*10
        datei(i)                       = datei(i) + date_int(i, 4)
        
        houri(i)                       =            date_int(i,12)*3600*10
        houri(i)                       = houri(i) + date_int(i,13)*3600
        houri(i)                       = houri(i) + date_int(i,15)*60*10
        houri(i)                       = houri(i) + date_int(i,16)*60
      end do

      hourf                            = houri * 1.

    else

      date_int                         = chartoint( stringtochar( data->SSTK_GDS4_SFC@initial_time ) ) - 48

      datei                            =         date_int(6 )*10000000
      datei                            = datei + date_int(7 )*1000000
      datei                            = datei + date_int(8 )*100000
      datei                            = datei + date_int(9 )*10000
      datei                            = datei + date_int(0 )*1000
      datei                            = datei + date_int(1 )*100
      datei                            = datei + date_int(3 )*10
      datei                            = datei + date_int(4 )

      houri                            =         date_int(12)*3600*10
      houri                            = houri + date_int(13)*3600
      houri                            = houri + date_int(15)*60*10
      houri                            = houri + date_int(16)*60
      hourf                            = houri * 1.

    end if

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    processed                        = True

  end if

  if(ftype .eq. "YOTC_sfc_fcst") then

    ; YOTC website   :  http://dss.ucar.edu/datasets/ds629.0/
    ; example dataset:  mss:/DSS/DS629.0/yt.oper.fc.sfc/2009/2009011512_000_021_003228_245128_400.yt.oper.fc.sfc

    native_dimname                   = "forecast_time0"
    att_array@time_dimname           = native_dimname
    att_array@time                   = data->$native_dimname$

    native_dimname                   = "g4_lat_1"
    att_array@lat_dimname            = native_dimname
    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)
    dim_lat                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    att_array@lat                    = new( (/dim_lat  /) , double)
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "g4_lon_2"
    att_array@lon_dimname            = native_dimname
    dim_lon                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    att_array@lon                    = new( (/dim_lon/) , double)
    att_array@lon                    = data->$native_dimname$

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    date_int                         = chartoint( stringtochar( data->TTR_GDS4_SFC@initial_time ) ) - 48
    date_hours                       = data->forecast_time0*3600
    dims                             = dimsizes(date_hours)

    datei                            = new( dims, integer)
    houri                            = new( dims, integer)

    do i                             = 0,dims-1
      datei(i)                       =            date_int( 6)*10000000
      datei(i)                       = datei(i) + date_int( 7)*1000000
      datei(i)                       = datei(i) + date_int( 8)*100000
      datei(i)                       = datei(i) + date_int( 9)*10000
      datei(i)                       = datei(i) + date_int( 0)*1000
      datei(i)                       = datei(i) + date_int( 1)*100
      datei(i)                       = datei(i) + date_int( 3)*10
      datei(i)                       = datei(i) + date_int( 4)
      
      houri(i)                       =            date_int(12)*3600*10
      houri(i)                       = houri(i) + date_int(13)*3600
      houri(i)                       = houri(i) + date_int(15)*60*10
      houri(i)                       = houri(i) + date_int(16)*60

      houri(i)                       = houri(i) + date_hours(i)
    end do
    hourf                            = houri * 1.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    delete(att_array@factor_var_names)
    att_array@factor                 = 1./21600.
    att_array@factor_var_names       = (/ "TP_GDS4_SFC", "TTR_GDS4_SFC" /)

    processed                        = True

  end if

  if(ftype .eq. "YOTC_sh") then

    ; YOTC website   :  http://dss.ucar.edu/datasets/ds629.0/
    ; example dataset:  mss:/DSS/DS629.0/yt.oper.an.ml/2008/2008050100_000_018_130128_132128_799.yt.oper.an.ml

    native_dimname                   = "initial_time0_hours"
    att_array@time_dimname           = native_dimname
    att_array@time                   = data->$native_dimname$

    native_dimname                   = "real_imaginary"
    att_array@real_imaginary_dimname = native_dimname

    native_dimname                   = "lv_HYBL1"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$

    native_dimname                   = "g50_lat_5"
    att_array@lat_dimname            = native_dimname
    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)
    dim_lat                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    gau_info                         = gaus(dim_lat/2)
    att_array@lat                    = new( (/dim_lat  /) , double)
    att_array@lat                    = gau_info(:,0)

    native_dimname                   = "g50_lon_6"
    att_array@lon_dimname            = native_dimname
    dim_lon                          = dim_lat*2
    del_lon                          = (360.D/dim_lon)
    att_array@lon                    = new( (/dim_lon/) , double)
    att_array@lon                    = ispan(0,dim_lon-1,1)*del_lon

    att_array@hyam                   = data->lv_HYBL1_a
    att_array@hybm                   = data->lv_HYBL1_b
    att_array@hyai                   = data->lv_HYBL_i2_a
    att_array@hybi                   = data->lv_HYBL_i2_b

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    date_int                         = chartoint( stringtochar( data->initial_time0 ) ) - 48
    dims                             = dimsizes(date_int)
    if(dims(0) .ne. 4) then
      print("ERROR:  Expected number of time slices is 4")
      print("        Got "+dims(0)+" instead")
      print("date_int = "+date_int)
      exit
    end if

    datei                            = new( 4, integer)
    houri                            = new( 4, integer)

    do i                             = 0,3
      datei(i)                       =            date_int(i, 6)*10000000
      datei(i)                       = datei(i) + date_int(i, 7)*1000000
      datei(i)                       = datei(i) + date_int(i, 8)*100000
      datei(i)                       = datei(i) + date_int(i, 9)*10000
      datei(i)                       = datei(i) + date_int(i, 0)*1000
      datei(i)                       = datei(i) + date_int(i, 1)*100
      datei(i)                       = datei(i) + date_int(i, 3)*10
      datei(i)                       = datei(i) + date_int(i, 4)
      
      houri(i)                       =            date_int(i,12)*3600*10
      houri(i)                       = houri(i) + date_int(i,13)*3600
      houri(i)                       = houri(i) + date_int(i,15)*60*10
      houri(i)                       = houri(i) + date_int(i,16)*60
    end do
    hourf                            = houri * 1.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    delete(att_array@cos_wgted_var_names)
    att_array@remove_cos_wgt         = True
    att_array@cos_wgted_var_names    = (/ "U_GDS50_HYBL", "V_GDS50_HYBL" /)

    processed                        = True

  end if

  if(ftype .eq. "ECMWF_gg") then

    ; YOTC website   :  http://dss.ucar.edu/datasets/ds629.0/
    ; example dataset:  mss:/DSS/DS629.0/yt.oper.an.ml/2008/2008050100_000_000_133128_248128_400.yt.oper.an.ml

    ; ECMWF VOCALS
    ; example dataset:  mss:/HANNAY/data/ECMWF/vocals/ECMWF_VOCALS_20061006_00_ml_gg.grb

    ; ECMWF e10
    ; example dataset:  mss:/HANNAY/data/ECMWF/e10/ECMWF_e10_GPCI_20030802_12_ml_gg.grb

    ; ECMWF EPIC (Martin Koehler data)
    ; example dataset:  mss:/HANNAY/data/ECMWF/EPIC/ECMWF_EPIC_20011016_18_ml_gg.grb

    ; ECMWF original
    ; example dataset:  mss:/HANNAY/data/ECMWF/analysis/ECMWF_GPCI_20030622_00_ml_gg.grb

    native_dimname                   = "lv_HYBL0"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$

    native_dimname                   = "g4_lat_3"
    att_array@lat_dimname            = native_dimname
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "g4_lon_4"
    att_array@lon_dimname            = native_dimname
    att_array@lon                    = data->$native_dimname$

    att_array@hyam                   = data->lv_HYBL0_a
    att_array@hybm                   = data->lv_HYBL0_b
    att_array@hyai                   = data->lv_HYBL_i1_a
    att_array@hybi                   = data->lv_HYBL_i1_b

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float
    ; --> Read time stamp from the Q_GDS4_HYBL variable directly

    date_int                         = chartoint( stringtochar( data->Q_GDS4_HYBL@initial_time ) ) - 48

    datei                            =         date_int(6 )*10000000
    datei                            = datei + date_int(7 )*1000000
    datei                            = datei + date_int(8 )*100000
    datei                            = datei + date_int(9 )*10000
    datei                            = datei + date_int(0 )*1000
    datei                            = datei + date_int(1 )*100
    datei                            = datei + date_int(3 )*10
    datei                            = datei + date_int(4 )

    houri                            =         date_int(12)*3600*10
    houri                            = houri + date_int(13)*3600
    houri                            = houri + date_int(15)*60*10
    houri                            = houri + date_int(16)*60
    hourf                            = houri * 1.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    processed                        = True

  end if

  if(ftype .eq. "ECMWF_sh") then

    ; ECMWF VOCALS
    ; example dataset:  mss:/HANNAY/data/ECMWF/vocals/ECMWF_VOCALS_20061006_00_ml_sh.grb

    ; ECMWF e10
    ; example dataset:  mss:/HANNAY/data/ECMWF/e10/ECMWF_e10_GPCI_20030802_12_ml_sh.grb

    ; ECMWF EPIC (Martin Koehler data)
    ; example dataset:  mss:/HANNAY/data/ECMWF/EPIC/ECMWF_EPIC_20011016_18_ml_sh.grb

    ; ECMWF original
    ; example dataset:  mss:/HANNAY/data/ECMWF/analysis/ECMWF_GPCI_20030622_00_ml_sh.grb

    native_dimname                   = "real_imaginary"
    att_array@real_imaginary_dimname = native_dimname

    native_dimname                   = "lv_HYBL3"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$

    native_dimname                   = "g50_lat_1"
    att_array@lat_dimname            = native_dimname
    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)
    dim_lat                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    gau_info                         = gaus(dim_lat/2)
    att_array@lat                    = new( (/dim_lat  /) , double)
    att_array@lat                    = gau_info(:,0)

    native_dimname                   = "g50_lon_2"
    att_array@lon_dimname            = native_dimname
    dim_lon                          = dim_lat*2
    del_lon                          = (360.D/dim_lon)
    att_array@lon                    = new( (/dim_lon/) , double)
    att_array@lon                    = ispan(0,dim_lon-1,1)*del_lon

    att_array@hyam                   = data->lv_HYBL3_a
    att_array@hybm                   = data->lv_HYBL3_b
    att_array@hyai                   = data->lv_HYBL_i4_a
    att_array@hybi                   = data->lv_HYBL_i4_b

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float
    ; --> Read time stamp from the Z_GDS50_HYBL variable directly

    date_int                         = chartoint( stringtochar( data->Z_GDS50_HYBL@initial_time ) ) - 48

    datei                            =         date_int(6 )*10000000
    datei                            = datei + date_int(7 )*1000000
    datei                            = datei + date_int(8 )*100000
    datei                            = datei + date_int(9 )*10000
    datei                            = datei + date_int(0 )*1000
    datei                            = datei + date_int(1 )*100
    datei                            = datei + date_int(3 )*10
    datei                            = datei + date_int(4 )

    houri                            =         date_int(12)*3600*10
    houri                            = houri + date_int(13)*3600
    houri                            = houri + date_int(15)*60*10
    houri                            = houri + date_int(16)*60
    hourf                            = houri * 1.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    att_array@logPs_to_Ps            = True
    att_array@logPs_var_names        = "LNSP_GDS50_HYBL"

    processed                        = True

  end if

  if(ftype .eq. "NASA_MERRA") then
 
    ; NASA MERRA website   :  http://disc.sci.gsfc.nasa.gov/daac-bin/DataHoldings.pl?LOOKUPID_List=MAI6NVANA
    ;
    ; example datasets     :      /project/cas/shea/MERRA/YOTC/MERRA300.prod.assim.const_2d_asm_Nx.00000000.hdf
    ;                             /project/cas/shea/MERRA/YOTC/MERRA300.prod.assim.inst6_3d_ana_Nv.20080515.hdf
    ;                         mss:/OLSON/MERRA/YOTC/MERRA300.prod.assim.inst6_3d_ana_Nv.20091103.hdf

    file_dim_names                   = getvardims     (data)

    native_dimname                   = "TIME_EOSGRID"
    att_array@time_dimname           = native_dimname
    att_array@time                   = data->$native_dimname$

    native_dimname                   = "Height_EOSGRID"
    if( any(file_dim_names .eq. native_dimname) ) then
      att_array@lev_dimname          = native_dimname
      att_array@lev                  = data->$native_dimname$
    end if

    native_dimname                   = "YDim_EOSGRID"
    att_array@lat_dimname            = native_dimname
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "XDim_EOSGRID"
    att_array@lon_dimname            = native_dimname
    att_array@lon                    = data->$native_dimname$

    if(isatt(att_array,"lev")) then

      att_array@hyai                 = data@ak/100000.D
      att_array@hybi                 = data@bk*1.D

      nlev                           = dimsizes(att_array@hyai) - 1

      att_array@hyam                 = new( nlev  , double)
      do k                           = 0,nlev-1
        att_array@hyam(k)            = ( att_array@hyai(k) + att_array@hyai(k+1) )*0.5D
      end do

      att_array@hybm                 = new( nlev  , double)
      do k                           = 0,nlev-1
        att_array@hybm(k)            = ( att_array@hybi(k) + att_array@hybi(k+1) )*0.5D
      end do

    end if

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    houri                            = doubletointeger(data->TIME_EOSGRID+0.1)
    hourf                            = houri*1.
    datei                            = houri
    datei                            = data->TIME_EOSGRID@begin_date
    dims                             = dimsizes(datei)
    if(max( houri - (houri/60)*60 ) .ne. 0) then
      print("ERROR:  'TIME_EOSGRID' does not contain an exact hour time stamp")
      print(data->TIME_EOSGRID)
      exit
    end if
    houri                            = (houri/60 )*3600
    hourf                            = (hourf/60.)*3600.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    processed                        = True

  end if

  if(ftype .eq. "NASA_MERRA_PREVOCA") then
 
    ; NASA MERRA Prevoca
    ;
    ; example datasets:  mss:/HANNAY/data/nasa/prevoca/d5_b10p19_jul06.ana.eta.20061020_12z.hdf

    native_dimname                   = "time"
    att_array@time_dimname           = native_dimname
    att_array@time                   = data->$native_dimname$

    native_dimname                   = "lev"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$

    native_dimname                   = "lat"
    att_array@lat_dimname            = native_dimname
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "lon"
    att_array@lon_dimname            = native_dimname
    att_array@lon                    = data->$native_dimname$

    att_array@hyai                   = data@ak/100000.D
    att_array@hybi                   = data@bk*1.D

    nlev                             = dimsizes(att_array@hyai) - 1

    att_array@hyam                   = new( nlev  , double)
    do k                             = 0,nlev-1
      att_array@hyam(k)              = ( att_array@hyai(k) + att_array@hyai(k+1) )*0.5D
    end do

    att_array@hybm                   = new( nlev  , double)
    do k                             = 0,nlev-1
      att_array@hybm(k)              = ( att_array@hybi(k) + att_array@hybi(k+1) )*0.5D
    end do

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    datei                            = data->time@begin_date
    houri                            = data->time@begin_time*1
    hourf                            = data->time@begin_time*1.
    dims                             = dimsizes(datei)
    if(max( houri - (houri/10000)*10000 ) .ne. 0) then
      print("ERROR:  'time@begin_time' does not contain an exact hour time stamp")
      print(data->time@begin_time)
      exit
    end if
    houri                            = (houri/10000 )*3600
    hourf                            = (hourf/10000.)*3600.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    delete(att_array@PREVOCA_u_v_var_names)
    att_array@PREVOCA_u_v_bugfix     = True
    att_array@PREVOCA_u_v_var_names  = (/ "uwnd", "vwnd" /)
    att_array@convert_theta_to_T     = True
    att_array@theta_var_names        = "theta"

    processed                        = True

  end if

  if(ftype .eq. "JRA_25") then

    ; Japanese 25-year Reanalysis Project 
    ; website         :  http://dss.ucar.edu/datasets/ds625.0/
    ; example datasets:  mss:/DSS/DS625.0/anl_mdl/1998/199808.anl_mdl.tar
    ;  (NOTE:  must un-tar monthly file before use)

    native_dimname                   = "lv_HYBL2"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$(dimsizes(data->$native_dimname$)-1:0)

    native_dimname                   = "g4_lat_0"
    att_array@lat_dimname            = native_dimname
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "g4_lon_1"
    att_array@lon_dimname            = native_dimname
    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)
    dim_lon                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    del_lon                          = (360.D/dim_lon)
    att_array@lon                    = new( (/dim_lon/) , double)
    att_array@lon                    = ispan(0,dim_lon-1,1)*del_lon

    att_array@hyai                   = (/ 0.0000000000000000D, 0.0008000000000000D, 0.0014600000000000D, 0.0025600000000000D, \
                                          0.0043300000000000D, 0.0071000000000000D, 0.0112000000000000D, 0.0170000000000000D, \
                                          0.0250000000000000D, 0.0352999990000000D, 0.0480999980000000D, 0.0626344300000000D, \
                                          0.0761050570000000D, 0.0883639980000000D, 0.0988765950000000D, 0.1072994920000000D, \
                                          0.1134470900000000D, 0.1172599790000000D, 0.1187773740000000D, 0.1181136090000000D, \
                                          0.1154385450000000D, 0.1109614490000000D, 0.1050948870000000D, 0.0981513060000000D, \
                                          0.0901928630000000D, 0.0814378200000000D, 0.0723235320000000D, 0.0630560150000000D, \
                                          0.0538116840000000D, 0.0447413480000000D, 0.0361580200000000D, 0.0281305770000000D, \
                                          0.0208627470000000D, 0.0144855000000000D, 0.0090642610000000D, 0.0046119540000000D, \
                                          0.0011056100000000D, 0.0000000000000000D, 0.0000000000000000D, 0.0000000000000000D, \
                                          0.0000000000000000D /)
    att_array@hybi                   = (/ 0.0000000000000000D, 0.0000000000000000D, 0.0000000000000000D, 0.0000000000000000D, \
                                          0.0000000000000000D, 0.0000000000000000D, 0.0000000000000000D, 0.0000000000000000D, \
                                          0.0000000000000000D, 0.0000000000000000D, 0.0000000000000000D, 0.0003660000000000D, \
                                          0.0038950000000000D, 0.0106360000000000D, 0.0211230000000000D, 0.0357010000000000D, \
                                          0.0545530000000000D, 0.0777400000000000D, 0.1052230000000000D, 0.1368860000000000D, \
                                          0.1725610000000000D, 0.2120390000000000D, 0.2539050000000000D, 0.2978490000000000D, \
                                          0.3448070000000000D, 0.3945620000000000D, 0.4456760000000000D, 0.4979440000000000D, \
                                          0.5511880000000000D, 0.6052590000000000D, 0.6588420000000000D, 0.7118690000000000D, \
                                          0.7631370000000000D, 0.8115140000000000D, 0.8559360000000000D, 0.8953880000000000D, \
                                          0.9288940000000000D, 0.9550000000000000D, 0.9750000000000000D, 0.9900000000000000D, \
                                          1.0000000000000000D /)

    nlev                             = dimsizes(att_array@hyai) - 1

    att_array@hyam                   = new( nlev  , double)
    do k                             = 0,nlev-1
      att_array@hyam(k)              = ( att_array@hyai(k) + att_array@hyai(k+1) )*0.5D
    end do

    att_array@hybm                   = new( nlev  , double)
    do k                             = 0,nlev-1
      att_array@hybm(k)              = ( att_array@hybi(k) + att_array@hybi(k+1) )*0.5D
    end do

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float
    ; --> Read time stamp from the PRES_GDS4_SFC variable directly

    date_int                         = chartoint( stringtochar( data->PRES_GDS4_SFC@initial_time ) ) - 48

    datei                            =         date_int(6 )*10000000
    datei                            = datei + date_int(7 )*1000000
    datei                            = datei + date_int(8 )*100000
    datei                            = datei + date_int(9 )*10000
    datei                            = datei + date_int(0 )*1000
    datei                            = datei + date_int(1 )*100
    datei                            = datei + date_int(3 )*10
    datei                            = datei + date_int(4 )

    houri                            =         date_int(12)*3600*10
    houri                            = houri + date_int(13)*3600
    houri                            = houri + date_int(15)*60*10
    houri                            = houri + date_int(16)*60
    hourf                            = houri * 1.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    att_array@Z_to_PHIS              = True
    att_array@Z_var_names            = "DIST_GDS4_SFC"

    processed                        = True

  end if

  if(ftype .eq. "Era_Interim_627.0_sc") then
 
    ; ERA-Interim data from ECMWF
    ; website         :  http://dss.ucar.edu/datasets/ds627.0/
    ; example datasets:  mss:/DSS/DS627.0/ei.oper.an.ml/200806/ei.oper.an.ml.regn128sc.2008061918

    native_dimname                   = "lv_HYBL2"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$

    native_dimname                   = "g4_lat_0"
    att_array@lat_dimname            = native_dimname
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "g4_lon_1"
    att_array@lon_dimname            = native_dimname
    att_array@lon                    = data->$native_dimname$

    att_array@hyam                   = data->lv_HYBL2_a
    att_array@hybm                   = data->lv_HYBL2_b
    att_array@hyai                   = data->lv_HYBL_i3_a
    att_array@hybi                   = data->lv_HYBL_i3_b

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float
    ; --> Read time stamp from the Z_GDS4_HYBL variable directly

    date_int                         = chartoint( stringtochar( data->Z_GDS4_HYBL@initial_time ) ) - 48

    datei                            =         date_int(6 )*10000000
    datei                            = datei + date_int(7 )*1000000
    datei                            = datei + date_int(8 )*100000
    datei                            = datei + date_int(9 )*10000
    datei                            = datei + date_int(0 )*1000
    datei                            = datei + date_int(1 )*100
    datei                            = datei + date_int(3 )*10
    datei                            = datei + date_int(4 )

    houri                            =         date_int(12)*3600*10
    houri                            = houri + date_int(13)*3600
    houri                            = houri + date_int(15)*60*10
    houri                            = houri + date_int(16)*60
    hourf                            = houri * 1.

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    att_array@logPs_to_Ps            = True
    att_array@logPs_var_names        = "LNSP_GDS4_HYBL"

    processed                        = True

  end if

  if(ftype .eq. "ERA40_ds117.2") then
 
    ; ERA-40 data from ECMWF
    ; website         :  http://dss.ucar.edu/datasets/ds117.2/
    ; example datasets:  mss:/DSS/U13003

    native_dimname                   = "initial_time0_hours"
    att_array@time_dimname           = native_dimname
    att_array@time                   = data->$native_dimname$

    native_dimname                   = "real_imaginary"
    att_array@real_imaginary_dimname = native_dimname

    native_dimname                   = "lv_HYBL4"
    att_array@lev_dimname            = native_dimname
    att_array@lev                    = data->$native_dimname$

    native_dimname                   = "g4_lat_7"
    att_array@lat_dimname            = native_dimname
    att_array@lat                    = data->$native_dimname$

    native_dimname                   = "g50_lat_2"
    att_array@lat_dimname_alt1       = native_dimname
    file_dim_names                   = getvardims     (data)
    file_dim_sizes                   = getfiledimsizes(data)
    dim_lat                          = file_dim_sizes( ind(file_dim_names .eq. native_dimname) )
    gau_info                         = gaus(dim_lat/2)
    att_array@lat_alt1               = new( (/dim_lat  /) , double)
    att_array@lat_alt1               = gau_info(:,0)

    native_dimname                   = "g4_lon_8"
    att_array@lon_dimname            = native_dimname
    att_array@lon                    = data->$native_dimname$

    native_dimname                   = "g50_lon_3"
    att_array@lon_dimname_alt1       = native_dimname
    dim_lon                          = dim_lat*2
    del_lon                          = (360.D/dim_lon)
    att_array@lon_alt1               = new( (/dim_lon/) , double)
    att_array@lon_alt1               = ispan(0,dim_lon-1,1)*del_lon

    att_array@hyam                   = data->lv_HYBL4_a
    att_array@hybm                   = data->lv_HYBL4_b
    att_array@hyai                   = data->lv_HYBL_i5_a
    att_array@hybi                   = data->lv_HYBL_i5_b

    ; Fill time stamp arrays based on time slices in file
    ;  datei = YYYYMMDD ; integer
    ;  houri = SSSSS    ; integer
    ;  hourf = SSSSS    ; float

    date_read                        = data->initial_time0_encoded
    datei                            = doubletointeger(date_read)/100
    datef                            = new( dimsizes(datei), double )
    datef                            = datei
    hourf                            = date_read - datef*100.
    hourf                            = hourf*3600.
    houri                            = doubletointeger(hourf)

    att_array@datei                  = datei
    att_array@houri                  = houri
    att_array@hourf                  = hourf

    ; Set processing flags

    att_array@logPs_to_Ps            = True
    att_array@logPs_var_names        = "LNSP_GDS50_HYBL"

    processed                        = True

  end if

  if(.not. processed) then
    print("Error in 'extract_analyses_info'")
    print("No recognized file type.  Add a block of code to extract dimension info from this new file type: '"+ftype+"'")
    exit
  end if

  return(att_array)

end

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
