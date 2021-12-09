undef("Create_field_Master_List")

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------

function Create_field_Master_List(spectral_trunc_setting:integer, field:string )

;
;********************************************************************************************;
;                                                                                            ;
; Create_field_Master_List                                                                   ;
;    spectral_trunc_setting : integer; Spectral truncation (Eulerian only)                   ;
;    "field" is a master variable carrying meta-data needed for file                         ;
;    processing                                                                              ;
;                                                                                            ;
;********************************************************************************************;

local  field_dim, ifield, Master_List

begin

;-------------
; Declarations
;-------------

  field_dim                   = 100
  Master_List                 = new( (/field_dim/), string )   ; CAM field name
  Master_List@dimensions      = new( (/field_dim/), string )   ; "2D" or "3D" designation
  Master_List@horz_interp     = new( (/field_dim/), string )   ; flag to indicate which type of horizontal interpolation to use
  Master_List@vert_interp     = new( (/field_dim/), string )   ; flag to indicate which type of vertical   interpolation to use
  Master_List@spec_trunc      = new( (/field_dim/), integer)   ; Spectral truncation to use (if Eulerian dycore)
  Master_List@bin_factor      = new( (/field_dim/), double )   ; bin factor to use if horizontal interpolation is "binning"
                                                               ; (typical value:  1.)
  Master_List@loutput         = new( (/field_dim/), logical)   ; Flag to output field (if False, field is used only for internal processing
                                                               ; only)
  Master_List@long_name       = new( (/field_dim/), string )   ; long name of field
  Master_List@units           = new( (/field_dim/), string )   ; units of field
  Master_List@min             = new( (/field_dim/), float  )   ; minimum possible value for this field (checked against final interpolated field)
  Master_List@max             = new( (/field_dim/), float  )   ; maximum possible value for this field (checked against final interpolated field)
  Master_List@post_process    = new( (/field_dim/), string )   ; post-processing flags (after interpolation):
                                                               ; "ps_adjust", "q_adjust", "cloud_fraction_adjust", "ts_adjust",
                                                               ; "cwat_adjust"

  do ifield = 0,field_dim-1
    Master_List@loutput     (ifield) = False
    Master_List@post_process(ifield) = ""
  end do

  post_process_flags = (/ "ps_adjust", "q_adjust", "cloud_fraction_adjust", "ts_adjust", "cwat_adjust", "" /)

  ;---------------------------------------------------------------------
  ;
  ; Interpolation flags
  ; 
  ;  Options for horizontal:   "no_interp"
  ;                            "cubic"    ; recommended only for non-tracer species
  ;                            "cubic_sp" ; should only be used for tracers and water products like Q, CLDLIQ, CLDICE, CLOUD, etc.
  ;                            "spectral" ; only for non-tracer species in the Eulerian dycore
  ;                            "binning"  ; "conservative remapping" --> recommended for all horizontal interpolations
  ;
  ;  Options for vertical:     "no_interp"
  ;                            "linear"
  ;                            "linear_logp"
  ;                            "linear_opt1"      ; customized for water products like Q, CLDLIQ, CLDICE, CLOUD, etc. (recommended)
  ;                            "linear_logp_opt1" ; customized for water products like Q, CLDLIQ, CLDICE, CLOUD, etc.
  ;                            "linear_opt2"      ; customized for U, V, W wind fields
  ;                            "linear_logp_opt2" ; customized for U, V, W wind fields (recommended)
  ;                            "quad_logp_opt1"   ; Can only be used on temperature fields ("T", "TCWAT", etc.) (recommended)
  ;                            "quad_opt1"        ; Can only be used on temperature fields ("T", "TCWAT", etc.)
  ;
  ;---------------------------------------------------------------------

  ifield = -1

  ifield                             = ifield + 1
  Master_List(ifield)                = "OMEGA"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Vertical wind"
  Master_List@units       (ifield)   = "Pa/s"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "U"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Zonal wind"
  Master_List@units       (ifield)   = "m/s"
  Master_List@min         (ifield)   = -200.
  Master_List@max         (ifield)   =  200.
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "US"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Zonal wind, staggered"
  Master_List@units       (ifield)   = "m/s"
  Master_List@min         (ifield)   = -200.
  Master_List@max         (ifield)   =  200.
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "VORT"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Relative Vorticity"
  Master_List@units       (ifield)   = "1/s"
  Master_List@min         (ifield)   = -1.e-1
  Master_List@max         (ifield)   =  1.e-1
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "V"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Meridional wind"
  Master_List@units       (ifield)   = "m/s"
  Master_List@min         (ifield)   = -150.
  Master_List@max         (ifield)   =  150.
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "VS"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Meridional wind, staggered"
  Master_List@units       (ifield)   = "m/s"
  Master_List@min         (ifield)   = -150.
  Master_List@max         (ifield)   =  150.
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "DIV"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Divergence"
  Master_List@units       (ifield)   = "1/s"
  Master_List@min         (ifield)   = -1.e-1
  Master_List@max         (ifield)   =  1.e-1
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "T"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "quad_logp_opt1"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  160.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "Q"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Specific humidity"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   =  0.
  Master_List@max         (ifield)   =  0.1
  Master_List@post_process(ifield)   = "q_adjust"              ; Adjust Q to eliminate RH < 0. or > 1.

  ifield                             = ifield + 1
  Master_List(ifield)                = "PS"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = spectral_trunc_setting
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Surface pressure"
  Master_List@units       (ifield)   = "Pa"
  Master_List@min         (ifield)   =   30000.
  Master_List@max         (ifield)   =  109000.
  Master_List@post_process(ifield)   = "ps_adjust"             ; Adjust PS for differences in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "PHIS"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "no_interp"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = field@output_phis
  Master_List@long_name   (ifield)   = "Surface geopotential"
  Master_List@units       (ifield)   = "m2/s2"
  Master_List@min         (ifield)   =  -2.e+3
  Master_List@max         (ifield)   =   1.e+5
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "SST_cpl"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "SST"
  Master_List@units       (ifield)   = "deg_C"
  Master_List@min         (ifield)   = -2.
  Master_List@max         (ifield)   = 40.
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "ice_cov"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Sea-ice concentration"
  Master_List@units       (ifield)   = "fraction"
  Master_List@min         (ifield)   =  0.
  Master_List@max         (ifield)   =  1.0001
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "SGH"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "no_interp"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "St dev of orography"
  Master_List@units       (ifield)   = "m"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "LANDFRAC"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "no_interp"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Fraction of sfc area covered by land"
  Master_List@units       (ifield)   = "fraction"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "LANDM_COSLAT"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "no_interp"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Landm"
  Master_List@units       (ifield)   = "unitless"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "TREFHT"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "2m temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = "ts_adjust"             ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "TS"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Surface temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = "ts_adjust"             ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "TSICE"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Ice temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = ""                      ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "SICTHK"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Sea ice thickness"
  Master_List@units       (ifield)   = "m"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "TS1"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "TS1      subsoil temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = ""                      ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "TS2"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "TS2      subsoil temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = ""                      ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "TS3"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "TS3      subsoil temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = ""                      ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "TS4"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "TS4      subsoil temperature"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = ""                      ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "TSOCN"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Ocean tempertare"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   = -40.
  Master_List@max         (ifield)   =  60.
  Master_List@post_process(ifield)   = ""                      ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "SNOWHICE"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Water equivalent snow depth"
  Master_List@units       (ifield)   = "m"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "ICEFRAC"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Fraction of sfc area covered by sea-ice"
  Master_List@units       (ifield)   = "fraction"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "PBLH"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "PBL height"
  Master_List@units       (ifield)   = "m"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "TPERT"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Perturbation temperature (eddies in PBL)"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "QPERT"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Perturbation specific humidity (eddies in PBL)"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "CLOUD"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Cloud fraction"
  Master_List@units       (ifield)   = "fraction"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = "cloud_fraction_adjust" ; Adjust CLOUD between 0. and 1.

  ifield                             = ifield + 1
  Master_List(ifield)                = "CWAT"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Total Grid box averaged Condensate Amount (liquid + ice)"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = "cwat_adjust"           ; Adjust Cloud Water species to eliminate negative values

  ifield                             = ifield + 1
  Master_List(ifield)                = "QCWAT"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "q associated with cloud water"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = "cwat_adjust"           ; Adjust Cloud Water species to eliminate negative values

  ifield                             = ifield + 1
  Master_List(ifield)                = "TCWAT"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "quad_logp_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "T associated with cloud water"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "LCWAT"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Total Grid box averaged Condensate Amount (liquid + ice)"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = "cwat_adjust"           ; Adjust Cloud Water species to eliminate negative values

  ifield                             = ifield + 1
  Master_List(ifield)                = "CLDLIQ"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Grid box averaged liquid condensate amount"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = "cwat_adjust"           ; Adjust Cloud Water species to eliminate negative values

  ifield                             = ifield + 1
  Master_List(ifield)                = "CLDICE"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Grid box averaged ice condensate amount"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = "cwat_adjust"           ; Adjust Cloud Water species to eliminate negative values

  ifield                             = ifield + 1
  Master_List(ifield)                = "TSICERAD"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Radiatively equivalent Surface T over seaice"
  Master_List@units       (ifield)   = "K"
  Master_List@min         (ifield)   =  200.
  Master_List@max         (ifield)   =  360.
  Master_List@post_process(ifield)   = ""                      ; Adjust TS fields for change in PHIS per ECMWF algorithm

  ifield                             = ifield + 1
  Master_List(ifield)                = "QBOT"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Lowest model level specific humidity"
  Master_List@units       (ifield)   = "kg/kg"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "RELHUM"
  Master_List@dimensions  (ifield)   = "3D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_opt1"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Relative Humidity"
  Master_List@units       (ifield)   = "percent"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "FLUT"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Upwelling longwave flux at top of model"
  Master_List@units       (ifield)   = "W/m2"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "PRECT"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = "Total precipitation rate "
  Master_List@units       (ifield)   = "m/s"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

;----------------------------------------------------
; If not in list above, assign field these attributes
;----------------------------------------------------

  ifield                             = ifield + 1
  Master_List(ifield)                = "generic"
  Master_List@dimensions  (ifield)   = ""
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "linear_logp_opt2"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = True
  Master_List@long_name   (ifield)   = ""
  Master_List@units       (ifield)   = ""
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ;---------------------------------------------------------------------------------
  ; Additional input fields which will not be output themselves ("loutput = False"), 
  ; but which will be used for processing other output fields.
  ;---------------------------------------------------------------------------------

  ifield                             = ifield + 1
  Master_List(ifield)                = "PHIS_input"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "binning"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = False
  Master_List@long_name   (ifield)   = "N/A"
  Master_List@units       (ifield)   = "N/A"
  Master_List@min         (ifield)   = -2.e+3
  Master_List@max         (ifield)   =  1.e+5
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "LANDFRAC_analysis"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "no_interp"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = False
  Master_List@long_name   (ifield)   = "N/A"
  Master_List@units       (ifield)   = "N/A"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  ifield                             = ifield + 1
  Master_List(ifield)                = "ICEFRAC_analysis"
  Master_List@dimensions  (ifield)   = "2D"
  Master_List@horz_interp (ifield)   = "no_interp"
  Master_List@vert_interp (ifield)   = "no_interp"
  Master_List@spec_trunc  (ifield)   = -1
  Master_List@bin_factor  (ifield)   = -1
  Master_List@loutput     (ifield)   = False
  Master_List@long_name   (ifield)   = "N/A"
  Master_List@units       (ifield)   = "N/A"
  Master_List@min         (ifield)   = -1.e30
  Master_List@max         (ifield)   =  1.e30
  Master_List@post_process(ifield)   = ""

  do ifield = 0,field_dim-1
    if( .not. any( Master_List@post_process(ifield) .eq. post_process_flags) ) then
      print("Error in 'Create_field_Master_List'")
      print("  Post-processing flag for field '"+Master_List(ifield)+"' is '"+Master_List@post_process(ifield)+"'")
      print("  Valid post-processing options are: '"+post_process_flags+"'")
      exit
    end if
    if( field@dycore .ne. "eul") then
      Master_List@spec_trunc(ifield) = -1
    end if
  end do

  delete(Master_List@_FillValue)

  return (Master_List)

end

; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
; ----------------------------------------------------------------------
