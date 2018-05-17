;************************************************
; This routine converts a reanalysis file into an IOP 
; file that can be used with SCAM
;
; John Truesdale, Andrew Gettelman, December 2017
;
; To use, run CAM with the following namelist parameters:
;
;       fincl2 = 'U','V','T','Q','OMEGA','TTEND_TOT','PTTEND','TAQ','TS','PS'
;	fincl2lonlat = '305e_62n'
;	nhtfrq = 0,-3
;	avgflag_pertape = 'A','I'
; 
; Where '305' is the longitude (0-360) & '62n' is latitude of the desired point
; Averaging can be instantaneous or averaged.
; MFILT can be set to mfilt = 1,100 or another suitably high number
;
;***********************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  ;===================================================================
  ; SET PATHS AND CASE NAME
  ;===================================================================

  case="cam175_sd_scamiop2"
  pth=""

; can put ouput in same directory, or another (default same)

  ipth=pth+case
  opth=ipth

  ;===================================================================
  ; Set location (in variable names)
  ;===================================================================

  locx="_262e"
  locy="_37n"

; pick a name for the iop (will go in filename)

  iopname = "sgptst"

  ;===================================================================
  ; User should not need to modify below here.
  ;===================================================================

; convert variable names....

  loc=locx+locy
  utmp="U"+loc
  vtmp="V"+loc
  ttmp="T"+loc
  qtmp="Q"+loc
  otmp="OMEGA"+loc
  taqtmp="TAQ"+loc
  tttmp="TTEND_TOT"+loc
  pttmp="PTTEND"+loc
  pstmp="PS"+loc
  tstmp="TS"+loc

; Find files

  hist_stream="*.h1.*"  ;mean

  flist=systemfunc("ls "+ipth+hist_stream)

  infile = addfiles(flist,"r")

;Set up output file

  outfilename = case+"."+iopname+".iop.nc"
  outfile = addfile(opth+outfilename,"c")

  ;===================================================================
  ; make time an UNLIMITED dimension; recommended  for most applications
  ;===================================================================
  filedimdef(outfile,"time",-1,True) 

  u = infile[:]->$utmp$
  v = infile[:]->$vtmp$
  t = infile[:]->$ttmp$
  q = infile[:]->$qtmp$
  omega = infile[:]->$otmp$

  p0= 100000
  p0!0 = "ncl_scalar"
  hyam= infile[0]->hyam
  hybm= infile[0]->hybm
  hyai= infile[0]->hyai
  hybi= infile[0]->hybi
  ps= infile[:]->$pstmp$

;Convert lev to Pa

  lev= infile[0]->lev*100
  lev!0 = "lev"
  lev&lev = lev
  ilev= infile[0]->ilev*100
  ilev!0 = "ilev"
  ilev&ilev = ilev

  divq = infile[:]->$taqtmp$

; Calculate divT

  divq!2="lat"
  divq!3="lon"
; printVarSummary(divq)
; exit

  divT = divq
  divT@long_name="T horz + vert + fixer tendency (TTEND_TOT-PTTEND)"
  divT@mdims=1
  divT@units="K/s"

  ttendtot = infile[:]->$tttmp$
  pttend   = infile[:]->$pttmp$
  divT=ttendtot-pttend
; 
; create tsec
;
  atime=infile[:]->time
  datestrs = str_split(atime@units," ")
  yyyymmddstr=str_split(datestrs(2),"-")
  yyyy=tointeger(yyyymmddstr(0))
  mm=tointeger(yyyymmddstr(1))
  dd=tointeger(yyyymmddstr(2))
  nbdate=yyyy*10000+mm*100+dd
  nbdate!0 = "ncl_scalar"
  if (dimsizes(datestrs).le.3) then
     nbsec=0
     tounits="seconds "+datestrs(1)+" "+datestrs(2)
  else
     hhmmssstrs=str_split(datestrs(3),":")
     nbsec=tointeger(hhmmssstrs(0))*3600+tointeger(hhmmssstrs(1))*60+tointeger(hhmmssstrs(2))
     tounits="seconds "+datestrs(1)+" "+datestrs(2)+" "+datestrs(3)
  end if
  tsec=(cd_convert(atime,tounits))
  tsec!0="time"
  tsec&time=atime
  tsec@units=tounits
  tsec@long_name=tounits
;
; output new variables directly
;
  outfile->omega = omega
  outfile->tsec  = tsec
  outfile->bdate  = nbdate
  outfile->nbdate  = infile[0]->nbdate
  outfile->nbsec  = infile[0]->nbsec
  outfile->u  = u
  outfile->v  = v
  outfile->t  = t
  outfile->q = infile[:]->$qtmp$
  outfile->Ps = infile[:]->$pstmp$
  outfile->Tsair = infile[:]->$tstmp$
  outfile->lev = lev
  outfile->ilev = ilev
  outfile->p0 = p0
  outfile->hyam= infile[0]->hyam
  outfile->hybm= infile[0]->hybm
  outfile->hyai= infile[0]->hyai
  outfile->hybi= infile[0]->hybi
  
;Rename dimensions in output file.

  system("ncrename -d  lat_physgrid"+locy+",lat "+opth+outfilename)
  system("ncrename -d  lon_physgrid"+locx+",lon "+opth+outfilename)
  system("ncrename -v  lat_physgrid"+locy+",lat "+opth+outfilename)
  system("ncrename -v  lon_physgrid"+locx+",lon "+opth+outfilename)

;add global atribute

  system("ncatted -a CAM_GENERATED_FORCING,global,c,c,'yes' "+opth+outfilename)

;now write the fields with correct dimensions

  delete(outfile)
  outfile2 = addfile(opth+outfilename,"w")

  outfile2->Q_dten = divq
  outfile2->divT3d = divT

end
