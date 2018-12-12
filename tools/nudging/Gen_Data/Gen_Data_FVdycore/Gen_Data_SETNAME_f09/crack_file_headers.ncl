;-----------------------------------------------------------------------
; Crack file headers
;-----------------------------------------------------------------------

begin

;------------------- USER-defined variables -------------------------

ncfpath    = (/"/DSS/DS627.0/ei.oper.an.ml/200806/"/)
mytmpdir   = (/"./"/)
;ncfname    = (/"ei.oper.an.ml.regn128sc.2008061918"/)
;ncfname    = (/"TMP/TMP_001.09/file0_ei.oper.an.ml.regn128sc.2008112900"/)
;ncfname    = (/"TMP/TMP_001.03/file0_ei.oper.an.ml.regn128sc.2009010100"/)
ncfname    = (/"TMP/TMP_001.08/file0_ei.oper.an.ml.regn128sc.2000021506"/)
;ftag       = (/""/)
ftag       = (/".grb"/)

;--------------------------------------------------------------------

print ("------------------------------------------------------------")
print ("msrcp -n mss:"+ncfpath+ncfname+" "+mytmpdir+ncfname)
print ("------------------------------------------------------------")
system("msrcp -n mss:"+ncfpath+ncfname+" "+mytmpdir+ncfname)
print ("------------------------------------------------------------")
print (" ")

;----------
; read file
;----------

data  = addfile(mytmpdir+ncfname+ftag, "r")
names = getfilevarnames(data)

;---------------------------------
; determine names of all variables
;---------------------------------

print(" variable name = "+names)

;------------------------
; Print global attributes
;------------------------

file_atts = getvaratts(data)
natts     = dimsizes(file_atts)
do i = 0,natts-1
  print(" ")
  print("-------------------------------------------------------------------------------")
  print(" ")
  print("file attribute = "+file_atts(i))
  print(data@$file_atts(i)$)
end do

;---------------------------------------
; determine dimensions of all variables
;---------------------------------------

nvars = dimsizes(names)
do i = 0,nvars-1
  field   = data->$names(i)$
  dims    = dimsizes(field)
  dimdims = dimsizes(dims)
  print(" ")
  print("-------------------------------------------------------------------------------")
  print(" ")
  print("variable name = "+names(i))
  printVarSummary(field)
  if(dimdims .eq. 1) then
    print("field = "+field)
  end if
  delete(dims)
  delete(field)
end do

exit
end
