;***************************************************************
; NCL script to copy or remap all variables of an FV netcdf file
; to a rectilinear Gaus grid.
; remapfv2eul.ncl
; John Truesdale, May 2018
;***************************************************************
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/esmf/ESMF_regridding.ncl"

;*************************************************************
;                        EX. USAGE
;*************************************************************
;env SRCFILE=/fs/cgd/csm/inputdata/atm/cam/inic/fv/cami-mam3_0000-01-01_0.9x1.25_L32_c141031.nc DSTLATLON=/fs/cgd/csm/inputdata/atm/cam/inic/gaus/cami_0000-01-01_128x256_L26_c030918.nc ncl < remapfv2eul.ncl
;************************************************************

begin
;-----The FV source files to be converted
srcFileName   = getenv("SRCFILE")

;-----The destination lat lon template file - any file that has lat lon coordinates
latlontemp    = getenv("DSTLATLON")

outpth="./"

;-----  shouldnt have to modify below here -------------

verbose=True   ;False
debug=False   ;False

;-----interplation method
map_method=getenv("MAPMETHOD")
if (ismissing(map_method)) then
delete(map_method)
map_method    = "bilinear"
end if

;-----directory to store weights for regridding
WGT_dir=getenv("MAPPATH")
if (ismissing(WGT_dir)) then
delete(WGT_dir)
WGT_dir       = systemfunc("pwd")
end if

src_file=addfile(srcFileName,"r")
dstlatlon_file=addfile(latlontemp,"r")

; Determine ingrid
srcdnames   = getvardims(src_file)       ; get level info from source file
srcdsizes   = getfiledimsizes(src_file)
srcnlat       = srcdsizes(ind(srcdnames.eq."lat" ));
srcnlon       = srcdsizes(ind(srcdnames.eq."lon" ));
ingrid="FV_"+tostring(srcnlat)+"x"+tostring(srcnlon)
haveslat=False
if (any(srcdnames.eq."slat")) then
   haveslat=True
   srcnslat       = srcdsizes(ind(srcdnames.eq."slat" ));
   ingrid_us="FV_"+tostring(srcnslat)+"x"+tostring(srcnlon)
else
   ingrid_us="SLAT_dimension_not_found"
end if
haveslon=False
if (any(srcdnames.eq."slon")) then
   haveslon=True
   srcnslon       = srcdsizes(ind(srcdnames.eq."slon" ));
   ingrid_vs="FV_"+tostring(srcnlat)+"x"+tostring(srcnslon)
else
   ingrid_vs="SLON_dimension_not_found"
end if

; Determine outgrid
dstdnames   = getvardims(dstlatlon_file)       ; get level info from source file
dstdsizes   = getfiledimsizes(dstlatlon_file)
dstnlat       = dstdsizes(ind(dstdnames.eq."lat" ));
dstnlon       = dstdsizes(ind(dstdnames.eq."lon" ));
outgrid="Gaus_"+tostring(dstnlat)+"x"+tostring(dstnlon)

srcFile=systemfunc("basename "+srcFileName)
suffix = get_file_suffix(srcFile,0)
srcbase= suffix@fBase
dstFileName   = outpth+srcbase+".regrid."+outgrid+".nc"
currdate=systemfunc("date +%y%m%d")

    ; 0. Set the source/destination file names, 
    ;    open the source file for reading, 
    ;    create a destination file for regridded data.
    ;------------------------------------------------

    print(" ")
    print("Regridding: ")
    print("SRC File:"+srcFileName)
    print("DST File:"+dstFileName)
    if(fileexists(dstFileName)) then
      system("rm "+dstFileName)
    end if

    setfileoption("nc","Format","LargeFile")
    dst_file=addfile(dstFileName,"c")


;
    ; 1. Generate a description file for the source grid (EUL).
    ;---------------------------------------------------
    srcGridName=WGT_dir+"/"+ingrid+"_SCRIP_desc.nc"
    srcGridName_us=WGT_dir+"/"+ingrid_us+"_US_SCRIP_desc.nc"
    srcGridName_vs=WGT_dir+"/"+ingrid_vs+"_VS_SCRIP_desc.nc"
    Opt               =True
    Opt@ForceOverWrite=True
    Opt@Title         ="FV Grid"  
    Opt@Debug         =False ; True
    Opt@PrintTimings  =False ; True
    if(isfilepresent(srcGridName)) then
      if (verbose) then print("Found srcGrid description "+srcGridName) end if
    else
      if (verbose) then print("Creating srcGrid description "+srcGridName) end if
      rectilinear_to_SCRIP(srcGridName,src_file->lat,src_file->lon,Opt)
    end if
    if (haveslat) then
       if(isfilepresent(srcGridName_us)) then
          if (verbose) then print("Found srcGrid description "+srcGridName_us) end if
       else
          if (verbose) then print("Creating srcGrid description "+srcGridName_us) end if
          rectilinear_to_SCRIP(srcGridName_us,src_file->slat,src_file->lon,Opt)
       end if
    end if
    if (haveslon) then
       if(isfilepresent(srcGridName_vs)) then
          if (verbose) then print("Found srcGrid description "+srcGridName_vs) end if
       else
         if (verbose) then print("Creating srcGrid description "+srcGridName_vs) end if
         rectilinear_to_SCRIP(srcGridName_vs,src_file->lat,src_file->slon,Opt)
       end if
    end if
    delete(Opt)

    ;
    ; 2. Generate a description file for the destination grid (EUL).
    ;-----------------------------------------------------
    dstGridName = WGT_dir+"/"+outgrid+"_SCRIP_desc.nc"

    if(isfilepresent(dstGridName)) then
      if (verbose) then print("Found dstGrid description "+dstGridName) end if 
    else
      if (verbose) then print("Creating dstGrid description "+dstGridName) end if
      Opt               =True
      Opt@ForceOverWrite=True
      Opt@Debug         =False ; True
      Opt@PrintTimings  =False ; True
      rectilinear_to_SCRIP(dstGridName,dstlatlon_file->lat,dstlatlon_file->lon,Opt)
      delete(Opt)
    end if

    ;
    ; 3. Generate the weights file, using the source and 
    ;    destination files created in #1 and #2.
    ;-----------------------------------------------------
    wgtFileName = WGT_dir+"/map_"+ingrid+"_to_"+outgrid+"_"+map_method+".nc"
    wgtFileName_us2u = WGT_dir+"/map_"+ingrid_us+"_US_to_"+outgrid+"_U_"+map_method+".nc"
    wgtFileName_vs2v = WGT_dir+"/map_"+ingrid_vs+"_VS_to_"+outgrid+"_V_"+map_method+".nc"

    Opt               = True
    Opt@InterpMethod  =map_method ;"bilinear"  "patch", "conserve"
    Opt@ForceOverWrite=True
    Opt@OverWrite=True
    Opt@SrcESMF       =False
    Opt@DstESMF       =False   
    Opt@Debug         =False ; True
    Opt@PrintTimings  =False ; True

    if(isfilepresent(wgtFileName)) then
       if (verbose) then print("Found WeightFiles "+wgtFileName) end if
    else
       if (verbose) then print("Creating WeightFiles "+wgtFileName) end if
       ESMF_regrid_gen_weights(srcGridName,dstGridName,wgtFileName,Opt)
    end if
    if (haveslat) then
    if(isfilepresent(wgtFileName_us2u)) then
      if (verbose) then print("Found WeightFiles "+wgtFileName_us2u) end if
    else
       if (verbose) then print("Creating WeightFile "+wgtFileName_us2u) end if
       ESMF_regrid_gen_weights(srcGridName_us,dstGridName,wgtFileName_us2u,Opt)
    end if
    end if
    if (haveslon) then
    if(isfilepresent(wgtFileName_vs2v)) then
      if (verbose) then print("Found WeightFiles "+wgtFileName_vs2v) end if
    else
       if (verbose) then print("Creating WeightFile "+wgtFileName_vs2v) end if
       ESMF_regrid_gen_weights(srcGridName_vs,dstGridName,wgtFileName_vs2v,Opt)
    end if
    end if
    delete(Opt)


    ;--- Specify a list of 1D variables on the *source file* that should NOT be copied
    var_in_exclude = (/"lat", "lon", "w_stag", "gw", "slat","slon", "area", "date_written", "time_written"/)  

    ;--- Specify a list of variables on the source file that should be directly copied 
    var_in_copy = (/"time_bnds"/)  

    ;--- Specify a list of variables to be regridded
    var_out        = "All_Variables"     ; to be regridded

    ;---Read from the weight file the method used to derive the remap weights
    wgt_file       = addfile(wgtFileName, "r")
    dst_grid_dims  = wgt_file->dst_grid_dims
    dst_mlon       = dst_grid_dims(0)
    dst_nlat       = dst_grid_dims(1)
    dst_lat        = wgt_file->yc_b(::dst_mlon)
    dst_lon        = wgt_file->xc_b(:dst_mlon-1)

    ;---Use the destination (EUL) grid info on the weight file to create lat/lon
    lat       = dst_lat                  ; get from weight file
    lat@long_name = "latitude"
    lat!0     = "lat"
    lat@units = "degrees_north"
    lat&lat   =  lat
    nlat      = dimsizes(lat)            ; same as dst_nlat

    lon       = dst_lon
    lon@long_name = "longitude"
    lon!0     = "lon"
    lon@units = "degrees_east"
    lon&lon   =  lon
    mlon      = dimsizes(lon)            ; same as dst_mlon


;---Get all variables on the FV file 
    var_in   = getfilevarnames( src_file )
    nvar_in  = dimsizes(var_in)
    Opt_RGRD  = True

    if (verbose) then print("creating new file by copying variable meta data from "+srcFile) end if

;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;                    PREDEFINE MODE
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    setfileoption(dst_file,"DefineMode",True)

;===================================================================
; create global attributes of the netCDF file containing regridded data
;===================================================================
    global    = True
    global@separator2= "---------------------------------------"
    copy_VarAtts(src_file, global)
    global@separator1= "------- SOURCE FILE ATTRIBUTES --------"

    if (isatt(src_file,"title")) then
        global@TITLE = "REMAPPED: "+src_file@title
    end if

    global@remap         = "NCL: ESMF_regrid_with_weights"
    global@map_method    = map_method
    global@creation_date = systemfunc("date")
    
    fileattdef( dst_file, global )        ; copy file attributes to netCDF file

;===================================================================
; predefine the coordinate dimension names and their sizes
;===================================================================
    dNames   = getvardims(src_file)       ; get level info from source file
    dSizes   = getfiledimsizes(src_file)
    ia       = ind(dNames.eq."lev" )   ; # of levels 
    ib       = ind(dNames.eq."ilev")
    klev     = dSizes(ia)
    klevi    = dSizes(ib)
                                       ; standard CAM dimensions
    dimNames = (/"time",    "lat",   "lon",    "lev",  "ilev", "nbnd", "chars", "scalar"/)  
    dimSizes = (/ -1   ,    nlat ,    mlon ,    klev,   klevi ,   2   ,    8   ,   1     /) 
    dimUnlim = (/ True ,    False,    False,    False , False ,  False,False   ,False    /)   
    filedimdef(dst_file,dimNames,dimSizes,dimUnlim)

;--- The following are explicitly added because they are 'special'

   filevardef(dst_file, "lat", typeof(lat), getvardims(lat)) 
   filevarattdef(dst_file,"lat" ,lat)           

   filevardef(dst_file, "lon", typeof(lon), getvardims(lon)) 
   filevarattdef(dst_file,"lon" ,lon)           

;   filevardef(dst_file, "gw", typeof(gw), getvardims(gw)) 
;   filevarattdef(dst_file,"gw" ,gw)           

;--- Loop over all variables and predfine meta data
;    do nv=0,nvar_in-1
    do nv=0,nvar_in-1
       if (.not.any(var_in(nv).eq.var_in_exclude)) then 
        if(var_out(0).eq."All_Variables"          .or. \
           any(var_in(nv).eq.var_out)             .or. \
           any(var_in(nv).eq.var_in_copy)       ) then

           rank_in = dimsizes(getfilevardimsizes(src_file, var_in(nv)))
           if (debug) then print(rank_in+var_in(nv)) end if
           if (rank_in .eq.1 .or. any(var_in(nv).eq.var_in_copy) ) then
              if (debug) then print(var_in(nv)+" rank 1") end if
              filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
                                             , getfilevardims(src_file,var_in(nv))  )
           end if
           if (rank_in .eq.2 .and. .not.any(var_in(nv).eq.var_in_copy)) then
              if (debug) then print(var_in(nv)+" rank 2 lat lon") end if
              filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
                                             , (/"lat", "lon" /)        )
           end if
           if (rank_in .eq.3 .and. .not.any(var_in(nv).eq.var_in_copy)) then
              if (debug) then print(var_in(nv)+" rank 3 time lat lon") end if
              filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
                                             , (/"time","lat", "lon" /)        )
           end if
           if (rank_in .eq.4.) then
              vdims = getfilevardims(src_file, var_in(nv))
              if (debug) then print(var_in(nv)+" rank 4 time lev lat lon") end if
              if (var_in(nv).eq."US") then
                 if (debug) then print("U rank 4 time lev lat lon") end if
                 filevardef(dst_file, "U", getfilevartypes(src_file,var_in(nv)) \
                     , (/ "time", vdims(1),  "lat", "lon" /)  )
              else if (var_in(nv).eq."VS") then
                 if (debug) then print("V rank 4 time lev lat lon") end if
                 filevardef(dst_file, "V", getfilevartypes(src_file,var_in(nv)) \
                     , (/ "time", vdims(1),  "lat", "lon" /)  )
              else
                 if (debug) then print(var_in(nv)+" rank 4 time lev lat lon") end if
                 filevardef(dst_file, var_in(nv), getfilevartypes(src_file,var_in(nv)) \
                     , (/ "time", vdims(1),  "lat", "lon" /)  )
	      end if
	      end if
	      delete(vdims)
	   end if
           dumAtts = new( 1, getfilevartypes(src_file,var_in(nv)))
           varAtts = getfilevaratts(src_file, var_in(nv))
           if (.not.ismissing(varAtts(0))) then
               nAtts = dimsizes(varAtts)
               do na=0,nAtts-1
                  dumAtts@$varAtts(na)$ = src_file->$var_in(nv)$@$varAtts(na)$
               end do
               if (var_in(nv).eq."US") then
                  filevarattdef(dst_file, "U" , dumAtts)
               else if (var_in(nv).eq."VS") then
                  filevarattdef(dst_file, "V" , dumAtts)
               else
                  filevarattdef(dst_file, var_in(nv) , dumAtts)
               end if
               end if
           end if
           delete([/varAtts, dumAtts/])   ; may change size next iteration
        end if
       end if
    end do       ; nv

;===================================================================
; explicitly exit file definition mode. **NOT REQUIRED in NCL**
;===================================================================
    setfileoption(dst_file,"DefineMode",False)
        
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;        Write the basic and regridded data values to the predefined structures
;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    if (verbose) then print("regridding all variables ... ") end if

    dst_file->lat = (/ lat /)
    dst_file->lon = (/ lon /)
;    dst_file->gw  = (/ gw  /)

    do nv=0,nvar_in-1
       
       if (.not.any(var_in(nv).eq.var_in_exclude)) then 
        if(var_out(0).eq."All_Variables"          .or. \
           any(var_in(nv).eq.var_out)             .or. \
           any(var_in(nv).eq.var_in_copy)       ) then

           rank_in = dimsizes(getfilevardimsizes(src_file, var_in(nv)))

           if (rank_in .eq.1 .or. any(var_in(nv).eq.var_in_copy) ) then
               dst_file->$var_in(nv)$ = (/ src_file->$var_in(nv)$ /)
           end if
           if (rank_in .ge.2 .and. .not.any(var_in(nv).eq.var_in_copy)) then
	                  if (debug) then print(var_in(nv)) end if
               if (var_in(nv).eq."VS") then
                  dst_file->V = (/ ESMF_regrid_with_weights(src_file->$var_in(nv)$,wgtFileName_vs2v,Opt_RGRD) /)
	       else if (var_in(nv).eq."US") then
		  dst_file->U = (/ ESMF_regrid_with_weights(src_file->$var_in(nv)$,wgtFileName_us2u,Opt_RGRD) /)
               else
                  dst_file->$var_in(nv)$ = (/ ESMF_regrid_with_weights(src_file->$var_in(nv)$,wgtFileName,Opt_RGRD) /)
               end if
               end if
           end if
        end if
       end if
       if (nv.ne.nvar_in-1) then system("echo -n .") else system("echo .") end if 
    end do       ; nv
    if (verbose) then print("Finished!!") end if
end
