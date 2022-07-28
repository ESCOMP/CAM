class h0scam:

    def __init__(self, xp, dir ):
        self.case    = xp
        self.archdir = dir

    def curtain(self,fld):
        import numpy as np
        import xarray as xr
        import glob
        import txtutil as tx

        

        xp=self.case
        dir=self.archdir
        #CplFreq=self.atm_ncpl
        #freqw = int( 86400./CplFreq )
        #freqs = str( freqw )
        #freqs = freqs.strip()+'S'
        #print('write interval='+freqs)
        

        fili=dir+'/atm_in'
        nhtfrq = tx.nmlread( fili, 'nhtfrq' )
        nhtfrq = nhtfrq.split(',')
        h0frq  = int( nhtfrq[0] )
        
        fili=dir+'/nuopc.runconfig'
        freqc = int( tx.nmlread( fili, 'atm_cpl_dt' ) )

        if (h0frq > 0):
            freqw=freqc*h0frq

        freqs = str( freqw )
        freqs = freqs.strip()+'S'
        print('write interval='+freqs)
 
        fl = sorted( glob.glob( dir +'/*cam.h0*') )
        nf = len( fl )

        #fl=fl[0:10]
        ird=0
        for f in fl:
            print(f)
            try:
                a=xr.open_dataset( f )
                print("Successfully opened"+f)
                nx=a.lon.size
                ny=a.lat.size
                nl=a.lev.size
                nli=a.ilev.size
                nt=a.time.size
  
                aa=a[fld] #.isel(time=0)
                print(aa.dims)

                if ('lon' in aa.dims) and ('lat' in aa.dims) and ('lev' not in aa.dims) and ('ilev' not in aa.dims):
                    print(fld+' is a surface var ' )
                    varType = 'surface'
                if ('lon' in aa.dims) and ('lat' in aa.dims) and ('lev' in aa.dims):
                    print(fld+' is a profile var ' )
                    varType = 'profile'
                if ('lon' in aa.dims) and ('lat' in aa.dims) and ('ilev' in aa.dims):
                    print(fld+' is an Iprofile var ' )
                    varType = 'iprofile'

                if (ird == 0):
                    #Cook up time array
                    timeData=a.time.data
                    #interval=a.time[1]-a.time[0]
                    #intersec=  ( interval.data.astype(int) / 10**9 ) # this is in NANOseconds 
                    bigTime = xr.cftime_range( timeData[0]  , periods=nt*nf, freq=freqs , calendar="noleap" )
                    if varType == 'surface':
                        dummy=np.zeros( nt*nf )
                        dummy=dummy.reshape( nt*nf ,1,1)
                        cu=xr.DataArray( dummy , coords=[bigTime,a.lat,a.lon], dims=['time','lat','lon'] , name=fld )                    
                    if varType == 'profile':
                        dummy=np.zeros( nt*nf*nl )
                        dummy=dummy.reshape( nt*nf,nl ,1,1)
                        cu=xr.DataArray( dummy , coords=[bigTime,a.lev,a.lat,a.lon], dims=['time','lev','lat','lon'], name=fld  )                    
                    if varType == 'iprofile':
                        dummy=np.zeros( nt*nf*nli )
                        dummy=dummy.reshape( nt*nf,nli ,1,1)
                        cu=xr.DataArray( dummy , coords=[bigTime,a.ilev,a.lat,a.lon], dims=['time','ilev','lat','lon'], name=fld  )                    
                        
                    print('Prepped XARRAY for data with time')

                if varType == 'surface':
                    #cu.values[  ird*nt :(ird+1)*nt-1, 0, 0] = aa[:,0,0] # why the heck is this wrong?????!!!!!
                    cu.values[  ird*nt :(ird+1)*nt , 0, 0] = aa[:,0,0]  
                if varType == 'profile' or varType == 'iprofile':
                    cu.values[  ird*nt :(ird+1)*nt , : , 0, 0] = aa[:,:,0,0]  

                ird=ird+1

            except ValueError:
                print('******** VALUE ERROR *************')
                print('File \n'+f+'\n probably no good')
                if (ird == 0):
                    exit('First dataset not valid')
                else:
                    ird=ird+1

            cu.to_netcdf('/project/amp/juliob/scam/'+xp+'_'+fld+'.nc')

        return cu


    def ncmerge(self):
        import numpy as np
        import xarray as xr
        import glob
        import txtutil as tx

        

        xp=self.case
        dir=self.archdir
             
        fl = sorted( glob.glob( dir +'/*cam.h0*') )
        nf = len( fl )

        flt = fl[0:10]
 
        #ds=xr.open_mfdataset( flt , concat_dim='time')
        #print("merged nc files")
        #ds.to_netcdf(path = dir +'/merged_h0.nc')
        #exit()

        ird=0
        for f in flt:
            print(f)
            try:
                a=xr.open_dataset( f )
                print("Successfully opened"+f)

                if (ird == 0):
                    b=a
                elif (ird>0):
                    b.merge(a,join='inner' )
                    b.to_netcdf( 'testmerge.nc' )

                ird=ird+1

            except ValueError:
                print('******** VALUE ERROR *************')
                print('File \n'+f+'\n probably no good')
                if (ird == 0):
                    exit('First dataset not valid')
                else:
                    ird=ird+1
 
        return
