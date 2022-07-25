class h0scam:

    def __init__(self, xp, dir ):
        self.case    = xp
        self.archdir = dir

    def curtain(self,fld):
        import numpy as np
        import xarray as xr
        import glob

        

        xp=self.case
        dir=self.archdir
        
        fl = sorted( glob.glob( dir +'/*cam.h0*') )
        nf = len( fl )

        ird=0
        for f in fl:
            print(f)
            try:
                a=xr.open_dataset( f )
                print("Successfully opened"+f)
                aa=a[fld] #.isel(time=0)
                nx=a.lon.size
                ny=a.lat.size
                nt=a.time.size


                if ('lon' in aa.dims) and ('lat' in aa.dims) and ('lev' not in aa.dims):
                    print(fld+' is a surface var ' )
                    varType = 'surface'
                elif ('lon' in aa.dims) and ('lat' in aa.dims) and ('lev' in aa.dims):
                    print(fld+' is a profile var ' )
                    varTtyp = 'profile'

                if (ird == 0):
                    #Cook up time array
                    timeData=a.time.data
                    interval=a.time[1]-a.time[0]
                    intersec=  ( interval.data.astype(int) / 10**9 ) # this mofo is in NANOseconds 
                    bigTime = xr.cftime_range( timeData[0]  , periods=nt*nf, freq='450S', calendar="noleap" )
                    if varType == 'surface':
                        dummy=np.zeros( nt*nf )
                        dummy=dummy.reshape( nt*nf ,1,1)
                        cu=xr.DataArray( dummy , coords=[bigTime,a.lat,a.lon], dims=['time','lat','lon'] )                    

                if varType == 'surface':
                    #cu.values[  ird*nt :(ird+1)*nt-1, 0, 0] = aa[:,0,0] # why the FUCK is this wrong?????!!!!!
                    cu.values[  ird*nt :(ird+1)*nt , 0, 0] = aa[:,0,0]  

                ird=ird+1

            except ValueError:
                print('******** VALUE ERROR *************')
                print('File \n'+f+'\n probably no good')
                if (ird == 0):
                    exit('First dataset not valid')
                else:
                    ird=ird+1
 
        return cu
