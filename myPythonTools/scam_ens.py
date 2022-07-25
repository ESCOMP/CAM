#!/usr/bin/env python

import scam_case as scm
import numpy as np

basecase = 'nCTOPb3_L58_080.0E_30.0N_2010-07-01'
base=scm.scam_case()
base=base.unpickle_base(basecase)

print(base.__dict__)

#exit()

lats = 32.+ np.arange(6)
lons = 85.+ np.arange(2)

Lons,Lats = np.meshgrid(lons,lats)


dims=Lons.shape


x=[]
n=0
for j in range( dims[1] ):
    for i in range( dims[0]):
        x.append( scm.scam_case() )
        n=n+1

Lons=Lons.reshape( dims[0]*dims[1] )
Lats=Lats.reshape( dims[0]*dims[1] )

N=n
for n in range(N):
    ee = 'x_E'+str(n).zfill(2)
    x[n].changeTag( ee )
    x[n].changeLon( Lons[n] )
    x[n].changeLat( Lats[n] )
    x[n].startdate = base.startdate

for n in range(N):
    x[n].spawn_case(basecase)

for n in range(N):
    x[n].ensemble_member_run()

