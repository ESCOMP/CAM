class scam_case:
    #import sys
    #import getopt as go
    #import os

    def __init__(self):
        import os

        self.scmlat=36.6
        self.scmlon=270.
        self.nlev=58
        self.name="base_case"
        self.tag="case_tag"
        self.startdate=[2010,4,1]
        self.atm_ncpl=192
        self.mfilt=1
        self.nsteps=31*self.atm_ncpl
        self.coupler="nuopc"
        self.compiler="intel"
        self.basecase="base_case"
        self.isbasecase=True
        self.NameByBuild=False
        self.cime_output_root="dir"
        self.ensemble_root="none"
        self.project="P93300642"  # Needed for Cheyenne

        host=os.environ['HOST']
 
        if ('izumi' in host):
            self.machine="izumi"
        elif ('cheyenne' in host):
            self.machine="cheyenne"

    def base_case(self):
        import subprocess as sp
        import os
        import pickle

        user=os.environ['USER']
        
        y   = self.startdate[0]
        m   = self.startdate[1]
        d   = self.startdate[2]
        lat = self.scmlat
        lon = self.scmlon
        lev = self.nlev
        tag = self.tag


        case_day = str(d).zfill(2)
        case_mon = str(m).zfill(2)
        case_yr  = str(y).zfill(4)

        case_date  = case_yr+case_mon+case_day
        case_sdate = case_yr+"-"+case_mon+"-"+case_day

        case_lat = str(abs(lat)).zfill(4)
        if lat < 0:
            case_lat = case_lat+'S'
        elif lat > 0:
            case_lat = case_lat+'N'
        if lon < 0:
            lon = lon+360.

        case_lon = str(abs(lon)).zfill(5)+'E'
        case_lev = 'L'+str(abs(lev))

        lonstr = str(lon)
        latstr = str(lat)
        levstr = str(lev)

        NstepStr  = str( self.nsteps )
        NcplStr   = str( self.atm_ncpl )
        CplStr    = str( self.coupler )
        CompilerStr = str( self.compiler )
        MachStr = str( self.machine )
        ProjStr = str( self.project )

        if (self.NameByBuild == True):
            case_tag = tag+'_'+MachStr+'_'+CompilerStr+'_'+CplStr
        else:
            case_tag = tag+'_'+case_lev+'_'+case_lon+'_'+case_lat+'_'+case_yr+'-'+case_mon+'-'+case_day

        if ( self.coupler=="nuopc"):
            COMPSET="FSCAM"
        elif ( self.coupler=="mct"):
            COMPSET="2000_CAM60%SCAM_CLM50%SP_CICE5%PRES_DOCN%DOM_SROF_SGLC_SWAV"


        #This is the actual case name for the 'base' case
        print(case_tag)
        self.name=case_tag
        self.basecase=case_tag
        self.isbasecase=True

        #-----------------------------------------------------
        #  This script assumes that you will create new cases in 
        #  in a directory '../../cases' from 
        #  ${CESMROOT}/cime/scripts
        #-----------------------------------------------------
        cmd1="mkdir -p ../../cases/"+case_tag

        if (MachStr == 'cheyenne'):
            cmd2="./create_newcase --case  ../../cases/"+case_tag+ " --compset "+ COMPSET + " --res T42_T42 --driver " + CplStr + " --user-mods-dir ../../cime_config/usermods_dirs/scam_STUB --walltime 01:00:00 --mach " + MachStr + " --pecount 1 --compiler "+ CompilerStr + " --project "+ ProjStr + " --run-unsupported"
        else:
            cmd2="./create_newcase --case  ../../cases/"+case_tag+ " --compset "+ COMPSET + " --res T42_T42 --driver " + CplStr + " --user-mods-dir ../../cime_config/usermods_dirs/scam_STUB --walltime 01:00:00 --mach " + MachStr + " --pecount 1 --compiler "+ CompilerStr + " --run-unsupported" 

        cd0 = 'cd ../../cases/'+case_tag +';'

        #sp.run('date')
        sp.run(cmd2,shell=True)

        cmd = ( "./case.setup" )
        sp.run(cd0 + cmd   ,    shell=True )

        #cmd = ( "cp ../../cime_config/usermods_dirs/scam_STUB/scripts/STUB_iop.nc ./")
        #sp.run(cd0 + cmd   ,    shell=True )

        cmd = ( "cp ../../myPythonTools/STUB_iop.nc ./")
        sp.run(cd0 + cmd   ,    shell=True )
        cmd = ( "cp ../../myPythonTools/user_nl_cice ./")
        sp.run(cd0 + cmd   ,    shell=True )

        if (self.nlev==58) and (self.machine=='cheyenne'):
            cmd = ( "cp ../../myPythonTools/user_nl_cam_L58_cheyenne ./user_nl_cam")
        elif (self.nlev==93) and (self.machine=='cheyenne'):
            cmd = ( "cp ../../myPythonTools/user_nl_cam_L93_cheyenne ./user_nl_cam")
        else:
            cmd = ( "cp ../../myPythonTools/user_nl_cam ./")

        sp.run(cd0 + cmd   ,    shell=True )

        cmd = ( "./xmlchange DOUT_S_ROOT='/project/amp/"+user+"/scam/archive/"+case_tag+"'" + ";" +
                "./xmlchange CAM_CONFIG_OPTS='-dyn eul -scam -phys cam_dev -nlev "+levstr+"'"+ ";" +
                "./xmlchange ATM_NCPL="+NcplStr+";"+ 
                "./xmlchange STOP_N="+NstepStr+";"+
                "./xmlchange START_TOD=00000"+";"+
                "./xmlchange STOP_OPTION=nsteps"+";"+
                "./xmlchange PTS_LAT="+latstr+";"+ 
                "./xmlchange PTS_LON="+lonstr+";"+
                "./xmlchange RUN_STARTDATE="+case_sdate+";" 
        )
        sp.run(cd0 + cmd   ,    shell=True )
     

        cmd = ( 
            "ncap2 --overwrite -s bdate="+case_date+" STUB_iop.nc STUB_iop.nc"+";"+
            "ncap2 --overwrite -s lat[lat]="+latstr+" STUB_iop.nc STUB_iop.nc"+";"+
            "ncap2 --overwrite -s lon[lon]="+lonstr+" STUB_iop.nc STUB_iop.nc"+";"
        )
        sp.run(cd0 + cmd   ,    shell=True )


        print("Machine  = "+self.machine)
        print("Compiler = "+self.compiler)
        print("created and setup case=")
        print("   ../../cases/"+case_tag )
        print("Should be ready to build and submit")
        fname = '../../cases/'+case_tag+'/'+'env_build.xml'

        # find CIME_OUTPUT_ROOT
        fob=open(fname,"r")
        linin = fob.readlines()
        for line in linin:
            if (line.find('entry id="CIME_OUTPUT_ROOT"') !=-1):
                linspl1 = line.split('value=')
                linspl2 = linspl1[1].split('"')
        self.cime_output_root = linspl2[1]
        fob.close()

        # write 'self' to pickle file in casedir
        fname = '../../cases/'+case_tag+'/'+'BaseCaseSelf.pkL'
        with open( fname, 'wb') as fob:
            pickle.dump( self, fob )
        fob.close()


    def spawn_case(self,basecase):
        #---------------------------------------
        #  This function is still under development
        #---------------------------------------
        import subprocess as sp
        import os
        import pickle
        import txtutil as tx

        user=os.environ['USER']
        
        fname = '../../cases/'+basecase+'/'+'BaseCaseSelf.pkL'
        with open( fname, 'rb') as fob:
            base=pickle.load( fob )
        fob.close()

        #---------------------------------
        # These sould be inherited from 
        # base case or hardwried here
        #---------------------------------
        self.basecase   = base.name
        self.isbasecase = False
        self.nlev       = base.nlev
        self.coupler    = base.coupler
        self.compiler   = base.compiler
        self.atm_ncpl   = base.atm_ncpl
        self.cime_output_root = base.cime_output_root
        lev = base.nlev

        #-----------------------------------
        # These are specified in scam_drv.py
        # or scam_ens.py
        #-----------------------------------
        y   = self.startdate[0]
        m   = self.startdate[1]
        d   = self.startdate[2]
        lat = self.scmlat
        lon = self.scmlon
        tag = self.tag
        

        case_day = str(d).zfill(2)
        case_mon = str(m).zfill(2)
        case_yr  = str(y).zfill(4)

        case_date  = case_yr+case_mon+case_day
        case_sdate = case_yr+"-"+case_mon+"-"+case_day

        case_lat = str(abs(lat)).zfill(4)
        if lat < 0:
            case_lat = case_lat+'S'
        elif lat > 0:
            case_lat = case_lat+'N'
        if lon < 0:
            lon = lon+360.

        case_lon = str(abs(lon)).zfill(5)+'E'
        case_lev = 'L'+str(abs(lev))

        lonstr = str(lon)
        latstr = str(lat)
        levstr = str(lev)

        case_tag  = tag+'_'+case_lev+'_'+case_lon+'_'+case_lat+'_'+case_yr+'-'+case_mon+'-'+case_day
        self.name = case_tag

        ensemble_root = self.cime_output_root + '/' + self.basecase +'_ENS'
        self.ensemble_root = ensemble_root


        # --------------
        # Clean before making new directory
        #----------------
        cmd = ("rm -rf "+ ensemble_root + "/"+case_tag)
        sp.run( cmd , shell=True )

        cmd = ("mkdir -p "+ ensemble_root + "/"+case_tag+"/bld")
        sp.run( cmd , shell=True )

        cmd = ("cp -r "+ self.cime_output_root + "/" +  base.name +"/run" + " "
            + ensemble_root + "/"+case_tag+"/run")
        sp.run( cmd , shell=True )
 
        cmd = ("cp "+ self.cime_output_root + "/" +  base.name +"/bld/cesm.exe" + " "
            + ensemble_root + "/"+case_tag+"/bld/")
        sp.run( cmd , shell=True )

        cmd = ( "cp ../../myPythonTools/STUB_iop.nc"  + " "
            + ensemble_root + "/"+case_tag+"/run/")
        sp.run(cmd   ,    shell=True )
        cmd = ( "cp ../../myPythonTools/ens_run.sh"  + " "
            + ensemble_root + "/"+case_tag+"/run/")
        sp.run(cmd   ,    shell=True )

        cd0 ="cd "+ ensemble_root + "/" +  case_tag +"/run ;" 

        cmd = ( 
            "ncap2 --overwrite -s bdate="+case_date+" STUB_iop.nc STUB_iop.nc"+";"+
            "ncap2 --overwrite -s lat[lat]="+latstr+" STUB_iop.nc STUB_iop.nc"+";"+
            "ncap2 --overwrite -s lon[lon]="+lonstr+" STUB_iop.nc STUB_iop.nc"+";"
        )
        sp.run(cd0 + cmd   ,    shell=True )

        fili= ensemble_root + "/" +  case_tag +"/run/atm_in"
        tx.nmled(fili,'iopfile','"STUB_iop.nc"')
        # Set history to make one file per day
        tx.nmled(fili,'mfilt',str(self.mfilt) )

        if (base.coupler=='nuopc'):
            fili= ensemble_root + "/" +  case_tag +"/run/nuopc.runconfig"
            tx.nmled(fili,'case_name',case_tag)
            tx.nmled(fili,'start_ymd',case_date)
            tx.nmled(fili,'scol_lat',latstr)
            tx.nmled(fili,'scol_lon',lonstr)

        if (base.coupler=='mct'):
            fili= ensemble_root + "/" +  case_tag +"/run/drv_in"
            tx.nmled(fili,'case_name',case_tag)
            tx.nmled(fili,'start_ymd',case_date)
            tx.nmled(fili,'scmlat',latstr)
            tx.nmled(fili,'scmlon',lonstr)

    def ensemble_member_run(self):
        import subprocess as sp
        import os

        cd0 ="cd "+ self.ensemble_root + "/" +  self.name +"/run ;" 
        cmd ="/usr/local/torque/bin/qsub ens_run.sh"

        sp.run(cd0 + cmd   ,    shell=True )

    def unpickle_base(self,basecase):
        #---------------------------------------
        #  This function is still under development
        #---------------------------------------
        import pickle

        fname = '../../cases/'+basecase+'/'+'BaseCaseSelf.pkL'
        with open( fname, 'rb') as fob:
            base=pickle.load( fob )
        fob.close()

        return base

    def changeTag(self,newtag):
        #---------------------------------------
        #  This function is still under development
        #---------------------------------------

        self.tag = newtag

    def changeLon(self,newlon):
        #---------------------------------------
        #  This function is still under development
        #---------------------------------------

        self.scmlon = newlon

    def changeLat(self,newlat):
        #---------------------------------------
        #  This function is still under development
        #---------------------------------------

        self.scmlat = newlat

