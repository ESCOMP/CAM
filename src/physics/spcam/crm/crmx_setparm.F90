module crmx_setparm_mod

contains

subroutine setparm
	
!       initialize parameters:

use crmx_vars
!use micro_params
use crmx_params
use crmx_microphysics, only: micro_setparm
use crmx_sgs, only: sgs_setparm

implicit none
	
integer icondavg, ierr

!NAMELIST /PARAMETERS/ dodamping, doupperbound, docloud, doprecip, &
!                dolongwave, doshortwave, dosgs, &
!                docoriolis, dosurface, dolargescale, doradforcing, &
!		nadams,fluxt0,fluxq0,tau0,tabs_s,z0,tauls,nelapse, &
!		dt, dx, dy, fcor, ug, vg, nstop, caseid, &
!		nstat, nstatfrq, nprint, nrestart, doradsimple, &
!		nsave3D, nsave3Dstart, nsave3Dend, dosfcforcing, &
!		donudging_uv, donudging_tq, dosmagor, doscalar,  &
!		timelargescale, longitude0, latitude0, day0, nrad, &
!		CEM,LES,OCEAN,LAND,SFC_FLX_FXD,SFC_TAU_FXD, soil_wetness, &
!                doensemble, nensemble, doxy, dowallx, dowally, &
!                nsave2D, nsave2Dstart, nsave2Dend, qnsave3D, & 
!                docolumn, save2Dbin, save2Davg, save3Dbin, &
!                save2Dsep, save3Dsep, dogzip2D, dogzip3D, restart_sep, &
!	        doseasons, doperpetual, doradhomo, dosfchomo, doisccp, &
!	        dodynamicocean, ocean_type, &
!		dosolarconstant, solar_constant, zenith_angle, rundatadir, &
!                dotracers, output_sep, perturb_type, &
!                doSAMconditionals, dosatupdnconditionals, &
!                doscamiopdata, iopfile, dozero_out_day0, &
!                nstatmom, nstatmomstart, nstatmomend, savemomsep, savemombin, &
!                nmovie, nmoviestart, nmovieend, nrestart_skip, &
!                bubble_x0,bubble_y0,bubble_z0,bubble_radius_hor, &
!               bubble_radius_ver,bubble_dtemp,bubble_dq, dosmoke, &
!               doclubb, doclubbnoninter, doclubb_sfc_fluxes, & ! added by dschanen UWM
!               docam_sfc_fluxes           ! added by mhwang 
	


!----------------------------------
!  Read namelist variables from the standard input:
!------------

!open(55,file='./'//trim(case)//'/prm', status='old',form='formatted') 
!read (55,PARAMETERS)
!close(55)

	doprecip        = .true.
	dosgs		= .true.
	dosurface	= .true.
	dodamping 	= .true.
	dt		= CRM_DT
	dx		= CRM_DX
	dy		= CRM_DY
        CEM             = .true.
#ifndef CLUBB_CRM
        doclubb         = .false.   ! then docloud must be .true.
        docloud         = .true.
#else
        doclubb         = .true.    ! then docloud must be .false.
        docloud         = .false.
        doclubbnoninter = .false.
        doclubb_sfc_fluxes = .false.
        docam_sfc_fluxes = .true.   ! update variables in cam, neither in sam nor in clubb +++mhwang
        nclubb          = 3 

#ifdef sam1mom
! for sam1mom, nclubb needs to be 1. 
! see comments in ./MICRO_SAM1MOM/microphysics.F90
        nclubb          = 3 
#endif

#endif
        rank            = 0   ! in MMF model, rank = 0
!------------------------------------
!  Set parameters 


        ! Allow only special cases for separate output:

        output_sep = output_sep.and.RUN3D
        if(output_sep)  save2Dsep = .true.

	if(RUN2D) dy=dx

	if(RUN2D.and.YES3D.eq.1) then
	  print*,'Error: 2D run and YES3D is set to 1. Exitting...'
	  call task_abort()
	endif
	if(RUN3D.and.YES3D.eq.0) then
	  print*,'Error: 3D run and YES3D is set to 0. Exitting...'
	  call task_abort()
	endif
#ifdef CLUBB_CRM
        if ( dx >= 1000. .and. LES ) then
          print*,'Error: Horizonatal grid spacing is >= 1000. meters'
          print*,'but LES is true.  Use CEM mode for coarse resolutions.'
          call task_abort()
        end if
#endif

	if(fcor.eq.-999.) fcor= 4*pi/86400.*sin(latitude0*pi/180.)
	fcorz = sqrt(4.*(2*pi/(3600.*24.))**2-fcor**2)	  
	
	if(ny.eq.1) dy=dx
	dtn = dt

	notopened2D = .true.
	notopened3D = .true.

!        call zero_instr_diag() ! initialize instruments output 
        call sgs_setparm() ! read in SGS options from prm file.
        call micro_setparm() ! read in microphysical options from prm file.

        if(dosmoke) then
           epsv=0.
        else    
           epsv=0.61
        endif   

        if(navgmom_x.lt.0.or.navgmom_y.lt.0) then  
            nstatmom        = 1
            nstatmomstart    = 99999999
            nstatmomend      = 999999999
        end if

        if(tautqls.eq.99999999.) tautqls = tauls

        masterproc = rank.eq.0
          
end subroutine setparm
end module crmx_setparm_mod
