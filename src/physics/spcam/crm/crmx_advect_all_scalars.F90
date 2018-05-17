subroutine advect_all_scalars()

  use crmx_vars
  use crmx_microphysics
  use crmx_sgs
  use crmx_crmtracers
#ifdef CLUBB_CRM
  use crmx_params, only: dotracers, doclubb, doclubbnoninter
#else
  use crmx_params, only: dotracers
#endif
  implicit none
  real dummy(nz)
  integer k


!---------------------------------------------------------
!      advection of scalars :

     call advect_scalar(t,tadv,twle,t2leadv,t2legrad,twleadv,.true.)
    
!
!    Advection of microphysics prognostics:
!

     do k = 1,nmicro_fields
        if(   k.eq.index_water_vapor             &! transport water-vapor variable no metter what
#ifdef CLUBB_CRM
!Added preprocessor directives. - nielsenb UWM 30 July 2008
        .or. ( docloud .or. doclubb .or. doclubbnoninter ) .and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#else
         .or. docloud.and.flag_precip(k).ne.1    & ! transport non-precipitation vars
#endif
         .or. doprecip.and.flag_precip(k).eq.1 ) &
           call advect_scalar(micro_field(:,:,:,k),mkadv(:,k),mkwle(:,k),dummy,dummy,dummy,.false.)
     end do

!
!    Advection of sgs prognostics:
!

     if(dosgs.and.advect_sgs) then
       do k = 1,nsgs_fields
           call advect_scalar(sgs_field(:,:,:,k),sgsadv(:,k),sgswle(:,k),dummy,dummy,dummy,.false.)
       end do
     end if


!
!   Precipitation fallout:
!
    if(doprecip) then

       total_water_prec = total_water_prec + total_water()

       call micro_precip_fall()

       total_water_prec = total_water_prec - total_water()


    end if

 ! advection of tracers:

     if(dotracers) then

        do k = 1,ntracers
         call advect_scalar(tracer(:,:,:,k),tradv(:,k),trwle(:,k),dummy,dummy,dummy,.false.)
        end do

     end if

end subroutine advect_all_scalars
