!-----------------------------------------------------------------------
!
! Manages the adjustment of ClOy, BrOy and IOy family components in response
! to conservation issues resulting from advection.
!
! Created as clybry_fam by: Francis Vitt
! Date: 21 May 2008
! Modified by Stacy Walters
! Date: 13 August 2008
!
! Extended for IOy family by: Carlos Ordonez
! (Implementation of IOy and update of both ClOy and BrOy
! Date: July 2012
!
! Updated by Rafa Fernandez
! Now it includes a logical condition in case iodine is not present
! so no need to maintain two different routines (clybry_fam and clybryiy_fam) 
! Date: Aug 2020
!-----------------------------------------------------------------------

module clybryiy_fam

  use shr_kind_mod,  only : r8 => shr_kind_r8
  use ppgrid,        only : pcols, pver
  use chem_mods,     only : gas_pcnst, adv_mass
  use constituents,  only : pcnst
  use short_lived_species,only: set_short_lived_species,get_short_lived_species

  implicit none

  save

  private
  public :: clybryiy_fam_set
  public :: clybryiy_fam_adj
  public :: clybryiy_fam_init

  integer :: id_cly,id_bry,id_iy

!rpf_CESM2_SLH
  integer :: id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2,id_clno2,id_cocl2,id_chcl2o2,id_cofcl
!rpf_CESM2_SLH

  integer :: id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr,id_br2,id_brno2
  integer :: id_i,id_i2,id_io,id_oio,id_hi,id_hoi,id_ino,id_ino2,id_iono2,id_ibr,id_icl,id_i2o2,id_i2o3,id_i2o4

  logical :: has_clybryiy

contains

  !------------------------------------------
  !------------------------------------------
  subroutine clybryiy_fam_init

    use mo_chem_utls, only : get_spc_ndx
    implicit none

!rpf_CESM2_SLH
    integer :: ids_clybry(16)
    integer :: ids_clybryiy(34)
!rpf_CESM2_SLH

    id_cly = get_spc_ndx('CLY')
    id_bry = get_spc_ndx('BRY')
    id_iy  = get_spc_ndx('IY')

    id_cl     = get_spc_ndx('CL')
    id_clo    = get_spc_ndx('CLO')
    id_hocl   = get_spc_ndx('HOCL')
    id_cl2    = get_spc_ndx('CL2')
    id_cl2o2  = get_spc_ndx('CL2O2')
    id_oclo   = get_spc_ndx('OCLO')
    id_hcl    = get_spc_ndx('HCL')
    id_clono2 = get_spc_ndx('CLONO2')
    id_clno2  = get_spc_ndx('CLNO2')

    id_br     = get_spc_ndx('BR')
    id_bro    = get_spc_ndx('BRO')
    id_hbr    = get_spc_ndx('HBR')
    id_brono2 = get_spc_ndx('BRONO2')
    id_brcl   = get_spc_ndx('BRCL')
    id_hobr   = get_spc_ndx('HOBR')
    id_br2    = get_spc_ndx('BR2')
    id_brno2  = get_spc_ndx('BRNO2')

    id_i      = get_spc_ndx('I')
    id_i2     = get_spc_ndx('I2')
    id_io     = get_spc_ndx('IO')
    id_oio    = get_spc_ndx('OIO')
    id_hi     = get_spc_ndx('HI')
    id_hoi    = get_spc_ndx('HOI')
    id_ino    = get_spc_ndx('INO')
    id_ino2   = get_spc_ndx('INO2')
    id_iono2  = get_spc_ndx('IONO2')
    id_ibr    = get_spc_ndx('IBR')
    id_icl    = get_spc_ndx('ICL')
    id_i2o2   = get_spc_ndx('I2O2')
    id_i2o3   = get_spc_ndx('I2O3')
    id_i2o4   = get_spc_ndx('I2O4')

!rpf_CESM2_SLH
    ids_clybry   = (/ id_cly,id_bry, &
             id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2, &
             id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr /)

    ids_clybryiy = (/ id_cly,id_bry,id_iy, &
             id_cl,id_clo,id_hocl,id_cl2,id_cl2o2,id_oclo,id_hcl,id_clono2,id_clno2, &
             id_br,id_bro,id_hbr,id_brono2,id_brcl,id_hobr,id_br2,id_brno2, &
             id_i,id_i2,id_io,id_oio,id_hi,id_hoi,id_ino,id_ino2,id_iono2,id_ibr,id_icl,id_i2o2,id_i2o3,id_i2o4 /)

    if ( id_iy>0 ) then
      has_clybryiy = all( ids_clybryiy(:) > 0 )
    else
      has_clybryiy = all( ids_clybry(:) > 0 )
    endif
!rpf_CESM2_SLH

  endsubroutine clybryiy_fam_init

!--------------------------------------------------------------
! set the ClOy, BrOy and IOy mass mixing ratios
!  - this is call before advection
!--------------------------------------------------------------
  subroutine clybryiy_fam_set( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : get_nstep
    use physics_buffer, only : physics_buffer_desc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

    real(r8) :: wrk(ncol,pver,3)
    real(r8) :: mmr(pcols,pver,gas_pcnst)
    integer  :: n, m

    if (.not. has_clybryiy) return

    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          mmr(:ncol,:,m) = q(:ncol,:, n)
       endif
    enddo
    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

!--------------------------------------------------------------
!       ... form updated chlorine, bromine atom mass mixing ratios
!--------------------------------------------------------------
    wrk(:,:,1) = cloy( mmr, pcols, ncol )
    wrk(:,:,2) = broy( mmr, pcols, ncol )
!rpf_CESM2_SLH
    if ( id_iy>0 ) then   
       wrk(:,:,3) = ioy ( mmr, pcols, ncol )
    endif
!rpf_CESM2_SLH

    mmr(:ncol,:,id_cly) = wrk(:,:,1)
    mmr(:ncol,:,id_bry) = wrk(:,:,2)
!rpf_CESM2_SLH
    if ( id_iy>0 ) then   
       mmr(:ncol,:,id_iy)  = wrk(:,:,3)
    endif
!rpf_CESM2_SLH

    call set_short_lived_species( mmr, lchnk, ncol, pbuf )
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          q(:ncol,:, n) = mmr(:ncol,:,m)
       endif
    enddo

  end subroutine clybryiy_fam_set

!--------------------------------------------------------------
! adjust the ClOy, BrOy and IOy individual family members 
!  - this is call after advection
!--------------------------------------------------------------
  subroutine clybryiy_fam_adj( ncol, lchnk, map2chm, q, pbuf )

    use time_manager,  only : is_first_step
    use physics_buffer, only : physics_buffer_desc

    implicit none

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in)    :: ncol, lchnk
    integer,  intent(in)    :: map2chm(pcnst)
    real(r8), intent(inout) :: q(pcols,pver,pcnst)
    type(physics_buffer_desc), pointer :: pbuf(:)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: factor(ncol,pver)
    real(r8) :: wrk(ncol,pver)
    real(r8) :: mmr(pcols,pver,gas_pcnst)

    integer  :: n, m

    if (.not. has_clybryiy) return

!--------------------------------------------------------------
!       ... CLY,BRY,IY are not adjusted until the end of the first timestep
!--------------------------------------------------------------
    if (is_first_step()) return

    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          mmr(:ncol,:,m) = q(:ncol,:, n)
       endif
    enddo
    call get_short_lived_species( mmr, lchnk, ncol, pbuf )

    !--------------------------------------------------------------
    !       ... form updated chlorine atom mass mixing ratio
    !--------------------------------------------------------------
    if ( id_cly>0 ) then
       wrk(:,:) = cloy( mmr, pcols, ncol )

       factor(:ncol,:) = mmr(:ncol,:,id_cly) / wrk(:ncol,:)
       !--------------------------------------------------------------
       !       ... adjust "group" members
       !--------------------------------------------------------------
       mmr(:ncol,:,id_cl)     = factor(:ncol,:)*mmr(:ncol,:,id_cl)
       mmr(:ncol,:,id_clo)    = factor(:ncol,:)*mmr(:ncol,:,id_clo)
       mmr(:ncol,:,id_hocl)   = factor(:ncol,:)*mmr(:ncol,:,id_hocl)
       mmr(:ncol,:,id_cl2)    = factor(:ncol,:)*mmr(:ncol,:,id_cl2)
       mmr(:ncol,:,id_cl2o2)  = factor(:ncol,:)*mmr(:ncol,:,id_cl2o2)
       mmr(:ncol,:,id_oclo)   = factor(:ncol,:)*mmr(:ncol,:,id_oclo)
       mmr(:ncol,:,id_hcl)    = factor(:ncol,:)*mmr(:ncol,:,id_hcl)
       mmr(:ncol,:,id_clono2) = factor(:ncol,:)*mmr(:ncol,:,id_clono2)
!rpf_CESM2_SLH
       if (id_clno2>0)   mmr(:ncol,:,id_clno2)    = factor(:ncol,:)*mmr(:ncol,:,id_clno2)
       if (id_chcl2o2>0) mmr(:ncol,:,id_chcl2o2)  = factor(:ncol,:)*mmr(:ncol,:,id_chcl2o2)
       if (id_cocl2>0)   mmr(:ncol,:,id_cocl2)    = factor(:ncol,:)*mmr(:ncol,:,id_cocl2)
       if (id_cofcl>0)   mmr(:ncol,:,id_cofcl)    = factor(:ncol,:)*mmr(:ncol,:,id_cofcl)
!rpf_CESM2_SLH
    endif

    !--------------------------------------------------------------
    !        ... form updated bromine atom mass mixing ratio
    !--------------------------------------------------------------
    if ( id_bry>0 ) then
       wrk(:,:) = broy( mmr, pcols, ncol )

       factor(:ncol,:) = mmr(:ncol,:,id_bry) / wrk(:ncol,:)
       !--------------------------------------------------------------
       !       ... adjust "group" members
       !--------------------------------------------------------------
       mmr(:ncol,:,id_br)     = factor(:ncol,:)*mmr(:ncol,:,id_br)
       mmr(:ncol,:,id_bro)    = factor(:ncol,:)*mmr(:ncol,:,id_bro)
       mmr(:ncol,:,id_hbr)    = factor(:ncol,:)*mmr(:ncol,:,id_hbr)
       mmr(:ncol,:,id_brono2) = factor(:ncol,:)*mmr(:ncol,:,id_brono2)
       mmr(:ncol,:,id_brcl)   = factor(:ncol,:)*mmr(:ncol,:,id_brcl)
       mmr(:ncol,:,id_hobr)   = factor(:ncol,:)*mmr(:ncol,:,id_hobr)
       if (id_br2>0)   mmr(:ncol,:,id_br2)    = factor(:ncol,:)*mmr(:ncol,:,id_br2)
       if (id_brno2>0) mmr(:ncol,:,id_brno2)  = factor(:ncol,:)*mmr(:ncol,:,id_brno2)
    endif

!rpf_CESM2_SLH
    !--------------------------------------------------------------
    !        ... form updated iodine atom mass mixing ratio
    !--------------------------------------------------------------
    if ( id_iy>0 ) then
       wrk(:,:) = ioy( mmr, pcols, ncol )

       factor(:ncol,:) = mmr(:ncol,:,id_iy) / wrk(:ncol,:)
       !--------------------------------------------------------------
       !       ... adjust "group" members
       !--------------------------------------------------------------
       mmr(:ncol,:,id_i)      = factor(:ncol,:)*mmr(:ncol,:,id_i)
       mmr(:ncol,:,id_i2)     = factor(:ncol,:)*mmr(:ncol,:,id_i2)
       mmr(:ncol,:,id_io)     = factor(:ncol,:)*mmr(:ncol,:,id_io)
       mmr(:ncol,:,id_oio)    = factor(:ncol,:)*mmr(:ncol,:,id_oio)
       mmr(:ncol,:,id_hi)     = factor(:ncol,:)*mmr(:ncol,:,id_hi)
       mmr(:ncol,:,id_hoi)    = factor(:ncol,:)*mmr(:ncol,:,id_hoi)
       mmr(:ncol,:,id_ino)    = factor(:ncol,:)*mmr(:ncol,:,id_ino)
       mmr(:ncol,:,id_ino2)   = factor(:ncol,:)*mmr(:ncol,:,id_ino2)
       mmr(:ncol,:,id_iono2)  = factor(:ncol,:)*mmr(:ncol,:,id_iono2)
       mmr(:ncol,:,id_ibr)    = factor(:ncol,:)*mmr(:ncol,:,id_ibr)
       mmr(:ncol,:,id_icl)    = factor(:ncol,:)*mmr(:ncol,:,id_icl)
       mmr(:ncol,:,id_i2o2)   = factor(:ncol,:)*mmr(:ncol,:,id_i2o2)
       mmr(:ncol,:,id_i2o3)   = factor(:ncol,:)*mmr(:ncol,:,id_i2o3)
       mmr(:ncol,:,id_i2o4)   = factor(:ncol,:)*mmr(:ncol,:,id_i2o4)
    endif
!rpf_CESM2_SLH

    call set_short_lived_species( mmr, lchnk, ncol, pbuf )
    do n = 1,pcnst
       m = map2chm(n)
       if( m > 0 ) then
          q(:ncol,:, n) = mmr(:ncol,:,m)
       endif
    enddo

  end subroutine clybryiy_fam_adj

!--------------------------------------------------------------
! private methods
!--------------------------------------------------------------

!--------------------------------------------------------------
! compute the mass mixing ratio of ClOy
!--------------------------------------------------------------
  function cloy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: cloy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    do k = 1,pver
       wrk(:) = q(:ncol,k,id_cl)           /adv_mass(id_cl)      &
              + q(:ncol,k,id_clo)          /adv_mass(id_clo)     &
              + q(:ncol,k,id_hocl)         /adv_mass(id_hocl)    &
              + 2._r8*( q(:ncol,k,id_cl2)  /adv_mass(id_cl2)     &
                      + q(:ncol,k,id_cl2o2)/adv_mass(id_cl2o2) ) &
              + q(:ncol,k,id_oclo)         /adv_mass(id_oclo)    &
              + q(:ncol,k,id_hcl)          /adv_mass(id_hcl)     &
              + q(:ncol,k,id_clono2)       /adv_mass(id_clono2)

       if (id_clno2>0) then
          wrk(:) =  wrk(:) &
              + q(:ncol,k,id_clno2)        /adv_mass(id_clno2)
       endif
!rpf_CESM2_SLH
       if (id_chcl2o2>0) then
          wrk(:) =  wrk(:) &
              + 2._r8*q(:ncol,k,id_chcl2o2)/adv_mass(id_chcl2o2)
       endif
       if (id_cocl2>0) then
          wrk(:) =  wrk(:) &
              + 2._r8*q(:ncol,k,id_cocl2)  /adv_mass(id_cocl2)
       endif
       if (id_cofcl>0) then
          wrk(:) =  wrk(:) &
              + q(:ncol,k,id_cofcl)        /adv_mass(id_cofcl)
       endif
!rpf_CESM2_SLH

       cloy(:,k) = adv_mass(id_cl) * wrk(:)
    end do

  end function cloy

!--------------------------------------------------------------
! compute the mass mixing ratio of BrOy
!--------------------------------------------------------------
  function broy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: broy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    do k = 1,pver
       wrk(:) = q(:ncol,k,id_br)         /adv_mass(id_br)      &
              + q(:ncol,k,id_bro)        /adv_mass(id_bro)     &
              + q(:ncol,k,id_hbr)        /adv_mass(id_hbr)     &
              + q(:ncol,k,id_brono2)     /adv_mass(id_brono2)  &
              + q(:ncol,k,id_brcl)       /adv_mass(id_brcl)    &
              + q(:ncol,k,id_hobr)       /adv_mass(id_hobr) 
       if (id_br2>0) then
          wrk(:) = wrk(:) &
              + 2._r8*q(:ncol,k,id_br2)  /adv_mass(id_br2)
       endif

       if (id_brno2>0) then
          wrk(:) = wrk(:) &
              + q(:ncol,k,id_brno2)      /adv_mass(id_brno2)   
       endif

       broy(:,k) = adv_mass(id_br) * wrk(:)
    end do

  end function broy

!rpf_CESM2_SLH
!--------------------------------------------------------------
! compute the mass mixing ratio of IOy
!--------------------------------------------------------------
  function ioy( q, pcols, ncol )

!--------------------------------------------------------------
!       ... dummy arguments
!--------------------------------------------------------------
    integer,  intent(in) :: pcols
    integer,  intent(in) :: ncol
    real(r8), intent(in) :: q(pcols,pver,gas_pcnst)

!--------------------------------------------------------------
!       ... function declaration
!--------------------------------------------------------------
    real(r8) :: ioy(ncol,pver)

!--------------------------------------------------------------
!       ... local variables
!--------------------------------------------------------------
    real(r8) :: wrk(ncol)
    integer  :: k

    do k = 1,pver
       wrk(:) = q(:ncol,k,id_i)             /adv_mass(id_i)      &
              + 2._r8*q(:ncol,k,id_i2)      /adv_mass(id_i2)     &
              + q(:ncol,k,id_io)            /adv_mass(id_io)     &
              + q(:ncol,k,id_oio)           /adv_mass(id_oio)    &
              + q(:ncol,k,id_hi)            /adv_mass(id_hi)     &
              + q(:ncol,k,id_hoi)           /adv_mass(id_hoi)    &
              + q(:ncol,k,id_ino)           /adv_mass(id_ino)    &
              + q(:ncol,k,id_ino2)          /adv_mass(id_ino2)   &
              + q(:ncol,k,id_iono2)         /adv_mass(id_iono2)  &
              + q(:ncol,k,id_ibr)           /adv_mass(id_ibr)    &
              + q(:ncol,k,id_icl)           /adv_mass(id_icl)    &
              + 2._r8*( q(:ncol,k,id_i2o2)  /adv_mass(id_i2o2)   &
                      + q(:ncol,k,id_i2o3)  /adv_mass(id_i2o3)   &
                      + q(:ncol,k,id_i2o4)  /adv_mass(id_i2o4) )

       ioy(:,k) = adv_mass(id_i) * wrk(:)
    end do

  end function ioy
!rpf_CESM2_SLH

end module clybryiy_fam
