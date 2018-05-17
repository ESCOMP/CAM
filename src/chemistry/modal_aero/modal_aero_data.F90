      module modal_aero_data

!--------------------------------------------------------------
! ... Basic aerosol mode parameters and arrays
!--------------------------------------------------------------
      use shr_kind_mod,    only: r8 => shr_kind_r8
      use constituents,    only: pcnst, cnst_mw, cnst_name, cnst_get_ind, cnst_set_convtran2, &
                                 cnst_set_spec_class, cnst_spec_class_aerosol, cnst_spec_class_undefined, &
                                 cnst_species_class, cnst_spec_class_gas
      use physics_buffer,  only: pbuf_add_field, dtype_r8
      use time_manager,    only: is_first_step
      use phys_control,    only: phys_getopts
      use infnan,          only: nan, assignment(=)
      use cam_logfile,     only: iulog
      use cam_abortutils,  only: endrun
      use spmd_utils,      only: masterproc
      use ppgrid,          only: pcols, pver, begchunk, endchunk
      use mo_tracname,     only: solsym
      use chem_mods,       only: gas_pcnst  
      use radconstants,    only: nswbands, nlwbands
      use shr_const_mod,   only: pi => shr_const_pi
      use rad_constituents,only: rad_cnst_get_info, rad_cnst_get_aer_props, rad_cnst_get_mode_props
      use physics_buffer,  only: physics_buffer_desc, pbuf_get_chunk

      implicit none
      private

      public :: modal_aero_data_init
      public :: modal_aero_data_reg
      public :: qqcw_get_field

      integer, public, protected :: nsoa = 0
      integer, public, protected :: npoa = 0
      integer, public, protected :: nbc = 0
      integer, public, protected :: nspec_max = 0
      integer, public, protected :: ntot_amode = 0
      integer, public, protected :: nSeaSalt=0, nDust=0
      integer, public, protected :: nSO4=0, nNH4=0

      !
      ! definitions for aerosol chemical components
      !
      
      real(r8), public, protected, allocatable :: specmw_amode(:,:)
      character(len=16), public, protected, allocatable :: modename_amode(:)

      integer, public, protected, allocatable :: nspec_amode(:)

      character(len=20), public, protected :: cnst_name_cw( pcnst )

      !   input mprognum_amode, mdiagnum_amode, mprogsfc_amode, mcalcwater_amode
      integer, public, protected, allocatable :: mprognum_amode(:)
      integer, public, protected, allocatable :: mdiagnum_amode(:)
      integer, public, protected, allocatable :: mprogsfc_amode(:)
      integer, public, protected, allocatable :: mcalcwater_amode(:)

      !   input dgnum_amode, dgnumlo_amode, dgnumhi_amode (units = m)
      real(r8), public, protected, allocatable :: dgnum_amode(:)
      real(r8), public, protected, allocatable :: dgnumlo_amode(:)
      real(r8), public, protected, allocatable :: dgnumhi_amode(:)

      !   input sigmag_amode
      real(r8), public, protected, allocatable :: sigmag_amode(:)

      !   input crystalization and deliquescence points
      real(r8), allocatable :: rhcrystal_amode(:)
      real(r8), allocatable :: rhdeliques_amode(:)


      integer, public, protected, allocatable :: &
           lmassptr_amode( :, : ),   &
           lmassptrcw_amode( :, : ), &
           numptr_amode( : ),        &
           numptrcw_amode( : )

      real(r8), public, protected, allocatable :: &
           alnsg_amode( : ),              &
           voltonumb_amode( : ),          &
           voltonumblo_amode( : ),        &
           voltonumbhi_amode( : ),        &
           alnv2n_amode( : ),             &
           alnv2nlo_amode( : ),           &
           alnv2nhi_amode( : ),           &
           specdens_amode(:,:),   &
           spechygro(:,:)

      integer, public, protected, allocatable ::  &
           lptr_so4_a_amode(:),     lptr_so4_cw_amode(:), &
           lptr_msa_a_amode(:),     lptr_msa_cw_amode(:), &
           lptr_nh4_a_amode(:),     lptr_nh4_cw_amode(:), &
           lptr_no3_a_amode(:),     lptr_no3_cw_amode(:), &
           lptr_nacl_a_amode(:),    lptr_nacl_cw_amode(:),&
           lptr_dust_a_amode(:),    lptr_dust_cw_amode(:)

      integer, public, protected :: &
           modeptr_accum,  modeptr_aitken,                               &
           modeptr_ufine,  modeptr_coarse,                               &
           modeptr_pcarbon,                                              &
           modeptr_finedust,  modeptr_fineseas,                          &
           modeptr_coardust,  modeptr_coarseas

      !2D lptr variables added by RCE to access speciated species
      integer, public, protected, allocatable :: &
           lptr2_bc_a_amode(:,:),   lptr2_bc_cw_amode(:,:), &
           lptr2_pom_a_amode(:,:),  lptr2_pom_cw_amode(:,:), &
           lptr2_soa_a_amode(:,:),  lptr2_soa_cw_amode(:,:), &
           lptr2_soa_g_amode(:)

      real(r8), public, protected :: specmw_so4_amode

      logical, public, protected :: soa_multi_species = .false.

      character(len=16), allocatable :: xname_massptr(:,:)     ! names of species in each mode
      character(len=16), allocatable :: xname_massptrcw(:,:)   ! names of cloud-borne species in each mode

      complex(r8), allocatable :: &
           specrefndxsw( :,:,: ), &
           specrefndxlw( :,:,: )

      character(len=8), allocatable :: &
           aodvisname(: ), &
           ssavisname(: )
      character(len=48), allocatable :: &
           aodvislongname(: ), &
           ssavislongname(: )

      character(len=8), allocatable :: &
           fnactname(: ), &
           fmactname(: ), &
           nactname(: )

      character(len=48), allocatable :: &
           fnactlongname(: ), &
           fmactlongname(: ), &
           nactlongname(: )


      !   threshold for reporting negatives from subr qneg3
      real(r8) :: qneg3_worst_thresh_amode(pcnst)

      integer :: qqcw(pcnst)=-1 ! Remaps modal_aero indices into pbuf
      
      logical :: convproc_do_aer
      logical :: cam_do_aero_conv = .true.
    contains

!--------------------------------------------------------------
!--------------------------------------------------------------
  subroutine modal_aero_data_reg

    character(len=6) :: xname_numptr, xname_numptrcw
    character(len=1) :: modechr
    integer :: m, l, iptr,i, n, tot_spec, idx
    character(len=3) :: trnum       ! used to hold mode number (as characters)

    character(len=20) :: dumStr1, specNameMode
    character(len=1000) :: msg
    character(len=32) :: spec_name, mode_type
    character(len=1) :: modestr

    call rad_cnst_get_info( 0, nmodes=ntot_amode )
    allocate( nspec_amode(ntot_amode) )
    allocate( numptr_amode(ntot_amode) )
    allocate( numptrcw_amode(ntot_amode) )
    allocate(modename_amode(ntot_amode))
    allocate(mprognum_amode(ntot_amode))
    allocate(mdiagnum_amode(ntot_amode))
    allocate(mprogsfc_amode(ntot_amode))
    allocate(mcalcwater_amode(ntot_amode))
    mprognum_amode(:) = 1
    mdiagnum_amode(:) = 0
    mprogsfc_amode(:) = 0
    if (ntot_amode==7) then
       mcalcwater_amode(:) = 1
    else
       mcalcwater_amode(:) = 0
    endif
    allocate(dgnum_amode(ntot_amode))
    allocate(dgnumlo_amode(ntot_amode))
    allocate(dgnumhi_amode(ntot_amode))
    allocate(sigmag_amode(ntot_amode))
    allocate(rhcrystal_amode(ntot_amode))
    allocate(rhdeliques_amode(ntot_amode))
    allocate( &
         alnsg_amode( ntot_amode ),              &   !
         voltonumb_amode( ntot_amode ),          &   !
         voltonumblo_amode( ntot_amode ),        &   !
         voltonumbhi_amode( ntot_amode ),        &   !
         alnv2n_amode( ntot_amode ),             &   !
         alnv2nlo_amode( ntot_amode ),           &   !
         alnv2nhi_amode( ntot_amode ),           &   !
         aodvisname(ntot_amode ), &
         ssavisname(ntot_amode ), &
         fnactname(ntot_amode ), &
         fmactname(ntot_amode ), &
         nactname(ntot_amode ), &
         fnactlongname(ntot_amode ), &
         fmactlongname(ntot_amode ), &
         nactlongname(ntot_amode ), &
         lptr_so4_a_amode(ntot_amode), lptr_so4_cw_amode(ntot_amode), &
         lptr_msa_a_amode(ntot_amode), lptr_msa_cw_amode(ntot_amode), &
         lptr_nh4_a_amode(ntot_amode), lptr_nh4_cw_amode(ntot_amode), &
         lptr_nacl_a_amode(ntot_amode), lptr_nacl_cw_amode(ntot_amode), &
         lptr_dust_a_amode(ntot_amode), lptr_dust_cw_amode(ntot_amode), &
         lptr_no3_a_amode(ntot_amode), lptr_no3_cw_amode(ntot_amode) &
    )

    allocate( &
         aodvislongname(ntot_amode ), &
         ssavislongname(ntot_amode ) &
    )

    do m = 1, ntot_amode
       call rad_cnst_get_info(0, m, mode_type=mode_type, nspec=nspec_amode(m))
       modename_amode(m) = mode_type
       ! count number of soa, poa, and bc bins in mode 1
       if (m==1) then
          do l = 1, nspec_amode(m)
             call rad_cnst_get_info(0, m, l, spec_name=spec_name )
             if (spec_name(:3) == 'soa') nsoa=nsoa+1
             if (spec_name(:3) == 'pom') npoa=npoa+1
             if (spec_name(:2) == 'bc' ) nbc =nbc +1
          enddo
       endif
    enddo

    soa_multi_species = nsoa > 1

    nspec_max = maxval( nspec_amode )

    allocate ( specdens_amode(nspec_max,ntot_amode) )
    allocate ( spechygro(nspec_max,ntot_amode) )
    allocate ( specmw_amode(nspec_max,ntot_amode) )
    allocate ( xname_massptr(nspec_max,ntot_amode) )
    allocate ( xname_massptrcw(nspec_max,ntot_amode) )
    specmw_amode = nan
    xname_massptr(:,:) = ' '
    xname_massptrcw(:,:) = ' '

    do m = 1, ntot_amode
       do l = 1, nspec_amode(m)
          call rad_cnst_get_info(0, m, l, spec_name=spec_name )
          xname_massptr(l,m) = spec_name
          write(modestr,'(I1)') m
          idx = index( xname_massptr(l,m), '_' )
          xname_massptrcw(l,m) = xname_massptr(l,m)(:idx-1)//'_c'//modestr
          if (xname_massptr(l,m)(:3) == 'dst') nDust=nDust+1
          if (xname_massptr(l,m)(:3) == 'ncl') nSeaSalt=nSeaSalt+1
          if (xname_massptr(l,m)(:3) == 'nh4') nNH4=nNH4+1
          if (xname_massptr(l,m)(:3) == 'so4') nSO4=nSO4+1
       enddo
    enddo
    
    allocate( &
         lmassptr_amode( nspec_max, ntot_amode ),&
         lmassptrcw_amode( nspec_max, ntot_amode ),&
         lptr2_pom_a_amode(ntot_amode,npoa),  lptr2_pom_cw_amode(ntot_amode,npoa), &
         lptr2_soa_a_amode(ntot_amode,nsoa),  lptr2_soa_cw_amode(ntot_amode,nsoa), &
         lptr2_bc_a_amode(ntot_amode,nbc),   lptr2_bc_cw_amode(ntot_amode,nbc), &
         lptr2_soa_g_amode(nsoa) &
         )
    lptr2_soa_g_amode = -999999

    allocate( specrefndxsw(nswbands,nspec_max,ntot_amode ) )
    allocate( specrefndxlw(nlwbands,nspec_max,ntot_amode) )

    do m = 1, ntot_amode
       if(nspec_amode(m).gt.nspec_max)then
          write(iulog,*)'modal_aero_data_reg: nspec_amode(m).gt.nspec_max '
          write(iulog,*)'modal_aero_data_reg: m,nspec_amode(m), nspec_max=',m, nspec_amode(m), nspec_max
          call endrun('modal_aero_data_reg: nspec_amode(m).gt.nspec_max ')
       end if
    end do

    call phys_getopts(convproc_do_aer_out = convproc_do_aer)
    if (convproc_do_aer) cam_do_aero_conv = .false.

    do m = 1, ntot_amode
       write(modechr,fmt='(I1)') m
       xname_numptr = 'num_a'//modechr
       xname_numptrcw  = 'num_c'//modechr

       if (masterproc) then
          write(iulog,9231) m, modename_amode(m)
          write(iulog,9232)                                          &
               'nspec                       ',                         &
               nspec_amode(m)
          write(iulog,9232)                                          &
               'mprognum, mdiagnum, mprogsfc',                         &
               mprognum_amode(m), mdiagnum_amode(m), mprogsfc_amode(m)
          write(iulog,9232)                                          &
               'mcalcwater                  ',                         &
               mcalcwater_amode(m)
       endif

       !    define species to hold interstitial & activated number
       call search_list_of_names(                                      &
            xname_numptr, numptr_amode(m), cnst_name, pcnst )
       if (numptr_amode(m) .le. 0) then
          write(iulog,9061) 'xname_numptr', xname_numptr, m
          call endrun('modal_aero_data_reg: numptr_amode(m) .le. 0')
       end if
       if (numptr_amode(m) .gt. pcnst) then
          write(iulog,9061) 'numptr_amode', numptr_amode(m), m
          write(iulog,9061) 'xname_numptr', xname_numptr, m
          call endrun('modal_aero_data_reg: numptr_amode(m) .gt. pcnst')
       end if

       call cnst_set_spec_class(numptr_amode(m), cnst_spec_class_aerosol)
       call cnst_set_convtran2(numptr_amode(m), cam_do_aero_conv)

       numptrcw_amode(m) = numptr_amode(m)  !use the same index for Q and QQCW arrays
       if (numptrcw_amode(m) .le. 0) then
          write(iulog,9061) 'xname_numptrcw', xname_numptrcw, m
          call endrun('modal_aero_data_reg: numptrcw_amode(m) .le. 0')
       end if
       if (numptrcw_amode(m) .gt. pcnst) then
          write(iulog,9061) 'numptrcw_amode', numptrcw_amode(m), m
          write(iulog,9061) 'xname_numptrcw', xname_numptrcw, m
          call endrun('modal_aero_data_reg: numptrcw_amode(m) .gt. pcnst')
       end if

       call pbuf_add_field(xname_numptrcw,'global',dtype_r8,(/pcols,pver/),iptr)
       call qqcw_set_ptr(numptrcw_amode(m),iptr)

       !   output mode information
       if ( masterproc ) then
          write(iulog,9233) 'numptr         ',                           &
               numptr_amode(m), xname_numptr
          write(iulog,9233) 'numptrcw       ',                           &
               numptrcw_amode(m), xname_numptrcw
       end if

       !   define the chemical species for the mode
       do l = 1, nspec_amode(m)

          call search_list_of_names(                                  &
               xname_massptr(l,m), lmassptr_amode(l,m), cnst_name, pcnst )
          if (lmassptr_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_massptr', xname_massptr(l,m), l, m
             write(iulog,'(10(a8,1x))')(cnst_name(i),i=1,pcnst)
             call endrun('modal_aero_data_reg: lmassptr_amode(l,m) .le. 0')
          end if
          call cnst_set_spec_class(lmassptr_amode(l,m), cnst_spec_class_aerosol)
          call cnst_set_convtran2(lmassptr_amode(l,m), cam_do_aero_conv)

          lmassptrcw_amode(l,m) = lmassptr_amode(l,m)  !use the same index for Q and QQCW arrays
          if (lmassptrcw_amode(l,m) .le. 0) then
             write(iulog,9062) 'xname_massptrcw', xname_massptrcw(l,m), l, m
             call endrun('modal_aero_data_reg: lmassptrcw_amode(l,m) .le. 0')
          end if
          call pbuf_add_field(xname_massptrcw(l,m),'global',dtype_r8,(/pcols,pver/),iptr)
          call qqcw_set_ptr(lmassptrcw_amode(l,m), iptr)

          if ( masterproc ) then
             write(iulog,9236) 'spec, massptr  ', l,                    &
                  lmassptr_amode(l,m), xname_massptr(l,m)
             write(iulog,9236) 'spec, massptrcw', l,                    &
                  lmassptrcw_amode(l,m), xname_massptrcw(l,m)
          end if

       end do

       !   set names for aodvis and ssavis
       write(unit=trnum,fmt='(i3)') m+100

       aodvisname(m) = 'AODVIS'//trnum(2:3)
       aodvislongname(m) = 'Aerosol optical depth for mode '//trnum(2:3)
       ssavisname(m) = 'SSAVIS'//trnum(2:3)
       ssavislongname(m) = 'Single-scatter albedo for mode '//trnum(2:3)
       fnactname(m) = 'FNACT'//trnum(2:3)
       fnactlongname(m) = 'Number faction activated for mode '//trnum(2:3)
       fmactname(m) = 'FMACT'//trnum(2:3)
       fmactlongname(m) = 'Fraction mass activated for mode'//trnum(2:3)
    end do

    ! At this point, species_class is either undefined or aerosol.
    ! For the "chemistry species" set the undefined ones to gas,
    ! and leave the aerosol ones as is
    do i = 1, gas_pcnst
       call cnst_get_ind(solsym(i), idx, abort=.false.)
       if (idx > 0) then
          if (cnst_species_class(idx) == cnst_spec_class_undefined) then
             call cnst_set_spec_class(idx, cnst_spec_class_gas)
          end if
       end if
    end do

       if (masterproc) write(iulog,9230)
9230   format( // '*** init_aer_modes mode definitions' )
9231   format( 'mode = ', i4, ' = "', a, '"' )
9232   format( 4x, a, 4(1x, i5 ) )
9233   format( 4x, a15, 4x, i7, '="', a, '"' )
9236   format( 4x, a15, i4, i7, '="', a, '"' )
9061   format( '*** subr modesmodal_aero_data_reg - bad ', a /                   &
            5x, 'name, m =  ', a, 5x, i5 )
9062   format( '*** subr modal_aero_data_reg - bad ', a /                       &
            5x, 'name, l, m =  ', a, 5x, 2i5 )
 end subroutine modal_aero_data_reg

!--------------------------------------------------------------
!--------------------------------------------------------------
    subroutine modal_aero_data_init(pbuf2d)

       type(physics_buffer_desc), pointer :: pbuf2d(:,:)

       !--------------------------------------------------------------
       ! ... local variables
       !--------------------------------------------------------------
       integer :: l, m, i, lchnk
       integer :: m_idx, s_idx, ndx

       character(len=3) :: trnum       ! used to hold mode number (as characters)
       integer :: qArrIndex
       integer  :: numaerosols     ! number of bulk aerosols in climate list
       character(len=20) :: bulkname
       complex(r8), pointer  :: refindex_aer_sw(:), &
            refindex_aer_lw(:)
       real(r8), pointer :: qqcw(:,:)
       real(r8), parameter :: huge_r8 = huge(1._r8)
       character(len=*), parameter :: routine='modal_aero_initialize'
       character(len=32) :: spec_type
       character(len=32) :: spec_name
       character(len=1) :: modestr
       integer :: soa_ndx

       !-----------------------------------------------------------------------

       ! Mode specific properties.
       do m = 1, ntot_amode
          call rad_cnst_get_mode_props(0, m, &
             sigmag=sigmag_amode(m), dgnum=dgnum_amode(m), dgnumlo=dgnumlo_amode(m), &
             dgnumhi=dgnumhi_amode(m), rhcrystal=rhcrystal_amode(m), rhdeliques=rhdeliques_amode(m))

          !   compute frequently used parameters: ln(sigmag),
          !   volume-to-number and volume-to-surface conversions, ...
          alnsg_amode(m) = log( sigmag_amode(m) )

          voltonumb_amode(m) = 1._r8 / ( (pi/6._r8)*                            &
             (dgnum_amode(m)**3._r8)*exp(4.5_r8*alnsg_amode(m)**2._r8) )
          voltonumblo_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
             (dgnumlo_amode(m)**3._r8)*exp(4.5_r8*alnsg_amode(m)**2._r8) )
          voltonumbhi_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
             (dgnumhi_amode(m)**3._r8)*exp(4.5_r8*alnsg_amode(m)**2._r8) )

          alnv2n_amode(m)   = log( voltonumb_amode(m) )
          alnv2nlo_amode(m) = log( voltonumblo_amode(m) )
          alnv2nhi_amode(m) = log( voltonumbhi_amode(m) )

       end do
       lptr2_soa_g_amode(:) = -1
       soa_ndx = 0
       do i = 1, pcnst
          if (cnst_name(i)(:4) == 'SOAG') then
             soa_ndx = soa_ndx+1
             lptr2_soa_g_amode(soa_ndx) = i
          endif
       enddo
       if (.not.any(lptr2_soa_g_amode>0)) call endrun('modal_aero_data_init: lptr2_soa_g_amode is not set properly')
       ! Properties of mode specie types.

       !     values from Koepke, Hess, Schult and Shettle, Global Aerosol Data Set 
       !     Report #243, Max-Planck Institute for Meteorology, 1997a
       !     See also Hess, Koepke and Schult, Optical Properties of Aerosols and Clouds (OPAC)
       !     BAMS, 1998.

       !      specrefndxsw(:ntot_aspectype)     = (/ (1.53,  0.01),   (1.53,  0.01),  (1.53,  0.01), &
       !                                           (1.55,  0.01),   (1.55,  0.01),  (1.90, 0.60), &
       !                                           (1.50, 1.0e-8), (1.50, 0.005) /)
       !      specrefndxlw(:ntot_aspectype)   = (/ (2.0, 0.5),   (2.0, 0.5), (2.0, 0.5), &
       !                                           (1.7, 0.5),   (1.7, 0.5), (2.22, 0.73), &
       !                                           (1.50, 0.02), (2.6, 0.6) /)
       !     get refractive indices from phys_prop files

       ! The following use of the rad_constituent interfaces makes the assumption that the
       ! prognostic modes are used in the mode climate (index 0) list.

       if (masterproc) write(iulog,9210)
       do m = 1, ntot_amode
          do l = 1, nspec_amode(m)
             qArrIndex =  lmassptr_amode(l,m)     !index of the species in the state%q array
             call rad_cnst_get_aer_props(0, m, l , &
                  refindex_aer_sw=refindex_aer_sw, &
                  refindex_aer_lw=refindex_aer_lw, &
                  density_aer=specdens_amode(l,m), &
                  hygro_aer=spechygro(l,m)         )

             if ( soa_multi_species ) then
                ! Molecular weight for the species
                specmw_amode(l,m) = cnst_mw(qArrIndex)
             else ! the follow preserves the molecular weights historically used in MAM
                call rad_cnst_get_info(0, m, l, spec_type=spec_type )
                select case( spec_type )
                case('sulfate')
                   if (ntot_amode==7) then
                      specmw_amode(l,m) = 96._r8
                   else
                      specmw_amode(l,m) = 115._r8
                   endif
                case('ammonium')
                   specmw_amode(l,m) = 18._r8
                case('p-organic','s-organic','black-c')
                   specmw_amode(l,m) = 12._r8
                case('seasalt')
                   specmw_amode(l,m) = 58.5_r8
                case('dust')
                   specmw_amode(l,m) = 135._r8
                case default
                   call endrun('modal_aero_data_init: species type not recognized: '//trim(spec_type))
                end select
             endif
 
             if(masterproc) then
                write(iulog,9212) '        name : ', cnst_name(qArrIndex)
                write(iulog,9213) ' density, MW : ', specdens_amode(l,m), specmw_amode(l,m)
                write(iulog,9213) '       hygro : ', spechygro(l,m)
             endif
             do i=1,nswbands
                specrefndxsw(i,l,m)=refindex_aer_sw(i)
                if(masterproc) write(iulog,9213) 'ref index sw    ', (specrefndxsw(i,l,m))
             end do
             do i=1,nlwbands
                specrefndxlw(i,l,m)=refindex_aer_lw(i)
                if(masterproc) write(iulog,9213) 'ref index ir    ', (specrefndxlw(i,l,m))
             end do

          enddo
       enddo

9210   format( // '*** init_aer_modes aerosol species-types' )
9211   format( 'spectype =', i4)
9212   format( 4x, a, 3x, '"', a, '"' )
9213   format( 4x, a, 5(1pe14.5) )

       !   set cnst_name_cw
       call initaermodes_set_cnstnamecw()


       !
       !   set the lptr_so4_a_amode(m), lptr_so4_cw_amode(m), ...
       !
       call initaermodes_setspecptrs

       !
       !   set threshold for reporting negatives from subr qneg3
       !   for aerosol number species set this to
       !      1e3 #/kg ~= 1e-3 #/cm3 for accum, aitken, pcarbon, ufine modes
       !      3e1 #/kg ~= 3e-5 #/cm3 for fineseas and finedust modes 
       !      1e0 #/kg ~= 1e-6 #/cm3 for other modes which are coarse
       !   for other species, set this to zero so that it will be ignored
       !      by qneg3
       !
       if ( masterproc ) write(iulog,'(/a)') &
            'mode, modename_amode, qneg3_worst_thresh_amode'
       qneg3_worst_thresh_amode(:) = 0.0_r8
       do m = 1, ntot_amode
          l = numptr_amode(m)
          if ((l <= 0) .or. (l > pcnst)) cycle

          if      (m == modeptr_accum) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_aitken) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_pcarbon) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8
          else if (m == modeptr_ufine) then
             qneg3_worst_thresh_amode(l) = 1.0e3_r8

          else if (m == modeptr_fineseas) then
             qneg3_worst_thresh_amode(l) = 3.0e1_r8
          else if (m == modeptr_finedust) then
             qneg3_worst_thresh_amode(l) = 3.0e1_r8

          else
             qneg3_worst_thresh_amode(l) = 1.0e0_r8
          end if

          if ( masterproc ) write(iulog,'(i3,2x,a,1p,e12.3)') &
               m, modename_amode(m), qneg3_worst_thresh_amode(l)
       end do

       if (is_first_step()) then
          ! initialize cloud bourne constituents in physics buffer

          do i = 1, pcnst
             do lchnk = begchunk, endchunk
                qqcw => qqcw_get_field(pbuf_get_chunk(pbuf2d,lchnk), i, lchnk, .true.)
                if (associated(qqcw)) then
                   qqcw = 1.e-38_r8
                end if
             end do
          end do
       end if

    end subroutine modal_aero_data_init

!--------------------------------------------------------------
!--------------------------------------------------------------
        subroutine qqcw_set_ptr(index, iptr)
          use cam_abortutils, only : endrun
          

          integer, intent(in) :: index, iptr

          if(index>0 .and. index <= pcnst ) then
             qqcw(index)=iptr
          else
             call endrun('qqcw_set_ptr: attempting to set qqcw pointer already defined')
          end if
        end subroutine qqcw_set_ptr

!--------------------------------------------------------------
!--------------------------------------------------------------
        function qqcw_get_field(pbuf, index, lchnk, errorhandle)
          use cam_abortutils, only : endrun
          use physics_buffer, only : physics_buffer_desc, pbuf_get_field

          integer, intent(in) :: index, lchnk
          real(r8), pointer :: qqcw_get_field(:,:)
          logical, optional :: errorhandle
          type(physics_buffer_desc), pointer :: pbuf(:)

          logical :: error

          nullify(qqcw_get_field)
          error = .false.
          if (index>0 .and. index <= pcnst) then
             if (qqcw(index)>0) then 
                call pbuf_get_field(pbuf, qqcw(index), qqcw_get_field)
             else
                error = .true.
             endif
          else
             error = .true.             
          end if

          if (error .and. .not. present(errorhandle)) then
             call endrun('qqcw_get_field: attempt to access undefined qqcw')
          end if

        end function qqcw_get_field

!----------------------------------------------------------------
!
!   nspec_max = maximum allowable number of chemical species
!       in each aerosol mode
!
!   ntot_amode = number of aerosol modes
!   ( ntot_amode_gchm = number of aerosol modes in gchm
!     ntot_amode_ccm2 = number of aerosol modes to be made known to ccm2
!       These are temporary until multi-mode activation scavenging is going.
!       Until then, ntot_amode is set to either ntot_amode_gchm or
!       ntot_amode_ccm2 depending on which code is active )
!
!   msectional - if positive, moving-center sectional code is utilized,
!       and each mode is actually a section.
!   msectional_concinit - if positive, special code is used to initialize
!       the mixing ratios of all the sections.
!
!   nspec_amode(m) = number of chemical species in aerosol mode m
!   nspec_amode_ccm2(m) = . . .  while in ccm2 code
!   nspec_amode_gchm(m) = . . .  while in gchm code
!   nspec_amode_nontracer(m) = number of "non-tracer" chemical
!       species while in gchm code
!   lspectype_amode(l,m) = species type/i.d. for chemical species l
!       in aerosol mode m.  (1=sulfate, others to be defined)
!   lmassptr_amode(l,m) = gchm r-array index for the mixing ratio
!       (moles-x/mole-air) for chemical species l in aerosol mode m
!       that is in clear air or interstitial air (but not in cloud water)
!   lmassptrcw_amode(l,m) = gchm r-array index for the mixing ratio
!       (moles-x/mole-air) for chemical species l in aerosol mode m
!       that is currently bound/dissolved in cloud water
!   lwaterptr_amode(m) = gchm r-array index for the mixing ratio
!       (moles-water/mole-air) for water associated with aerosol mode m
!       that is in clear air or interstitial air
!   lkohlercptr_amode(m) = gchm r-array index for the kohler "c" parameter
!       for aerosol mode m.  This is defined on a per-dry-particle-mass basis:
!           c = r(i,j,k,lkohlercptr_amode) * [rhodry * (4*pi/3) * rdry^3]
!   numptr_amode(m) = gchm r-array index for the number mixing ratio
!       (particles/mole-air) for aerosol mode m that is in clear air or
!       interstitial are (but not in cloud water).  If zero or negative,
!       then number is not being simulated.
!   ( numptr_amode_gchm(m) = same thing but for within gchm
!     numptr_amode_ccm2(m) = same thing but for within ccm2
!       These are temporary, to allow testing number in gchm before ccm2 )
!   numptrcw_amode(m) = gchm r-array index for the number mixing ratio
!       (particles/mole-air) for aerosol mode m
!       that is currently bound/dissolved in cloud water
!   lsfcptr_amode(m) = gchm r-array index for the surface area mixing ratio
!       (cm^2/mole-air) for aerosol mode m that is in clear air or
!       interstitial are (but not in cloud water).  If zero or negative,
!       then surface area is not being simulated.
!   lsfcptrcw_amode(m) = gchm r-array index for the surface area mixing ratio
!       (cm^2/mole-air) for aerosol mode m that is currently
!       bound/dissolved in cloud water.
!   lsigptr_amode(m) = gchm r-array index for sigmag for aerosol mode m
!       that is in clear air or interstitial are (but not in cloud water).
!       If zero or negative, then the constant sigmag_amode(m) is used.
!   lsigptrcw_amode(m) = gchm r-array index for sigmag for aerosol mode m
!       that is currently bound/dissolved in cloud water.
!       If zero or negative, then the constant sigmag_amode(m) is used.
!   lsigptrac_amode(m) = gchm r-array index for sigmag for aerosol mode m
!       for combined clear-air/interstial plus bound/dissolved in cloud water.
!       If zero or negative, then the constant sigmag_amode(m) is used.
!
!   dgnum_amode(m) = geometric dry mean diameter (m) of the number
!       distribution for aerosol mode m.
!       (Only used when numptr_amode(m) is zero or negative.)
!   dgnumlo_amode(m), dgnumhi_amode(m) = lower and upper limits on the
!       geometric dry mean diameter (m) of the number distribution
!       (Used when mprognum_amode>0, to limit dgnum to reasonable values)
!   sigmag_amode(m) = geometric standard deviation for aerosol mode m
!   sigmaglo_amode(m), sigmaghi_amode(m) = lower and upper limits on the
!       geometric standard deviation of the number distribution
!       (Used when mprogsfc_amode>0, to limit sigmag to reasonable values)
!   alnsg_amode(m) = alog( sigmag_amode(m) )
!   alnsglo_amode(m), alnsghi_amode(m) = alog( sigmaglo/hi_amode(m) )
!   voltonumb_amode(m) = ratio of number to volume for mode m
!   voltonumblo_amode(m), voltonumbhi_amode(m) = ratio of number to volume
!       when dgnum = dgnumlo_amode or dgnumhi_amode, respectively
!   voltosfc_amode(m), voltosfclo_amode(m), voltosfchi_amode(m) - ratio of
!       surface to volume for mode m (like the voltonumb_amode's)
!   alnv2n_amode(m), alnv2nlo_amode(m), alnv2nhi_amode(m) -
!       alnv2n_amode(m) = alog( voltonumblo_amode(m) ), ...
!   alnv2s_amode(m), alnv2slo_amode(m), alnv2shi_amode(m) -
!       alnv2s_amode(m) = alog( voltosfclo_amode(m) ), ...
!   rhcrystal_amode(m) = crystalization r.h. for mode m
!   rhdeliques_amode(m) = deliquescence r.h. for mode m
!   (*** these r.h. values are 0-1 fractions, not 0-100 percentages)
!
!   mcalcwater_amode(m) - if positive, water content for mode m will be
!       calculated and stored in rclm(k,lwaterptr_amode(m)).  Otherwise, no.
!   mprognum_amode(m) - if positive, number mixing-ratio for mode m will
!       be prognosed.  Otherwise, no.
!   mdiagnum_amode(m) - if positive, number mixing-ratio for mode m will
!       be diagnosed and put into rclm(k,numptr_amode(m)).  Otherwise, no.
!   mprogsfc_amode(m) - if positive, surface area mixing-ratio for mode m will
!       be prognosed, and sigmag will vary temporally and spatially.
!       Otherwise, sigmag is constant.
!       *** currently surface area is not prognosed when msectional>0 ***
!
!   ntot_aspectype = overall number of aerosol chemical species defined (over all modes)
!   specdens_amode(l) = dry density (kg/m^3) of aerosol chemical species type l
!   specmw_amode(l) = molecular weight (kg/kmol) of aerosol chemical species type l
!   specname_amode(l) = name of aerosol chemical species type l
!   specrefndxsw(l) = complex refractive index (visible wavelengths)
!                   of aerosol chemical species type l
!   specrefndxlw(l) = complex refractive index (infrared wavelengths)
!                   of aerosol chemical species type l
!   spechygro(l) = hygroscopicity of aerosol chemical species type l
!
!   lptr_so4_a_amode(m), lptr_so4_cw_amode(m) = gchm r-array index for the
!       mixing ratio for sulfate associated with aerosol mode m
!       ("a" and "cw" phases)
!   (similar for msa, oc, bc, nacl, dust)
!
!   modename_amode(m) = character-variable name for mode m,
!       read from mirage2.inp
!   modeptr_accum - mode index for the main accumulation mode
!       if modeptr_accum = 1, then mode 1 is the main accumulation mode,
!       and modename_amode(1) = "accum"
!   modeptr_aitken - mode index for the main aitken mode
!       if modeptr_aitken = 2, then mode 2 is the main aitken mode,
!       and modename_amode(2) = "aitken"
!   modeptr_ufine - mode index for the ultrafine mode
!       if modeptr_ufine = 3, then mode 3 is the ultrafine mode,
!       and modename_amode(3) = "ufine"
!   modeptr_coarseas - mode index for the coarse sea-salt mode
!       if modeptr_coarseas = 4, then mode 4 is the coarse sea-salt mode,
!       and modename_amode(4) = "coarse_seasalt"
!   modeptr_coardust - mode index for the coarse dust mode
!       if modeptr_coardust = 5, then mode 5 is the coarse dust mode,
!       and modename_amode(5) = "coarse_dust"
!
!   specdens_XX_amode = dry density (kg/m^3) of aerosol chemical species type XX
!       where XX is so4, om, bc, dust, seasalt
!       contains same values as the specdens_amode array
!       allows values to be referenced differently
!   specmw_XX_amode = molecular weight (kg/kmol) of aerosol chemical species type XX
!       contains same values as the specmw_amode array
!
!-----------------------------------------------------------------------


!--------------------------------------------------------------
!
! ... aerosol size information for the current chunk
!
!--------------------------------------------------------------
!
!  dgncur = current geometric mean diameters (cm) for number distributions
!  dgncur_a - for unactivated particles, dry
!             (in physics buffer as DGNUM)
!  dgncur_awet - for unactivated particles, wet at grid-cell ambient RH
!             (in physics buffer as DGNUMWET)
!
!  the dgncur are computed from current mass and number
!  mixing ratios in the grid cell, BUT are then adjusted to be within
!  the bounds defined by dgnumlo/hi_amode
!
!  v2ncur = current (number/volume) ratio based on dgncur and sgcur
!              (volume in cm^3/whatever, number in particles/whatever)
!         == 1.0 / ( pi/6 * dgncur**3 * exp(4.5*((log(sgcur))**2)) )
!  v2ncur_a - for unactivated particles
!             (currently just defined locally)
!

     !==============================================================
     subroutine search_list_of_names(                                &
          name_to_find, name_id, list_of_names, list_length )
       !
       !   searches for a name in a list of names
       !
       !   name_to_find - the name to be found in the list  [input]
       !   name_id - the position of "name_to_find" in the "list_of_names".
       !       If the name is not found in the list, then name_id=0.  [output]
       !   list_of_names - the list of names to be searched  [input]
       !   list_length - the number of names in the list  [input]
       !
       character(len=*), intent(in):: name_to_find, list_of_names(:)
       integer, intent(in) :: list_length
       integer, intent(out) :: name_id
       
       integer :: i
       name_id = -999888777
       if (name_to_find .ne. ' ') then
          do i = 1, list_length
             if (name_to_find .eq. list_of_names(i)) then
                name_id = i
                exit
             end if
          end do
       end if
     end subroutine search_list_of_names


     !==============================================================
     subroutine initaermodes_setspecptrs
       !
       !   sets the lptr_so4_a_amode(m), lptr_so4_cw_amode(m), ...
       !       and writes them to iulog
       !   ALSO sets the mode-pointers:  modeptr_accum, modeptr_aitken, ...
       !       and writes them to iulog
       !   ALSO sets values of specdens_XX_amode and specmw_XX_amode
       !       (XX = so4, om, bc, dust, seasalt)
       !
       implicit none

       !   local variables
       integer :: i, l, l2, lmassa, lmassc, m
       character(len=1000) :: msg
       character*8 :: dumname
       character*3 :: tmpch3
       integer, parameter :: init_val=-999888777
       integer :: bc_ndx, soa_ndx, pom_ndx

       !   all processes set the pointers

       modeptr_accum = init_val
       modeptr_aitken = init_val
       modeptr_ufine = init_val
       modeptr_coarse = init_val
       modeptr_pcarbon = init_val
       modeptr_fineseas = init_val
       modeptr_finedust = init_val
       modeptr_coarseas = init_val
       modeptr_coardust = init_val
       do m = 1, ntot_amode
          if (modename_amode(m) .eq. 'accum') then
             modeptr_accum = m
          else if (modename_amode(m) .eq. 'aitken') then
             modeptr_aitken = m
          else if (modename_amode(m) .eq. 'ufine') then
             modeptr_ufine = m
          else if (modename_amode(m) .eq. 'coarse') then
             modeptr_coarse = m
          else if (modename_amode(m) .eq. 'primary_carbon') then
             modeptr_pcarbon = m
          else if (modename_amode(m) .eq. 'fine_seasalt') then
             modeptr_fineseas = m
          else if (modename_amode(m) .eq. 'fine_dust') then
             modeptr_finedust = m
          else if (modename_amode(m) .eq. 'coarse_seasalt') then
             modeptr_coarseas = m
          else if (modename_amode(m) .eq. 'coarse_dust') then
             modeptr_coardust = m
          end if
       end do

       lptr2_pom_a_amode   = init_val
       lptr2_pom_cw_amode  = init_val
       lptr2_soa_a_amode   = init_val
       lptr2_soa_cw_amode  = init_val
       lptr2_bc_a_amode    = init_val
       lptr2_bc_cw_amode   = init_val

       do m = 1, ntot_amode

          lptr_so4_a_amode(m)    = init_val
          lptr_so4_cw_amode(m)   = init_val
          lptr_msa_a_amode(m)    = init_val
          lptr_msa_cw_amode(m)   = init_val
          lptr_nh4_a_amode(m)    = init_val
          lptr_nh4_cw_amode(m)   = init_val
          lptr_no3_a_amode(m)    = init_val
          lptr_no3_cw_amode(m)   = init_val
          lptr_nacl_a_amode(m)   = init_val
          lptr_nacl_cw_amode(m)  = init_val
          lptr_dust_a_amode(m)   = init_val
          lptr_dust_cw_amode(m)  = init_val

          pom_ndx = 0
          soa_ndx = 0
          bc_ndx = 0

          do l = 1, nspec_amode(m)
             lmassa  = lmassptr_amode(l,m)
             lmassc = lmassptrcw_amode(l,m)

             if (lmassa > 0 .and. lmassa <= pcnst) then
                write( msg, '(2a,3(1x,i12),2x,a)' ) &
                     'subr initaermodes_setspecptrs error setting lptr_', &
                     ' - m, l, lmassa, cnst_name = ', m, l, lmassa, cnst_name(lmassa)
             else
                write( msg, '(2a,3(1x,i12),2x,a)' ) &
                     'subr initaermodes_setspecptrs error setting lptr_', &
                     ' - m, l, lmassa, cnst_name = ', m, l, lmassa, 'UNDEF '
                call endrun( trim(msg) )
             end if

             tmpch3 = cnst_name(lmassa)(:3)
             select case (tmpch3)
             case('so4')
                lptr_so4_a_amode(m)  = lmassa
                lptr_so4_cw_amode(m) = lmassc
             case('msa')
                lptr_msa_a_amode(m)  = lmassa
                lptr_msa_cw_amode(m) = lmassc
             case('nh4')
                lptr_nh4_a_amode(m)  = lmassa
                lptr_nh4_cw_amode(m) = lmassc
             case('no3')
                lptr_no3_a_amode(m)  = lmassa
                lptr_no3_cw_amode(m) = lmassc
             case('dst')
                lptr_dust_a_amode(m)  = lmassa
                lptr_dust_cw_amode(m) = lmassc
             case('ncl')
                lptr_nacl_a_amode(m)  = lmassa
                lptr_nacl_cw_amode(m) = lmassc
             case('pom')
                pom_ndx = pom_ndx+1
                lptr2_pom_a_amode(m,pom_ndx)  = lmassa
                lptr2_pom_cw_amode(m,pom_ndx) = lmassc
             case('soa')
                soa_ndx = soa_ndx+1
                lptr2_soa_a_amode(m,soa_ndx)  = lmassa
                lptr2_soa_cw_amode(m,soa_ndx) = lmassc
             case('bc_','bcf','bcb')
                bc_ndx = bc_ndx+1
                lptr2_bc_a_amode(m,bc_ndx)  = lmassa
                lptr2_bc_cw_amode(m,bc_ndx) = lmassc
             case default
                call endrun( trim(msg) )
             end select
          end do ! l
       end do ! m

       specmw_so4_amode = 1.0_r8

       do m = 1, ntot_amode
          do l = 1, nspec_amode(m)
             dumname = trim(adjustl(xname_massptr(l,m)))
             tmpch3  = trim(adjustl(dumname(:3)))
             if(trim(adjustl(tmpch3)) == 'so4' .or. trim(adjustl(tmpch3)) == 'SO4') then
                specmw_so4_amode = specmw_amode(l,m) 
             endif
          enddo
       enddo


       !   masterproc writes out the pointers
       if ( .not. ( masterproc ) ) return

       write(iulog,9230)
       write(iulog,*) 'modeptr_accum    =', modeptr_accum
       write(iulog,*) 'modeptr_aitken   =', modeptr_aitken
       write(iulog,*) 'modeptr_ufine    =', modeptr_ufine
       write(iulog,*) 'modeptr_coarse   =', modeptr_coarse
       write(iulog,*) 'modeptr_pcarbon  =', modeptr_pcarbon
       write(iulog,*) 'modeptr_fineseas =', modeptr_fineseas
       write(iulog,*) 'modeptr_finedust =', modeptr_finedust
       write(iulog,*) 'modeptr_coarseas =', modeptr_coarseas
       write(iulog,*) 'modeptr_coardust =', modeptr_coardust
       
       dumname = 'none'
       write(iulog,9240)
       write(iulog,9000) 'sulfate    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_so4_a_amode(m), lptr_so4_cw_amode(m),  'so4' )
       end do

       write(iulog,9000) 'msa        '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_msa_a_amode(m), lptr_msa_cw_amode(m),  'msa' )
       end do

       write(iulog,9000) 'ammonium   '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_nh4_a_amode(m), lptr_nh4_cw_amode(m),  'nh4' )
       end do

       write(iulog,9000) 'nitrate    '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_no3_a_amode(m), lptr_no3_cw_amode(m),  'no3' )
       end do

       write(iulog,9000) 'p-organic  '
       do m = 1, ntot_amode
          do i = 1, npoa
             write(dumname,'(a,i2.2)') 'pom', i
             call initaermodes_setspecptrs_write2b( m, &
                  lptr2_pom_a_amode(m,i), lptr2_pom_cw_amode(m,i),  dumname(1:5) )
          end do
       end do

       write(iulog,9000) 's-organic  '
       do m = 1, ntot_amode
          do i = 1, nsoa
             write(dumname,'(a,i2.2)') 'soa', i
             call initaermodes_setspecptrs_write2b( m, &
                  lptr2_soa_a_amode(m,i), lptr2_soa_cw_amode(m,i), dumname(1:5) )
          end do
       end do
       do i = 1, nsoa
          l = lptr2_soa_g_amode(i)
          write(iulog,'(i4,2x,i12,2x,a,20x,a,i2.2,a)') i, l, cnst_name(l), 'lptr2_soa', i, '_g'
       end do

       write(iulog,9000) 'black-c    '    
       do m = 1, ntot_amode
          do i = 1, nbc
             write(dumname,'(a,i2.2)') 'bc', i
             call initaermodes_setspecptrs_write2b( m, &
                  lptr2_bc_a_amode(m,i), lptr2_bc_cw_amode(m,i), dumname(1:5) )
          end do
       end do

       write(iulog,9000) 'seasalt   '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_nacl_a_amode(m), lptr_nacl_cw_amode(m),  'nacl' )
       end do

       write(iulog,9000) 'dust       '
       do m = 1, ntot_amode
          call initaermodes_setspecptrs_write2( m,                    &
               lptr_dust_a_amode(m), lptr_dust_cw_amode(m),  'dust' )
       end do

9000   format( a )
9230   format(                                                         &
            / 'mode-pointer output from subr initaermodes_setspecptrs' )
9240   format(                                                         &
            / 'species-pointer output from subr initaermodes_setspecptrs' / &
            'mode', 12x, 'id  name_a  ', 12x, 'id  name_cw' )

       return
     end subroutine initaermodes_setspecptrs


     !==============================================================
     subroutine initaermodes_setspecptrs_write2(                     &
          m, laptr, lcptr, txtdum )
       !
       !   does some output for initaermodes_setspecptrs

       use constituents, only: pcnst, cnst_name

       implicit none

       !   subr arguments
       integer, intent(in) :: m, laptr, lcptr
       character*(*), intent(in) ::  txtdum

       !   local variables
       character*8 dumnamea, dumnamec

       dumnamea = 'none'
       dumnamec = 'none'
       if (laptr > pcnst .or. lcptr > pcnst ) then
          call endrun('initaermodes_setspecptrs_write2: ERROR')
       endif
       if (laptr .gt. 0) dumnamea = cnst_name(laptr)
       if (lcptr .gt. 0) dumnamec = cnst_name(lcptr)
       write(iulog,9241) m, laptr, dumnamea, lcptr, dumnamec, txtdum

9241   format( i4, 2( 2x, i12, 2x, a ),                                &
            4x, 'lptr_', a, '_a/cw_amode' )

       return
     end subroutine initaermodes_setspecptrs_write2


     !==============================================================
     subroutine initaermodes_setspecptrs_write2b( &
          m, laptr, lcptr, txtdum )
       !
       !   does some output for initaermodes_setspecptrs

       use constituents, only: pcnst, cnst_name

       implicit none

       !   subr arguments
       integer, intent(in) :: m, laptr, lcptr
       character*(*), intent(in) :: txtdum

       !   local variables
       character*8 dumnamea, dumnamec

       dumnamea = 'none'
       dumnamec = 'none'
       if (laptr .gt. 0) dumnamea = cnst_name(laptr)
       if (lcptr .gt. 0) dumnamec = cnst_name(lcptr)
       write(iulog,9241) m, laptr, dumnamea, lcptr, dumnamec, txtdum

9241   format( i4, 2( 2x, i12, 2x, a ),                                &
            4x, 'lptr2_', a, '_a/cw_amode' )

       return
     end subroutine initaermodes_setspecptrs_write2b

     !==============================================================
     subroutine initaermodes_set_cnstnamecw
       !
       !   sets the cnst_name_cw
       !
       use constituents, only: pcnst, cnst_name
       implicit none

       !   subr arguments (none)

       !   local variables
       integer j, l, la, lc, ll, m

       !   set cnst_name_cw
       cnst_name_cw = ' '
       do m = 1, ntot_amode
          do ll = 0, nspec_amode(m)
             if (ll == 0) then
                la = numptr_amode(m)
                lc = numptrcw_amode(m)
             else
                la = lmassptr_amode(ll,m)
                lc = lmassptrcw_amode(ll,m)
             end if
             if ((la < 1) .or. (la > pcnst) .or.   &
                  (lc < 1) .or. (lc > pcnst)) then
                write(*,'(/2a/a,5(1x,i10))')   &
                     '*** initaermodes_set_cnstnamecw error',   &
                     ' -- bad la or lc',   &
                     '    m, ll, la, lc, pcnst =', m, ll, la, lc, pcnst
                call endrun( '*** initaermodes_set_cnstnamecw error' )
             end if
             do j = 2, len( cnst_name(la) ) - 1
                if (cnst_name(la)(j:j+1) == '_a') then
                   cnst_name_cw(lc) = cnst_name(la)
                   cnst_name_cw(lc)(j:j+1) = '_c'
                   exit
                else if (cnst_name(la)(j:j+1) == '_A') then
                   cnst_name_cw(lc) = cnst_name(la)
                   cnst_name_cw(lc)(j:j+1) = '_C'
                   exit
                end if
             end do
             if (cnst_name_cw(lc) == ' ') then
                write(*,'(/2a/a,3(1x,i10),2x,a)')   &
                     '*** initaermodes_set_cnstnamecw error',   &
                     ' -- bad cnst_name(la)',   &
                     '    m, ll, la, cnst_name(la) =',   &
                     m, ll, la, cnst_name(la)
                call endrun( '*** initaermodes_set_cnstnamecw error' )
             end if
          end do   ! ll = 0, nspec_amode(m)
       end do   ! m = 1, ntot_amode

       if ( masterproc ) then
          write(*,'(/a)') 'l, cnst_name(l), cnst_name_cw(l)'
          do l = 1, pcnst
             write(*,'(i4,2(2x,a))') l, cnst_name(l), cnst_name_cw(l)
          end do
       end if

       return
     end subroutine initaermodes_set_cnstnamecw

  end module modal_aero_data


