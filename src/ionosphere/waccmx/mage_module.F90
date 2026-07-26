module mage_module
  !
  !   Module used to exchange data back and forth with MAGE
  !   WACCM-X Receives (POT, mean energy, energy flux)
  !   MAGE Receives (Pedersen Conductance, Hall Conductance)
  !

  use shr_kind_mod,   only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,    only: iulog
  use spmd_utils,     only: masterproc
  use edyn_maggrid,   only: nmlat, nmlonp1
  use edyn_maggrid,   only: gmlon     ! magnetic latitudes (nmlat) (radians)
  use edyn_maggrid,   only: gmlat     ! magnetic longtitudes (nmlonp1) (radians)
  use edyn_mpi,       only: ntask, mytid
  use edyn_params,    only: pi, dtr, rtd
  use edynamo,        only: azigm1, azigm2 ! Hall and Ped conductances
  use cam_abortutils, only: endrun

  implicit none

  private
  public :: mage_init
  public :: mage_advance

  integer, parameter :: &! For MPI-based coupling tag
       myAppId = 67,       &! waccmxID
       voltId  = 116,      &! voltronID
       gamId   = 45,       &! gameraID
       rcmId   = 34,       &! rcmID
       mageId  = 26,       &! mageID
       hidraId = 40,       &! hidraID
       hidraNId = 54,      &! hidraNID
       hidraSId = 59,      &! hidraSID
       tiegcmId = 57        ! tiegcmID

  ! Change the following parameters if adding more variables
  integer, parameter :: &
       nmixinapex = 3, &! APEX Potential, Flux (for auroral bc)
       nmixingeo  = 0, &! GEO Flux, Energy
       nmixoutapex= 2, &! APEX Pedersen, Hall
       nmixoutgeo = 0, &! No GEO exports
       nhoutvar   = 9, &! TN, UN, VN, OMEGA, O2, O1, NO, Z, HE
       nhinvar    = 0   !

  integer,dimension(:),allocatable :: IAm

  integer :: &
       CplComm, CplCommSize, CplRank, & ! Global coupling communicator
       mixCplRank,hidraCplRank,       & ! Direct mix coupling communicator
       hidraNCplRank,hidraSCplRank      ! Direct hidra coupling communicator

  integer :: nlatp2=0, nlonp1=0 ! Needed for geographic conversions later

contains

  !-----------------------------------------------------------------------

  subroutine mage_advance(iyear, imo, iday, iutsec, &
                          phihm, mage_efxm, mage_kevm)

    use cam_history_support, only: fillvalue

    !
    !    Read MAGE outputs from mage_ncfile file, returning electric potential,
    !    auroral mean energy and energy flux at current date and time,
    !    and the data is linearly interpolated to the model time
    !
    !
    !    Args:

    integer,  intent(in)    :: iyear
    integer,  intent(in)    :: imo
    integer,  intent(in)    :: iday
    integer,  intent(in)    :: iutsec
    real(r8), intent(out)   :: phihm(nmlonp1,nmlat)
    real(r8), intent(out)   :: mage_efxm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: mage_kevm(nmlonp1,nmlat) ! on geomag grid

    if (mytid .eq. 0) then
       write(iulog,'(A,I4,1X,I2.2,1X,I2.2,1X,I2,":",I2.2,":",I2.2)') "WCMX Coupling at ", &
             iyear, imo, iday, iutsec/3600, mod(iutsec,3600)/60, mod(iutsec,60)
    endif

    phihm = 0._r8
    mage_efxm = 0._r8
    mage_kevm = 0._r8

    call update_mage(phihm, mage_efxm, mage_kevm)

  end subroutine mage_advance

  !-----------------------------------------------------------------------

  subroutine mage_init()
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_ERROR

    integer :: i,ierr
    character(len=*), parameter :: preface = 'mage_module::mage_init : WACCMX send to REMIX : '

    ! Initialize the MPI Coupling interface

    call mp_coupling()
    if (mytid .eq. 0) then
       i = 0
       if (masterproc) write(iulog,*) "WCMX: Start init_mpi_remix"

       i = i + 1
       call mpi_send(nmlat, 1, MPI_INTEGER, mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
       if (ierr == MPI_ERROR) then
          call endrun(preface//'nmlat')
       end if

       i = i + 1
       call mpi_send(nmlonp1, 1, MPI_INTEGER, mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
       if (ierr == MPI_ERROR) then
          call endrun(preface//'nmlonp1')
       end if

       i = i + 1
       call mpi_send(gmlat, nmlat, MPI_DOUBLE_PRECISION, mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
       if (ierr == MPI_ERROR) then
          call endrun(preface//'gmlat')
       end if

       i = i + 1
       ! Making a note here. For whatever reason, WACCM-X's conductance is flipped 180 degrees compared to TIEGCM.
       !I'm going to hack it here so that it's straight up 180 degrees shifted but I do not know why and what the consequences are.
       call mpi_send(gmlon, nmlonp1, MPI_DOUBLE_PRECISION, mixCplRank, (myAppId+voltId)*100+i, CplComm, ierr)
       if (ierr == MPI_ERROR) then
          call endrun(preface//'gmlon, nmlonp1')
       end if

       if (masterproc) write(iulog,*)  "WCMX: Done init_mpi_remix"
    endif

  end subroutine mage_init

  !-----------------------------------------------------------------------

  subroutine mp_coupling()
    use mpi, only: mpi_comm_size, mpi_comm_free, MPI_COMM_WORLD, MPI_ERROR, MPI_COMM_NULL, MPI_INTEGER

    integer :: ierr,color,i
    integer :: tmpComm
    character(len=*), parameter :: preface = 'mage_module::mp_coupling : '

    ! Create a second communicator to transfer data between TIEGCM and MAGE
    ! This communicator only includes the root processes
    ! TIEGCM root sends/receives data from/to MAGE root

    color = mageId
    tmpComm = MPI_COMM_WORLD

    if (masterproc) write(iulog,*) "WACCMX: Get Coupling Comm ",ntask

    call mp_get_coupling_comm(tmpComm,color, 0, CplComm)

    call mpi_comm_size(CplComm,CplCommSize,ierr)
    if (ierr == MPI_ERROR) then
       write(6,"('>>> Error from mpi_comm_size: ierr=',i4)") ierr
       call endrun(preface//'mpi_comm_size')
    endif

    ! At most two processes will register in CplComm (TIEGCM root, REMIX root)
    ! Therefore CplCommSize is either 1 or 2
    if (CplCommSize == 1) then

       ! Only one process registered in CplComm (TIEGCM root)
       ! No coupling will take place, free up the resources
       call mpi_comm_free(CplComm,ierr)
       if (ierr == MPI_ERROR) then
          write(6,"('>>> Error from mpi_comm_free: ierr=',i4)") ierr
          call endrun(preface//'mpi_comm_free')
       endif
       CplComm = MPI_COMM_NULL
    else

       ! There is another process registered in CplComm (REMIX root)
       call mpi_comm_rank(CplComm,CplRank,ierr)
       !write(*,*) "W cplrank: ",CplRank
       if (ierr == MPI_ERROR) then
          write(6,"('>>> Error from mpi_comm_rank: ierr=',i4)") ierr
          call endrun(preface//'mpi_comm_rank')
       endif

       if (.not.allocated(IAm)) allocate(IAm(CplCommSize))

       if (mytid .eq. 0) then
          IAm(CplRank+1) = myAppId*100+1
       else
          IAm(CplRank+1) = myAppId*100
       endif

       do i=1,CplCommSize
          call mpi_bcast(IAm(i), 1, MPI_INTEGER, i-1, CplComm, ierr)
          if (ierr == MPI_ERROR) then
             write(6,"('>>> Error from mpi_bcast: ierr=',i4)") ierr
             call endrun(preface//'mpi_bcast')
          endif
       enddo

       do i=1,CplCommSize
          ! Assign rank if match
          select case (IAm(i))
          case (voltId)
             mixCplRank = i-1
             if (mytid .eq. 0) write(iulog,*) "W coupling to remix",mixCplRank
          case (gamId)
             if (mytid .eq. 0) write(iulog,*) "W not coupling to Gam yet"
          case (rcmId)
             if (mytid .eq. 0) write(iulog,*) "W not coupling to RCM yet"
          case (hidraNId)
             hidraNCplRank = i-1
             if (mytid .eq. 0) write(iulog,*) "W coupling to hidraN"
          case (hidraSId)
             hidraSCplRank = i-1
             if (mytid .eq. 0) write(iulog,*) "W coupling to hidraS"
          case (hidraId)
             hidraCplRank = i-1
             if (mytid .eq. 0) write(iulog,*) "W coupling to hidra"
          case (myAppId,myAppId*100,myAppId*100+1)
             if (mytid .eq. 0 .and. IAm(i) .eq. myAppId*100+1) write(iulog,*) "W is W", i-1
             IAm(i) = myAppId
          case default
             if (mytid .eq. 0) write(iulog,*) "W does not know about", &
                  " this Coupling ID: ", IAm(i)
          end select
       enddo
    endif

    if (masterproc) write(iulog,'(A,I0,A,I0,A,I0,A)') "W COUPLING to ",CplCommSize, &
                          " Models on ",CplRank," Rank on ",CplComm," Comm"

  end subroutine mp_coupling

  !-----------------------------------------------------------------------

  subroutine mp_get_coupling_comm(couplingPool,appId,key,coupledComm)
    use mpi,only: MPI_comm_rank, MPI_comm_split, MPI_ERROR, MPI_IN_PLACE, MPI_INTEGER, MPI_MAX

    integer, intent(inout) :: couplingPool
    integer, intent(in) :: appId, key
    integer, intent(inout) :: coupledComm

    integer :: ierr, myRank, appIdCpy
    character(len=*), parameter :: preface = 'mage_module::mp_get_coupling_comm : '

    appIdCpy = appId

    ! mpi_bcast doesn't interact well with intent(in)

    ! tell everyone I'm the broadcasting root
    ! broadcast which app I'm creating a communicator with, split with it, and then
    ! create a smaller pool that excludes that app
    call MPI_comm_rank(couplingPool, myRank, ierr)
    if (ierr == MPI_ERROR) then
       write(6,"('>>> Error from MPI_comm_rank: ierr=',i4)") ierr
       call endrun(preface//'MPI_comm_rank')
    endif

    call MPI_Allreduce(MPI_IN_PLACE, myRank, 1, MPI_INTEGER, MPI_MAX, couplingPool, ierr)
    if (ierr == MPI_ERROR) then
       write(6,"('>>> Error from MPI_Allreduce: ierr=',i4)") ierr
       call endrun(preface//'MPI_Allreduce')
    endif

    ! This Bcast is causing a lot of issues. I don't know if this is needed or
    ! if it will cause problems for voltron and other models. The behavior here is odd.
    call MPI_Bcast(appIdCpy, 1, MPI_INTEGER, myRank, couplingPool, ierr)
    if (ierr == MPI_ERROR) then
       write(6,"('>>> Error from MPI_Bcast: ierr=',i4)") ierr
       call endrun(preface//'MPI_Bcast')
    endif

    call MPI_comm_split(couplingPool, appIdCpy, key, coupledComm, ierr)
    if (ierr == MPI_ERROR) then
       write(6,"('>>> Error from MPI_comm_split: ierr=',i4)") ierr
       call endrun(preface//'MPI_comm_split')
    endif

    ! key is never used when making the exclusion pool, 0 is used to preserve order
    call MPI_comm_split(couplingPool, myAppId, 0, couplingPool, ierr)
    if (ierr == MPI_ERROR) then
       write(6,"('>>> Error from MPI_comm_split: ierr=',i4)") ierr
       call endrun(preface//'MPI_comm_split')
    endif

  end subroutine mp_get_coupling_comm

  !-----------------------------------------------------------------------

  subroutine update_mage(phihm, mage_efxm, mage_kevm)

    real(r8), intent(out)   :: phihm(nmlonp1,nmlat)
    real(r8), intent(out)   :: mage_efxm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: mage_kevm(nmlonp1,nmlat) ! on geomag grid

    call import_mage(phihm, mage_efxm, mage_kevm)
    call export_mage()

  end subroutine update_mage

  !-----------------------------------------------------------------------

  subroutine import_mage(phihm, mage_efxm, mage_kevm)

    real(r8), intent(out)   :: phihm(nmlonp1,nmlat)
    real(r8), intent(out)   :: mage_efxm(nmlonp1,nmlat) ! on geomag grid
    real(r8), intent(out)   :: mage_kevm(nmlonp1,nmlat) ! on geomag grid

    ! 1. Receive arrays (pot, eng, flx) from REMIX
    ! 2. Broadcast to all MPI tasks
    ! 3. Set periodic boundaries for gpot, geng, gflx, gpotm
    ! 4. Clean up pole values for the dynamo solver

    integer :: i,j
    real(r8),dimension(nlatp2,nlonp1,nmixingeo) :: gvar2d
    real(r8),dimension(nmlat,nmlonp1,nmixinapex) :: avar2d

    avar2d = 0.0_r8
    gvar2d = 0.0_r8

    call import_remix(avar2d,gvar2d)

    ! Unit conversion
    avar2d(:,:,1) = avar2d(:,:,1)*1e3_r8 ! potential: kV -> V

    ! Process the imported data

    if (nmixoutapex .ne. 0) then
       do j=1,nmlat
          do i=1,nmlonp1
             phihm(i,j) = avar2d(j,i,1)
             !if (avar2d(j,i,2) .ne. avar2d(j,i,2)) write(*,*) "AVAR2?? ",avar2d(j,i,2),j,i
             !if (avar2d(j,i,3) .ne. avar2d(j,i,3)) write(*,*) "AVAR3?? ",avar2d(j,i,3),j,i
             !if (avar2d(j,i,2) .lt. 0) write(*,*) "AVAR4?? ",avar2d(j,i,2),j,i
             !if (avar2d(j,i,3) .lt. 0) write(*,*) "AVAR5?? ",avar2d(j,i,3),j,i
             mage_kevm(i,j) = max(avar2d(j,i,2),0.5_r8) ! keV
             mage_efxm(i,j) = max(avar2d(j,i,2) * avar2d(j,i,3),0.01_r8)*1.602e-9_r8 ! Convert from keV/cm^s to mW/m^2
          enddo
       enddo
    endif

  end subroutine import_mage

  !-------------------------------------------------------------------

  subroutine import_remix(avar2d,gvar2d)
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_ERROR

    integer :: ierr
    real(r8),dimension(:,:,:) :: avar2d,gvar2d
    character(len=*), parameter :: preface = 'mage_module::import_remix : '

    if (nmixingeo .ne. 0) then
       call MPI_BCAST(gvar2d, nlatp2*nlonp1*nmixingeo, MPI_DOUBLE_PRECISION, mixCplRank, CplComm, ierr)
       if (ierr == MPI_ERROR) then
          write(6,"('>>> Error from MPI_BCAST: ierr=',i4)") ierr
          call endrun(preface//'MPI_BCAST')
       endif
    endif
    if (nmixinapex .ne. 0) then
       call MPI_BCAST(avar2d, nmlat*nmlonp1*nmixinapex, MPI_DOUBLE_PRECISION, mixCplRank, CplComm, ierr)
       if (ierr == MPI_ERROR) then
          write(6,"('>>> Error from MPI_BCAST: ierr=',i4)") ierr
          call endrun(preface//'MPI_BCAST')
       endif
       if (masterproc) then
          write(iulog,*) "WCMX Done Import2: ",mixCplRank,nmlat*nmlonp1*nmixinapex
          write(iulog,*) "WCMX Import Check1: ",minval(avar2d(:,:,1)),maxval(avar2d(:,:,1))
          write(iulog,*) "WCMX Import Check2: ",minval(avar2d(:,:,2)),maxval(avar2d(:,:,2))
          write(iulog,*) "WCMX Import Check3: ",minval(avar2d(:,:,3)),maxval(avar2d(:,:,3))
       endif
    endif

  end subroutine import_remix

  !-----------------------------------------------------------------------

  subroutine export_mage()

    logical :: hidra_prep
    real(r8),dimension(:,:,:), allocatable :: mixgeoout,mixapexout
    integer :: nreq,i
    character(len=*), parameter :: preface = 'mage_module::export_mage : '

    ! Prepare the export data
    hidra_prep = .false.
    nreq = 0
    !if (mytid == 0) write(iulog,*) "W Starting Export Prep"

    do i=1,CplCommSize
       ! Skip Self
       if (i == CplRank+1) continue
       ! Assign rank if match
       select case (IAm(i))
       case (voltId)
          call prep_export_remix(mixapexout,mixgeoout)
          nreq = nreq + 1
       case (gamId)
          !write(*,*) "T not coupling to Gam yet"
       case (rcmId)
          !write(*,*) "T not coupling to RCM yet"
       case (hidraNId)
          if (.not.hidra_prep) then
             !write(*,*) "Prep H Export"
             !call prep_export_hidra(hidraout)
             hidra_prep = .true.
          endif
          nreq = nreq + 1
       case (hidraSId)
          if (.not.hidra_prep) then
             !write(*,*) "Prep H Export"
             !call prep_export_hidra(hidraout)
             hidra_prep = .true.
          endif
          nreq = nreq + 1
       case (hidraId)
          if (.not.hidra_prep) then
             !call prep_export_hidra(hidraout)
             hidra_prep = .true.
          endif
          nreq = nreq + 1
       case (myAppId)
          !write(*,*) "T is T"
       case default
          if (IAm(i) .eq. 0) cycle
          if (mytid == 0) &
               write(iulog,*) "W does not know about this Coupling ID: ", IAm(i)
       end select
    enddo

    ! Allocate request array here

    ! Send the export data
    do i=1,CplCommSize
       ! Assign rank if match
       select case (IAm(i))
       case (voltId)
          call export_remix(mixapexout,mixgeoout)
       case (gamId)
          !write(*,*) "T not coupling to Gam yet"
       case (rcmId)
          !write(*,*) "T not coupling to RCM yet"
       case (hidraNId)
          !write(*,*) "Export to HidraN"
          !call export_hidra(hidraout,myAppId+hidraNId,hidraNCplRank)
       case (hidraSId)
          !write(*,*) "Export to HidraS"
          !call export_hidra(hidraout,myAppId+hidraSId,hidraSCplRank)
       case (hidraId)
          !call export_hidra(hidraout,myAppId+hidraId,hidraCplRank)
       case (myAppId)
          !write(*,*) "T is T"
       case default
          if (IAm(i) .eq. 0) cycle
          if (mytid == 0) &
               write(iulog,*) "W does not know about this Coupling ID: ", IAm(i)
       end select
    enddo
  end subroutine export_mage

  !-----------------------------------------------------------------------

  subroutine export_remix(avar2d,gvar2d)
    use mpi, only: MPI_DOUBLE_PRECISION, MPI_ERROR

    real(r8),dimension(:,:,:) :: avar2d,gvar2d
    integer :: ierr
    character(len=*), parameter :: preface = 'mage_module::export_remix : '

    ! Export the data
    if (mytid .eq. 0) then
       !write(iulog,*) "WCMX Waiting to export1:",mixCplRank,(myAppId+voltId)*100,nlatp2*nlonp1*nmixoutgeo
       if ( nmixoutgeo .ne. 0) then
          call mpi_send(gvar2d, nlatp2*nlonp1*nmixoutgeo, MPI_DOUBLE_PRECISION, mixCplRank, &
               (myAppId+voltId)*100, CplComm, ierr)
          if (ierr == MPI_ERROR) then
             write(6,"('>>> Error from mpi_send: ierr=',i4)") ierr
             call endrun(preface//'mpi_comm_size')
          endif
       endif

       !write(iulog,*) "WCMX Waiting to export2:",mixCplRank,(myAppId+voltId)*100,nmlat*nmlonp1*nmixoutapex
       !write(iulog,*) "WCMX export3: ",shape(avar2d)
       if ( nmixoutapex .ne. 0) then
          call mpi_send(avar2d, nmlat*nmlonp1*nmixoutapex, MPI_DOUBLE_PRECISION, mixCplRank, &
               (myAppId+voltId)*100, CplComm, ierr)
          if (ierr == MPI_ERROR) then
             write(6,"('>>> Error from mpi_send: ierr=',i4)") ierr
             call endrun(preface//'mpi_comm_size')
          endif
       endif
       !write(iulog,*) "WCMX Done to export2:"

    endif

  end subroutine export_remix
  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  subroutine prep_export_remix(avar2d,gvar2d)
    use edyn_mpi,only: mlat0,mlat1,mlon0,mlon1, mp_gather_edyn

    ! Local variables ..................................................
    integer :: i,j,v
    real(r8),dimension(:,:,:),allocatable :: avar2d,gvar2d

    real(r8) :: amsub(mlon0:mlon1,mlat0:mlat1,nmixoutapex)
    real(r8) :: amglb(nmlonp1,nmlat,nmixoutapex)

    ! Begin ............................................................

    ! Prepare data for export
    amsub(:,:,1) = azigm2(mlon0:mlon1,mlat0:mlat1)  ! Ped
    amsub(:,:,2) = azigm1(mlon0:mlon1,mlat0:mlat1) ! Hall

    where(amsub < 0.2_r8) amsub = 0.2_r8

    call mp_gather_edyn(amsub,mlon0,mlon1,mlat0,mlat1, &
         amglb,nmlonp1,nmlat,nmixoutapex)


    if (mytid .eq. 0) then
       !write(iulog,*) "WCMX: Start Preparing Export"

       if (nmixoutapex .ne. 0) then
          if (.not. allocated(avar2d)) then
             allocate(avar2d(nmlat,nmlonp1,nmixoutapex))
          endif

          do v=1,nmixoutapex
             do j = 1,nmlat
                do i = 1,nmlonp1
                   avar2d(j,i,v) = amglb(i,j,v)
                enddo
             enddo
          enddo

       endif

       if (nmixoutgeo .ne. 0) then
          ! Prepare geographic output arrays
          if (.not. allocated(gvar2d)) then
             allocate(gvar2d(nlatp2,nlonp1,nmixoutgeo))
          endif

          gvar2d(:,:,1) = 0._r8 ! put stuff here

       endif
    endif

    !if (mytid .eq. 0) write(iulog,*) "WCMX: Done Preparing Export"

  end subroutine prep_export_remix

end module mage_module
