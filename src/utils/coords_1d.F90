module coords_1d

! This module defines the Coords1D type, which is intended to to cache
! commonly used information derived from a collection of sets of 1-D
! coordinates.

use shr_kind_mod, only: r8 => shr_kind_r8

implicit none
private
save

public :: Coords1D

type :: Coords1D
   ! Number of sets of coordinates in the object.
   integer :: n = 0
   ! Number of coordinates in each set.
   integer :: d = 0

   ! All fields below will be allocated with first dimension "n".
   ! The second dimension is d+1 for ifc, d for mid, del, and rdel, and
   ! d-1 for dst and rdst.

   ! Cell interface coordinates.
   real(r8), allocatable :: ifc(:,:)
   ! Coordinates at cell mid-points.
   real(r8), allocatable :: mid(:,:)
   ! Width of cells.
   real(r8), allocatable :: del(:,:)
   ! Distance between cell midpoints.
   real(r8), allocatable :: dst(:,:)
   ! Reciprocals: 1/del and 1/dst.
   real(r8), allocatable :: rdel(:,:)
   real(r8), allocatable :: rdst(:,:)
 contains
   procedure :: section
   procedure :: finalize
end type Coords1D

interface Coords1D
   module procedure new_Coords1D_from_fields
   module procedure new_Coords1D_from_int
end interface

contains

! Constructor to create an object from existing data.
function new_Coords1D_from_fields(ifc, mid, del, dst, &
     rdel, rdst) result(coords)
  real(r8), USE_CONTIGUOUS intent(in) :: ifc(:,:)
  real(r8), USE_CONTIGUOUS intent(in) :: mid(:,:)
  real(r8), USE_CONTIGUOUS intent(in) :: del(:,:)
  real(r8), USE_CONTIGUOUS intent(in) :: dst(:,:)
  real(r8), USE_CONTIGUOUS intent(in) :: rdel(:,:)
  real(r8), USE_CONTIGUOUS intent(in) :: rdst(:,:)
  type(Coords1D) :: coords

  coords = allocate_coords(size(ifc, 1), size(ifc, 2) - 1)

  coords%ifc = ifc
  coords%mid = mid
  coords%del = del
  coords%dst = dst
  coords%rdel = rdel
  coords%rdst = rdst

end function new_Coords1D_from_fields

! Constructor if you only have interface coordinates; derives all the other
! fields.
function new_Coords1D_from_int(ifc) result(coords)
  real(r8), USE_CONTIGUOUS intent(in) :: ifc(:,:)
  type(Coords1D) :: coords

  coords = allocate_coords(size(ifc, 1), size(ifc, 2) - 1)

  coords%ifc = ifc
  coords%mid = 0.5_r8 * (ifc(:,:coords%d)+ifc(:,2:))
  coords%del = coords%ifc(:,2:) - coords%ifc(:,:coords%d)
  coords%dst = coords%mid(:,2:) - coords%mid(:,:coords%d-1)
  coords%rdel = 1._r8/coords%del
  coords%rdst = 1._r8/coords%dst

end function new_Coords1D_from_int

! Create a new Coords1D object that is a subsection of some other object,
! e.g. if you want only the first m coordinates, use d_bnds=[1, m].
!
! Originally this used pointers, but it was found to actually be cheaper
! in practice just to make a copy, especially since pointers can impede
! optimization.
function section(self, n_bnds, d_bnds)
  class(Coords1D), intent(in) :: self
  integer, intent(in) :: n_bnds(2), d_bnds(2)
  type(Coords1D) :: section

  section = allocate_coords(n_bnds(2)-n_bnds(1)+1, d_bnds(2)-d_bnds(1)+1)

  section%ifc = self%ifc(n_bnds(1):n_bnds(2),d_bnds(1):d_bnds(2)+1)
  section%mid = self%mid(n_bnds(1):n_bnds(2),d_bnds(1):d_bnds(2))
  section%del = self%del(n_bnds(1):n_bnds(2),d_bnds(1):d_bnds(2))
  section%dst = self%dst(n_bnds(1):n_bnds(2),d_bnds(1):d_bnds(2)-1)
  section%rdel = self%rdel(n_bnds(1):n_bnds(2),d_bnds(1):d_bnds(2))
  section%rdst = self%rdst(n_bnds(1):n_bnds(2),d_bnds(1):d_bnds(2)-1)

end function section

! Quick utility to get allocate each array with the correct size.
function allocate_coords(n, d) result(coords)
  integer, intent(in) :: n, d
  type(Coords1D) :: coords

  coords%n = n
  coords%d = d

  allocate(coords%ifc(coords%n,coords%d+1))
  allocate(coords%mid(coords%n,coords%d))
  allocate(coords%del(coords%n,coords%d))
  allocate(coords%dst(coords%n,coords%d-1))
  allocate(coords%rdel(coords%n,coords%d))
  allocate(coords%rdst(coords%n,coords%d-1))

end function allocate_coords

! Deallocate and reset to initial state.
subroutine finalize(self)
  class(Coords1D), intent(inout) :: self

  self%n = 0
  self%d = 0

  call guarded_deallocate(self%ifc)
  call guarded_deallocate(self%mid)
  call guarded_deallocate(self%del)
  call guarded_deallocate(self%dst)
  call guarded_deallocate(self%rdel)
  call guarded_deallocate(self%rdst)

contains

  subroutine guarded_deallocate(array)
    real(r8), allocatable :: array(:,:)

    if (allocated(array)) deallocate(array)

  end subroutine guarded_deallocate

end subroutine finalize

end module coords_1d
