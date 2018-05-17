!-------------------------------------------------------------------
! shared utilities for diagnostics summations
!-------------------------------------------------------------------
module sums_utils
  use shr_kind_mod, only : r8 => shr_kind_r8
  use shr_kind_mod, only : CL => SHR_KIND_CL
  use shr_kind_mod, only : CXX => SHR_KIND_CXX

  implicit none

  !-------------------------------------------------------------------
  ! object which holds the terms of a summation
  !-------------------------------------------------------------------  
  type sums_grp_t
    character(len=64) :: name
    integer :: nmembers = 0
    character(len=64), allocatable :: term(:)
    real(r8), allocatable :: multipler(:)
  endtype sums_grp_t

contains

  !-------------------------------------------------------------------
  ! parses summation strings
  !-------------------------------------------------------------------
  subroutine parse_sums(sums, ngrps, grps)

    character(len=CL),   intent(in ) :: sums(:)
    integer,             intent(out) :: ngrps
    type(sums_grp_t), allocatable, intent(out) :: grps(:)

    integer :: ndxs(512)
    integer :: nelem, i,j,k
    character(len=CXX) :: tmp_str, tmp_name

    character(len=8) :: xchr ! multipler
    real(r8) :: xdbl

    logical :: more_to_come
    integer, parameter :: maxgrps = 100
    character(len=CXX) :: sums_grps(maxgrps)

    character(len=CXX) :: sum_string

    sums_grps(:) = ' '

    ! combine lines that have a trailing "+" with the next line
    i=1
    j=1
    do while( len_trim(sums(i)) > 0 )

       k = scan(sums(i), '+', back=.true. )
       more_to_come = k == len_trim(sums(i)) ! line ends with "+"

       if ( more_to_come ) then
          sums_grps(j) = trim(sums_grps(j)) // trim(adjustl(sums(i)))
       else
          sums_grps(j) = trim(sums_grps(j)) // trim(adjustl(sums(i)))
          j = j+1
       endif
       i = i+1

    end do

    ngrps = j-1

    ! a group is  a summation of terms

    ! parse the individual sum strings...  and form the groupings
    has_grps: if (ngrps>0) then

       allocate( grps(ngrps) )

       ! from shr_megan_mod ... should be generalized and shared...
       grploop: do i = 1,ngrps

          ! parse out the term names
          ! from first parsing out the terms in the summation equation ("+" separates the terms)

          sum_string = sums_grps(i)
          j = scan( sum_string, '=' )
          nelem = 1
          ndxs(nelem) = j ! ndxs stores the index of each term of the equation

          ! find indices of all the terms in the equation
          tmp_str = trim( sum_string(j+1:) )
          j = scan( tmp_str, '+' )
          do while(j>0)
             nelem = nelem+1
             ndxs(nelem) = ndxs(nelem-1) + j
             tmp_str = tmp_str(j+1:)
             j = scan( tmp_str, '+' )
          enddo
          ndxs(nelem+1) = len(sum_string)+1

          grps(i)%nmembers = nelem ! number of terms 
          grps(i)%name =  trim(adjustl( sum_string(:ndxs(1)-1))) ! thing to the left of the "=" is used as the name of the group

          ! now that we have the number of terms in the summation allocate memory for the terms
          allocate(grps(i)%term(nelem)) 
          allocate(grps(i)%multipler(nelem))

          ! now parse out the multiplier from the terms 
          elmloop: do k = 1,nelem

             grps(i)%multipler(k) = 1._r8
             ! get the term name which follows the '*' operator if the is one
             tmp_name = adjustl(sum_string(ndxs(k)+1:ndxs(k+1)-1))
             j = scan( tmp_name, '*' )
             if (j>0) then
                xchr = tmp_name(1:j-1) ! get the multipler (left of the '*')
                read( xchr, * ) xdbl   ! convert the string to a real
                grps(i)%multipler(k) = xdbl ! store the multiplier
                tmp_name = adjustl(tmp_name(j+1:)) ! get the term name (right of the '*')
             endif
             grps(i)%term(k) = trim(tmp_name)

          enddo elmloop
       enddo grploop
    endif has_grps

  end subroutine parse_sums


end module sums_utils
