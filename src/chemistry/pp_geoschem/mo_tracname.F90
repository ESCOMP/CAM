
      module mo_tracname
!-----------------------------------------------------------
! 	... List of advected and non-advected trace species, and
!           surface fluxes for the advected species.
!-----------------------------------------------------------

      use chem_mods, only : nTracersMax

      implicit none

! modified to an arbitrary high #, was gas_pcnst. this would cause a memory
! overflow overwrite in mo_sim_dat, which allocates :273 larger than
! the default specified gas_pcnst (hplin, 5/16/20)
      character(len=16) :: solsym(284)   ! species names

      end module mo_tracname
