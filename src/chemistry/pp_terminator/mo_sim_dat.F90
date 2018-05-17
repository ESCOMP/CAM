
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,     only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use chem_mods,     only : diag_map
      use chem_mods,     only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,     only : pht_alias_lst, pht_alias_mult
      use chem_mods,     only : extfrc_lst, inv_lst, slvd_lst
      use chem_mods,     only : enthalpy_cnt, cph_enthalpy, cph_rid, num_rnts, rxntot
      use cam_abortutils,only : endrun
      use mo_tracname,   only : solsym
      use chem_mods,     only : frc_from_dataset
      use chem_mods,     only : is_scalar, is_vector
      use shr_kind_mod,  only : r8 => shr_kind_r8
      use cam_logfile,   only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      is_scalar = .true.
      is_vector = .false.

      clscnt(:) = (/      0,     0,     0,     3,     0 /)

      cls_rxt_cnt(:,4) = (/      0,     1,     1,     3 /)

      solsym(:  3) = (/ 'CL              ','CL2             ','RHO             ' /)

      adv_mass(:  3) = (/    35.452700_r8,    70.905400_r8,     1.007400_r8 /)

      crb_mass(:  3) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8 /)

      fix_mass(:  1) = (/ 0.00000000_r8 /)

      clsmap(:  3,4) = (/    1,   2,   3 /)

      permute(:  3,4) = (/    3,   2,   1 /)

      diag_map(:  3) = (/    1,   2,   5 /)

      inv_lst(:  1) = (/ 'M               ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(     1:     2) = (/ 'toy_k1                          ', 'toy_k2                          ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2 /)
      allocate( num_rnts(rxntot-phtcnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate num_rnts; error = ',ios
         call endrun
      end if
      num_rnts(:) = (/      1,     2 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
