!-------------------------------------------------------------------------------
!> Call-back routine for deallocate_obs
!!
!! This routine calls the routine PDAFomi_deallocate_obs
!! for each observation type
!!
SUBROUTINE deallocate_obs_pdafomi(step)
!$AGRIF_DO_NOT_TREAT

  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_deallocate_obs
  ! Include observation types (rename generic name)
  USE mod_obs_ssh_NEMO_pdafomi, ONLY: obs_ssh_NEMO => thisobs
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs
  USE mod_assimilation_pdaf, ONLY: indx_dom_l_par, indx_dom_l_child

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(obs_ssh_NEMO)
  CALL PDAFomi_deallocate_obs(obs_B)

  ! Tidy-up from init_n_domains_pdaf
  IF(ALLOCATED(indx_dom_l_par)) DEALLOCATE(indx_dom_l_par)

#if defined key_agrif
  ! Tidy-up from init_n_domains_pdaf
  IF(ALLOCATED(indx_dom_l_child)) DEALLOCATE(indx_dom_l_child)
#endif

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE deallocate_obs_pdafomi
