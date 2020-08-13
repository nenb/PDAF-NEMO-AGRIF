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
  USE mod_obs_ssh_par_pdafomi, ONLY: obs_ssh_par => thisobs
  USE mod_obs_ssh_child_pdafomi, ONLY: obs_ssh_child => thisobs
  USE mod_obs_fake_ssh_par_pdafomi, ONLY: obs_fake_ssh_par => thisobs
  USE mod_obs_fake_ssh_child_pdafomi, ONLY: obs_fake_ssh_child => thisobs
  USE mod_assimilation_pdaf, ONLY: indx_dom_l_par, indx_dom_l_child

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(obs_ssh_par)
  CALL PDAFomi_deallocate_obs(obs_ssh_child)
  CALL PDAFomi_deallocate_obs(obs_fake_ssh_par)
  CALL PDAFomi_deallocate_obs(obs_fake_ssh_child)

  ! Tidy-up from init_n_domains_pdaf
  IF(ALLOCATED(indx_dom_l_par)) DEALLOCATE(indx_dom_l_par)

#if defined key_agrif
  ! Tidy-up from init_n_domains_pdaf
  IF(ALLOCATED(indx_dom_l_child)) DEALLOCATE(indx_dom_l_child)
#endif

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE deallocate_obs_pdafomi
