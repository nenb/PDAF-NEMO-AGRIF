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
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step   !< Current time step


! *************************************
! *** Deallocate observation arrays ***
! *************************************

  CALL PDAFomi_deallocate_obs(obs_A)
  CALL PDAFomi_deallocate_obs(obs_B)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE deallocate_obs_pdafomi
