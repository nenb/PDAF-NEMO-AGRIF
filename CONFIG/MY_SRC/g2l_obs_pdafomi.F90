!-------------------------------------------------------------------------------
!> Call-back routine for g2l_obs
!!
!! This routine calls the routine PDAFomi_g2l_obs
!! for each observation type
!!
SUBROUTINE g2l_obs_pdafomi(domain_p, step, dim_obs_f, dim_obs_l, ostate_f, &
     ostate_l)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_g2l_obs
  ! Include observation types (rename generic name)
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f  !< Dimension of full PE-local observation vector
  INTEGER, INTENT(in) :: dim_obs_l  !< Dimension of local observation vector
  REAL, INTENT(in)    :: ostate_f(dim_obs_f)   !< Full PE-local obs.ervation vector
  REAL, INTENT(out)   :: ostate_l(dim_obs_l)   !< Observation vector on local domain


! *******************************************************
! *** Perform localization of some observation vector ***
! *** to the current local analysis domain.           ***
! *******************************************************

  CALL PDAFomi_g2l_obs(obs_A_l, obs_A, ostate_f, ostate_l)
  CALL PDAFomi_g2l_obs(obs_B_l, obs_B, ostate_f, ostate_l)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE g2l_obs_pdafomi
