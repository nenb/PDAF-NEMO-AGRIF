!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_l
!!
!! This routine calls the routine PDAFomi_init_obs_l
!! for each observation type
!!
SUBROUTINE init_obs_l_pdafomi(domain_p, step, dim_obs_l, observation_l)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obs_l
  ! Include observation types (rename generic name)
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p   !< Index of current local analysis domain index
  INTEGER, INTENT(in) :: step       !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l  !< Local dimension of observation vector
  REAL, INTENT(out)   :: observation_l(dim_obs_l) !< Local observation vector


! *******************************************
! *** Initialize local observation vector ***
! *******************************************

  CALL PDAFomi_init_obs_l(obs_A_l, obs_A, observation_l)
  CALL PDAFomi_init_obs_l(obs_B_l, obs_B, observation_l)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_obs_l_pdafomi
