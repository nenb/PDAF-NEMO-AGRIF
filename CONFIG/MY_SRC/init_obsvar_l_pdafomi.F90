!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar_l
!!
!! This routine calls the routine PDAFomi_init_obsvar_l
!! for each observation type
!!
SUBROUTINE init_obsvar_l_pdafomi(domain_p, step, dim_obs_l, obs_l, meanvar_l)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obsvar_l
  ! Include observation types (rename generic name)
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: domain_p      !< Index of current local analysis domain
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_l     !< Local dimension of observation vector
  REAL, INTENT(in) :: obs_l(dim_obs_l) !< Local observation vector
  REAL, INTENT(out)   :: meanvar_l     !< Mean local observation error variance

! *** Local variables ***
  INTEGER :: cnt_obs_l


! ***********************************
! *** Compute local mean variance ***
! ***********************************

  ! Initialize observation counter (it will be incremented in PDAFomi_init_obsvar_f)
  cnt_obs_l = 0

  CALL PDAFomi_init_obsvar_l(obs_A_l, obs_A, meanvar_l, cnt_obs_l)
  CALL PDAFomi_init_obsvar_l(obs_B_l, obs_B, meanvar_l, cnt_obs_l)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_obsvar_l_pdafomi
