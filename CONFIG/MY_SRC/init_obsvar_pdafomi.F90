!-------------------------------------------------------------------------------
!> Call-back routine for init_obsvar
!!
!! This routine calls the routine PDAFomi_init_obsvar_f
!! for each observation type
!!
SUBROUTINE init_obsvar_pdafomi(step, dim_obs_p, obs_p, meanvar)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obsvar_f
  ! Include observation types (rename generic name)
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step          !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p     !< PE-local dimension of observation vector
  REAL, INTENT(in) :: obs_p(dim_obs_p) !< PE-local observation vector
  REAL, INTENT(out)   :: meanvar       !< Mean observation error variance

! *** Local variables ***
  INTEGER :: cnt_obs_f


! *****************************
! *** Compute mean variance ***
! *****************************

  ! Initialize observation counter (it will be incremented in PDAFomi_init_obsvar_f)
  cnt_obs_f = 0

  CALL PDAFomi_init_obsvar_f(obs_A, meanvar, cnt_obs_f)
  CALL PDAFomi_init_obsvar_f(obs_B, meanvar, cnt_obs_f)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_obsvar_pdafomi
