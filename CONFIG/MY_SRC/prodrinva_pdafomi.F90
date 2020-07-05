!-------------------------------------------------------------------------------
!> Call-back routine for prodRinvA
!!
!! This routine calls the routine PDAFomi_prodRinvA
!! for each observation type
!!
SUBROUTINE prodrinva_pdafomi(step, dim_obs_p, ncol, obs_p, A_p, C_p)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_prodRinvA
  ! Include observation types (rename generic name)
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step              !< Current time step
  INTEGER, INTENT(in) :: dim_obs_p         !< Dimension of PE-local observation vector
  INTEGER, INTENT(in) :: ncol              !< Number of columns in A_p and C_p
  REAL, INTENT(in)    :: obs_p(dim_obs_p)  !< PE-local vector of observations
  REAL, INTENT(in)    :: A_p(dim_obs_p, ncol) !< Input matrix
  REAL, INTENT(out)   :: C_p(dim_obs_p, ncol) !< Output matrix


! *************************************
! *** Compute                       ***
! ***                -1             ***
! ***           C = R   A           ***
! *************************************

  CALL PDAFomi_prodRinvA(obs_A, ncol, A_p, C_p)
  CALL PDAFomi_prodRinvA(obs_B, ncol, A_p, C_p)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE prodrinva_pdafomi
