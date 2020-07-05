!-------------------------------------------------------------------------------
!> Call-back routine for init_obs_f
!!
!! This routine calls the routine PDAFomi_init_obs_f
!! for each observation type
!!
SUBROUTINE init_obs_f_pdafomi(step, dim_obs_f, observation_f)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_obs_f
  ! Include observation types (rename generic name)
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step        !< Current time step
  INTEGER, INTENT(in) :: dim_obs_f   !< Dimension of full observation vector
  REAL, INTENT(out)   :: observation_f(dim_obs_f) !< Full observation vector

! *** local variables ***
  INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector


! ******************************************
! *** Initialize full observation vector ***
! ******************************************

  ! Initialize offset (it will be incremented in PDAFomi_init_obs_f)
  offset_obs_f = 0

  ! The order of the calls has to be consistent with those in obs_op_f_pdafomi
  CALL PDAFomi_init_obs_f(obs_A, dim_obs_f, observation_f, offset_obs_f)
  CALL PDAFomi_init_obs_f(obs_B, dim_obs_f, observation_f, offset_obs_f)
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_obs_f_pdafomi
