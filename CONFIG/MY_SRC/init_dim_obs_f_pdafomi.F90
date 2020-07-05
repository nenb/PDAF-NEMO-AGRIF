!> Call-back routine for init_dim_obs_f
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_f_TYPE.
!!
SUBROUTINE init_dim_obs_f_pdafomi(step, dim_obs_f)
!$AGRIF_DO_NOT_TREAT
  ! Include functions for different observations
  USE mod_obs_A_pdafomi, ONLY: assim_A, init_dim_obs_f_A
  USE mod_obs_B_pdafomi, ONLY: assim_B, init_dim_obs_f_B

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs_f !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_f_A ! Observation dimensions
  INTEGER :: dim_obs_f_B ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_f_A = 0
  dim_obs_f_B = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_A) CALL init_dim_obs_f_A(step, dim_obs_f_A)
  IF (assim_B) CALL init_dim_obs_f_B(step, dim_obs_f_B)

  dim_obs_f = dim_obs_f_A + dim_obs_f_B
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_dim_obs_f_pdafomi
