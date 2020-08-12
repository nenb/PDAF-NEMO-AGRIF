!> Call-back routine for init_dim_obs_f
!!
!! This routine calls the observation-specific
!! routines init_dim_obs_f_TYPE.
!!
SUBROUTINE init_dim_obs_f_pdafomi(step, dim_obs_f)
!$AGRIF_DO_NOT_TREAT

  ! Include functions for different observations
  USE mod_obs_ssh_par_pdafomi, ONLY: assim_ssh_par, init_dim_obs_f_ssh_par
  USE mod_obs_ssh_child_pdafomi, ONLY: assim_ssh_child, init_dim_obs_f_ssh_child

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step      !< Current time step
  INTEGER, INTENT(out) :: dim_obs_f !< Dimension of full observation vector

! *** Local variables ***
  INTEGER :: dim_obs_f_ssh_par ! Observation dimensions
  INTEGER :: dim_obs_f_ssh_child ! Observation dimensions


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

  ! Initialize number of observations
  dim_obs_f_ssh_par = 0
  dim_obs_f_ssh_child = 0

  ! Call observation-specific routines
  ! The routines are independent, so it is not relevant
  ! in which order they are called
  IF (assim_ssh_par) CALL init_dim_obs_f_ssh_par(step, dim_obs_f_ssh_par)
  IF (assim_ssh_child) CALL init_dim_obs_f_ssh_child(step, dim_obs_f_ssh_child)

  dim_obs_f = dim_obs_f_ssh_par + dim_obs_f_ssh_child

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_dim_obs_f_pdafomi
