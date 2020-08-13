!-------------------------------------------------------------------------------
!> Call-back routine for obs_op_f
!!
!! This routine calls the observation-specific
!! routines obs_op_f_TYPE.
!!
SUBROUTINE obs_op_f_pdafomi(step, dim_p, dim_obs_f, state_p, ostate_f)
!$AGRIF_DO_NOT_TREAT

  ! Include functions for different observations
  USE mod_obs_ssh_par_pdafomi, ONLY: obs_op_f_ssh_par
  USE mod_obs_ssh_child_pdafomi, ONLY: obs_op_f_ssh_child
  USE mod_obs_fake_ssh_par_pdafomi, ONLY: obs_op_f_fake_ssh_par
  USE mod_obs_fake_ssh_child_pdafomi, ONLY: obs_op_f_fake_ssh_child

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step                 !< Current time step
  INTEGER, INTENT(in) :: dim_p                !< PE-local state dimension
  INTEGER, INTENT(in) :: dim_obs_f            !< Dimension of full observed state
  REAL, INTENT(in)    :: state_p(dim_p)       !< PE-local model state
  REAL, INTENT(inout) :: ostate_f(dim_obs_f)  !< PE-local full observed state

! *** local variables
  INTEGER :: offset_obs_f     ! Count offset of an observation type in full obs. vector


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

  ! Initialize offset
  offset_obs_f = 0

  ! The order of the calls determines how the different observations
  ! are ordered in the full state vector
  CALL obs_op_f_ssh_par(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)
  CALL obs_op_f_ssh_child(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)
  CALL obs_op_f_fake_ssh_par(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)
  CALL obs_op_f_fake_ssh_child(dim_p, dim_obs_f, state_p, ostate_f, offset_obs_f)

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE obs_op_f_pdafomi
