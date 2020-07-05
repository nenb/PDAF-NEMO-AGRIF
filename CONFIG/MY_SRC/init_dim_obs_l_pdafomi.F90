!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs_f, dim_obs_l)
!$AGRIF_DO_NOT_TREAT
  ! Include PDAFomi function
  USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l
  ! Include observation types
  USE mod_obs_A_pdafomi, ONLY: obs_A => thisobs, obs_A_l => thisobs_l
  USE mod_obs_B_pdafomi, ONLY: obs_B => thisobs, obs_B_l => thisobs_l

  ! Include localization radius and local coordinates
  USE mod_assimilation_pdaf, &
       ONLY: local_range, coords_l

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
  INTEGER :: dim_obs_l_A ! Dimension of observation type A
  INTEGER :: dim_obs_l_B ! Dimension of observation type B
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Initialize offsets (they are incremented in PDAFomi_init_dim_obs_l)
  offset_obs_l = 0
  offset_obs_f = 0

  ! Call init_dim_obs_l specific for each observation
  ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
  CALL PDAFomi_init_dim_obs_l(obs_A_l, obs_A, coords_l, local_range, dim_obs_l_A, &
       offset_obs_l, offset_obs_f)
  CALL PDAFomi_init_dim_obs_l(obs_B_l, obs_B, coords_l, local_range, dim_obs_l_B, &
       offset_obs_l, offset_obs_f)

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_A + dim_obs_l_B
!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_dim_obs_l_pdafomi
