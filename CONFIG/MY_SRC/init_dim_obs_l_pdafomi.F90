!-------------------------------------------------------------------------------
!> Call-back routine for init_dim_obs_l
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type
!!
SUBROUTINE init_dim_obs_l_pdafomi(domain_p, step, dim_obs_f, dim_obs_l)
!$AGRIF_DO_NOT_TREAT

  ! Include functions for different observations
  USE mod_obs_ssh_par_pdafomi, ONLY: init_dim_obs_l_ssh_par
  USE mod_obs_ssh_child_pdafomi, ONLY: init_dim_obs_l_ssh_child
  USE mod_kind_pdaf
  USE mod_assimilation_pdaf, &
       ONLY: indx_dom_l_par, indx_dom_l_child, num_domains_par
  USE mod_statevector_pdaf, &
       ONLY: mpi_subd_lat_par, mpi_subd_lon_child, mpi_subd_lon_par
  USE mod_parallel_pdaf, &
       ONLY: nldj_child, nldj_par, nldi_child, nldi_par
  USE mod_agrif_pdaf, &
       ONLY: tmask_par

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: domain_p   !< Index of current local analysis domain
  INTEGER, INTENT(in)  :: step       !< Current time step
  INTEGER, INTENT(in)  :: dim_obs_f  !< Full dimension of observation vector
  INTEGER, INTENT(out) :: dim_obs_l  !< Local dimension of observation vector

! *** local variables ***
  INTEGER :: dim_obs_l_ssh_par ! Dimension of observation ssh_par
  INTEGER :: dim_obs_l_ssh_child ! Dimension of observation type B
  INTEGER :: offset_obs_l, offset_obs_f  ! local and full offsets
  INTEGER :: domain_p_child ! Counter for local analysis domain on child grid
  INTEGER :: i, j           ! Coordinates for local analysis domain
  INTEGER :: i0, j0         ! Halo offset for local PE
  INTEGER :: i_par, j_par   ! Grid coordinates for local analysis domain
  LOGICAL :: land           ! Indicate whether surface gridpoint land


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

  ! Initialize offsets (they are incremented in PDAFomi_init_dim_obs_l)
  offset_obs_l = 0
  offset_obs_f = 0

  ! Hack for dealing with case when no valid local domains on PE. See
  ! init_n_domains for further details.
  land = .FALSE.
  IF (domain_p == 1 .AND. num_domains_par > 0) THEN
     ! Compute halo offset
     i0 = nldi_par - 1
     j0 = nldj_par - 1
     ! Compute i,j indices for parent grid.
     i_par = indx_dom_l_par(1, domain_p)
     j_par = indx_dom_l_par(2, domain_p)
     ! Check whether the local domain is actually a land point (and
     ! hence not a valid local domain).
     IF(tmask_par(i_par+i0, j_par+j0, 1) == 0.0_pwp) land = .TRUE.
  END IF

  ! Set local observation dimension to 0 if local domain is on land. This is a hack
  ! so that no update will be computed for these (invalid) local domains. Otherwise
  ! call PDAF-OMI routines.
  IF(land) THEN
     dim_obs_l_ssh_par = 0
     dim_obs_l_ssh_child = 0
  ELSE
     ! Call init_dim_obs_l specific for each observation
     ! The order of the calls has to be consistent with that in obs_op_f_pdafomi
     CALL init_dim_obs_l_ssh_par(domain_p, step, dim_obs_f, dim_obs_l_ssh_par, offset_obs_l, offset_obs_f)
     CALL init_dim_obs_l_ssh_child(domain_p, step, dim_obs_f, dim_obs_l_ssh_child, offset_obs_l, offset_obs_f)
  END IF

  ! Compute overall local observation dimension
  dim_obs_l = dim_obs_l_ssh_par + dim_obs_l_ssh_child

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_dim_obs_l_pdafomi
