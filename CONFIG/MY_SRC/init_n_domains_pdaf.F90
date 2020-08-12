!$Id: init_n_domains_pdaf.F90 343 2020-01-21 14:21:42Z lnerger $
!> Set number of local analysis domains
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF
!!
!! The routine is called in PDAF_X_update 
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains. 
!! It has to set the number of local analysis 
!! domains for the PE-local domain.
!!
!! Implementation for the 2D online example
!! without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE init_n_domains_pdaf(step, n_domains_p)
!$AGRIF_DO_NOT_TREAT

  USE mod_kind_pdaf
  USE mod_assimilation_pdaf, &
       ONLY: indx_dom_l_par, indx_dom_l_child, num_domains_child, &
       num_domains_par
  USE mod_statevector_pdaf, &
       ONLY: mpi_subd_lon_child, mpi_subd_lon_par, mpi_subd_lat_child, &
       mpi_subd_lat_par
  USE mod_parallel_pdaf, &
       ONLY: nldi_par, nldj_par, nldi_child, nldj_child, nimpp_par, &
       njmpp_par, mype_filter
  USE mod_agrif_pdaf, &
       ONLY: tmask_par, tmask_child

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step        !< Current time step
  INTEGER, INTENT(out) :: n_domains_p !< PE-local number of analysis domains

  ! Local variables
  INTEGER :: i, j, cnt      ! Counters
  INTEGER :: i_glbl, j_glbl ! Indices on global grid
  INTEGER :: i0, j0         ! Halo offset


! ************************************
! *** Initialize number of domains ***
! ************************************

  ! *******************************************
  ! The number of local domains is defined as
  ! the number of grid points at the surface
  ! where tmask is 1 ie horizontal localization
  ! is used, and land points are ignored.
  ! *******************************************

! ***********
! Parent grid
! ***********

  ! Initialisation.
  n_domains_p = 0

  ! Compute halo offset on parent grid.
  i0 = nldi_par - 1
  j0 = nldj_par - 1

  ! Use horizontal localization.
  DO j = 1, mpi_subd_lat_par
     DO i = 1, mpi_subd_lon_par
        ! Only include as local domain if ocean at surface.
        mask: IF(tmask_par(i+i0, j+j0, 1) == 1.0_pwp) THEN
#if defined key_agrif
!!$           ! **********************************
!!$           ! THIS METHOD IS CURRENTLY NOT USED!
!!$           ! **********************************
!!$           !
!!$           ! *******************************************
!!$           ! In the case of NEMO+AGRIF, the calculation
!!$           ! is slightly more complicated. The same
!!$           ! method of horizontal localization is used
!!$           ! as before, except now any points in the
!!$           ! NEMO grid that are also contained in the
!!$           ! AGRIF grid are ignored. This is so that an
!!$           ! update is only performed on the AGRIF grid,
!!$           ! which is necessary for theoretical reasons.
!!$           ! These updated values on the AGRIF grid are
!!$           ! then used to update the respective values
!!$           ! on the NEMO grid at the end of the PDAF
!!$           ! analysis step.
!!$           ! *******************************************
!!$           !
!!$           ! Only include NEMO surface ocean points outside
!!$           ! or on boundaries of AGRIF grid.
!!$           i_glbl = i - 1 + nimpp_par + i0
!!$           j_glbl = j - 1 + njmpp_par + j0
!!$           IF(j_glbl > 134 .AND. j_glbl < 268) THEN
!!$              IF(i_glbl > 435 .AND. i_glbl < 715) THEN
!!$                 CYCLE
!!$              END IF
!!$           END IF
!!$           ! *******************************************
           n_domains_p = n_domains_p + 1
#else
           ! Include all NEMO surface ocean points.
           n_domains_p = n_domains_p + 1
#endif
        END IF mask
     END DO
  END DO

  ! Store number of domains on parent grid for later.
  num_domains_par = n_domains_p

  ! ***********************************
  ! Store local domain i,j index values
  ! ***********************************

  ! Treat case where no valid local domains on parent grid.
  IF(n_domains_p > 0) THEN
     ! Store i,j indices of local domain for use later.
     ! (Deallocated in deallocate_obs_pdafomi.)
     ALLOCATE(indx_dom_l_par(2,n_domains_p))
     ! Counter for local domain
     cnt=0
     DO j = 1, mpi_subd_lat_par
        DO i = 1, mpi_subd_lon_par
           indx_mask: IF(tmask_par(i+i0, j+j0, 1) == 1.0_pwp) THEN
#if defined key_agrif
!!$              ! **********************************
!!$              ! THIS METHOD IS CURRENTLY NOT USED!
!!$              ! **********************************
!!$              !
!!$              ! *******************************************
!!$              ! In the case of NEMO+AGRIF, the calculation
!!$              ! is slightly more complicated. The same
!!$              ! method of horizontal localization is used
!!$              ! as before, except now any points in the
!!$              ! NEMO grid that are also contained in the
!!$              ! AGRIF grid are ignored. This is so that an
!!$              ! update is only performed on the AGRIF grid,
!!$              ! which is necessary for theoretical reasons.
!!$              ! These updated values on the AGRIF grid are
!!$              ! then used to update the respective values
!!$              ! on the NEMO grid at the end of the PDAF
!!$              ! analysis step.
!!$              ! *******************************************
!!$              !
!!$              ! Only include NEMO surface ocean points outside
!!$              ! or on boundaries of AGRIF grid.
!!$              i_glbl = i - 1 + nimpp_par + i0
!!$              j_glbl = j - 1 + njmpp_par + j0
!!$              IF(j_glbl > 134 .AND. j_glbl < 268) THEN
!!$                 IF(i_glbl > 435 .AND. i_glbl < 715) THEN
!!$                    CYCLE
!!$                 END IF
!!$              END IF
!!$              ! *******************************************
              cnt=cnt+1
              indx_dom_l_par(1,cnt) = i
              indx_dom_l_par(2,cnt) = j
#else
              ! Include all NEMO surface ocean points.
              cnt=cnt+1
              indx_dom_l_par(1,cnt) = i
              indx_dom_l_par(2,cnt) = j
#endif
           END IF indx_mask
        END DO
     END DO
  END IF


! **********
! Child grid
! **********

#if defined key_agrif
  ! Compute halo offset on child grid.
  i0 = nldi_child - 1
  j0 = nldj_child - 1

  ! Counter for number of local domains on *child grid*.
  cnt=0

  ! Use horizontal localization, hence number of domains corresponds to
  ! number of (x,y) points on parent grid (plus child grid if AGRIF used).
  DO j = 1, mpi_subd_lat_child
     DO i = 1, mpi_subd_lon_child
        ! Only include as local domain if ocean at surface.
        IF(tmask_child(i+i0, j+j0, 1) == 1.0_pwp) THEN
           n_domains_p = n_domains_p + 1
           cnt = cnt + 1
        END IF
     END DO
  END DO

  ! Store number of domains on child grid for later.
  num_domains_child = cnt

  ! ***********************************
  ! Store local domain i,j index values
  ! ***********************************

  ! Treat case where no valid local domains on child grid.
  IF(cnt > 0 ) THEN
     ! Store i,j indices of local domain for use later.
     ! (Deallocated in deallocate_obs_pdafomi.)
     ALLOCATE(indx_dom_l_child(2,cnt))
     ! Counter for local domain.
     cnt=0
     DO j = 1, mpi_subd_lat_child
        DO i = 1, mpi_subd_lon_child
           ! Only include as local domain if ocean at surface.
           IF(tmask_child(i+i0, j+j0, 1) == 1.0_pwp) THEN
              cnt=cnt+1
              indx_dom_l_child(1,cnt) = i
              indx_dom_l_child(2,cnt) = j
           END IF
        END DO
     END DO
  END IF
#endif

! ****************
! No local domains
! ****************

  ! Hack for when no valid local domains on PE. We set the dimension of the
  ! local domain to 1, and ensure (see init_dim_obs_l) that the dimension of the
  ! observations for this local domain is zero. Hence no update will be performed
  ! on this local domain.
  IF(n_domains_p == 0) THEN
     WRITE(*,'(8x,a, i3)') 'WARNING: No valid local domains on PE:', mype_filter
     n_domains_p = 1
     num_domains_par = 1
     ALLOCATE(indx_dom_l_par(2,1))
     indx_dom_l_par(1,1) = 1
     indx_dom_l_par(2,1) = 1
  END IF

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_n_domains_pdaf
