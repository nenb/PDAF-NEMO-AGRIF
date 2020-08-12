MODULE mod_obs_op_pdaf
!$AGRIF_DO_NOT_TREAT

!  !DESCRIPTION:
!  This module contains the observation operators for NEMO/AGRIF

  !  !USES:
  USE mod_kind_pdaf
  USE PDAFomi_obs_f, &
       ONLY: obs_f, PDAFomi_gather_obsstate_f

  IMPLICIT NONE

  INTEGER :: mem_par = 0  ! Current ensemble member - parent grid
  INTEGER :: mem_child = 0  ! Current ensemble member - child_grid

CONTAINS

  SUBROUTINE ssh_par_obs_op_gcirc(thisobs, nrows, state_p, obs_f_all, offset_obs)

    !  !DESCRIPTION:
    !  Use great-circle linear interpolation for ssh field on parent grid.

    USE mod_statevector_pdaf, &
         ONLY: halo_2d_par
    USE mod_parallel_pdaf, &
         ONLY: mype_filter, n_modeltasks

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL(pwp), INTENT(in)    :: state_p(:)     !< PE-local model state (dim_p)
    REAL(pwp), INTENT(inout) :: obs_f_all(:)   !< Full observed state for all observation types (nobs_f_all)
    INTEGER, INTENT(inout) :: offset_obs   !< Offset of current observation in overall observation vector

! *** Local variables ***
    INTEGER :: i, row                      ! Counters
    REAL(pwp), ALLOCATABLE :: ostate_p(:)  ! Local observed part of state vector
    REAL(pwp) :: rrows                     ! Real-value for nrows
    REAL(pwp) :: mstate                    ! Local state vector value


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Check if required arrays are allocated (assuming that they are initialzed in this case)
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: ssh_par_obs_op_gcirc - thisobs%id_obs_p is not allocated'
       END IF
       IF (.NOT.ALLOCATED(thisobs%icoeff_p)) THEN
          WRITE (*,*) 'ERROR: ssh_par_obs_op_gcirc - thisobs%icoeff_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by weighted averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF

       rrows = REAL(nrows)

       ! Determine current ensemble member
       mem_par=MOD(mem_par,n_modeltasks)
       mem_par=mem_par+1

       DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = 0.0
          DO row = 1, nrows
             ! If index is negative, this means use halo region
             IF (thisobs%id_obs_p(row,i) < 0.0_pwp) THEN
                ! Use index value of 1 for SSH.
                mstate = halo_2d_par(-thisobs%id_obs_p(row,i),mem_par,1)
             ELSE
                mstate = state_p(thisobs%id_obs_p(row,i))
             END IF
             ostate_p(i) = ostate_p(i) + thisobs%icoeff_p(row,i)*mstate
          END DO
       ENDDO

       ! *** Store offset (mandatory!)
       thisobs%off_obs_f = offset_obs

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate_f(thisobs, ostate_p, obs_f_all, offset_obs)

       ! *** Clean up
       DEALLOCATE(ostate_p)

    END IF doassim

  END SUBROUTINE ssh_par_obs_op_gcirc

  SUBROUTINE ssh_child_obs_op_gcirc(thisobs, nrows, state_p, obs_f_all, offset_obs)

    !  !DESCRIPTION:
    !  Use great-circle linear interpolation for ssh field on child grid.

    USE mod_statevector_pdaf, &
         ONLY: halo_2d_child
    USE mod_parallel_pdaf, &
         ONLY: mype_filter, n_modeltasks

    IMPLICIT NONE

! *** Arguments ***
    TYPE(obs_f), INTENT(inout) :: thisobs  !< Data type with full observation
    INTEGER, INTENT(in) :: nrows           !< Number of values to be averaged
    REAL(pwp), INTENT(in)    :: state_p(:)     !< PE-local model state (dim_p)
    REAL(pwp), INTENT(inout) :: obs_f_all(:)   !< Full observed state for all observation types (nobs_f_all)
    INTEGER, INTENT(inout) :: offset_obs   !< Offset of current observation in overall observation vector

! *** Local variables ***
    INTEGER :: i, row                      ! Counters
    REAL(pwp), ALLOCATABLE :: ostate_p(:)  ! Local observed part of state vector
    REAL(pwp) :: rrows                     ! Real-value for nrows
    REAL(pwp) :: mstate                    ! Local state vector value


! *********************************************
! *** Perform application of measurement    ***
! *** operator H on vector or matrix column ***
! *********************************************

    doassim: IF (thisobs%doassim == 1) THEN

       ! Check if required arrays are allocated (assuming that they are initialzed in this case)
       IF (.NOT.ALLOCATED(thisobs%id_obs_p)) THEN
          WRITE (*,*) 'ERROR: ssh_child_obs_op_gcirc - thisobs%id_obs_p is not allocated'
       END IF
       IF (.NOT.ALLOCATED(thisobs%icoeff_p)) THEN
          WRITE (*,*) 'ERROR: ssh_child_obs_op_gcirc - thisobs%icoeff_p is not allocated'
       END IF

       ! *** PE-local: Initialize observed part state vector by weighted averaging

       IF (thisobs%dim_obs_p>0) THEN
          ALLOCATE(ostate_p(thisobs%dim_obs_p))
       ELSE
          ALLOCATE(ostate_p(1))
       END IF

       rrows = REAL(nrows)

       ! Determine current ensemble member
       mem_child=MOD(mem_child,n_modeltasks)
       mem_child=mem_child+1

       DO i = 1, thisobs%dim_obs_p
          ostate_p(i) = 0.0
          DO row = 1, nrows
             ! If index is negative, this means use halo region
             IF (thisobs%id_obs_p(row,i) < 0.0_pwp) THEN
                ! Use index value of 1 for SSH.
                mstate = halo_2d_child(-thisobs%id_obs_p(row,i),mem_child,1)
             ELSE
                mstate = state_p(thisobs%id_obs_p(row,i))
             END IF
             ostate_p(i) = ostate_p(i) + thisobs%icoeff_p(row,i)*mstate
          END DO
       ENDDO

       ! *** Store offset (mandatory!)
       thisobs%off_obs_f = offset_obs

       ! *** Global: Gather full observed state vector
       CALL PDAFomi_gather_obsstate_f(thisobs, ostate_p, obs_f_all, offset_obs)

       ! *** Clean up
       DEALLOCATE(ostate_p)

    END IF doassim

  END SUBROUTINE ssh_child_obs_op_gcirc

!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_obs_op_pdaf
