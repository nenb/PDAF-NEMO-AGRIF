MODULE mod_agrif_pdaf
!$AGRIF_DO_NOT_TREAT

!  !DESCRIPTION:
!  This module contains auxiliary functions for exchanging data between
!  NEMO/AGRIF and PDAF.

  !  !USES:
  USE mod_kind_pdaf

  IMPLICIT NONE

  ! Indices for location of S and T fields in 4D array
  INTEGER :: jp_tem_child
  INTEGER :: jp_tem_par
  INTEGER :: jp_sal_child
  INTEGER :: jp_sal_par

  ! Time constants for different grids
  INTEGER :: nitend_par
  INTEGER :: nit000_par
  INTEGER :: nitend_child
  INTEGER :: nit000_child
  REAL(pwp) :: rdt_par
  REAL(pwp) :: rdt_child

!$AGRIF_END_DO_NOT_TREAT

CONTAINS

  SUBROUTINE calc_grid_cnst()

    ! !DESCRIPTION:
    ! Compute value of NEMO, AGRIF grid constants and store in PDAF.
    !
    ! NOTE: When AGRIF is used, this subroutine is called twice, once
    ! on each grid.

    ! !USES:
    USE in_out_manager, ONLY: nitend, nit000
    USE par_oce, ONLY: jpiglo, jpjglo, jpk, jp_tem, jp_sal
    USE dom_oce, ONLY: nldi, nldj, nlei, nlej, nimpp, njmpp, rdt
    USE mod_parallel_pdaf, ONLY: jpiglo_par, jpjglo_par, jpk_par, &
         jpiglo_child, jpjglo_child, jpk_child, nldi_par, &
         nldj_par, nlei_par, nlej_par, nimpp_par, njmpp_par, &
         nldi_child, nldj_child, nlei_child, nlej_child, nimpp_child, &
         njmpp_child
    IMPLICIT NONE


#if defined key_agrif
    IF (Agrif_Root()) THEN
       jpiglo_par = jpiglo
       jpjglo_par = jpjglo
       jpk_par = jpk
       jp_tem_par = jp_tem
       jp_sal_par = jp_sal
       nit000_par = nit000
       nitend_par = nitend
       nldi_par = nldi
       nldj_par = nldj
       nlei_par = nlei
       nlej_par = nlej
       nimpp_par = nimpp
       njmpp_par = njmpp
       rdt_par = rdt
    ELSE
       jpiglo_child = jpiglo
       jpjglo_child = jpjglo
       jpk_child = jpk
       jp_tem_child = jp_tem
       jp_sal_child = jp_sal
       nit000_child = nit000
       nitend_child = nitend
       nldi_child = nldi
       nldj_child = nldj
       nlei_child = nlei
       nlej_child = nlej
       nimpp_child = nimpp
       njmpp_child = njmpp
       rdt_child = rdt
    ENDIF
#else
    jpiglo_par = jpiglo
    jpjglo_par = jpjglo
    jpk_par = jpk
    jp_tem_par = jp_tem
    jp_sal_par = jp_sal
    nit000_par = nit000
    nitend_par = nitend
    nldi_par = nldi
    nldj_par = nldj
    nlei_par = nlei
    nlej_par = nlej
    nimpp_par = nimpp
    njmpp_par = njmpp
    rdt_par = rdt
#endif

  END SUBROUTINE calc_grid_cnst

  SUBROUTINE fill2d_statevector(state_p, grid)

    ! !DESCRIPTION:
    ! Fill state vector on a grid with 2d state variable.

    ! !USES:
    USE par_kind
    USE oce, ONLY: sshb
    USE mod_parallel_pdaf, &
         ONLY: nldj_child, nldj_par, nldi_child, nldi_par
    USE mod_statevector_pdaf, &
         ONLY: ssh_p_offset_par, ssh_p_offset_child, mpi_subd_lat_child, &
         mpi_subd_lat_par, mpi_subd_lon_child, mpi_subd_lon_par
    IMPLICIT NONE

    ! !ARGUMENTS
    REAL(pwp), POINTER, INTENT(inout) :: state_p(:,:)  ! PE-local model state
    CHARACTER(len=*), INTENT(in) :: grid      ! Flag for parent/child grid

    ! *** local variables ***
    INTEGER :: i, j        ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


! ***********
! Parent grid
! ***********
    IF (grid=='par') THEN
       ! Set the starting position after the halo region
       j0 = nldj_par - 1
       i0 = nldi_par - 1

       ! Fill state vector with 2d variables
       ! SSH
       DO j = 1, mpi_subd_lat_par
          DO i = 1, mpi_subd_lon_par
             state_p(i+(j-1)*mpi_subd_lon_par + ssh_p_offset_par, 1) = sshb(i+i0,j+j0)
          END DO
       END DO
! **********
! Child grid
! **********
    ELSE IF (grid=='child') THEN
       ! Set the starting position after the halo region
       j0 = nldj_child - 1
       i0 = nldi_child - 1

       ! Fill state vector with 2d variables
       ! SSH
       DO j = 1, mpi_subd_lat_child
          DO i = 1, mpi_subd_lon_child
             state_p(i+(j-1)*mpi_subd_lon_child + ssh_p_offset_child, 1) = sshb(i+i0,j+j0)
          END DO
       END DO
    END IF

  END SUBROUTINE fill2d_statevector

  SUBROUTINE distrib2d_statevector(state_p, grid)

    ! !DESCRIPTION:
    ! Distribute state vector 2d state variable.

    ! !USES:
    USE mod_statevector_pdaf, &
         ONLY: ssh_p_offset_par, ssh_p_offset_child, mpi_subd_lat_child, &
         mpi_subd_lat_par, mpi_subd_lon_child, mpi_subd_lon_par
    USE mod_parallel_pdaf, &
         ONLY: nldj_child, nldj_par, nldi_child, nldi_par
    USE par_kind
    USE oce, ONLY: sshb

    IMPLICIT NONE

    ! !ARGUMENTS
    REAL(pwp), POINTER, INTENT(inout) :: state_p(:,:)  ! PE-local model state
    CHARACTER(len=*), INTENT(in) :: grid      ! Flag for parent/child grid

    ! *** local variables ***
    INTEGER :: i, j        ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


! ***********
! Parent grid
! ***********
    IF (grid=='par') THEN
       ! Set the starting position after the halo region
       j0 = nldj_par - 1
       i0 = nldi_par - 1

       ! Fill state vector with 2d variables
       ! SSH
       DO j = 1, mpi_subd_lat_par
          DO i = 1, mpi_subd_lon_par
             sshb(i+i0,j+j0) = state_p(i+(j-1)*mpi_subd_lon_par + ssh_p_offset_par, 1)
          END DO
       END DO
! **********
! Child grid
! **********
    ELSE IF (grid=='child') THEN
       ! Set the starting position after the halo region
       j0 = nldj_child - 1
       i0 = nldi_child - 1

       ! Fill state vector with 2d variables
       ! SSH
       DO j = 1, mpi_subd_lat_child
          DO i = 1, mpi_subd_lon_child
             sshb(i+i0,j+j0) = state_p(i+(j-1)*mpi_subd_lon_child + ssh_p_offset_child, 1)
          END DO
       END DO
    END IF

  END SUBROUTINE distrib2d_statevector

  SUBROUTINE fill3d_statevector(state_p, grid)

    ! !DESCRIPTION:
    ! Fill state vector on a grid with 3d state variable.

    ! !USES:
    USE mod_statevector_pdaf, &
         ONLY: t_p_offset_par, t_p_offset_child, &
         s_p_offset_par, s_p_offset_child, u_p_offset_par, u_p_offset_child, &
         v_p_offset_par, v_p_offset_child, mpi_subd_lat_child, mpi_subd_lat_par, &
         mpi_subd_lon_child, mpi_subd_lon_par, mpi_subd_vert_child, &
         mpi_subd_vert_par
    USE mod_parallel_pdaf, &
         ONLY: nldj_child, nldj_par, nldi_par, nldi_child
    USE par_kind
    USE oce, ONLY: tsb, ub, vb
    IMPLICIT NONE

    ! !ARGUMENTS
    REAL(pwp), POINTER, INTENT(inout) :: state_p(:,:)   ! PE-local model state
    CHARACTER(len=*), INTENT(in) :: grid       ! Flag for parent/child grid

    ! *** local variables ***
    INTEGER :: i, j, k     ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


! ***********
! Parent grid
! ***********
    IF (grid=='par') THEN
       ! Set the starting position after the halo region
       j0 = nldj_par - 1
       i0 = nldi_par - 1

       ! Fill state vector with 3d variables
       ! T
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                state_p(i+(j-1)*mpi_subd_lon_par + (k-1)*mpi_subd_lat_par*mpi_subd_lon_par &
                     + t_p_offset_par, 1) =  tsb(i+i0,j+j0,k,jp_tem_par)
             END DO
          END DO
       END DO
       ! S
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                state_p(i+(j-1)*mpi_subd_lon_par + (k-1)*mpi_subd_lat_par*mpi_subd_lon_par &
                     + s_p_offset_par, 1) = tsb(i+i0,j+j0,k,jp_sal_par)
             END DO
          END DO
       END DO
       ! U
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                state_p(i+(j-1)*mpi_subd_lon_par + (k-1)*mpi_subd_lat_par*mpi_subd_lon_par &
                     + u_p_offset_par, 1) = ub(i+i0,j+j0,k)
             END DO
          END DO
       END DO
       ! V
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                state_p(i+(j-1)*mpi_subd_lon_par + (k-1)*mpi_subd_lat_par*mpi_subd_lon_par &
                     + v_p_offset_par, 1) = vb(i+i0,j+j0,k)
             END DO
          END DO
       END DO
! **********
! Child grid
! **********
    ELSE IF (grid=='child') THEN
       ! Set the starting position after the halo region
       j0 = nldj_child - 1
       i0 = nldi_child - 1

       ! Fill state vector with 3d variables
       ! T
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child &
                     + t_p_offset_child, 1) =  tsb(i+i0,j+j0,k,jp_tem_child)
             END DO
          END DO
       END DO
       ! S
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child &
                     + s_p_offset_child, 1) = tsb(i+i0,j+j0,k,jp_sal_child)
             END DO
          END DO
       END DO
       ! U
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child &
                     + u_p_offset_child, 1) = ub(i+i0,j+j0,k)
             END DO
          END DO
       END DO
       ! V
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child &
                     + v_p_offset_child, 1) = vb(i+i0,j+j0,k)
             END DO
          END DO
       END DO
    END IF

  END SUBROUTINE fill3d_statevector

  SUBROUTINE distrib3d_statevector(state_p, grid)

    ! !DESCRIPTION:
    ! Distribute state vector 3d state variable.

    ! !USES:
    USE mod_statevector_pdaf, &
         ONLY: t_p_offset_par, t_p_offset_child, &
         s_p_offset_par, s_p_offset_child, u_p_offset_par, u_p_offset_child, &
         v_p_offset_par, v_p_offset_child, mpi_subd_lat_child, mpi_subd_lat_par, &
         mpi_subd_lon_child, mpi_subd_lon_par, mpi_subd_vert_child, &
         mpi_subd_vert_par
    USE mod_parallel_pdaf, &
         ONLY: nldj_child, nldj_par, nldi_par, nldi_child
    USE par_kind
    USE oce, ONLY: tsb, ub, vb
    IMPLICIT NONE

    ! !ARGUMENTS
    REAL(pwp), POINTER, INTENT(inout) :: state_p(:,:)   ! PE-local model state
    CHARACTER(len=*), INTENT(in) :: grid       ! Flag for parent/child grid

    ! *** local variables ***
    INTEGER :: i, j, k     ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


! ***********
! Parent grid
! ***********
    IF (grid=='par') THEN
       ! Set the starting position after the halo region
       j0 = nldj_par - 1
       i0 = nldi_par - 1

       ! Distribute state vector 3d state variable.
       ! T
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                tsb(i+i0,j+j0,k,jp_tem_par) = &
                     state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + t_p_offset_par, 1)
             END DO
          END DO
       END DO
       ! S
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                tsb(i+i0,j+j0,k,jp_sal_par) = &
                     state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + s_p_offset_par, 1)
             END DO
          END DO
       END DO
       ! U
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                ub(i+i0,j+j0,k) = &
                     state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + u_p_offset_par, 1)
             END DO
          END DO
       END DO
       ! V
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                vb(i+i0,j+j0,k) = &
                     state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + v_p_offset_par, 1)
             END DO
          END DO
       END DO
! **********
! Child grid
! **********
    ELSE IF (grid=='child') THEN
       ! Set the starting position after the halo region
       j0 = nldj_child - 1
       i0 = nldi_child - 1

       ! Distribute state vector 3d state variable.
       ! T
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                tsb(i+i0,j+j0,k,jp_tem_child) = &
                     state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + t_p_offset_child, 1)
             END DO
          END DO
       END DO
       ! S
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                tsb(i+i0,j+j0,k,jp_sal_child) = &
                     state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + s_p_offset_child, 1)
             END DO
          END DO
       END DO
       ! U
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                ub(i+i0,j+j0,k) = &
                     state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + u_p_offset_child, 1)
             END DO
          END DO
       END DO
       ! V
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                vb(i+i0,j+j0,k) = &
                     state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + v_p_offset_child, 1)
             END DO
          END DO
       END DO
    END IF
  END SUBROUTINE distrib3d_statevector

  SUBROUTINE update_halo(grid)

    ! !DESCRIPTION:
    ! Update halo regions for state variables.

    ! !USES:
    USE par_kind
    USE par_oce, ONLY: jpi, jpj, jpk
    USE oce, ONLY: sshb, tsb, ub, vb
    USE lbclnk, ONLY: lbc_lnk

    IMPLICIT NONE

    ! !ARGUMENTS
    CHARACTER(len=*), INTENT(in) :: grid      ! Flag for parent/child grid

    ! *** local variables ***
    INTEGER :: i, j, k                        ! Counters
    REAL(pwp), ALLOCATABLE :: tmp(:,:,:), tmp2(:,:,:)      ! Temporary 3D array for 4D arrays


! ***********
! Parent grid
! ***********
    IF (grid=='par') THEN
       ! Fill halo regions
       CALL lbc_lnk(sshb, 'T', 1.)
       CALL lbc_lnk(ub, 'U', -1.)
       CALL lbc_lnk(vb, 'V', -1.)

       ! Fill halo regions of 4D arrays, use original
       ! NEMO/AGRIF variables for dimensions
       ALLOCATE(tmp(jpi,jpj,jpk))
       ! T
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tmp(i,j,k) = tsb(i,j,k,jp_tem_par)
             END DO
          END DO
       END DO
       CALL lbc_lnk(tmp, 'T', 1.)
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tsb(i,j,k,jp_tem_par) = tmp(i,j,k)
             END DO
          END DO
       END DO
       ! S
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tmp(i,j,k) = tsb(i,j,k,jp_sal_par)
             END DO
          END DO
       END DO
       CALL lbc_lnk(tmp, 'T', 1.)
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tsb(i,j,k,jp_sal_par) = tmp(i,j,k)
             END DO
          END DO
       END DO
       DEALLOCATE(tmp)
! **********
! Child grid
! **********
    ELSE IF (grid=='child') THEN
       ! Fill halo regions
       CALL lbc_lnk(sshb, 'T', 1.)
       CALL lbc_lnk(ub, 'U', -1.)
       CALL lbc_lnk(vb, 'V', -1.)

       ! Fill halo regions of 4D arrays, use original
       ! NEMO/AGRIF variables for dimensions
       ALLOCATE(tmp(jpi,jpj,jpk))
       ! T
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tmp(i,j,k) = tsb(i,j,k,jp_tem_child)
             END DO
          END DO
       END DO
       CALL lbc_lnk(tmp, 'T', 1.)
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tsb(i,j,k,jp_tem_child) = tmp(i,j,k)
             END DO
          END DO
       END DO
       ! S
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tmp(i,j,k) = tsb(i,j,k,jp_sal_child)
             END DO
          END DO
       END DO
       CALL lbc_lnk(tmp, 'T', 1.)
       DO k = 1, jpk
          DO j = 1, jpj
             DO i = 1, jpi
                tsb(i,j,k,jp_sal_child) = tmp(i,j,k)
             END DO
          END DO
       END DO
       DEALLOCATE(tmp)
    END IF

  END SUBROUTINE update_halo

END MODULE mod_agrif_pdaf
