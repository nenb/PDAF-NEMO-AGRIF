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

  ! i-, j- indexes for each PE
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nimppt_par
  INTEGER, ALLOCATABLE, DIMENSION(:) :: nimppt_child
  INTEGER, ALLOCATABLE, DIMENSION(:) :: njmppt_par
  INTEGER, ALLOCATABLE, DIMENSION(:) :: njmppt_child

  ! Longitude, lattitude indexes for each PE
  REAL(pwp), ALLOCATABLE, DIMENSION(:,:) :: glamt_par
  REAL(pwp), ALLOCATABLE, DIMENSION(:,:) :: glamt_child
  REAL(pwp), ALLOCATABLE, DIMENSION(:,:) :: gphit_par
  REAL(pwp), ALLOCATABLE, DIMENSION(:,:) :: gphit_child

  ! Masks for each PE
  REAL(pwp), ALLOCATABLE, DIMENSION(:,:,:) :: tmask_par
  REAL(pwp), ALLOCATABLE, DIMENSION(:,:,:) :: tmask_child

  ! Thresholds for rejecting observations/updates
  REAL(pwp) :: lowlim_ssh_NEMO
  REAL(pwp) :: upplim_ssh_NEMO
  REAL(pwp) :: lowlim_temp_NEMO
  REAL(pwp) :: upplim_temp_NEMO
  REAL(pwp) :: lowlim_sal_NEMO
  REAL(pwp) :: upplim_sal_NEMO
  REAL(pwp) :: lowlim_uvel_NEMO
  REAL(pwp) :: upplim_uvel_NEMO
  REAL(pwp) :: lowlim_vvel_NEMO
  REAL(pwp) :: upplim_vvel_NEMO

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
    USE par_oce, ONLY: jpiglo, jpjglo, jpk, jpi, jpj, jpnij, &
         jp_tem, jp_sal, jpni, jpnj
    USE dom_oce, ONLY: nldi, nldj, nlei, nlej, nimpp, njmpp, rdt, &
         nimppt, njmppt, gphit, glamt, narea, tmask
    USE mod_parallel_pdaf, ONLY: jpiglo_par, jpjglo_par, jpk_par, &
         jpi_par, jpj_par, jpnij_par, jpni_par, jpnj_par, jpiglo_child, &
         jpjglo_child, jpk_child, nldi_par, narea_par, nldj_par, nlei_par, &
         nlej_par, nimpp_par, njmpp_par, nldi_child, nldj_child, &
         nlei_child, nlej_child, nimpp_child, njmpp_child, narea_child, &
         jpi_child, jpj_child, jpnij_child, jpni_child, jpnj_child
    IMPLICIT NONE


#if defined key_agrif
    IF (Agrif_Root()) THEN
       jpi_par = jpi
       jpj_par = jpj
       jpnij_par = jpnij
       jpnj_par = jpnj
       jpni_par = jpni
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
       narea_par = narea

       ALLOCATE(nimppt_par(jpnij_par))
       ALLOCATE(njmppt_par(jpnij_par))
       ALLOCATE(gphit_par(jpi_par,jpj_par))
       ALLOCATE(glamt_par(jpi_par,jpj_par))
       ALLOCATE(tmask_par(jpi_par,jpj_par,jpk_par))
       nimppt_par(:) = nimppt(:)
       njmppt_par(:) = njmppt(:)
       gphit_par(:,:) = gphit(:,:)
       glamt_par(:,:) = glamt(:,:)
       tmask_par(:,:,:) = tmask(:,:,:)
    ELSE
       jpi_child = jpi
       jpj_child = jpj
       jpnij_child = jpnij
       jpnj_child = jpnj
       jpni_child = jpni
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
       narea_child = narea

       ALLOCATE(nimppt_child(jpnij_child))
       ALLOCATE(njmppt_child(jpnij_child))
       ALLOCATE(gphit_child(jpi_child,jpj_child))
       ALLOCATE(glamt_child(jpi_child,jpj_child))
       ALLOCATE(tmask_child(jpi_child,jpj_child,jpk_child))
       nimppt_child(:) = nimppt(:)
       njmppt_child(:) = njmppt(:)
       gphit_child(:,:) = gphit(:,:)
       glamt_child(:,:) = glamt(:,:)
       tmask_child(:,:,:) = tmask(:,:,:)
    ENDIF
#else
    jpi_par = jpi
    jpj_par = jpj
    jpnij_par = jpnij
    jpnj_par = jpnj
    jpni_par = jpni
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
    narea_par = narea
    ALLOCATE(nimppt_par(jpnij_par))
    ALLOCATE(njmppt_par(jpnij_par))
    ALLOCATE(gphit_par(jpi_par,jpj_par))
    ALLOCATE(glamt_par(jpi_par,jpj_par))
    ALLOCATE(tmask_par(jpi_par,jpj_par,jpk_par))
    nimppt_par(:) = nimppt(:)
    njmppt_par(:) = njmppt(:)
    gphit_par(:,:) = gphit(:,:)
    glamt_par(:,:) = glamt(:,:)
    tmask_par(:,:,:) = tmask(:,:,:)
#endif

  END SUBROUTINE calc_grid_cnst

  SUBROUTINE fill2d_statevector(state_p, halo_p, grid)

    ! !DESCRIPTION:
    ! Fill state vector on a grid with 2d state variable.

    ! !USES:
    USE par_kind
    USE oce, ONLY: sshb
    USE mod_parallel_pdaf, &
         ONLY: nldj_child, nldj_par, nldi_child, nldi_par, task_id, &
         n_modeltasks
    USE mod_statevector_pdaf, &
         ONLY: ssh_p_offset_par, ssh_p_offset_child, mpi_subd_lat_child, &
         mpi_subd_lat_par, mpi_subd_lon_child, mpi_subd_lon_par

    IMPLICIT NONE

    ! !ARGUMENTS
    REAL(pwp), POINTER, INTENT(inout) :: state_p(:,:)  ! PE-local model state
    REAL(pwp), INTENT(inout) :: halo_p(:,:,:)          ! PE-local model halo
    CHARACTER(len=*), INTENT(in) :: grid      ! Flag for parent/child grid

    ! *** local variables ***
    INTEGER :: i, j        ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain
    INTEGER :: mem         ! Ensemble member


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

       ! Fill halo array with 2d variables. (Needed for observation operator.)
       !
       ! ******************************
       ! Explanation of indexing method
       ! ******************************
       !
       ! We wish to include halo regions right of/above of the PE. We start counting
       ! in the bottom right counter and proceed for the entire column (there are
       ! mpi_subd_lat_par such points). We then move to the top left corner and
       ! continue counting from here (the index for counting will begin at mpi_subd_
       ! lat_par + 1). We proceed for the entire row, resulting in mpi_subd_lat_par +
       ! mpi_subd_lon_par points. Finally, we include the point in the top right hand
       ! corner as the final point.
       !
       ! ******************************
       !
       ! SSH : 1

       ! Determine ensemble member
       mem=task_id

       ! Begin in bottom right corner
       DO j = 1, mpi_subd_lat_par
          halo_p(j,mem,1) = sshb(mpi_subd_lon_par+i0+1,j+j0)
       END DO
       ! Continue from top left corner
       DO i = 1, mpi_subd_lon_par
          halo_p(i+mpi_subd_lat_par,mem,1) = sshb(i+i0,mpi_subd_lat_par+j0+1)
       END DO
       ! Finish in top right corner
       halo_p(mpi_subd_lon_par+mpi_subd_lat_par+1,mem,1) = &
            sshb(mpi_subd_lon_par+i0+1, mpi_subd_lat_par+j0+1)

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

       ! Fill halo array with 2d variables. (Needed for observation operator.)
       !
       ! ******************************
       ! Explanation of indexing method
       ! ******************************
       !
       ! We wish to include halo regions right of/above of the PE. We start counting
       ! in the bottom right counter and proceed for the entire column (there are
       ! mpi_subd_lat_child such points). We then move to the top left corner and
       ! continue counting from here (the index for counting will begin at mpi_subd_
       ! lat_child + 1). We proceed for the entire row, resulting in mpi_subd_lat_child +
       ! mpi_subd_lon_child points. Finally, we include the point in the top right hand
       ! corner as the final point.
       !
       ! ******************************
       !
       ! SSH : 1

       ! Determine ensemble member
       mem=task_id

       ! Begin in bottom right corner
       DO j = 1, mpi_subd_lat_child
          halo_p(j,mem,1) = sshb(mpi_subd_lon_child+i0+1,j+j0)
       END DO
       ! Continue from top left corner
       DO i = 1, mpi_subd_lon_child
          halo_p(i+mpi_subd_lat_child,mem,1) = sshb(i+i0,mpi_subd_lat_child+j0+1)
       END DO
       ! Finish in top right corner
       halo_p(mpi_subd_lon_child+mpi_subd_lat_child+1,mem,1) = &
            sshb(mpi_subd_lon_child+i0+1, mpi_subd_lat_child+j0+1)
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
             ! Check to return sensible values
             IF( (state_p(i+(j-1)*mpi_subd_lon_par + ssh_p_offset_par, 1) > &
                  lowlim_ssh_NEMO) .AND.  (state_p(i+(j-1)*mpi_subd_lon_par + &
                  ssh_p_offset_par, 1) < upplim_ssh_NEMO) ) THEN
                sshb(i+i0,j+j0) = state_p(i+(j-1)*mpi_subd_lon_par &
                     + ssh_p_offset_par, 1)
             END IF
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
             ! Check to return sensible values
             IF( (state_p(i+(j-1)*mpi_subd_lon_child + ssh_p_offset_child, 1) > &
                  lowlim_ssh_NEMO) .AND.  (state_p(i+(j-1)*mpi_subd_lon_child + &
                  ssh_p_offset_child, 1) < upplim_ssh_NEMO) ) THEN
                sshb(i+i0,j+j0) = state_p(i+(j-1)*mpi_subd_lon_child &
                     + ssh_p_offset_child, 1)
             END IF
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
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + t_p_offset_par, 1) &
                     > lowlim_temp_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                        (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + t_p_offset_par, 1) &
                        < upplim_temp_NEMO) THEN
                      tsb(i+i0,j+j0,k,jp_tem_par) = state_p(i+(j-1)*mpi_subd_lon_par &
                           + (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + t_p_offset_par, 1)
                   END IF
                END IF
             END DO
          END DO
       END DO
       ! S
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + s_p_offset_par, 1) &
                     > lowlim_sal_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                        (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + s_p_offset_par, 1) &
                        < upplim_sal_NEMO) THEN
                     tsb(i+i0,j+j0,k,jp_sal_par) = &
                          state_p(i+(j-1)*mpi_subd_lon_par + &
                          (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + s_p_offset_par, 1)
                  END IF
               END IF
             END DO
          END DO
       END DO
       ! U
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + u_p_offset_par, 1) &
                     > lowlim_uvel_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                        (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + u_p_offset_par, 1) &
                        < upplim_uvel_NEMO) THEN
                      ub(i+i0,j+j0,k) = &
                           state_p(i+(j-1)*mpi_subd_lon_par + &
                           (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + u_p_offset_par, 1)
                   END IF
                END IF
             END DO
          END DO
       END DO
       ! V
       DO k = 1, mpi_subd_vert_par
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                     (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + v_p_offset_par, 1) &
                     > lowlim_vvel_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_par + &
                        (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + v_p_offset_par, 1) &
                        < upplim_vvel_NEMO) THEN
                      vb(i+i0,j+j0,k) = &
                           state_p(i+(j-1)*mpi_subd_lon_par + &
                           (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + v_p_offset_par, 1)
                   END IF
                END IF
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
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + t_p_offset_child, 1) &
                     > lowlim_temp_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                        (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + t_p_offset_child, 1) &
                        < upplim_temp_NEMO) THEN
                      tsb(i+i0,j+j0,k,jp_tem_child) = &
                           state_p(i+(j-1)*mpi_subd_lon_child + &
                           (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + t_p_offset_child, 1)
                   END IF
                END IF
             END DO
          END DO
       END DO
       ! S
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + s_p_offset_child, 1) &
                     > lowlim_sal_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                        (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + s_p_offset_child, 1) &
                        < upplim_sal_NEMO) THEN
                      tsb(i+i0,j+j0,k,jp_sal_child) = &
                           state_p(i+(j-1)*mpi_subd_lon_child + &
                           (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + s_p_offset_child, 1)
                   END IF
                END IF
             END DO
          END DO
       END DO
       ! U
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + u_p_offset_child, 1) &
                     > lowlim_uvel_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                        (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + u_p_offset_child, 1) &
                        < upplim_uvel_NEMO) THEN
                      ub(i+i0,j+j0,k) = &
                           state_p(i+(j-1)*mpi_subd_lon_child + &
                           (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + u_p_offset_child, 1)
                   END IF
                END IF
             END DO
          END DO
       END DO
       ! V
       DO k = 1, mpi_subd_vert_child
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                ! Check to return sensible values
                IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                     (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + v_p_offset_child, 1) &
                     > lowlim_vvel_NEMO) THEN
                   IF( state_p(i+(j-1)*mpi_subd_lon_child + &
                        (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + v_p_offset_child, 1) &
                        < upplim_vvel_NEMO) THEN
                      vb(i+i0,j+j0,k) = &
                           state_p(i+(j-1)*mpi_subd_lon_child + &
                           (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + v_p_offset_child, 1)
                   END IF
                END IF
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
