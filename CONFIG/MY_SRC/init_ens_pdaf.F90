!$Id: init_ens.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: init_ens --- Initialize ensemble
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)
!$AGRIF_DO_NOT_TREAT

! !DESCRIPTION:
! User-supplied routine for PDAF.
! Used in the filters: SEIK/LSEIK/ETKF/LETKF/ESTKF/LESTKF
!
! The routine is called when the filter is
! initialized in PDAF\_filter\_init.  It has
! to initialize an ensemble of dim\_ens states.
! Typically, the ensemble will be directly read from files.
!
! The routine is called by all filter processes and
! initializes the ensemble for the PE-local domain.
!
! Implementation for the 2D online example
! without parallelization.
!
! !REVISION HISTORY:
! 2013-02 - Lars Nerger - Initial code based on offline_1D
! Later revisions - see svn log
!
! !USES:
  USE mod_kind_pdaf
  USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens, jpiglo_par, jpjglo_par, &
       jpiglo_child, jpjglo_child, jpk_child, jpk_par
  USE mod_assimilation_pdaf, ONLY: istate_t_par, istate_u_par, istate_v_par, &
       screen, wght, istate_t_child, istate_u_child, &
       istate_v_child
  USE mod_statevector_pdaf, ONLY: fill2d_par_ensarray, fill3d_par_ensarray, &
       fill2d_child_ensarray, fill3d_child_ensarray
  USE netcdf

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL(pwp), INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK.
  ! It is available here only for convenience and can be used freely.
  REAL(pwp), INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL(pwp), INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

  ! *** local variables ***
  INTEGER :: s, i, j, var                 ! Counters
  INTEGER :: stat(20000)                  ! Status flag for NetCDF commands
  INTEGER :: ncid_in                      ! ID for NetCDF file
  INTEGER :: id_dimx,id_dimy,id_dimz      ! IDs for dimensions
  INTEGER :: dim_lon, dim_lat, dim_vert   ! Initial state dimensions
  CHARACTER(len=100) :: istate_ncfile(3)      ! Files holding initial state estimate
  CHARACTER(len=20), DIMENSION(3) :: zdim_list ! z dimension labels for ic files
  REAL(pwp) :: dimens_mean               ! Constant used in generating ens weights

  ! List of z dimension labels for initial state files
  DATA zdim_list / 'deptht', 'depthu', 'depthv' /


  ! **********************************************************************************
  ! ****************  Method for generating ensemble members  ************************
  ! **********************************************************************************
  !
  !  We define X0 as our initial state estimate. We perform a control run from this
  !  initial state estimate and save snapshots of model fields every 6h. We denote
  !  by Xi_ctrl the snapshot at time (6*i) hours of our control run.
  !
  !  For an ensemble of size N, we consider the set of Xj_ctrl for 0<=j<N. We define
  !  X_ensmean as mean(Xj_ctrl) for this set. We define the constant Wght_j as
  !  Wght_j = 1 + (N - 1 - dimens_mean) - abs(j - dimens_mean), where
  !  dimens_mean=(N-1)/2.
  !
  !  We then define our initial state for ensemble member Xj as
  !  Xj = X_ensmean + ( Wght_j * (Xj_ctrl - X_ensmean) ).
  !
  ! **********************************************************************************

  ALLOCATE(wght(dim_ens))

  dimens_mean=(REAL(dim_ens) - 1.0)/2.0

  ! Compute weighting factor for ensemble perturbations
  DO i = 1, dim_ens
     !wght(i)= (REAL(dim_ens) - dimens_mean - ABS(REAL(i) - 1.0 - dimens_mean))
     ! *WARNING:* Variable weights currently leads to NEMO crashing. Use constant
     ! weight until/if solution found.
     wght(i) = 1.0_pwp
  END DO

  IF(screen > 1 ) THEN
     IF(mype_ens == 0) WRITE (*, *) 'Ensemble weights:', wght(:)
  END IF


  ! **********************************************************************************
  !
  ! Initial state files are first checked to confirm that file dimensions consistent
  ! with dimensions for model run. Then ensemble array of state vectors is filled
  ! via call to routines in state vector module.
  !
  ! **********************************************************************************

! ***********
! Parent grid
! ***********

  ! Files holding initial state estimate
  istate_ncfile(1) = TRIM(istate_t_par)
  istate_ncfile(2) = TRIM(istate_u_par)
  istate_ncfile(3) = TRIM(istate_v_par)

  IF(screen > 1) THEN
     IF(mype_ens == 0) WRITE (*,'(/9x, a, 3x, a, 3x, a, 3x, a)')&
          "Initial state estimate files:", istate_ncfile(1),&
          istate_ncfile(2), istate_ncfile(3)
  END IF

  ! Check dimensions for each state variable file
  parent_dimcheck:DO var = 1, size(istate_ncfile)

     ! *******************************************
     ! *** Open file containing initial state ***
     ! *******************************************
     s = 1
     stat(s) = NF90_OPEN(istate_ncfile(var) , NF90_NOWRITE, ncid_in)

     DO j = 1, s
        IF (stat(j) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in opening initial state file:',&
                istate_ncfile(var)
           CALL abort_parallel()
        END IF
     END DO

     ! ************************
     ! *** Check dimensions ***
     ! ************************

     s=1
     stat(s) = NF90_INQ_DIMID(ncid_in, 'x', id_dimx)
     s = s + 1
     stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimx, len=dim_lon)
     s = s + 1
     stat(s) = NF90_INQ_DIMID(ncid_in, 'y', id_dimy)
     s = s + 1
     stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimy, len=dim_lat)
     s = s + 1
     stat(s) = NF90_INQ_DIMID(ncid_in, TRIM(zdim_list(var)), id_dimz)
     s = s + 1
     stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimz, len=dim_vert)

     DO j = 1, s
        IF (stat(j) .NE. NF90_NOERR) THEN
           WRITE(*, '(/9x, a, 3x, i1, 3x, a, 3x, a)') &
                'NetCDF error in reading dimension:', j,&
                ' from initial state file:', istate_ncfile(var)
           CALL abort_parallel()
        END IF
     END DO

     ! Compare NEMO global dimensions to dimensions in file.
     IF (dim_lon == jpiglo_par .AND. dim_lat == jpjglo_par &
          .AND. dim_vert == jpk_par) THEN
        IF (screen > 1 .AND. mype_ens == 0) WRITE (*,'(/9x, a)') &
             'Dimensions in initial state file valid.'
     ELSE
        WRITE (*,'(/9x, a, 3x, a, 3x, a)') 'ERROR: Initial state file:', &
             istate_ncfile(var), 'has invalid dimensions.'
        CALL abort_parallel()
     END IF

     ! *******************************************
     ! *** Close file containing initial state ***
     ! *******************************************

     s = 1
     stat(s) = NF90_CLOSE(ncid_in)

     DO j = 1, s
        IF (stat(j) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in closing initial state file:',&
                istate_ncfile(var)
           CALL abort_parallel()
        END IF
     END DO

  END DO parent_dimcheck


  ! ***************************************************
  ! *** Initialize ensemble array of state vectors  ***
  ! ***************************************************

  IF (mype_ens == 0) THEN
     WRITE (*,'(/1x,a)') '------- Reading Initial State -------------'
     WRITE (*,'(/1x,a)') 'Calling NEMO 2dfill_ensarray'
     WRITE (*,'(/1x,a)') 'Calling NEMO 3dfill_ensarray'
  END IF

  CALL fill2d_par_ensarray(dim_p, dim_ens, wght, istate_ncfile(1), ens_p)
  ! 'T' also fills salinity ('S') ensemble.
  CALL fill3d_par_ensarray(dim_p, dim_ens, wght, istate_ncfile(1), 'T', ens_p)
  CALL fill3d_par_ensarray(dim_p, dim_ens, wght, istate_ncfile(2), 'U', ens_p)
  CALL fill3d_par_ensarray(dim_p, dim_ens, wght, istate_ncfile(3), 'V', ens_p)


#if defined key_agrif
! ***********
! Child grid
! ***********

  ! Files holding initial state estimate
  istate_ncfile(1) = TRIM(istate_t_child)
  istate_ncfile(2) = TRIM(istate_u_child)
  istate_ncfile(3) = TRIM(istate_v_child)

  IF(screen > 1) THEN
     IF(mype_ens == 0) WRITE (*,'(/9x, a, 3x, a, 3x, a, 3x, a)')&
          "Initial state estimate files:", istate_ncfile(1),&
          istate_ncfile(2), istate_ncfile(3)
  END IF

  ! Check dimensions for each state variable file
  child_dimcheck:DO var = 1, size(istate_ncfile)

     ! *******************************************
     ! *** Open file containing initial state ***
     ! *******************************************
     s = 1
     stat(s) = NF90_OPEN(istate_ncfile(var) , NF90_NOWRITE, ncid_in)

     DO j = 1, s
        IF (stat(j) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in opening initial state file:',&
                istate_ncfile(var)
           CALL abort_parallel()
        END IF
     END DO

     ! ************************
     ! *** Check dimensions ***
     ! ************************

     s=1
     stat(s) = NF90_INQ_DIMID(ncid_in, 'x', id_dimx)
     s = s + 1
     stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimx, len=dim_lon)
     s = s + 1
     stat(s) = NF90_INQ_DIMID(ncid_in, 'y', id_dimy)
     s = s + 1
     stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimy, len=dim_lat)
     s = s + 1
     stat(s) = NF90_INQ_DIMID(ncid_in, TRIM(zdim_list(var)), id_dimz)
     s = s + 1
     stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimz, len=dim_vert)

     DO j = 1, s
        IF (stat(j) .NE. NF90_NOERR) THEN
           WRITE(*, '(/9x, a, 3x, i1, 3x, a, 3x, a)') &
                'NetCDF error in reading dimension:', j,&
                ' from initial state file:', istate_ncfile(var)
           CALL abort_parallel()
        END IF
     END DO

     ! Compare NEMO global dimensions to dimensions in file.
     IF (dim_lon == jpiglo_child .AND. dim_lat == jpjglo_child &
          .AND. dim_vert == jpk_child) THEN
        IF (screen > 1 .AND. mype_ens == 0) WRITE (*,'(/9x, a)') &
             'Dimensions in initial state file valid.'
     ELSE
        WRITE (*,'(/9x, a, 3x, a, 3x, a)') 'ERROR: Initial state file:', &
             istate_ncfile(var), 'has invalid dimensions.'
        CALL abort_parallel()
     END IF

     ! *******************************************
     ! *** Close file containing initial state ***
     ! *******************************************

     s = 1
     stat(s) = NF90_CLOSE(ncid_in)

     DO j = 1, s
        IF (stat(j) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in closing initial state file:',&
                istate_ncfile(var)
           CALL abort_parallel()
        END IF
     END DO

  END DO child_dimcheck


  ! ***************************************************
  ! *** Initialize ensemble array of state vectors  ***
  ! ***************************************************

  IF (mype_ens == 0) THEN
     WRITE (*,'(/1x,a)') '------- Reading Initial State -------------'
     WRITE (*,'(/1x,a)') 'Calling AGRIF 2dfill_ensarray'
     WRITE (*,'(/1x,a)') 'Calling AGRIF 3dfill_ensarray'
  END IF

  CALL fill2d_child_ensarray(dim_p, dim_ens, wght, istate_ncfile(1), ens_p)
  ! 'T' also fills salinity ('S') ensemble.
  CALL fill3d_child_ensarray(dim_p, dim_ens, wght, istate_ncfile(1), 'T', ens_p)
  CALL fill3d_child_ensarray(dim_p, dim_ens, wght, istate_ncfile(2), 'U', ens_p)
  CALL fill3d_child_ensarray(dim_p, dim_ens, wght, istate_ncfile(3), 'V', ens_p)
#endif

  ! ********
  ! Clean up
  ! ********

  DEALLOCATE(wght)

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_ens_pdaf
