!$Id: init_ens.F90 1589 2015-06-12 11:57:58Z lnerger $
!BOP
!
! !ROUTINE: init_ens --- Initialize ensemble
!
! !INTERFACE:
SUBROUTINE init_ens_pdaf(filtertype, dim_p, dim_ens, state_p, Uinv, &
     ens_p, flag)

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

  USE netcdf
  USE mod_assimilation, ONLY: istate_fname_t, istate_fname_u, istate_fname_v, &
       screen, wght
  USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens
  USE mod_statevector, ONLY: fill2d_ensarray, fill3d_ensarray
  USE par_oce, ONLY: jpiglo,jpjglo,jpk

  IMPLICIT NONE

! !ARGUMENTS:
  INTEGER, INTENT(in) :: filtertype              ! Type of filter to initialize
  INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens                 ! Size of ensemble
  REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state
  ! It is not necessary to initialize the array 'state_p' for SEIK.
  ! It is available here only for convenience and can be used freely.
  REAL, INTENT(inout) :: Uinv(dim_ens-1,dim_ens-1) ! Array not referenced for SEIK
  REAL, INTENT(out)   :: ens_p(dim_p, dim_ens)   ! PE-local state ensemble
  INTEGER, INTENT(inout) :: flag                 ! PDAF status flag

  ! *** local variables ***
  REAL :: dimens_mean
  INTEGER :: s, i, j                         ! Counters
  INTEGER :: stat(20000)                  ! Status flag for NetCDF commands
  INTEGER :: ncid_in                      ! ID for NetCDF file
  INTEGER :: id_dimx,id_dimy,id_dimz   ! IDs for dimensions
  INTEGER :: dim_lon, dim_lat, dim_vert ! Initial state dimensions
  CHARACTER(len=100) :: istate_ncfile(3)     ! Files holding initial state estimate


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
     wght(i)= (REAL(dim_ens) - dimens_mean - ABS(REAL(i) - 1.0 - dimens_mean))
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

  ! Files holding initial state estimate
  istate_ncfile(1) = Trim(istate_fname_t)
  istate_ncfile(2) = Trim(istate_fname_u)
  istate_ncfile(3) = Trim(istate_fname_v)

  IF(screen > 1) THEN
     IF(mype_ens == 0) WRITE (*,'(/9x, a, 3x, a, 3x, a, 3x, a)')&
          "Initial state estimate files:", istate_ncfile(1),&
          istate_ncfile(2), istate_ncfile(3)
  END IF

  ! Check dimensions for each state variable file
  file_dimcheck:DO i=1,3

     ! *******************************************
     ! *** Open file containing initial state ***
     ! *******************************************
     s = 1
     stat(s) = NF90_OPEN(istate_ncfile(i) , NF90_NOWRITE, ncid_in)

     DO j = 1, s
        IF (stat(j) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in opening initial state file:',&
                istate_ncfile(i)
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
     stat(s) = NF90_INQ_DIMID(ncid_in, 'deptht', id_dimz)
     s = s + 1
     stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimz, len=dim_vert)

     DO j = 1, s
        IF (stat(i) .NE. NF90_NOERR) THEN
           WRITE(*, '(/9x, a, 3x, i1, 3x, a, 3x, a)') &
                'NetCDF error in reading dimension:', j,&
                ' from initial state file:', istate_ncfile(i)
           CALL abort_parallel()
        END IF
     END DO

     ! Compare NEMO global dimensions to dimensions in file.
     IF (dim_lon == jpiglo .AND. dim_lat == jpjglo &
          .AND. dim_vert == jpk) THEN
        IF (screen > 1 .AND. mype_ens == 0) WRITE (*,'(/9x, a)') &
             'Dimensions in initial state file valid.'
     ELSE
        WRITE (*,'(/9x, a, 3x, a, 3x, a)') 'ERROR: Initial state file:', &
             istate_ncfile(i), 'has invalid dimensions.'
        CALL abort_parallel()
     END IF

     ! *******************************************
     ! *** Close file containing initial state ***
     ! *******************************************

     s = 1
     stat(s) = NF90_CLOSE(ncid_in)

     DO j = 1, s
        IF (stat(i) .NE. NF90_NOERR) THEN
           WRITE(*,'(/9x, a, 3x, a)') &
                'NetCDF error in closing initial state file:',&
                istate_ncfile(i)
           CALL abort_parallel()
        END IF
     END DO

  END DO file_dimcheck


  ! ***************************************************
  ! *** Initialize ensemble array of state vectors  ***
  ! ***************************************************

  IF (mype_ens == 0) THEN
     WRITE (*,'(/1x,a)') '------- Reading Initial State -------------'
     WRITE (*,'(/1x,a)') 'Calling 2dfill_ensarray'
     WRITE (*,'(/1x,a)') 'Calling 3dfill_ensarray'
  END IF

  CALL fill2d_ensarray(dim_p, dim_ens, wght, istate_ncfile(1), ens_p)
  !CALL fill3d_ensarray(dim_p, dim_ens, wght, istate_ncfile(1), ens_p)

  DEALLOCATE(wght)

END SUBROUTINE init_ens_pdaf
