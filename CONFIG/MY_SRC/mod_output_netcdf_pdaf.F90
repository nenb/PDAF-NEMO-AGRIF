!$Id: mod_output_netcdf_pdaf.F90 61 2019-02-01 08:49:36Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_output_netcdf_pdaf
!$AGRIF_DO_NOT_TREAT
! !DESCRIPTION:
! This module provides routines to initialize
! NetCDF output files for the assimilation
! and to write output into the files.
!
! !REVISION HISTORY:
! 2010-01 - Lars Nerger - Initial code
! Later revisions - see SVN log
!
! !USES:
  USE mod_kind_pdaf
  USE mod_parallel_pdaf, &
       ONLY: abort_parallel

  IMPLICIT NONE
  SAVE
  PUBLIC

! !PUBLIC DATA MEMBERS:
! Name of the NetCDF output file
  CHARACTER(len=100) :: file_asml = 'asml.nc'
  LOGICAL :: write_stats     = .FALSE.       ! Whether to write ensemble statistics
  LOGICAL :: write_ens       = .TRUE.        ! Whether to write full ensemble
  INTEGER :: output_lev=10                   ! Number of vertical levels to ouput
                                             ! in order 1 -> max(output_lev)
!EOP

! Private variables
  INTEGER, PRIVATE :: fileid     ! Id of netcdf file

CONTAINS
!BOP
!
! !ROUTINE: init_netcdf_asml  --- initialize netcdf output
!
! !INTERFACE:
  SUBROUTINE init_netcdf_asml(step, dt, dimx, dimy, dimz, output_var, &
       output_dim, filtertype, subtype, dim_ens, forget, local_range, &
       locweight, srange, rms_obs, delt_obs, total_steps, dim_lag)

! !DESCRIPTION:
! This routine initializes the netcdf file

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN) :: step       ! Initial time step
    REAL(pwp), INTENT(IN)    :: dt         ! Size of time step
    INTEGER, INTENT(IN) :: dimx       ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy       ! Dimension: latitude
    INTEGER, INTENT(IN) :: dimz       ! Dimension: vertical
    CHARACTER(len=*), INTENT(IN) :: output_var  ! State variable to be written to file
    CHARACTER(len=2), INTENT(IN) :: output_dim  ! Dimension of state variable (2D/3D)
    INTEGER, INTENT(IN) :: filtertype ! Type of filter
    INTEGER, INTENT(IN) :: subtype    ! Sub-type of filter
    INTEGER, INTENT(IN) :: dim_ens    ! ensemble_size
    REAL(pwp), INTENT(IN)    :: forget     ! forgetting factor
    REAL(pwp), INTENT(IN)    :: local_range  ! Localization radius
    INTEGER, INTENT(IN) :: locweight    ! Type of localization
    REAL(pwp), INTENT(IN)    :: rms_obs      ! RMS error of observations
    REAL(pwp), INTENT(IN)    :: srange       ! Support range for 5th order polynomial
                                        !   and range for 1/e for exponential weighting
    INTEGER, INTENT(IN) :: delt_obs     ! Number of time steps between two analysis steps
    INTEGER, INTENT(IN) :: total_steps  ! Total number of time steps in experiment
    INTEGER, INTENT(IN) :: dim_lag             ! Size of lag for smoothing
!EOP

! Local variables
    INTEGER :: i, s                   ! Counters
    INTEGER :: dimid_1                ! Dimension IDs
    INTEGER :: dimid_x, dimid_y       ! Dimension IDs
    INTEGER :: dimid_z                ! Dimension IDs
    INTEGER :: dimid_step             ! Dimension ID
    INTEGER :: dimid_ens, dimid_ensp1 ! Dimension IDs
    INTEGER :: dimid_lag              ! Dimension ID
    INTEGER :: id_tmp                 ! Variable IDs
    INTEGER :: stat(250)              ! Array for status flag
    INTEGER :: dimarray(2)            ! Array for dimensions
    CHARACTER(len=100) :: fname       ! Name of netcdf file
    CHARACTER(len=100) :: attstr      ! String to write attributes

! *** Initialize file ***

! Print screen information
    WRITE (*, '(/1x, a, 3x, a)') 'Initialize NetCDF output for', trim(output_var)
    IF (write_ens) WRITE (*,'(5x,a,3x,a)') '--> Write ensemble states for',&
         trim(output_var)
    IF (write_stats) WRITE (*,'(5x,a,3x,a)') '--> Write higher-order ensemble&
         &statistics for', trim(output_var)

! Initialize name for netcdf file
    WRITE(fname,'(a)') TRIM(output_var)//'_'//TRIM(file_asml)

! Overwrite existing files if present (CLOBBER)
    s = 1
    stat(s) = NF90_CREATE(TRIM(fname), NF90_CLOBBER, fileid)
    s = s + 1

    attstr  = 'Assimilation performed using PDAF-NEMO-AGRIF model'
    stat(s) = NF90_PUT_ATT(fileid, NF90_GLOBAL, 'title', &
         TRIM(attstr))
    s = s + 1

! Define Dimensions
    stat(s) = NF90_DEF_DIM(fileid, 'x', dimx, dimid_x)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'y', dimy, dimid_y)
    s = s + 1
    IF (output_dim == '3D') THEN
       stat(s) = NF90_DEF_DIM(fileid, 'depth', dimz, dimid_z)
       s = s + 1
    END IF
    stat(s) = NF90_DEF_DIM(fileid, 'mem', dim_ens, dimid_ens)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'iteration', NF90_UNLIMITED, dimid_step)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'one', 1, dimid_1)
    s = s + 1
    stat(s) = NF90_DEF_DIM(fileid, 'dim_ensp1', dim_ens+1, dimid_ensp1)
    s = s + 1
    dimlag_switch: IF (dim_lag > 0) THEN
       stat(s) = NF90_DEF_DIM(fileid, 'dim_lag', dim_lag, dimid_lag)
       s = s + 1
    END IF dimlag_switch

! Define variables characterizing the experiment
    stat(s) = NF90_DEF_VAR(fileid, 'filtertype', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'subtype', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'dim_ens', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'forget', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'step_null', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'total_steps', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'local_range', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'locweight', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'srange', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'rms_obs', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'delt_obs', NF90_INT, DimId_1, Id_tmp)
    s = s + 1

! Define variables
    stat(s) = NF90_DEF_VAR(fileid, 'mrmse_for_null', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'mrmse_ana_null', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'mtrmse_for_null', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'mtrmse_ana_null', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'mrmse_for_step', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'mrmse_ana_step', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'mtrmse_for_step', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'mtrmse_ana_step', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1

    smootherA: IF (dim_lag > 0) THEN
       stat(s) = NF90_DEF_VAR(fileid, 'mrmse_smoother_null', NF90_DOUBLE, DimId_lag, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'mtrmse_smoother_null', NF90_DOUBLE, DimId_lag, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'mrmse_smoother_step', NF90_DOUBLE, DimId_lag, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'mtrmse_smoother_step', NF90_DOUBLE, DimId_lag, Id_tmp)
       s = s + 1
    END IF smootherA

    stat(s) = NF90_DEF_VAR(fileid, 'step', NF90_INT, DimId_step, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'step_ini', NF90_INT, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'time', NF90_DOUBLE, DimId_step, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'time_ini', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'rmse_ini', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'rmse_for', NF90_DOUBLE, DimId_step, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'rmse_ana', NF90_DOUBLE, DimId_step, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'trmse_ini', NF90_DOUBLE, DimId_1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'trmse_for', NF90_DOUBLE, DimId_step, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'trmse_ana', NF90_DOUBLE, DimId_step, Id_tmp)
    s = s + 1

    smootherB: IF (dim_lag > 0) THEN
       dimarray(1) = dimid_lag
       dimarray(2) = dimid_step
       stat(s) = NF90_DEF_VAR(fileid, 'rmse_smoother', NF90_DOUBLE, dimarray, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'trmse_smoother', NF90_DOUBLE, dimarray, Id_tmp)
       s = s + 1
    END IF smootherB

    stat(s) = NF90_DEF_VAR(fileid, 'hist_true_null', NF90_INT, Dimid_ensp1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'hist_mean_null', NF90_INT, Dimid_ensp1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'hist_true_step', NF90_INT, Dimid_ensp1, Id_tmp)
    s = s + 1
    stat(s) = NF90_DEF_VAR(fileid, 'hist_mean_step', NF90_INT, Dimid_ensp1, Id_tmp)
    s = s + 1

    writestats: IF (write_stats) THEN
       stat(s) = NF90_DEF_VAR(fileid, 'skewness_ini', NF90_DOUBLE, DimId_1, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'kurtosis_ini', NF90_DOUBLE, DimId_1, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'skewness_for', NF90_DOUBLE, DimId_step, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'kurtosis_for', NF90_DOUBLE, DimId_step, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'skewness_ana', NF90_DOUBLE, DimId_step, Id_tmp)
       s = s + 1
       stat(s) = NF90_DEF_VAR(fileid, 'kurtosis_ana', NF90_DOUBLE, DimId_step, Id_tmp)
       s = s + 1
    END IF writestats

    writeens: IF (write_ens) THEN

       IF (output_dim == '2D') THEN
          ! Initialise the ensemble for NEMO 2D state variable
          CALL init_2dens(dimid_x, dimid_y, dimid_ens, dimid_1, dimid_step, &
               output_var)
       ELSE IF (output_dim == '3D') THEN
          ! Initialise the ensemble for NEMO 3D state variable
          CALL init_3dens(dimid_x, dimid_y, dimid_z, dimid_ens, dimid_1, &
               dimid_step, output_var)
       END IF

    END IF writeens

    stat(s) = NF90_ENDDEF(fileid)
    s = s + 1

! Write variables characterizing the experiment
    stat(s) = NF90_INQ_VARID(fileid, 'filtertype', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, filtertype)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'subtype', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, subtype)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'dim_ens', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, dim_ens)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'forget', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, forget)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'step_null', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, step)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'total_steps', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, total_steps)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'local_range', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, local_range)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'locweight', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, locweight)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'srange', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, srange)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'rms_obs', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, rms_obs)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, 'delt_obs', Id_tmp)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_tmp, delt_obs)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in file initialization, no.', i
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE init_netcdf_asml

  SUBROUTINE init_2dens(id_x, id_y,id_ens, id_1, id_step, statevar)

! !DESCRIPTION:
! ! This routine initializes the ensemble for different 2D NEMO state variables

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN)    :: id_x        ! ID for x dimension
    INTEGER, INTENT(IN)    :: id_y        ! ID for y dimension
    INTEGER, INTENT(IN)    :: id_ens      ! ID for ensemble
    INTEGER, INTENT(IN)    :: id_1        ! ID for initial state (dim=1)
    INTEGER, INTENT(IN)    :: id_step     ! ID for timestep
    CHARACTER(len=*), INTENT(IN) :: statevar     ! Name of state variable

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: dimarray(4)             ! Array for dimensions
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: id_statevar             ! ID for state variable
    CHARACTER(len=100) :: statevar_ana ! Label for analysis
    CHARACTER(len=8)   :: anastr       ! Dummy string for label subscript


    anastr='_ens_ana'
    WRITE(statevar_ana,'(a)') TRIM(statevar)//TRIM(anastr)

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_ens
    dimarray(4) = id_step
    s=1
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ana), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in NEMO ensemble variable initialization:', statevar
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE init_2dens

  SUBROUTINE init_3dens(id_x, id_y, id_z, id_ens, id_1, id_step, statevar)

! !DESCRIPTION:
! ! This routine initializes the ensemble for different 3D NEMO state variables

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    INTEGER, INTENT(IN)    :: id_x        ! ID for x dimension
    INTEGER, INTENT(IN)    :: id_y        ! ID for y dimension
    INTEGER, INTENT(IN)    :: id_z        ! ID for vertical dimension
    INTEGER, INTENT(IN)    :: id_ens      ! ID for ensemble
    INTEGER, INTENT(IN)    :: id_1        ! ID for initial state (dim=1)
    INTEGER, INTENT(IN)    :: id_step     ! ID for timestep
    CHARACTER(len=*), INTENT(IN) :: statevar     ! Name of state variable

! Local variables
    INTEGER :: i, s                    ! Counters
    INTEGER :: dimarray(5)             ! Array for dimensions
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: id_statevar             ! ID for state variable
    CHARACTER(len=100) :: statevar_ana ! Label for analysis
    CHARACTER(len=8)   :: anastr       ! Dummy string for label subscript


    anastr='_ens_ana'
    WRITE(statevar_ana,'(a)') TRIM(statevar)//TRIM(anastr)

    dimarray(1) = id_x
    dimarray(2) = id_y
    dimarray(3) = id_z
    dimarray(4) = id_ens
    dimarray(5) = id_step
    s=1
    stat(s) = NF90_DEF_VAR(fileid, TRIM(statevar_ana), NF90_DOUBLE, dimarray, id_statevar)
    s = s + 1

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in NEMO ensemble variable initialization:', statevar
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE init_3dens

!BOP
!
! !ROUTINE: write_netcdf_asml  --- write netcdf output during assimilation
!
! !INTERFACE:
  SUBROUTINE write_netcdf_asml(calltype, dimx, dimy, dimz, file_pos, lev_tot, &
       output_var, output_dim, rmse, trmse, mrmse_null, mtrmse_null, &
       mrmse_step, mtrmse_step, dim_ens, ens, hist_true, hist_mean, skewness, &
       kurtosis, dim_lag, rmse_s, trmse_s, mrmse_s_null, mtrmse_s_null, &
       mrmse_s_step, mtrmse_s_step)

! !DESCRIPTION:
! This routine writes the netcdf file.

! !USES:
    USE netcdf

    IMPLICIT NONE

! !ARGUMENTS:
    CHARACTER(len=3)    :: calltype    ! Type of output call
    INTEGER, INTENT(IN) :: dimx        ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy        ! Dimension: latitude
    INTEGER, INTENT(IN) :: dimz        ! Dimension: vertical
    INTEGER, INTENT(IN) :: file_pos    ! Dimension: time
    INTEGER, INTENT(IN) :: lev_tot     ! Total number of vertical levels to output
    CHARACTER(len=*), INTENT(IN) :: output_var  ! State variable to be written to file
    CHARACTER(len=2), INTENT(IN) :: output_dim  ! Dimension of state variable (2D/3D)
    REAL(pwp), INTENT(IN)    :: rmse        ! Estimated RMS error
    REAL(pwp), INTENT(IN)    :: trmse       ! True RMS error
    REAL(pwp), INTENT(IN)    :: mrmse_null  ! Time-mean estimated RMS error from step 0
    REAL(pwp), INTENT(IN)    :: mtrmse_null ! Time-mean true RMS error from step 0
    REAL(pwp), INTENT(IN)    :: mrmse_step  ! Time-mean estimated RMS error from stepnull_means
    REAL(pwp), INTENT(IN)    :: mtrmse_step ! Time-mean true RMS error from stepnull_means
    INTEGER, INTENT(IN) :: dim_ens     ! Ensemble size
    REAL(pwp), INTENT(IN)    :: ens(dimx, dimy, dimz, dim_ens) ! Ensemble
    INTEGER, INTENT(IN) :: hist_true(dim_ens+1, 2) ! Rank histogram about true state
    INTEGER, INTENT(IN) :: hist_mean(dim_ens+1, 2) ! Rank histogram about ensemble mean
    REAL(pwp), INTENT(IN)    :: skewness                ! Skewness of ensemble
    REAL(pwp), INTENT(IN)    :: kurtosis                ! Kurtosis of ensemble
    ! RMS errors for smoother
    INTEGER, INTENT(IN) :: dim_lag             ! Size of lag for smoothing
    REAL(pwp), INTENT(IN) :: rmse_s(dim_lag)        ! Estimated RMS error
    REAL(pwp), INTENT(IN) :: trmse_s(dim_lag)       ! True RMS error
    REAL(pwp), INTENT(IN) :: mrmse_s_null(dim_lag)  ! Time-mean estimated RMS error from step 0
    REAL(pwp), INTENT(IN) :: mtrmse_s_null(dim_lag) ! Time-mean true RMS error from step 0
    REAL(pwp), INTENT(IN) :: mrmse_s_step(dim_lag)  ! Time-mean estimated RMS error from stepnull_means
    REAL(pwp), INTENT(IN) :: mtrmse_s_step(dim_lag) ! Time-mean true RMS error from stepnull_means
!EOP

! Local variables
    INTEGER :: i, s, idx               ! Counters
    INTEGER :: ID_time, ID_step        ! Variable IDs
    INTEGER :: ID_rmse, ID_trmse       ! Variable IDs
    INTEGER :: ID_mrmseN, ID_mtrmseN   ! Variable IDs
    INTEGER :: ID_mrmseS, ID_mtrmseS   ! Variable IDs
    INTEGER :: ID_hist_true_null, ID_hist_mean_null ! Variable IDs
    INTEGER :: ID_hist_true_step, ID_hist_mean_step ! Variable IDs
    INTEGER :: ID_skew, ID_kurt        ! Variable IDs
    INTEGER :: Id_rmse_s, Id_trmse_s, Id_mrmseN_s      ! Variable IDs for smoother output
    INTEGER :: Id_mtrmseN_s, Id_mrmseS_s, Id_mtrmseS_s ! Variable IDs for smoother output
    INTEGER :: stat(250)               ! Array for status flag
    INTEGER :: pos(1)                  ! Position index for writing
    INTEGER :: pos2(2)                 ! Position index for writing
    INTEGER :: cnt2(2)                 ! Count index for writing
    INTEGER :: id_2dens                ! 2D variable ID for ensemble
    INTEGER :: id_3dens                ! 3D variable ID for ensemble
    CHARACTER(len=100)   :: fname           ! Name of netcdf file
    CHARACTER(len=100)   :: varstr          ! State variable name
    CHARACTER(len=8)     :: anastr          ! Dummy string for variable name subscript


    ! Initialize name for netcdf file
    WRITE(fname,'(a)') TRIM(output_var)//'_'//TRIM(file_asml)

    s = 1
    stat(s) = NF90_OPEN(TRIM(fname), NF90_WRITE, fileid)

    WRITE (*,'(/9x, a, 3x, a)') 'Writing to file:', fname

    IF (stat(1) .NE. NF90_NOERR) THEN
       WRITE(*,'(/9x, a, 3x, a)') &
            'NetCDF error in opening assimilation file:', fname
       CALL abort_parallel()
    END IF

! Inquire variable Ids
    s = 1

    IF (write_ens) THEN
       IF (output_dim == '2D') THEN
          ! Inquire 2D variable ID
          anastr='_ens_ana'
          WRITE(varstr,'(a)') TRIM(output_var)//TRIM(anastr)

          stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_2dens)
          s = s + 1
       ELSE IF (output_dim == '3D') THEN
          ! Inquire 3D variable ID
          anastr='_ens_ana'
          WRITE(varstr,'(a)') TRIM(output_var)//TRIM(anastr)

          stat(s) = NF90_INQ_VARID(fileid, TRIM(varstr), Id_3dens)
          s = s + 1
       END IF
    END IF

    stat(s) = NF90_INQ_VARID(fileid, "rmse_ana", Id_rmse)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "trmse_ana", Id_trmse)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "mrmse_ana_null", Id_mrmseN)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "mtrmse_ana_null", Id_mtrmseN)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "mrmse_ana_step", Id_mrmseS)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "mtrmse_ana_step", Id_mtrmseS)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "hist_true_null", Id_hist_true_null)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "hist_mean_null", Id_hist_mean_null)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "hist_true_step", Id_hist_true_step)
    s = s + 1
    stat(s) = NF90_INQ_VARID(fileid, "hist_mean_step", Id_hist_mean_step)
    s = s + 1
    writestatsC: IF (write_stats) THEN
       stat(s) = NF90_INQ_VARID(fileid, "skewness_ana", Id_skew)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "kurtosis_ana", Id_kurt)
       s = s + 1
    END IF writestatsC

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in preparing output, no.', i
          CALL abort_parallel()
       END IF
    END DO
    s = 1

    smootherA: IF (dim_lag > 0 .AND. calltype=='ana') THEN
       stat(s) = NF90_INQ_VARID(fileid, "rmse_smoother", Id_rmse_s)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "trmse_smoother", Id_trmse_s)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "mrmse_smoother_null", Id_mrmseN_s)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "mtrmse_smoother_null", Id_mtrmseN_s)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "mrmse_smoother_step", Id_mrmseS_s)
       s = s + 1
       stat(s) = NF90_INQ_VARID(fileid, "mtrmse_smoother_step", Id_mtrmseS_s)
       s = s + 1
    END IF smootherA

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in preparing smoother output, no.', i
          CALL abort_parallel()
       END IF
    END DO
    s = 1

    ! Write ensemble only at specified intervals
    IF (write_ens) THEN
       IF (output_dim == '2D') THEN
          ! Write 2D ensemble
          CALL write_2dens(Id_2dens, dim_ens, dimx, dimy, ens, &
               file_pos)
       ELSE IF (output_dim == '3D') THEN
          ! Write 3D ensemble
          CALL write_3dens(Id_3dens, dim_ens, dimx, dimy, lev_tot, ens, &
               file_pos)
       END IF
    END IF

    pos(1) = file_pos
    stat(s) = NF90_PUT_VAR(fileid, Id_rmse, rmse, pos)
    s = s + 1

    pos(1) = file_pos
    stat(s) = NF90_PUT_VAR(fileid, Id_trmse, trmse, pos)
    s = s + 1

    ! Write RMS errors from smoothing
    smootherB: IF (dim_lag > 0 .AND. calltype=='ana') THEN
       pos2(1) = 1
       pos2(2) = file_pos
       cnt2(1) = dim_lag
       cnt2(2) = 1
       stat(s) = NF90_PUT_VAR(fileid, Id_rmse_s, rmse_s, pos2, cnt2)
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_trmse_s, trmse_s, pos2, cnt2)
       s = s + 1
    END IF smootherB

    IF (write_stats) THEN
       ! Ensemble statistics
       pos(1) = file_pos
       stat(s) = NF90_PUT_VAR(fileid, Id_skew, skewness, pos)
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_kurt, kurtosis, pos)
       s = s + 1
    END IF


    IF (calltype=='ana') THEN
       ! Histogram information - written at each time
       stat(s) = NF90_PUT_VAR(fileid, Id_hist_true_null, hist_true(:, 1))
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_hist_mean_null, hist_mean(:, 1))
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_hist_true_step, hist_true(:, 2))
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_hist_mean_step, hist_mean(:, 2))
       s = s + 1
    END IF

    ! Mean errors are written at each time
    stat(s) = NF90_PUT_VAR(fileid, Id_mrmseN, mrmse_null)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_mtrmseN, mtrmse_null)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_mrmseS, mrmse_step)
    s = s + 1
    stat(s) = NF90_PUT_VAR(fileid, Id_mtrmseS, mtrmse_step)
    s = s + 1

    ! Write mean smoother errors at each time
    smootherC: IF (dim_lag > 0 .AND. calltype=='ana') THEN
       stat(s) = NF90_PUT_VAR(fileid, Id_mrmseN_s, mrmse_s_null)
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_mtrmseN_s, mtrmse_s_null)
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_mrmseS_s, mrmse_s_step)
       s = s + 1
       stat(s) = NF90_PUT_VAR(fileid, Id_mtrmseS_s, mtrmse_s_step)
       s = s + 1
    END IF smootherC

    DO i = 1,  s - 1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in writing output, no.', i
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE write_netcdf_asml

  SUBROUTINE write_2dens(idvar, dim_ens, dimx, dimy, ens, iter)

! !DESCRIPTION:
! ! This routine writes 2D NEMO ensemble variables to file

! !USES:
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(IN) :: idvar        ! 2D variable ID
    INTEGER, INTENT(IN) :: dim_ens      ! Ensemble size
    INTEGER, INTENT(IN) :: dimx         ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy         ! Dimension: latitude
    REAL(pwp), INTENT(IN)    :: ens(:,:,:,:) ! Ensemble
    INTEGER, INTENT(IN) :: iter         ! Iteration number

    ! Local variables
    INTEGER :: i, j, k, s                     ! Counters
    INTEGER :: stat(20)                       ! Array for status flag
    INTEGER :: pos4(4)                        ! Position index for writing
    INTEGER :: cnt4(4)                        ! Count index for writing
    REAL(pwp), ALLOCATABLE :: tmp_array(:,:,:)     ! Array for 3D field select


    ! Select 3D field from 4D array
    ALLOCATE(tmp_array(dimx,dimy,dim_ens))
    DO k = 1, dim_ens
       DO j = 1, dimy
          DO i = 1, dimx
             tmp_array(i,j,k)=ens(i,j,1,k)
          END DO
       END DO
    END DO

    s=1
    pos4(1) = 1
    pos4(2) = 1
    pos4(3) = 1
    pos4(4) = iter
    cnt4(1) = dimx
    cnt4(2) = dimy
    cnt4(3) = dim_ens
    cnt4(4) = 1
    stat(s) = NF90_PUT_VAR(fileid, idvar, tmp_array, pos4, cnt4)
    s = s + 1

    DO i = 1,  s-1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in NEMO 2D ensemble write.'
          CALL abort_parallel()
       END IF
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE write_2dens

  SUBROUTINE write_3dens(idvar, dim_ens, dimx, dimy, lev_num, ens, iter)

! !DESCRIPTION:
! ! This routine writes 3D NEMO ensemble variables to file

! !USES:
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS:
    INTEGER, INTENT(IN) :: idvar        ! 3D variable ID
    INTEGER, INTENT(IN) :: dim_ens      ! Ensemble size
    INTEGER, INTENT(IN) :: dimx         ! Dimension: longitude
    INTEGER, INTENT(IN) :: dimy         ! Dimension: latitude
    INTEGER, INTENT(IN) :: lev_num      ! Total number of levels to output
    REAL(pwp), INTENT(IN)    :: ens(:,:,:,:) ! Ensemble
    INTEGER, INTENT(IN) :: iter         ! Iteration number

    ! Local variables
    INTEGER :: i, j, k, mem, s                  ! Counters
    INTEGER :: stat(20)                         ! Array for status flag
    INTEGER :: pos5(5)                          ! Position index for writing
    INTEGER :: cnt5(5)                          ! Count index for writing
    REAL(pwp), ALLOCATABLE :: tmp_array(:,:,:,:)     ! Array for 3D field select


    ALLOCATE(tmp_array(dimx,dimy,lev_num,dim_ens))
    DO mem = 1, dim_ens
       DO k = 1, lev_num
          DO j = 1, dimy
             DO i = 1, dimx
                tmp_array(i,j,k,mem)=ens(i,j,k,mem)
             END DO
          END DO
       END DO
    END DO

    s=1
    pos5(1) = 1
    pos5(2) = 1
    pos5(3) = 1
    pos5(4) = 1
    pos5(5) = iter
    cnt5(1) = dimx
    cnt5(2) = dimy
    cnt5(3) = lev_num
    cnt5(4) = dim_ens
    cnt5(5) = 1
    stat(s) = NF90_PUT_VAR(fileid, idvar, tmp_array, pos5, cnt5)
    s = s + 1

    DO i = 1,  s-1
       IF (stat(i) /= NF90_NOERR) THEN
          WRITE(*, *) 'NetCDF error in NEMO 3D ensemble write.'
          CALL abort_parallel()
       END IF
    END DO

    DEALLOCATE(tmp_array)

  END SUBROUTINE write_3dens
!BOP
!
! !ROUTINE: close_netcdf_asml  --- close netcdf file
!
! !INTERFACE:
  SUBROUTINE close_netcdf_asml()

! !DESCRIPTION:
! This routine closes the netcdf file.

! !USES:
    USE netcdf

    IMPLICIT NONE

!EOP

! Local variables
    INTEGER :: stat(50)                ! Array for status flag

! Close file

    stat(1) = NF90_CLOSE(fileid)
    IF (stat(1) /= NF90_NOERR) THEN
       WRITE(*, *) 'NetCDF error in closing assimilation file.'
       CALL abort_parallel()
    END IF

  END SUBROUTINE close_netcdf_asml
!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_output_netcdf_pdaf
