!$Id: prepoststep_ens_pdaf.F90 1864 2017-12-20 19:53:30Z lnerger $
!BOP
!
! !ROUTINE: prepoststep_ens_pdaf --- Used-defined Pre/Poststep routine for PDAF
!
! !INTERFACE:
SUBROUTINE prepoststep_ens_pdaf(step, dim_p, dim_ens, dim_ens_p, dim_obs_p, &
     state_p, Uinv, ens_p, flag)

  ! !DESCRIPTION:
  ! User-supplied routine for PDAF.
  ! Used in the filters: SEIK/EnKF/LSEIK/ETKF/LETKF/ESTKF/LESTKF
  ! 
  ! The routine is called for global filters (e.g. SEIK)
  ! before the analysis and after the ensemble transformation.
  ! For local filters (e.g. LSEIK) the routine is called
  ! before and after the loop over all local analysis
  ! domains.
  ! The routine provides full access to the state 
  ! estimate and the state ensemble to the user.
  ! Thus, user-controlled pre- and poststep 
  ! operations can be performed here. For example 
  ! the forecast and the analysis states and ensemble
  ! covariance matrix can be analyzed, e.g. by 
  ! computing the estimated variances. 
  ! For the offline mode, this routine is the place
  ! in which the writing of the analysis ensemble
  ! can be performed.
  !
  ! If a user considers to perform adjustments to the 
  ! estimates (e.g. for balances), this routine is 
  ! the right place for it.
  !
  ! Implementation for the 2D offline example
  ! with parallelization.
  !
  ! !REVISION HISTORY:
  ! 2013-02 - Lars Nerger - Initial code based on offline_1D
  ! Later revisions - see svn log
  !
  ! !USES:
  USE kind_pdaf
  USE mod_assimilation, &
       ONLY: screen, filtertype, subtype, forget, local_range, &
       locweight, srange, rms_obs, delt_obs, dim_lag, iter, &
       output_ssh, output_t, output_s, output_u, output_v
  USE mod_parallel_pdaf, &
       ONLY: mype_world, mype_filter, npes_filter, COMM_filter, &
       gather_ens
  USE mod_statevector, &
       ONLY: ssh_dim_state, ssh_p_dim_state, t_dim_state, &
       t_p_dim_state, s_dim_state, s_p_dim_state, u_dim_state, &
       u_p_dim_state, v_p_dim_state, v_dim_state, ssh_p_offset, &
       t_p_offset, s_p_offset, u_p_offset, v_p_offset, mpi_subd_lat, &
       mpi_subd_lon, mpi_subd_vert
  USE output_netcdf_asml, &
       ONLY: init_netcdf_asml, write_netcdf_asml, close_netcdf_asml, &
       output_lev
  USE in_out_manager, &
       ONLY: nit000, nitend
  USE dom_oce, ONLY: rdt
  USE par_oce, ONLY: jpiglo, jpjglo, jpk

  IMPLICIT NONE

  ! !ARGUMENTS:
  INTEGER, INTENT(in) :: step        ! Current time step (negative for call after forecast)
  INTEGER, INTENT(in) :: dim_p       ! PE-local state dimension
  INTEGER, INTENT(in) :: dim_ens     ! Size of state ensemble
  INTEGER, INTENT(in) :: dim_ens_p   ! PE-local size of ensemble
  INTEGER, INTENT(in) :: dim_obs_p   ! PE-local dimension of observation vector
  REAL(pwp), INTENT(inout) :: state_p(dim_p) ! PE-local forecast/analysis state
  ! The array 'state_p' is initialised and can be used freely here (not for SEEK!)
  REAL(pwp), INTENT(inout) :: Uinv(dim_ens-1, dim_ens-1) ! Inverse of matrix U
  REAL(pwp), INTENT(inout) :: ens_p(dim_p, dim_ens)      ! PE-local state ensemble
  INTEGER, INTENT(in) :: flag        ! PDAF status flag

  ! !CALLING SEQUENCE:
  ! Called by: PDAF_get_state      (as U_prepoststep)
  ! Called by: PDAF_X_update       (as U_prepoststep)
  !EOP

  ! *** local variables ***
  INTEGER :: i, member                ! Counters
  LOGICAL, SAVE :: firsttime = .TRUE. ! Routine is called for first time?
  CHARACTER(len=20) :: var            ! State variable to be written to file
  CHARACTER(len=3)  :: ensstr         ! String for ensemble member
  CHARACTER(len=2)  :: stepstr        ! String for time step
  CHARACTER(len=3)  :: anastr         ! String for call type
  INTEGER :: rank_ssh                 ! PE rank for ssh write to file
  INTEGER :: rank_T                   ! PE rank for T write to file
  INTEGER :: rank_S                   ! PE rank for S write to file
  INTEGER :: rank_U                   ! PE rank for U write to file
  INTEGER :: rank_V                   ! PE rank for V write to file
  REAL(pwp), ALLOCATABLE :: ens_ssh(:,:,:,:) ! Global state ensemble
  REAL(pwp), ALLOCATABLE :: ens_T(:,:,:,:)   ! Global state ensemble
  REAL(pwp), ALLOCATABLE :: ens_S(:,:,:,:)   ! Global state ensemble
  REAL(pwp), ALLOCATABLE :: ens_U(:,:,:,:)   ! Global state ensemble
  REAL(pwp), ALLOCATABLE :: ens_V(:,:,:,:)   ! Global state ensemble

  ! *** local variables for netcdf files - NOT CURRENTLY USED ***
  REAL(pwp) :: invdim_ens                   ! Inverse ensemble size
  REAL(pwp) :: invdim_ensm1                 ! Inverse of ensemble size minus 1
  REAL(pwp) :: rmse_est                     ! estimated RMS error
  REAL(pwp) :: rmse_true                    ! true RMS error
  REAL(pwp), ALLOCATABLE :: variance_p(:)  ! local variance
  INTEGER, SAVE, ALLOCATABLE :: hist_true(:,:) ! Array for rank histogram about true state
  INTEGER, SAVE, ALLOCATABLE :: hist_mean(:,:) ! Array for rank histogram about ensemble mean
  ! Variables for mean errors from step 0
  REAL(pwp), SAVE :: sum_rmse_est_null(2) = 0.0  ! RMS error estimate accumulated over time
  REAL(pwp), SAVE :: sum_rmse_true_null(2) = 0.0 ! True RMS error accumulated over time
  INTEGER, SAVE :: nsum_null(2) = 0         ! Length of sums over time
  REAL(pwp) :: mrmse_est_null = 0.0              ! Time-mean of estimated RMS error
  REAL(pwp) :: mrmse_true_null = 0.0             ! Time-mean of true RMS error
  ! Variables for sum from step stepnull_means
  REAL(pwp), SAVE :: sum_rmse_est_step(2) = 0.0  ! RMS error estimate accumulated over time
  REAL(pwp), SAVE :: sum_rmse_true_step(2) = 0.0 ! True RMS error accumulated over time
  INTEGER, SAVE :: nsum_step(2) = 0         ! Length of sums over time
  REAL(pwp) :: mrmse_est_step = 0.0              ! Time-mean of estimated RMS error
  REAL(pwp) :: mrmse_true_step = 0.0             ! Time-mean of true RMS error
  REAL(pwp) :: skewness                          ! Skewness of ensemble
  REAL(pwp) :: kurtosis                          ! Kurtosis of ensemble
  ! Variables for smoother erros
  REAL(pwp), Allocatable :: rmse_s(:)        ! estimated RMS error of smoothed states
  REAL(pwp), ALLOCATABLE :: trmse_s(:)       ! true RMS error of smoothed states
  REAL(pwp), ALLOCATABLE :: mrmse_s_null(:)  ! Time-mean of estimated smoother RMS error
  REAL(pwp), ALLOCATABLE :: mtrmse_s_null(:) ! Time-mean of true smoother RMS error
  REAL(pwp), ALLOCATABLE :: mrmse_s_step(:)  ! Time-mean of estimated smootherRMS error
  REAL(pwp), ALLOCATABLE :: mtrmse_s_step(:) ! Time-mean of true smoother RMS error


  ! **********************
  ! *** INITIALIZATION ***
  ! **********************

  ! Variables to be initialised only on PE of rank 0 on mype_filter communicator
  mype0:IF (mype_filter == 0) THEN
     ! Create all netcdf output files
     firststep:IF (firsttime) THEN
        ! SSH
        IF (output_ssh == .TRUE.) THEN
           WRITE (*, '(8x, a)') 'Initialize netcdf file for SSH analysis field.'
           var = 'sossheig'
           CALL init_netcdf_asml(nit000, rdt, jpiglo, jpjglo, output_lev, trim(var), &
                '2D', filtertype, subtype, dim_ens, forget, local_range, &
                locweight, srange, rms_obs, delt_obs, nitend, dim_lag)
           CALL close_netcdf_asml()
        END IF

        ! T
        IF (output_T == .TRUE.) THEN
           WRITE (*, '(8x, a)') 'Initialize netcdf file for T analysis field.'
           var = 'votemper'
           CALL init_netcdf_asml(nit000, rdt, jpiglo, jpjglo, output_lev, trim(var), &
                '3D', filtertype, subtype, dim_ens, forget, local_range, &
                locweight, srange, rms_obs, delt_obs, nitend, dim_lag)
           CALL close_netcdf_asml()
        END IF

        ! S
        IF (output_S == .TRUE.) THEN
           WRITE (*, '(8x, a)') 'Initialize netcdf file for S analysis field.'
           var = 'vosaline'
           CALL init_netcdf_asml(nit000, rdt, jpiglo, jpjglo, output_lev, trim(var), &
                '3D', filtertype, subtype, dim_ens, forget, local_range, &
                locweight, srange, rms_obs, delt_obs, nitend, dim_lag)
           CALL close_netcdf_asml()
        END IF

        ! U
        IF (output_U == .TRUE.) THEN
           WRITE (*, '(8x, a)') 'Initialize netcdf file for U analysis field.'
           var = 'vozocrtx'
           CALL init_netcdf_asml(nit000, rdt, jpiglo, jpjglo, output_lev, trim(var), &
                '3D', filtertype, subtype, dim_ens, forget, local_range, &
                locweight, srange, rms_obs, delt_obs, nitend, dim_lag)
           CALL close_netcdf_asml()
        END IF

        ! V
        IF (output_V == .TRUE.) THEN
           WRITE (*, '(8x, a)') 'Initialize netcdf file for V analysis field.'
           var = 'vomecrty'
           CALL init_netcdf_asml(nit000, rdt, jpiglo, jpjglo, output_lev, trim(var), &
                '3D', filtertype, subtype, dim_ens, forget, local_range, &
                locweight, srange, rms_obs, delt_obs, nitend, dim_lag)
           CALL close_netcdf_asml()
        END IF
     ELSE firststep
        IF (step<0) THEN
           WRITE (*, '(8x, a)') 'XIOS writes forecast fields, prepoststep does nothing.'
        ELSE
           WRITE (*, '(8x, a)') 'Begin prepoststep write of ensemble of analysis fields.'
        END IF
     END IF firststep
  END IF mype0

  ! Values to be initialised on all PEs on mype_filter communicator
  IF (firsttime) THEN
     anastr = 'ini'
     WRITE (stepstr, '(i2.2)') step
     iter = 1
  ELSE
     IF (step<0) THEN
        anastr = 'for'
        WRITE (stepstr, '(i2.2)') -step
     ELSE
        anastr = 'ana'
        WRITE (stepstr, '(i2.2)') step
     END IF
  END IF

  ! *******************************
  ! *** Compute variance vector ***
  ! *******************************
!!$  ! Allocate fields
!!$  ALLOCATE(variance_p(dim_p))
!!$  ALLOCATE(variance(dim_state))
!!$
!!$  ! Initialize numbers
!!$  rmserror_ssh_est  = 0.0
!!$  invdim_ens    = 1.0 / REAL(dim_ens)  
!!$  invdim_ensm1  = 1.0 / REAL(dim_ens - 1)
!!$
!!$  ! *** Compute mean state
!!$  IF (mype_filter == 0) WRITE (*, '(8x, a)') '--- compute ensemble mean'
!!$
!!$  state_p = 0.0
!!$  DO member = 1, dim_ens
!!$     DO i = 1, dim_p
!!$        state_p(i) = state_p(i) + ens_p(i, member)
!!$     END DO
!!$  END DO
!!$  state_p(:) = invdim_ens * state_p(:)
!!$
!!$  ! *** Compute sampled variances ***
!!$  variance_p(:) = 0.0
!!$  DO member = 1, dim_ens
!!$     DO j = 1, dim_p
!!$        variance_p(j) = variance_p(j) &
!!$             + (ens_p(j, member) - state_p(j)) &
!!$             * (ens_p(j, member) - state_p(j))
!!$     END DO
!!$  END DO
!!$  variance_p(:) = invdim_ensm1 * variance_p(:)


  ! ******************************************************
  ! *** Assemble global variance vector on filter PE 0 ***
  ! ******************************************************

!!$  PE0_a: IF (mype_filter /= 0) THEN
!!$
!!$     ! send sub-fields from PEs /=0
!!$     CALL MPI_send(variance_p(1 : dim_p), dim_p, &
!!$          MPI_DOUBLE_PRECISION,0, mype_filter, COMM_filter, MPIerr)
!!$
!!$  ELSE PE0_a
!!$     ! receive and assemble variance field
!!$
!!$     ! On PE 0 init variance directly
!!$     variance(1 : dim_p) = variance_p(1 : dim_p)
!!$
!!$     ! Receive part of variance field from PEs > 0 into 
!!$     ! correct part of global variance
!!$
!!$     off_p = 0
!!$
!!$     DO i = 2, npes_filter
!!$        ! Increment offset
!!$        off_p = off_p + ssh_p_dim_state
!!$
!!$        ! Receive variance part
!!$        CALL MPI_recv(variance(1 + off_p), ssh_p_dim_state, &
!!$             MPI_DOUBLE_PRECISION, i - 1, i - 1, COMM_filter, MPIstatus, MPIerr)
!!$     END DO
!!$
!!$  END IF PE0_a
!!$
!!$  DEALLOCATE(variance_p)


  ! ************************************************************
  ! *** Compute RMS errors according to sampled covar matrix ***
  ! ************************************************************

!!$  ! total estimated RMS error
!!$  DO i = 1, ssh_p_dim_state
!!$     rmserror_ssh_est = rmserror_ssh_est + variance(i)
!!$  ENDDO
!!$  rmserror_ssh_est = SQRT(rmserror_ssh_est / ssh_p_dim_state)

  ! DUMMY STATISTICS FOR NOW
  ALLOCATE(hist_true(dim_ens+1,2))
  ALLOCATE(hist_mean(dim_ens+1,2))
  hist_true=0.0
  hist_mean=0.0

  ! *******************
  ! *** File output ***
  ! *******************

  ! Only write analysis fields to file.
  ana: IF (anastr == 'ana') THEN
     ! Global ensembles for file output.
     ! *MUST* be allocated before entry into gather_X blocks.
     ! 'ens_ssh' allocated a jpk value as ugly hack for
     ! reusing the same routines for all state variables.
     ALLOCATE(ens_ssh(jpiglo, jpjglo, jpk, dim_ens))
     ALLOCATE(ens_t(jpiglo, jpjglo, jpk, dim_ens))
     ALLOCATE(ens_s(jpiglo, jpjglo, jpk, dim_ens))
     ALLOCATE(ens_u(jpiglo, jpjglo, jpk, dim_ens))
     ALLOCATE(ens_v(jpiglo, jpjglo, jpk, dim_ens))

     ! Strategy is to gather global ensemble on a single PE with rank = k and
     ! then to write to file. Each state variable is written on a PE with a
     ! different rank. The value of each such rank is *hardcoded*.
     rank_ssh = 0
     rank_T = 1
     rank_S = 2
     rank_U = 3
     rank_V = 4
     !
     ! SSH
     gather_ssh: IF (output_ssh == .TRUE.) THEN
        ! Gather ssh subset of ens_p on PE with *RANK 0* and store in ens_ssh.
        ! Set mpi_subd_vert=1 as SSH is a 2D field (ugly hack).
        CALL gather_ens(rank_ssh, mpi_subd_lon, mpi_subd_lat, 1, jpiglo, &
             jpjglo, jpk, ssh_p_dim_state, ssh_p_offset, dim_ens, ens_p, ens_ssh)
     END IF gather_ssh

     ! T
     gather_T: IF (output_T == .TRUE.) THEN
        ! Gather T subset of ens_p on PE with *RANK 1* and store in ens_t.
        CALL gather_ens(rank_T, mpi_subd_lon, mpi_subd_lat, mpi_subd_vert, jpiglo, &
             jpjglo, jpk, t_p_dim_state, t_p_offset, dim_ens, ens_p, ens_t)
     END IF gather_T

     ! S
     gather_S: IF (output_S == .TRUE.) THEN
        ! Gather S subset of ens_p on PE with *RANK 2* and store in ens_s.
        CALL gather_ens(rank_S, mpi_subd_lon, mpi_subd_lat, mpi_subd_vert, jpiglo, &
             jpjglo, jpk, s_p_dim_state, s_p_offset, dim_ens, ens_p, ens_s)
     END IF gather_S

     ! U
     gather_U: IF (output_U == .TRUE.) THEN
        ! Gather U subset of ens_p on PE with *RANK 3* and store in ens_u.
        CALL gather_ens(rank_U, mpi_subd_lon, mpi_subd_lat, mpi_subd_vert, jpiglo, &
             jpjglo, jpk, u_p_dim_state, u_p_offset, dim_ens, ens_p, ens_u)
     END IF gather_U

     ! V
     gather_v: IF (output_V == .TRUE.) THEN
        ! Gather V subset of ens_p on PE with *RANK 4* and store in ens_v.
        CALL gather_ens(rank_V, mpi_subd_lon, mpi_subd_lat, mpi_subd_vert, jpiglo, &
             jpjglo, jpk, v_p_dim_state, v_p_offset, dim_ens, ens_p, ens_v)
     END IF gather_v

     ! SSH
     write_ssh: IF (output_ssh == .TRUE. .AND. mype_filter == rank_ssh) THEN
        var='sossheig'
        WRITE (*, '(8x, a, 3x, a, 1x, a, 1x, i3)') '--- write ensemble for', var, &
             'on PE with rank:', rank_ssh
        CALL write_netcdf_asml(anastr, jpiglo, jpjglo, jpk, iter, &
             output_lev, trim(var), '2D', rmse_est, rmse_true, mrmse_est_null, &
             mrmse_true_null, mrmse_est_step, mrmse_true_step, dim_ens, ens_ssh, &
             hist_true, hist_mean, skewness, kurtosis, dim_lag, rmse_s, trmse_s, &
             mrmse_s_null, mtrmse_s_null, mrmse_s_step, mtrmse_s_step)
        CALL close_netcdf_asml()
     END IF write_ssh

     ! T
     write_T: IF (output_T == .TRUE. .AND. mype_filter == rank_T) THEN
        var='votemper'
        WRITE (*, '(8x, a, 3x, a, 1x, a, 1x, i3)') '--- write ensemble for', var, &
             'on PE with rank:', rank_T
        CALL write_netcdf_asml(anastr, jpiglo, jpjglo, jpk, iter, &
             output_lev, trim(var), '3D', rmse_est, rmse_true, mrmse_est_null, &
             mrmse_true_null, mrmse_est_step, mrmse_true_step, dim_ens, ens_t, &
             hist_true, hist_mean, skewness, kurtosis, dim_lag, rmse_s, trmse_s, &
             mrmse_s_null, mtrmse_s_null, mrmse_s_step, mtrmse_s_step)
        CALL close_netcdf_asml()
     END IF write_T

     ! S
     write_S: IF (output_S == .TRUE. .AND. mype_filter == rank_S) THEN
        var='vosaline'
        WRITE (*, '(8x, a, 3x, a, 1x, a, 1x, i3)') '--- write ensemble for', var, &
             'on PE with rank:', rank_S
        CALL write_netcdf_asml(anastr, jpiglo, jpjglo, jpk, iter, &
             output_lev, trim(var), '3D', rmse_est, rmse_true, mrmse_est_null, &
             mrmse_true_null, mrmse_est_step, mrmse_true_step, dim_ens, ens_s, &
             hist_true, hist_mean, skewness, kurtosis, dim_lag, rmse_s, trmse_s, &
             mrmse_s_null, mtrmse_s_null, mrmse_s_step, mtrmse_s_step)
        CALL close_netcdf_asml()
     END IF write_S

     ! U
     write_U: IF (output_U == .TRUE. .AND. mype_filter == rank_U) THEN
        var='vozocrtx'
        WRITE (*, '(8x, a, 3x, a, 1x, a, 1x, i3)') '--- write ensemble for', var, &
             'on PE with rank:', rank_U
        CALL write_netcdf_asml(anastr, jpiglo, jpjglo, jpk, iter, &
             output_lev, trim(var), '3D', rmse_est, rmse_true, mrmse_est_null, &
             mrmse_true_null, mrmse_est_step, mrmse_true_step, dim_ens, ens_u, &
             hist_true, hist_mean, skewness, kurtosis, dim_lag, rmse_s, trmse_s, &
             mrmse_s_null, mtrmse_s_null, mrmse_s_step, mtrmse_s_step)
        CALL close_netcdf_asml()
     END IF write_U

     ! V
     write_V: IF (output_V == .TRUE. .AND. mype_filter == rank_V) THEN
        var='vomecrty'
        WRITE (*, '(8x, a, 3x, a, 1x, a, 1x, i3)') '--- write ensemble for', var, &
             'on PE with rank:', rank_V
        CALL write_netcdf_asml(anastr, jpiglo, jpjglo, jpk, iter, &
             output_lev, trim(var), '3D', rmse_est, rmse_true, mrmse_est_null, &
             mrmse_true_null, mrmse_est_step, mrmse_true_step, dim_ens, ens_v, &
             hist_true, hist_mean, skewness, kurtosis, dim_lag, rmse_s, trmse_s, &
             mrmse_s_null, mtrmse_s_null, mrmse_s_step, mtrmse_s_step)
        CALL close_netcdf_asml()
     END IF write_V

     ! Increment netcdf counter for next analysis time step
     iter = iter + 1

     ! Tidy up
     DEALLOCATE(ens_ssh,ens_t,ens_s,ens_u,ens_v)

  END IF ana

  ! ********************
  ! *** finishing up ***
  ! ********************

  firsttime = .FALSE.

  ! Tidy up statistics
  DEALLOCATE(hist_true, hist_mean)

END SUBROUTINE prepoststep_ens_pdaf
