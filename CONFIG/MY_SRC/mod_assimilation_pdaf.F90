!$Id: mod_assimilation_pdaf.F90 1866 2017-12-21 09:05:27Z lnerger $
!BOP
!
! !MODULE:
MODULE mod_assimilation_pdaf
!$AGRIF_DO_NOT_TREAT
  ! !DESCRIPTION:
  ! This module provides variables needed for the 
  ! assimilation within the routines of the dummy model.
  ! For simplicity, all assimilation-related variables
  ! are stored here, even if they are only used in
  ! the main program for the filter initialization.
  ! Most variables can be specified as a command line 
  ! argument.
  !
  ! Implementation for the 2D online example
  ! with parallelization.
  !
  ! !REVISION HISTORY:
  ! 2013-09 - Lars Nerger - Initial code
  ! Later revisions - see svn log
  !
  ! !USES:
  USE mod_kind_pdaf
  IMPLICIT NONE
  SAVE
  !EOP

  ! *** Model- and data specific variables ***

  REAL(pwp), POINTER :: state_p_pointer(:,:) ! Pointer to PDAF state_p array
  INTEGER :: status_pointer                  ! PDAF state_p pointer status flag

  INTEGER :: dim_state           ! Global model state dimension
  INTEGER :: dim_state_p         ! Model state dimension for PE-local domain
  INTEGER :: dim_state_p_par     ! dim_state_p on parent grid
  INTEGER :: dim_state_p_child   ! dim_state_p on child grid

  INTEGER :: dim_obs_p                    ! Process-local number of observations
  REAL(pwp), ALLOCATABLE    :: obs_p(:)        ! Vector holding process-local observations
  INTEGER, ALLOCATABLE :: obs_index_p(:)  ! Vector holding state-vector indices of observations
  REAL(pwp), ALLOCATABLE    :: obs_f(:)        ! Vector holding full vector of observations
  REAL(pwp), ALLOCATABLE :: coords_obs_f(:,:)  ! Array for full observation coordinates
  INTEGER, ALLOCATABLE :: obs_index_l(:)  ! Vector holding local state-vector indices of observations
  REAL(pwp), ALLOCATABLE    :: distance_l(:)   ! Vector holding distances of local observations
  REAL(pwp), ALLOCATABLE    :: wght(:)         ! Vector holding weights for ensemble initialization


  ! *** Below are the generic variables used for configuring PDAF ***
  ! *** Their values are set in init_PDAF                         ***

  ! !PUBLIC MEMBER FUNCTIONS:
  ! ! Settings for time stepping - available as command line options
  LOGICAL :: model_error   ! Control application of model error
  REAL(pwp)    :: model_err_amp ! Amplitude for model error

  ! ! Settings for observations - available as command line options
  INTEGER   :: delt_obs       ! time step interval between assimilation steps
  INTEGER   :: child_dt_fac   ! time step factor between child and parent grid
  REAL(pwp) :: rms_obs        ! RMS error size for observation generation
  INTEGER   :: dim_obs        ! Number of observations

  ! ! General control of PDAF - available as command line options
  INTEGER :: screen       ! Control verbosity of PDAF
  ! (0) no outputs, (1) progess info, (2) add timings
  ! (3) debugging output
  INTEGER :: dim_ens      ! Size of ensemble for SEIK/LSEIK/EnKF/ETKF
  ! Number of EOFs to be used for SEEK
  INTEGER :: filtertype   ! Select filter algorithm:
  !   SEEK (0), SEIK (1), EnKF (2), LSEIK (3), ETKF (4)
  !   LETKF (5), ESTKF (6), LESTKF (7), NETF (9), LNETF (10)
  INTEGER :: subtype      ! Subtype of filter algorithm
  !   SEEK: 
  !     (0) evolve normalized modes
  !     (1) evolve scaled modes with unit U
  !     (2) fixed basis (V); variable U matrix
  !     (3) fixed covar matrix (V,U kept static)
  !   SEIK:
  !     (0) ensemble forecast; new formulation
  !     (1) ensemble forecast; old formulation
  !     (2) fixed error space basis
  !     (3) fixed state covariance matrix
  !     (4) SEIK with ensemble transformation
  !   EnKF:
  !     (0) analysis for large observation dimension
  !     (1) analysis for small observation dimension
  !   LSEIK:
  !     (0) ensemble forecast;
  !     (2) fixed error space basis
  !     (3) fixed state covariance matrix
  !     (4) LSEIK with ensemble transformation
  !   ETKF:
  !     (0) ETKF using T-matrix like SEIK
  !     (1) ETKF following Hunt et al. (2007)
  !       There are no fixed basis/covariance cases, as
  !       these are equivalent to SEIK subtypes 2/3
  !   LETKF:
  !     (0) LETKF using T-matrix like SEIK
  !     (1) LETKF following Hunt et al. (2007)
  !       There are no fixed basis/covariance cases, as
  !       these are equivalent to LSEIK subtypes 2/3
  !   ESTKF:
  !     (0) standard ESTKF 
  !       There are no fixed basis/covariance cases, as
  !       these are equivalent to SEIK subtypes 2/3
  !   LESTKF:
  !     (0) standard LESTKF 
  !       There are no fixed basis/covariance cases, as
  !       these are equivalent to LSEIK subtypes 2/3
  !   NETF:
  !     (0) standard NETF 
  !   LNETF:
  !     (0) standard LNETF 
  INTEGER :: incremental  ! Perform incremental updating in LSEIK
  INTEGER :: dim_lag      ! Number of time instances for smoother

  ! ! Filter settings - available as command line options
  !    ! General
  INTEGER :: type_forget  ! Type of forgetting factor
  REAL(pwp)    :: forget       ! Forgetting factor for filter analysis
  INTEGER :: dim_bias     ! dimension of bias vector
  !    ! SEEK
  INTEGER :: int_rediag   ! Interval to perform re-diagonalization in SEEK
  REAL(pwp)    :: epsilon      ! Epsilon for gradient approx. in SEEK forecast
  !    ! ENKF
  INTEGER :: rank_analysis_enkf  ! Rank to be considered for inversion of HPH
  !    ! SEIK/ETKF/ESTKF/LSEIK/LETKF/LESTKF
  INTEGER :: type_trans    ! Type of ensemble transformation
  ! SEIK/LSEIK:
  ! (0) use deterministic omega
  ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
  ! (2) use product of (0) with random orthonormal matrix with
  !     eigenvector (1,...,1)^T
  ! ETKF/LETKF with subtype=4:
  ! (0) use deterministic symmetric transformation
  ! (2) use product of (0) with random orthonormal matrix with
  !     eigenvector (1,...,1)^T
  ! ESTKF/LESTKF:
  ! (0) use deterministic omega
  ! (1) use random orthonormal omega orthogonal to (1,...,1)^T
  ! (2) use product of (0) with random orthonormal matrix with
  !     eigenvector (1,...,1)^T
  ! NETF/LNETF:
  ! (0) use random orthonormal transformation orthogonal to (1,...,1)^T
  ! (1) use identity transformation
  !    ! LSEIK/LETKF/LESTKF
  REAL(pwp)    :: local_range   ! Range for local observation domain
  INTEGER :: locweight     ! Type of localizing weighting of observations
  !   (0) constant weight of 1
  !   (1) exponentially decreasing with SRANGE
  !   (2) use 5th-order polynomial
  !   (3) regulated localization of R with mean error variance
  !   (4) regulated localization of R with single-point error variance
  REAL(pwp)    :: srange        ! Support range for 5th order polynomial
  !   or radius for 1/e for exponential weighting
  !    ! SEIK-subtype4/LSEIK-subtype4/ESTKF/LESTKF
  INTEGER :: type_sqrt     ! Type of the transform matrix square-root 
  !   (0) symmetric square root, (1) Cholesky decomposition

  !    ! File input - available as namelist option
  ! Parent grid
  CHARACTER (len=110) :: istate_fname_t  ! file for t initial state estimate
  CHARACTER (len=110) :: istate_fname_u  ! file for u initial state estimate
  CHARACTER (len=110) :: istate_fname_v  ! file for v initial state estimate
  ! Child grid
  CHARACTER (len=110) :: istate_fname_t_child  ! file for t initial state estimate
  CHARACTER (len=110) :: istate_fname_u_child  ! file for u initial state estimate
  CHARACTER (len=110) :: istate_fname_v_child  ! file for v initial state estimate
  !    ! Other variables - _NOT_ available as command line options!
  INTEGER :: covartype     ! For SEIK: Definition of ensemble covar matrix
  ! (0): Factor (r+1)^-1 (or N^-1)
  ! (1): Factor r^-1 (or (N-1)^-1) - real ensemble covar.
  ! This setting is only for the model part; The definition
  ! of P has also to be specified in PDAF_filter_init.
  ! Only for upward-compatibility of PDAF!
  REAL(pwp)    :: time          ! model time

  !    ! NEMO-AGRIF specific variables
  LOGICAL :: euler_flag = .FALSE.  ! Flag for using euler timestep in NEMO after assimilation

  !    ! Time counter for netcdf output
  INTEGER :: iter

  !    ! Control of output
  LOGICAL :: output_ssh = .TRUE.
  LOGICAL :: output_S   = .TRUE.
  LOGICAL :: output_T   = .TRUE.
  LOGICAL :: output_U   = .TRUE.
  LOGICAL :: output_V   = .TRUE.

CONTAINS

  !$Id: assimilate_pdaf.F90 1411 2013-09-25 14:04:41Z lnerger $
  !BOP
  !
  ! !ROUTINE: assimilate_pdaf - Routine to control perform analysis step
  !
  ! !INTERFACE:
  SUBROUTINE assimilate_pdaf()

    ! !DESCRIPTION:
    ! This routine is called during the model integrations at each time 
    ! step. It check whether the forecast phase is completed. If so, 
    ! PDAF_put_state_X is called to perform the analysis step.
    !
    ! !REVISION HISTORY:
    ! 2013-08 - Lars Nerger - Initial code for NEMO
    ! Later revisions - see svn log
    !
    ! !USES:
    USE mod_kind_pdaf
    USE mod_parallel_pdaf, &     ! Parallelization variables
         ONLY: mype_world, abort_parallel
    USE mod_agrif_pdaf, &
         ONLY: fill2d_statevector, fill3d_statevector

    IMPLICIT NONE

    ! !CALLING SEQUENCE:
    ! Called by: step
    ! CAlls: PDAF_assimilate_X
    !EOP

    ! Local variables
    INTEGER :: status_pdaf          ! PDAF status flag
    INTEGER, SAVE :: fill_cnt = 1   ! Counter to determine whether to fill statevector

    ! External subroutines
    EXTERNAL :: collect_state_pdaf, & ! Routine to collect a state vector from model fields
         init_dim_obs_pdaf, &         ! Initialize Dimension Of Observation Vector
         obs_op_pdaf, &               ! Implementation of the Observation operator
         init_obs_pdaf, &             ! Routine to provide vector of measurements
         prepoststep_ens_pdaf, &      ! User supplied pre/poststep routine
         prodRinvA_pdaf, &            ! Provide product R^-1 A for some matrix A
         init_obsvar_pdaf, &          ! Initialize mean observation error variance
         next_observation_pdaf, &     ! Provide time step, model time, &
                                ! and dimension of next observation
         distribute_state_pdaf        ! Routine to distribute a state vector to model fields
    EXTERNAL :: init_n_domains_pdaf, &   ! Provide number of local analysis domains
         init_dim_l_pdaf, &              ! Initialize state dimension for local ana. domain
         init_dim_obs_l_pdaf,&           ! Initialize dim. of obs. vector for local ana. domain
         g2l_state_pdaf, &               ! Get state on local ana. domain from global state
         l2g_state_pdaf, &               ! Init global state from state on local analysis domain
         g2l_obs_pdaf, &                 ! Restrict a global obs. vector to local analysis domain
         init_obs_l_pdaf, &              ! Provide vector of measurements for local ana. domain
         prodRinvA_l_pdaf, &             ! Provide product R^-1 A for some local matrix A
         init_obsvar_l_pdaf, &           ! Initialize local mean observation error variance
         init_obs_f_pdaf, &              ! Provide full vector of measurements for PE-local domain
         obs_op_f_pdaf, &                ! Obs. operator for full obs. vector for PE-local domain
         init_dim_obs_f_pdaf             ! Get dimension of full obs. vector for PE-local domain


    ! **************************************************************************
    ! ************************** Collect state variables ***********************
    ! **************************************************************************
    !
    ! Solution adopted is to collect variables in two stages: i) first collect
    ! on the child grid and ii) second collect on the parent grid. If no child
    ! grid exists then only the second step is done. Because collection is done
    ! here, the collect_state call-back routine from PDAF is now redundant.
    !
    ! **************************************************************************
    ! **************************************************************************

#if defined key_agrif
    IF ( fill_cnt == ((child_dt_fac+1)*delt_obs) - 1 ) THEN
       CALL fill2d_statevector(state_p_pointer, 'child')
       CALL fill3d_statevector(state_p_pointer, 'child')
    ENDIF
#endif
    IF ( fill_cnt == (child_dt_fac+1)*delt_obs ) THEN
       CALL fill2d_statevector(state_p_pointer, 'par')
       CALL fill3d_statevector(state_p_pointer, 'par')
       ! Reset counter for next assimilation step
       fill_cnt = 0
    END IF

    ! Increment counter
    fill_cnt = fill_cnt+1

    ! *********************************
    ! *** Call assimilation routine ***
    ! *********************************

    IF (filtertype == 4) THEN
       CALL PDAF_assimilate_etkf(collect_state_pdaf,&
            distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,&
            init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf,&
            init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
    ELSE IF (filtertype == 5) THEN
       CALL PDAF_assimilate_letkf(collect_state_pdaf,&
            distribute_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf,&
            init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf,&
            prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,&
            init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf,&
            g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf,&
            next_observation_pdaf, status_pdaf)
    ELSE IF (filtertype == 6) THEN
       CALL PDAF_assimilate_estkf(collect_state_pdaf,&
            distribute_state_pdaf, init_dim_obs_pdaf, obs_op_pdaf,&
            init_obs_pdaf, prepoststep_ens_pdaf, prodRinvA_pdaf,&
            init_obsvar_pdaf, next_observation_pdaf, status_pdaf)
    ELSEIF (filtertype == 7) THEN
       CALL PDAF_assimilate_lestkf(collect_state_pdaf,&
            distribute_state_pdaf, init_dim_obs_f_pdaf, obs_op_f_pdaf,&
            init_obs_f_pdaf, init_obs_l_pdaf, prepoststep_ens_pdaf,&
            prodRinvA_l_pdaf, init_n_domains_pdaf, init_dim_l_pdaf,&
            init_dim_obs_l_pdaf, g2l_state_pdaf, l2g_state_pdaf,&
            g2l_obs_pdaf, init_obsvar_pdaf, init_obsvar_l_pdaf,&
            next_observation_pdaf, status_pdaf)
    END IF


    ! Check for errors during execution of PDAF
    IF (status_pdaf /= 0) THEN
       WRITE (*,'(/1x,a6,i3,a43,i4,a1/)') &
            'ERROR ', status_pdaf, &
            ' in PDAF_put_state - stopping! (PE ', mype_world,')'
       CALL  abort_parallel()
    END IF

  END SUBROUTINE assimilate_pdaf
!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_assimilation_pdaf
