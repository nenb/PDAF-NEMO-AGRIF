!$Id: obs_A_pdafomi.F90 496 2020-06-09 15:26:17Z lnerger $
!> PDAF-OMI observation module for type A observations
!!
!! This module handles operations for one data type (called 'module-type' below):
!! TYPE = A
!!
!! __Observation type A:__
!! The observation type A in this tutorial is the full set of observations except
!! for three observations at the locations (8,5), (12,15), and (4,30) which are
!! removed and only used for observation type B.
!!
!! The subroutines in this module are for the particular handling of
!! a single observation type.
!! The routines are called by the different call-back routines of PDAF.
!! Most of the routines are generic so that in practice only 2 routines
!! need to be adapted for a particular data type. These are the routines
!! for the initialization of the observation information (init_dim_obs_f)
!! and for the observation operator (obs_op_f).
!!
!! The module and the routines are named according to the observation type.
!! This allows to distinguish the observation type and the routines in this
!! module from other observation types.
!!
!! The module uses two derived data type (obs_f and obs_l), which contain
!! all information about the full and local observations. Only variables
!! of the type obs_f need to be initialized in this module. The variables
!! in the type obs_l are initilized by the generic routines from PDAFomi.
!!
!! These 2 routines need to be adapted for the particular observation type:
!! * init_dim_obs_f_TYPE \n
!!           Count number of process-local and full observations; 
!!           initialize vector of observations and their inverse variances;
!!           initialize coordinate array and index array for indices of
!!           observed elements of the state vector.
!! * obs_op_f_TYPE \n
!!           observation operator to get full observation vector of this type. Here
!!           one has to choose a proper observation operator or implement one.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_obs_fake_ssh_child_pdafomi
!$AGRIF_DO_NOT_TREAT

  USE mod_kind_pdaf
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, abort_parallel
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
  USE netcdf
  USE mod_agrif_pdaf, &
       ONLY: lowlim_ssh_child, upplim_ssh_child

  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_fake_ssh_child    !< Whether to assimilate this data type
  REAL(pwp) :: rms_fake_ssh_child    !< Observation error standard deviation (for constant errors)
  LOGICAL :: twin_exp_fake_ssh_child ! Whether to perform an identical twin experiment
  REAL(pwp) :: noise_amp_fake_ssh_child  ! Standard deviation for Gaussian noise in twin experiment

  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.
  CHARACTER(len=200) :: file_fake_ssh_child ! netcdf file holding observations

! ***********************************************************************
! *** The following two data types are used in PDAFomi                ***
! *** They are declared in PDAFomi and only listed here for reference ***
! ***********************************************************************

! Data type to define the full observations by internally shared variables of the module
!   TYPE obs_f
!           Mandatory variables to be set in init_dim_obs_f
!      INTEGER :: doassim                   ! Whether to assimilate this observation type
!      INTEGER :: disttype                  ! Type of distance computation to use for localization
!      INTEGER :: ncoord                    ! Number of coordinates use for distance computation
!      LOGICAL :: use_global_obs=.true.     ! Whether to use (T) global full obs. 
!                                           ! or (F) obs. restricted to those relevant for a process domain
!      INTEGER, ALLOCATABLE :: id_obs_p(:,:) ! indices of observed field in state vector
!           
!           Optional variables - they can be set in init_dim_obs_f
!      REAL, ALLOCATABLE :: icoeff_p(:,:)   ! Interpolation coefficients for obs. operator
!      REAL, ALLOCATABLE :: domainsize(:)   ! Size of domain for periodicity (<=0 for no periodicity)
!      INTEGER :: obs_err_type=0            ! Type of observation error: (0) Gauss, (1) Laplace
!
!           The following variables are set in the routine PDAFomi_gather_obs_f
!      INTEGER :: dim_obs_p                 ! number of PE-local observations
!      INTEGER :: dim_obs_f                 ! number of full observations
!      REAL, ALLOCATABLE :: obs_f(:)        ! Full observed field
!      REAL, ALLOCATABLE :: ocoord_f(:,:)   ! Coordinates of full observation vector
!      REAL, ALLOCATABLE :: ivar_obs_f(:)   ! Inverse variance of full observations
!      INTEGER :: dim_obs_g                 ! global number of observations 
!                                           ! (only if full obs. are restricted to process domain))
!      INTEGER, ALLOCATABLE :: id_obs_f_lim(:) ! Indices of domain-relevant full obs. in global vector of obs.
!                                           ! (only if full obs. are restricted to process domain))
!
!           Mandatory variable to be set in obs_op_f
!      INTEGER :: off_obs_f                 ! Offset of this observation in overall full obs. vector
!   END TYPE obs_f

! Data type to define the local observations by internally shared variables of the module
!   TYPE obs_l
!      INTEGER :: dim_obs_l                 ! number of local observations
!      INTEGER :: off_obs_l                 ! Offset of this observation in overall local obs. vector
!      INTEGER, ALLOCATABLE :: id_obs_l(:)  ! Indices of local observations in full obs. vector 
!      REAL, ALLOCATABLE :: distance_l(:)   ! Distances of local observations
!      REAL, ALLOCATABLE :: ivar_obs_l(:)   ! Inverse variance of local observations
!   END TYPE obs_l
! ***********************************************************************

! Declare instances of observation data types used here
! We use generic names here, but one could renamed the variables
  TYPE(obs_f), TARGET, PUBLIC :: thisobs      ! full observation
  TYPE(obs_l), TARGET, PUBLIC :: thisobs_l    ! local observation

!$OMP THREADPRIVATE(thisobs_l)


!-------------------------------------------------------------------------------

CONTAINS

!> Initialize information on the module-type observation
!!
!! The routine is called by each filter process.
!! at the beginning of the analysis step before 
!! the loop through all local analysis domains.
!! 
!! It has to count the number of observations of the
!! observation type handled in this module according
!! to the current time step for all observations 
!! required for the analyses in the loop over all local 
!! analysis domains on the PE-local state domain.
!!
!! The following four variables have to be initialized in this routine
!! * thisobs\%doassim     - Whether to assimilate this type of observations
!! * thisobs\%disttype    - type of distance computation for localization with this observaton
!! * thisobs\%ncoord      - number of coordinates used for distance computation
!! * thisobs\%id_obs_p    - index of module-type observation in PE-local state vector
!!
!! Optional is the use of
!! * thisobs\%icoeff_p    - Interpolation coefficients for obs. operator (only if interpolation is used)
!! * thisobs\%domainsize  - Size of domain for periodicity for disttype=1 (<0 for no periodicity)
!! * thisobs\%obs_err_type - Type of observation errors for particle filter and NETF
!! * thisobs\%use_global obs - Whether to use global observations or restrict the observations to the relevant ones
!!                          (default: .true.: use global full observations)
!!
!! The following variables are set in the routine gather_obs_f
!! * thisobs\%dim_obs_p   - PE-local number of module-type observations
!! * thisobs\%dim_obs_f   - full number of module-type observations
!! * thisobs\%obs_f       - full vector of module-type observations
!! * thisobs\%ocoord_f    - coordinates of observations in OBS_MOD_F
!! * thisobs\%ivar_obs_f  - full vector of inverse obs. error variances of module-type
!! * thisobs\%dim_obs_g   - Number of global observations (only if if use_global_obs=.false)
!! * thisobs\%id_obs_f_lim - Ids of full observations in global observations (if use_global_obs=.false)
!!
  SUBROUTINE init_dim_obs_f_fake_ssh_child(step, dim_obs_f)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs_f
    USE mod_assimilation_pdaf, &
         ONLY: filtertype, local_range_child, delt_obs
    USE mod_statevector_pdaf, &
         ONLY: mpi_subd_lon_child, mpi_subd_lat_child, ssh_p_offset_child
    USE mod_parallel_pdaf, &
         ONLY: nldj_child, nldi_child, COMM_filter, jpiglo_child, jpjglo_child, &
         nimpp_child, njmpp_child
    USE mod_agrif_pdaf, &
         ONLY: glamt_child, gphit_child, ndastp_child

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)    :: step       !< Current time step
    INTEGER, INTENT(inout) :: dim_obs_f  !< Dimension of full observation vector

! *** Local variables ***
    ! Local variables for NetCDF file
    INTEGER :: i, j, s           ! Counters
    INTEGER :: step_count = 0    ! Counter for number of times mod_obs_a is called
    INTEGER :: stat(20000)       ! Status flag for NetCDF commands
    INTEGER :: ncid_in           ! ID for NetCDF file
    INTEGER :: id_var            ! IDs for fields
    INTEGER :: pos(3),cnt(3)     ! Vectors for 3D field
    REAL(pwp), ALLOCATABLE :: obs(:,:,:)            ! Global observation field
    ! Local variables for counting number of PE-local observations
    INTEGER :: dim_obs_p                       ! Number of process-local observations
    INTEGER :: cnt_p, cnt0_p                   ! Counters
    INTEGER :: i_obs, j_obs                   ! Coordinates for observation and gridbox
    INTEGER :: i0, j0                          ! Halo offset for local PE
    REAL(pwp), ALLOCATABLE :: obs_p(:)         ! PE-local observation vector
    REAL(pwp), ALLOCATABLE :: ivar_obs_p(:)    ! PE-local inverse observation error variance
    REAL(pwp), ALLOCATABLE :: ocoord_p(:,:)    ! PE-local observation coordinates
    REAL(pwp), PARAMETER :: rad_conv = 3.141592653589793/180.0   ! Degree to radian conversion


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs_fake_ssh_child'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_fake_ssh_child) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 3   ! 3=Haversine

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! Format of ndastp_child is YYYYMMDD
    step_count=MOD(ndastp_child,100)
    IF(mype_filter == 0) WRITE(*,'(/9x, a, 2i)') &
         'mod_obs_fake_ssh_child current date:', ndastp_child

    ! *****************************************
    ! *** Open file containing observations ***
    ! *****************************************

    s = 1
    stat(s) = NF90_OPEN(file_fake_ssh_child, NF90_NOWRITE, ncid_in)
    s = s + 1

    DO j = 1, s-1
       IF (stat(s) /= NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'ERROR: NetCDF error in opening obs file:', file_fake_ssh_child
          CALL abort_parallel()
       END IF
    END DO

    ! *************************
    ! *** Read observations ***
    ! *************************

    ALLOCATE (obs(jpiglo_child, jpjglo_child, 1))

    s = 1
    stat(s) = NF90_INQ_VARID(ncid_in, 'sossheig', id_var)
    s = s + 1

    pos = (/ 1, 1, step_count /)
    cnt = (/ jpiglo_child, jpjglo_child, 1 /)

    stat(s) = NF90_GET_VAR(ncid_in, id_var, obs, start=pos, count=cnt)
    s = s + 1

    DO j = 1, s-1
       IF (stat(j) .NE. NF90_NOERR) THEN
          WRITE(*, '(/9x, a, 3x, i1, 3x, a, 3x, a)') &
               'NetCDF error in reading dimension:', j,&
               ' from obs file:', file_fake_ssh_child
          CALL abort_parallel()
       END IF
    END DO

    ! *******************************************
    ! *** Close file containing initial state ***
    ! *******************************************

    s = 1
    stat(s) = NF90_CLOSE(ncid_in)
    s = s + 1

    DO j = 1, s - 1
       IF (stat(j) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in closing obs file:', file_fake_ssh_child
          CALL abort_parallel()
       END IF
    END DO


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************


    ! *** Count valid observations that lie within the process sub-domain ***
    cnt_p = 0

    ! Compute halo offset on child grid.
    i0 = nldi_child - 1
    j0 = nldj_child - 1

    ! DO loop to count valid observations on PE-local
    DO j = 1, mpi_subd_lat_child
       DO i = 1, mpi_subd_lon_child
          ! Convert to global coordinates
          i_obs = nimpp_child + i0 + i - 1
          j_obs = njmpp_child + j0 + j - 1
          ! Reject non-existent and unrealistic obs
          IF ( obs(i_obs, j_obs, 1) .GT. lowlim_ssh_child) THEN
             IF ( obs(i_obs, j_obs, 1) .LT. upplim_ssh_child) THEN
                cnt_p=cnt_p+1
             END IF
          END IF
       END DO
    END DO

    ! Set number of local observations
    dim_obs_p = cnt_p

    IF(cnt_p == 0) THEN
       WRITE(*,'(/9x, a, i3, 3x, a, i3)') 'WARNING: No fake_ssh_child observations on &
            & PE:', mype_filter, 'step_count:', step_count
    END IF

    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***
 
    obs_nonzero:IF(dim_obs_p>0) THEN
       ! Allocate process-local observation arrays
       ALLOCATE(obs_p(dim_obs_p))
       ALLOCATE(ivar_obs_p(dim_obs_p))
       ALLOCATE(ocoord_p(2, dim_obs_p))
       ! Allocate process-local index array
       ! This array has a many rows as required for the observation operator
       ! 1 if observations are at grid points; >1 if interpolation is required
       ALLOCATE(thisobs%id_obs_p(1, dim_obs_p))

       ! Initialise counters.
       cnt_p = 0
       cnt0_p = 0

       ! Compute halo offset on child grid.
       i0 = nldi_child - 1
       j0 = nldj_child - 1

       DO j = 1, mpi_subd_lat_child
          DO i = 1, mpi_subd_lon_child
             ! State vector index counter for observation operator.
             cnt0_p = cnt0_p + 1
             ! Convert to global coordinates.
             i_obs = nimpp_child + i0 + i - 1
             j_obs = njmpp_child + j0 + j - 1
             ! Reject non-existent and unrealistic obs.
             IF ( obs(i_obs, j_obs, 1) .GT. lowlim_ssh_child) THEN
                IF ( obs(i_obs, j_obs, 1) .LT. upplim_ssh_child) THEN
                   ! Total number of valid observations on local PE
                   cnt_p = cnt_p + 1
                   ! Observation
                   obs_p(cnt_p) = obs(i_obs, j_obs, 1)
                   ! Observation coordinates - must be in radians
                   ocoord_p(1, cnt_p) = glamt_child(i+i0,j+j0)*rad_conv
                   ocoord_p(2, cnt_p) = gphit_child(i+i0,j+j0)*rad_conv
                   ! Coordinates for observation operator (gridpoint)
                   thisobs%id_obs_p(1, cnt_p) = cnt0_p + ssh_p_offset_child
                END IF
             END IF
          END DO
       END DO
    ELSE
       ! No observations on PE, create dummy arrays to pass to PDAF-OMI
       ALLOCATE(obs_p(1))
       ALLOCATE(ivar_obs_p(1))
       ALLOCATE(ocoord_p(2, 1))
       ALLOCATE(thisobs%id_obs_p(1, 1))
       obs_p=-999999.0
       ivar_obs_p=-999999.0
       ocoord_p=-999999.0
       thisobs%id_obs_p=-999999.0
    END IF obs_nonzero


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ! *** Set inverse observation error variances ***

    ivar_obs_p(:) = 1.0 / (rms_fake_ssh_child*rms_fake_ssh_child)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs_f(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, local_range_child, dim_obs_f)


! *********************************************************
! *** For twin experiment: Read synthetic observations  ***
! *********************************************************

!     IF (twin_experiment .AND. filtertype/=11) THEN
!        CALL read_syn_obs(file_syntobs_TYPE, dim_obs_f, thisobs%obs_f, 0, 1-mype_filter)
!     END IF


! ********************
! *** Finishing up ***
! ********************

    ! Deallocate all local arrays
    DEALLOCATE(obs)
    DEALLOCATE(obs_p, ocoord_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_f_fake_ssh_child



!-------------------------------------------------------------------------------
!> Implementation of observation operator 
!!
!! This routine applies the full observation operator
!! for the type of observations handled in this module
!! It has to append the observations to ostate_f from
!! position OFFSET_OBS+1. For the return value OFFSET_OBS
!! has to be incremented by the number of added observations.
!!
!! One can choose a proper observation operator from
!! PDAFOMI_OBS_OP or add one to that module or 
!! implement another observation operator here.
!!
!! The order of the calls to this routine for different modules
!! is important because it influences the offset of the 
!! module-type observation in the overall full observation vector.
!!
!! Outputs for within the module are:
!! * thisobs\%off_obs_f - Offset of full module-type observation in overall full obs. vector
!!
!! The routine is called by all filter processes.
!!
  SUBROUTINE obs_op_f_fake_ssh_child(dim_p, dim_obs_f, state_p, ostate_f, offset_obs)

    USE PDAFomi, &
         ONLY: PDAFomi_obs_op_f_gridpoint

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in) :: dim_p                 !< PE-local state dimension
    INTEGER, INTENT(in) :: dim_obs_f             !< Dimension of full observed state (all observed fields)
    REAL(pwp), INTENT(in)    :: state_p(dim_p)        !< PE-local model state
    REAL(pwp), INTENT(inout) :: ostate_f(dim_obs_f)   !< Full observed state
    INTEGER, INTENT(inout) :: offset_obs         !< input: offset of module-type observations in ostate_f
                                                 !< output: input + number of added observations


! ******************************************************
! *** Apply observation operator H on a state vector ***
! ******************************************************

    IF (thisobs%doassim==1) THEN
       ! observation operator for observed grid point values
       CALL PDAFomi_obs_op_f_gridpoint(thisobs, state_p, ostate_f, offset_obs)
    END IF

  END SUBROUTINE obs_op_f_fake_ssh_child



!-------------------------------------------------------------------------------
!> Initialize local information on the module-type observation
!!
!! The routine is called during the loop over all local
!! analysis domains. It has to initialize the information
!! about local observations of the module type. It returns
!! number of local observations of the module type for the
!! current local analysis domain in DIM_OBS_L and the full
!! and local offsets of the observation in the overall
!! observation vector.
!!
!! This routine calls the routine PDAFomi_init_dim_obs_l
!! for each observation type. The call allows to specify a
!! different localization radius and localization functions
!! for each observation type and  local analysis domain.
!!
  SUBROUTINE init_dim_obs_l_fake_ssh_child(domain_p, step, dim_obs_f, dim_obs_l, &
       off_obs_l, off_obs_f)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l, PDAFomi_set_debug_flag

    ! Include localization radius and local coordinates
    USE mod_assimilation_pdaf, &   
         ONLY: coords_l, local_range_child, locweight, srange_child

    IMPLICIT NONE

! *** Arguments ***
    INTEGER, INTENT(in)  :: domain_p     !< Index of current local analysis domain
    INTEGER, INTENT(in)  :: step         !< Current time step
    INTEGER, INTENT(in)  :: dim_obs_f    !< Full dimension of observation vector
    INTEGER, INTENT(out) :: dim_obs_l    !< Local dimension of observation vector
    INTEGER, INTENT(inout) :: off_obs_l  !< Offset in local observation vector
    INTEGER, INTENT(inout) :: off_obs_f  !< Offset in full observation vector


! **********************************************
! *** Initialize local observation dimension ***
! **********************************************

    ! Gridpoint debugging
!!$    IF (domain_p == 1019 .AND. mype_filter==20) THEN
!!$       CALL PDAFomi_set_debug_flag(1019)
!!$    ELSE IF (domain_p == 959 .AND. mype_filter==9) THEN
!!$       CALL PDAFomi_set_debug_flag(959)
!!$    ELSE IF (domain_p == 899 .AND. mype_filter==9) THEN
!!$       CALL PDAFomi_set_debug_flag(899)
!!$    ELSE
!!$       CALL PDAFomi_set_debug_flag(0)
!!$    ENDIF

    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, local_range_child, srange_child, &
         dim_obs_l, off_obs_l, off_obs_f)

  END SUBROUTINE init_dim_obs_l_fake_ssh_child

!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_obs_fake_ssh_child_pdafomi
