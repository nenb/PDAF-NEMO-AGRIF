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
!!           observation operator to get full observation vector of this type.
!!           one has to choose a proper observation operator or implement one.
!!
!! __Revision history:__
!! * 2019-06 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
MODULE mod_obs_ssh_NEMO_pdafomi
!$AGRIF_DO_NOT_TREAT

  USE mod_kind_pdaf
  USE mod_parallel_pdaf, &
       ONLY: mype_filter, abort_parallel
  USE PDAFomi, &
       ONLY: obs_f, obs_l   ! Declaration of observation data types
  USE netcdf
  USE mod_agrif_pdaf, &
       ONLY: lowlim_ssh_NEMO, upplim_ssh_NEMO
 
  IMPLICIT NONE
  SAVE

  ! Variables which are inputs to the module (usually set in init_pdaf)
  LOGICAL :: assim_ssh_NEMO    !< Whether to assimilate this data type
  REAL(pwp) :: rms_ssh_NEMO    !< Observation error standard deviation (for constant errors)
  LOGICAL :: twin_exp_ssh_NEMO ! Whether to perform an identical twin experiment
  REAL(pwp) :: noise_amp_ssh_NEMO  ! Standard deviation for Gaussian noise in twin experiment


  ! One can declare further variables, e.g. for file names which can
  ! be use-included in init_pdaf() and initialized there.
  CHARACTER(len=200) :: file_ssh_NEMO ! netcdf file holding observations
  REAL(pwp), ALLOCATABLE :: coord2grid(:,:)  ! Array to convert PE-local observation
                                             ! coordinates to PE-local gridbox coordinates.
  REAL(pwp), PARAMETER :: rad_conv = 3.141592653589793/180.0   ! Degree to radian conversion

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
  SUBROUTINE init_dim_obs_f_ssh_NEMO(step, dim_obs_f)

    USE PDAFomi, &
         ONLY: PDAFomi_gather_obs_f
    USE mod_assimilation_pdaf, &
         ONLY: filtertype, local_range, delt_obs
    USE mod_statevector_pdaf, &
         ONLY: mpi_subd_lon_par, mpi_subd_lat_par, ssh_p_offset_par
    USE mod_parallel_pdaf, &
         ONLY: nldj_child, nldj_par, nldi_child, nldi_par, COMM_filter
    USE mod_agrif_pdaf, &
         ONLY: glamt_par, gphit_par
    USE dom_oce, &
         ONLY: ndastp ! Current date

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
    INTEGER :: id_dimx,id_dimy   ! IDs for dimensions
    INTEGER :: dim_lon, dim_lat  ! NetCDF dimension values
    INTEGER :: id_var            ! IDs for fields
    INTEGER :: pos(3),cnt(3)     ! Vectors for 3D field
    REAL(pwp), ALLOCATABLE :: longitude_ssh_NEMO(:) ! Global observation longitude values
    REAL(pwp), ALLOCATABLE :: latitude_ssh_NEMO(:)  ! Global observation latitude values
    REAL(pwp), ALLOCATABLE :: obs(:,:,:)            ! Global observation field
    REAL(pwp) :: scale_fac                          ! Scale factor for observation field
    ! Local variables for counting number of PE-local observations
    INTEGER :: dim_obs_p                       ! Number of process-local observations
    INTEGER :: cnt_p, ind_obs                  ! Counters
    INTEGER :: i_obs, j_obs, i_pe, j_pe        ! Coordinates for observation and gridbox
    INTEGER :: i0, j0                          ! Halo offset for local PE
    REAL(pwp), ALLOCATABLE :: obs_p(:)         ! PE-local observation vector
    REAL(pwp), ALLOCATABLE :: ivar_obs_p(:)    ! PE-local inverse observation error variance
    REAL(pwp), ALLOCATABLE :: ocoord_p(:,:)    ! PE-local observation coordinates
    REAL(pwp), ALLOCATABLE :: ogrid_p(:,:)     ! PE-local observation gridbox coordinates
    LOGICAL :: firsttime = .TRUE.  ! Flag for routines that need only be called once
    ! Local variables for observation operator
    REAL(pwp) :: lon_pe(4), lat_pe(4)      ! Longitude/latitude values of gridbox
    REAL(pwp) :: lon_obs, lat_obs          ! Longitude/latitude values of observation
    REAL(pwp) :: dist_obs_a, dist_obs_b    ! Distances for great-circle weighted interpolation
    REAL(pwp) :: dist_obs_c, dist_obs_d    ! Distances for great-circle weighted interpolation
    REAL(pwp) :: w_a, w_b, w_c, w_d, norm  ! Weights for great-circle weighted interpolation


! *********************************************
! *** Initialize full observation dimension ***
! *********************************************

    IF (mype_filter==0) &
         WRITE (*,'(8x,a)') 'Assimilate observations - obs ssh_NEMO'

    ! Store whether to assimilate this observation type (used in routines below)
    IF (assim_ssh_NEMO) thisobs%doassim = 1

    ! Specify type of distance computation
    thisobs%disttype = 3   ! 3=Haversine

    ! Number of coordinates used for distance computation
    ! The distance compution starts from the first row
    thisobs%ncoord = 2


! **********************************
! *** Read PE-local observations ***
! **********************************

    ! *************************************************************

    ! Logic here is to initially just work with a single file name
    ! which will contain a sequence of times.

    ! Eventually will need to create some sort of archive system
    ! for the different obs files, and introduce a counter to
    ! know which obs file to open. Useful error handling is to
    ! include condition to test whether there are 0 obs at this
    ! timestep.

    ! Once obs file is open, initial logic is to read in a single
    ! time. If necessary, can start introducing composite values at
    ! a later time step.

    ! *************************************************************

    ! Format of ndastp is YYYYMMDD
    step_count=MOD(ndastp,100)
    IF(mype_filter == 0) WRITE(*,'(/9x, a, 2i)') &
         'mod_obs_ssh_NEMO current date:', ndastp

    ! *****************************************
    ! *** Open file containing observations ***
    ! *****************************************

    s = 1
    stat(s) = NF90_OPEN(file_ssh_NEMO, NF90_NOWRITE, ncid_in)
    s = s + 1

    DO j = 1, s-1
       IF (stat(s) /= NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'ERROR: NetCDF error in opening obs file:', file_ssh_NEMO
          CALL abort_parallel()
       END IF
    END DO

    ! ***********************
    ! *** Read dimensions ***
    ! ***********************

    s = 1
    stat(s) = NF90_INQ_DIMID(ncid_in, 'longitude', id_dimx)
    s = s + 1
    stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimx, len=dim_lon)
    s = s + 1
    stat(s) = NF90_INQ_DIMID(ncid_in, 'latitude', id_dimy)
    s = s + 1
    stat(s) = NF90_INQUIRE_DIMENSION(ncid_in, id_dimy, len=dim_lat)
    s = s + 1

    ALLOCATE(longitude_ssh_NEMO(dim_lon))
    stat(s) = NF90_INQ_VARID(ncid_in, 'longitude', id_dimx)
    s = s + 1
    stat(s) = NF90_GET_VAR(ncid_in, id_dimx, longitude_ssh_NEMO)
    s = s + 1

    ! Work in longitude format between 0 and 360 degrees
    DO i = 1 , dim_lon
       IF(longitude_ssh_NEMO(i) < 0.0_pwp) longitude_ssh_NEMO(i) = &
            longitude_ssh_NEMO(i) + 360.0_pwp
    END DO

    IF ( (MAXVAL(longitude_ssh_NEMO) >= 360.0_pwp) .OR. &
         (MINVAL(longitude_ssh_NEMO) < 0.0_pwp) ) THEN
       WRITE(*,*) 'ERROR: ssh_NEMO longitude has non-standard format.'
       CALL abort_parallel()
    END IF

    ALLOCATE(latitude_ssh_NEMO(dim_lat))
    stat(s) = NF90_INQ_VARID(ncid_in, 'latitude', id_dimy)
    s = s + 1
    stat(s) = NF90_GET_VAR(ncid_in, id_dimy, latitude_ssh_NEMO)
    s = s + 1

    DO j = 1, s-1
       IF (stat(j) .NE. NF90_NOERR) THEN
          WRITE(*, '(/9x, a, 3x, i3, 3x, a, 3x, i1, 3x, a, 3x, a)') &
               'ERROR: NetCDF error', stat(j), 'in reading dimension:', j, &
               ' from obs file:', file_ssh_NEMO
          CALL abort_parallel()
       END IF
    END DO

    ! *************************
    ! *** Read observations ***
    ! *************************

    ALLOCATE (obs(dim_lon, dim_lat, 1))

    s = 1
    stat(s) = NF90_INQ_VARID(ncid_in, 'sla', id_var)
    s = s + 1

    pos = (/ 1, 1, step_count /)
    cnt = (/ dim_lon, dim_lat, 1 /)

    stat(s) = NF90_GET_VAR(ncid_in, id_var, obs, start=pos, count=cnt)
    s = s + 1
    stat(s) = NF90_GET_ATT(ncid_in, id_var, 'scale_factor', scale_fac)
    s = s + 1

    DO j = 1 , dim_lat
       DO i = 1 , dim_lon
             obs(i,j,1) = obs(i,j,1) * scale_fac
       END DO
    END DO

    DO j = 1, s-1
       IF (stat(j) .NE. NF90_NOERR) THEN
          WRITE(*, '(/9x, a, 3x, i1, 3x, a, 3x, a)') &
               'NetCDF error in reading dimension:', j,&
               ' from obs file:', file_ssh_NEMO
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
               'NetCDF error in closing obs file:', file_ssh_NEMO
          CALL abort_parallel()
       END IF
    END DO


! *********************************
! *** Determine obs on local PE ***
! *********************************

    ! NEMO has a curvilinear grid. This makes determining if obs on
    ! local PE non-trivial. See notes in following routine for details.
    ! *WARNING*: Longitude needs format 0-360 degrees for this routine.
    IF (firsttime) THEN
       CALL obs_in_grid_lcl(longitude_ssh_NEMO, latitude_ssh_NEMO, &
            mpi_subd_lon_par, mpi_subd_lat_par, coord2grid)
       firsttime = .FALSE.
    ELSE
       ! Do not call this routine again during later assimilation timesteps.
       ! *IMPORTANT* This assumes that observation coordinates are time
       ! invariant.
       IF (mype_filter ==0) THEN
          WRITE(*,'(/9x, a)') 'WARNING: Reusing coord2grid from first &
               & assimilation timestep. This assumes that longitude/latitude &
               & observation coordinates are the same for all time steps.'
       END IF
    END IF


! ***********************************************************
! *** Count available observations for the process domain ***
! *** and initialize index and coordinate arrays.         ***
! ***********************************************************

    ! *** Count valid observations that lie within the process sub-domain ***
    cnt_p = 0
    ! Use SIZE(coord2grid(:,1)) as upper limit for obs on local PE
    DO ind_obs = 1 , SIZE(coord2grid(:,1))
       i_obs = INT(coord2grid(ind_obs,1))
       j_obs = INT(coord2grid(ind_obs,2))
       ! Reject non-existent and unrealistic obs
       IF ( obs(i_obs, j_obs, 1) .GT. lowlim_ssh_NEMO) THEN
          IF ( obs(i_obs, j_obs, 1) .LT. upplim_ssh_NEMO) THEN
             i_pe = INT(coord2grid(ind_obs,3))
             j_pe = INT(coord2grid(ind_obs,4))
             ! Reject obs in gridbox with one land vertex.
             ! Current interpolation algorithm possibly unstable for land vertex.
             IF(.NOT. land_vertex_2d(i_pe,j_pe) ) THEN
                cnt_p=cnt_p+1
             END IF
          END IF
       END IF
    END DO

    ! Set number of local observations
    dim_obs_p = cnt_p

    IF(cnt_p == 0) THEN
       WRITE(*,'(/9x, a, i3, 3x, a, i3)') 'WARNING: No ssh_NEMO observations on &
            & local domain:', mype_filter, 'step_count:', step_count
    END IF

    ! *** Initialize vector of observations on the process sub-domain ***
    ! *** Initialize coordinate array of observations on the process sub-domain ***

    obs_nonzero:IF(dim_obs_p>0) THEN
       ! Allocate process-local observation arrays
       ALLOCATE(obs_p(dim_obs_p))
       ALLOCATE(ivar_obs_p(dim_obs_p))
       ALLOCATE(ocoord_p(2, dim_obs_p))
       ALLOCATE(ogrid_p(4, dim_obs_p))
       ! Allocate process-local index array
       ! This array has a many rows as required for the observation operator
       ! 1 if observations are at grid points; >1 if interpolation is required
       ALLOCATE(thisobs%id_obs_p(4, dim_obs_p))
       cnt_p = 0
       ! Use SIZE(coord2grid(:,1)) as upper limit for obs on local PE
       DO ind_obs = 1 , SIZE(coord2grid(:,1))
          i_obs = INT(coord2grid(ind_obs,1))
          j_obs = INT(coord2grid(ind_obs,2))
          ! Reject non-existent and unrealistic obs
          IF ( obs(i_obs, j_obs, 1) .GT. lowlim_ssh_NEMO) THEN
             IF ( obs(i_obs, j_obs, 1) .LT. upplim_ssh_NEMO) THEN
                i_pe = INT(coord2grid(ind_obs,3))
                j_pe = INT(coord2grid(ind_obs,4))
                ! Reject obs in gridbox with one land vertex.
                ! Current interpolation algorithm possibly unstable for land vertex.
                IF(.NOT. land_vertex_2d(i_pe,j_pe) ) THEN
                   ! Total number of valid observations on local PE
                   cnt_p = cnt_p + 1
                   ! Observation
                   obs_p(cnt_p) = obs(i_obs, j_obs, 1)
                   ! Observation coordinates - must be in radians
                   ocoord_p(1, cnt_p) = longitude_ssh_NEMO(i_obs)*rad_conv
                   ocoord_p(2, cnt_p) = latitude_ssh_NEMO(j_obs)*rad_conv
                   ! PE coordinates for observation gridbox
                   ogrid_p(1, cnt_p) = i_pe
                   ogrid_p(2, cnt_p) = j_pe
                   ! Index coordinates for observation
                   ogrid_p(3, cnt_p) = i_obs
                   ogrid_p(4, cnt_p) = j_obs
                   !
                   ! ***************************************************************************
                   ! Determine gridpoints for obs operator (great-circle weighted interpolation)
                   ! ***************************************************************************
                   !
                   ! The issue here is when gridpoints not contained in the state vector
                   ! are required for interpolation i.e. when halo regions are required. For the
                   ! interpolation, only halo regions to the East and North of a grid box are
                   ! required. (This is because of the convention that all grid boxes are counted
                   ! by their bottom left coordinate.) These East and North halo regions are
                   ! contained in a separate vector (see mod_statevector) and are indexed as
                   ! follows. The bottom right corner is given an index of 1, the grid box
                   ! directly above the bottom right corner is given an index of 2, and so on
                   ! until we reach the grid box directly below the top right corner (this is a
                   ! total of mpi_subd_lat_par points). We then move to the top left corner and
                   ! continue indexing from mpi_subd_lat_par + 1. We continue in a horizontal
                   ! fashion along the North boundary until eventually we reach the top right
                   ! corner. The top right corner is then indexed as (mpi_subd_lat_par +
                   ! mpi_subd_lon_par + 1).
                   !
                   ! A further refinement to the above indexing scheme is introduced here: all
                   ! halo region indices are assigned a negative value. This is necessary so
                   ! that the observation operator routine (see mod_obs_op) can determine
                   ! whether an index belongs to the state vector or the halo vector.
                   !
                   ! To determine how many halo values are required for each grid box, there
                   ! are 4 cases to consider. The first case is for the gridbox in the top right
                   ! corner. The second case is for all gridboxes on the Eastern boundary and the
                   ! third case is for all gridboxes on the Northern boundary. The fourth case is
                   ! when no halo values are required.
                   !
                   ! Case 1: Halo regions required for 3 gridpoints
                   IF (i_pe == mpi_subd_lon_par .AND. j_pe == mpi_subd_lat_par) THEN
                      ! Bottom left
                      thisobs%id_obs_p(1, cnt_p) = i_pe + (j_pe-1)*mpi_subd_lon_par &
                           + ssh_p_offset_par
                      ! Bottom right
                      thisobs%id_obs_p(2, cnt_p) = -mpi_subd_lat_par
                      ! Top left
                      thisobs%id_obs_p(3, cnt_p) = -mpi_subd_lat_par - mpi_subd_lon_par
                      ! Top right
                      thisobs%id_obs_p(4, cnt_p) = -mpi_subd_lat_par - mpi_subd_lon_par - 1
                   ! Case 2: Halo regions required for 2 gridpoints - Eastern boundary
                   ELSE IF (i_pe == mpi_subd_lon_par) THEN
                      ! Bottom left
                      thisobs%id_obs_p(1, cnt_p) = i_pe + (j_pe-1)*mpi_subd_lon_par &
                           + ssh_p_offset_par
                      ! Bottom right
                      thisobs%id_obs_p(2, cnt_p) = -j_pe
                      ! Top left
                      thisobs%id_obs_p(3, cnt_p) = i_pe + j_pe*mpi_subd_lon_par &
                           + ssh_p_offset_par
                      ! Top right
                      thisobs%id_obs_p(4, cnt_p) = -j_pe - 1
                   ! Case 3: Halo regions required for 2 gridpoints - Northern boundary
                   ELSE IF (j_pe == mpi_subd_lat_par) THEN
                      ! Bottom left
                      thisobs%id_obs_p(1, cnt_p) = i_pe + (j_pe-1)*mpi_subd_lon_par &
                           + ssh_p_offset_par
                      ! Bottom right
                      thisobs%id_obs_p(2, cnt_p) = (i_pe+1) + &
                           (j_pe-1)* mpi_subd_lon_par + ssh_p_offset_par
                      ! Top left
                      thisobs%id_obs_p(3, cnt_p) = -mpi_subd_lat_par - i_pe
                      ! Top right
                      thisobs%id_obs_p(4, cnt_p) = -mpi_subd_lat_par - i_pe - 1
                   ! Case 4: No halo regions required
                   ELSE
                      ! Bottom left
                      thisobs%id_obs_p(1, cnt_p) = i_pe + (j_pe-1)*mpi_subd_lon_par &
                           + ssh_p_offset_par
                      ! Bottom right
                      thisobs%id_obs_p(2, cnt_p) = (i_pe+1) + &
                           (j_pe-1)* mpi_subd_lon_par + ssh_p_offset_par
                      ! Top left
                      thisobs%id_obs_p(3, cnt_p) = i_pe + j_pe*mpi_subd_lon_par &
                           + ssh_p_offset_par
                      ! Top right
                      thisobs%id_obs_p(4, cnt_p) = (i_pe+1) + &
                           j_pe*mpi_subd_lon_par + ssh_p_offset_par
                   END IF
                   !
                   ! ****************************************************************
                END IF
             END IF
          END IF
       END DO
    ELSE
       ! No observations on PE, create dummy arrays to pass to PDAF-OMI
       ALLOCATE(obs_p(1))
       ALLOCATE(ivar_obs_p(1))
       ALLOCATE(ocoord_p(2, 1))
       ALLOCATE(ogrid_p(2, 1))
       ALLOCATE(thisobs%id_obs_p(4, 1))
       obs_p=-999999.0
       ivar_obs_p=-999999.0
       ocoord_p=-999999.0
       ogrid_p=-999999.0
       thisobs%id_obs_p=-999999.0
    END IF obs_nonzero


! **********************************************************************
! *** Initialize interpolation coefficients for observation operator ***
! **********************************************************************

    interp_nonzero:IF(dim_obs_p>0) THEN
       ! Allocate array of interpolation coefficients. As ID_OBS_P, the number
       ! of rows corresponds to the number of grid points using the the interpolation
       ALLOCATE(thisobs%icoeff_p(4, dim_obs_p))

       DO i = 1, dim_obs_p
          ! i,j coordinates for *bottom left* of grid box.
          i_pe = ogrid_p(1,i)
          j_pe = ogrid_p(2,i)

          ! Compute halo offset
          i0 = nldi_par - 1
          j0 = nldj_par - 1

          ! Convert from gridpoint to lon,lat
          ! bottom left
          lon_pe(1) = glamt_par( i_pe+i0, j_pe+j0 )
          ! Convert to 0-360 degree format
          IF(lon_pe(1) < 0.0_pwp) &
               lon_pe(1) = lon_pe(1) + 360.0_pwp
          lat_pe(1) = gphit_par( i_pe+i0, j_pe+j0 )

          ! bottom right
          lon_pe(2) = glamt_par( i_pe+i0+1, j_pe+j0 )
          ! Convert to 0-360 degree format
          IF(lon_pe(2) < 0.0_pwp) &
               lon_pe(2) = lon_pe(2) + 360.0_pwp
          lat_pe(2) = gphit_par( i_pe+i0+1, j_pe+j0 )

          ! top left
          lon_pe(3) = glamt_par( i_pe+i0, j_pe+j0+1 )
          ! Convert to 0-360 degree format
          IF(lon_pe(3) < 0.0_pwp) &
               lon_pe(3) = lon_pe(3) + 360.0_pwp
          lat_pe(3) = gphit_par( i_pe+i0, j_pe+j0+1 )

          ! top right
          lon_pe(4) = glamt_par( i_pe+i0+1, j_pe+j0+1 )
          ! Convert to 0-360 degree format
          IF(lon_pe(4) < 0.0_pwp) &
               lon_pe(4) = lon_pe(4) + 360.0_pwp
          lat_pe(4) = gphit_par( i_pe+i0+1, j_pe+j0+1 )

          ! i,j coordinates for obs.
          i_obs = ogrid_p(3,i)
          j_obs = ogrid_p(4,i)

          ! Obs coordinates - *degrees*, not radians
          lon_obs = longitude_ssh_NEMO(i_obs)
          lat_obs = latitude_ssh_NEMO(j_obs)

          ! Compute interpolation weights for case obs on vertex
          IF ( ABS( lon_pe(1)-lon_obs ) .LT. 1e-6_pwp .AND. &
               ABS( lat_pe(1)-lat_obs ) .LT. 1e-6_pwp ) THEN
             w_a=1.0_pwp
             w_b=0.0_pwp
             w_c=0.0_pwp
             w_d=0.0_pwp
             norm=1.0_pwp
          ELSE ! Otherwise use great circle-weighted interpolation
             dist_obs_a=grt_cir_dis(lon_obs, lat_obs, &
                  lon_pe(1), lat_pe(1))
             dist_obs_b=grt_cir_dis(lon_obs, lat_obs, &
                  lon_pe(2), lat_pe(2))
             dist_obs_c=grt_cir_dis(lon_obs, lat_obs, &
                  lon_pe(3), lat_pe(3))
             dist_obs_d=grt_cir_dis(lon_obs, lat_obs, &
                  lon_pe(4), lat_pe(4))
             ! This method follows NEMO docs...
             w_a=dist_obs_b*dist_obs_c*dist_obs_d
             w_b=dist_obs_a*dist_obs_c*dist_obs_d
             w_c=dist_obs_a*dist_obs_b*dist_obs_d
             w_d=dist_obs_a*dist_obs_b*dist_obs_c
             norm= w_a + w_b + w_c + w_d
          END IF

          ! Compute normalised weight
          thisobs%icoeff_p(1, i) = w_a / norm
          thisobs%icoeff_p(2, i) = w_b / norm
          thisobs%icoeff_p(3, i) = w_c / norm
          thisobs%icoeff_p(4, i) = w_d / norm
       END DO
    ELSE
       ALLOCATE(thisobs%icoeff_p(4, 1))
       thisobs%icoeff_p=-999999.0
    END IF interp_nonzero


! ****************************************************************
! *** Define observation errors for process-local observations ***
! ****************************************************************

    ! *** Set inverse observation error variances ***

    ivar_obs_p(:) = 1.0 / (rms_ssh_NEMO*rms_ssh_NEMO)


! ****************************************
! *** Gather global observation arrays ***
! ****************************************

    CALL PDAFomi_gather_obs_f(thisobs, dim_obs_p, obs_p, ivar_obs_p, ocoord_p, &
         thisobs%ncoord, local_range, dim_obs_f)


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
    DEALLOCATE(longitude_ssh_NEMO, latitude_ssh_NEMO)
    DEALLOCATE(obs_p, ocoord_p, ogrid_p, ivar_obs_p)

    ! Arrays in THISOBS have to be deallocated after the analysis step
    ! by a call to deallocate_obs() in prepoststep_pdaf.

  END SUBROUTINE init_dim_obs_f_ssh_NEMO


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
  SUBROUTINE obs_op_f_ssh_NEMO(dim_p, dim_obs_f, state_p, ostate_f, offset_obs)

    USE mod_obs_op_pdaf, &
         ONLY: ssh_par_obs_op_gcirc

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
       ! Great-circle interpolation obs operator for ssh_NEMO
       CALL ssh_par_obs_op_gcirc(thisobs, 4, state_p, ostate_f, offset_obs)
    END IF

  END SUBROUTINE obs_op_f_ssh_NEMO


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
  SUBROUTINE init_dim_obs_l_ssh_NEMO(domain_p, step, dim_obs_f, dim_obs_l, &
       off_obs_l, off_obs_f)

    ! Include PDAFomi function
    USE PDAFomi, ONLY: PDAFomi_init_dim_obs_l, PDAFomi_set_debug_flag

    ! Include localization radius and local coordinates
    USE mod_assimilation_pdaf, &   
         ONLY: coords_l, local_range, locweight, srange

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

!!$    ! Gridpoint debugging
!!$    IF (domain_p == 1019 .AND. mype_filter==9) THEN
!!$       CALL PDAFomi_set_debug_flag(1019)
!!$    ELSE IF (domain_p == 959 .AND. mype_filter==9) THEN
!!$       CALL PDAFomi_set_debug_flag(959)
!!$    ELSE IF (domain_p == 899 .AND. mype_filter==9) THEN
!!$       CALL PDAFomi_set_debug_flag(899)
!!$    ELSE
!!$       CALL PDAFomi_set_debug_flag(0)
!!$    ENDIF

    CALL PDAFomi_init_dim_obs_l(thisobs_l, thisobs, coords_l, &
         locweight, local_range, srange, &
         dim_obs_l, off_obs_l, off_obs_f)

  END SUBROUTINE init_dim_obs_l_ssh_NEMO


!-------------------------------------------------------------------------------
  SUBROUTINE obs_in_grid_lcl (longitude, latitude, idim_pe, jdim_pe, obs_gp_lcl)
    ! **********************************************************
    !                     ROUTINE obs_in_grid_lcl
    !
    ! Purpose : Determine number of observations on local PE,
    !           and gridbox surrounding each observation.
    !
    ! **********************************************************

    USE mod_parallel_pdaf, &
         ONLY: nldj_par, nldi_par, mype_ens, COMM_filter
    USE mod_agrif_pdaf, &
         ONLY: glamt_par, gphit_par
    USE mod_parallel_pdaf, &
         ONLY: jpni_par, jpnj_par

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    INTEGER, INTENT(in)          :: idim_pe      ! Size of PE along i-axis
    INTEGER, INTENT(in)          :: jdim_pe      ! Size of PE along j-axis
    REAL(pwp), INTENT(in)        :: longitude(:) ! Longitude coordinates for obs
    REAL(pwp), INTENT(in)        :: latitude(:)  ! Latitude coordiantes for obs
    REAL(pwp), ALLOCATABLE, INTENT(out) :: obs_gp_lcl(:,:) ! Coordinates of obs
    ! on local PE, as well as coordinates of local PE gridbox surrounding obs

    ! Local variables
    INTEGER :: i, j, ii, jj, s     ! Counters
    INTEGER :: col, row, mid_col   ! Values to compute modified local PE dimension
    INTEGER :: cnt, wrap_cnt       ! Counters for wraparound case
    INTEGER :: i_pe, j_pe          ! Local PE coordinates
    INTEGER :: i0, j0              ! Halo offset for local PE
    INTEGER :: ilim, jlim          ! Limits for local PE
    REAL(pwp) :: lat_obs, lon_obs ! Longitude/latitude values for obs
    REAL(pwp) :: lonmax, lonmin   ! Bounds for dealing with wraparound case
    REAL(pwp) :: c_prod_13, c_prod_21, c_prod_34, c_prod_42 ! Cross products
    REAL(pwp), ALLOCATABLE :: grid_lon(:,:,:), grid_lat(:,:,:) ! Gridbox coordinates
    REAL(pwp), ALLOCATABLE :: grid_lon_wrap(:,:), grid_lat_wrap(:,:) ! Gridbox
    ! coordinates for wraparound case
    REAL(pwp), ALLOCATABLE   :: grid_lon_max(:,:),grid_lon_min(:,:) ! Useful limits
    REAL(pwp), ALLOCATABLE   :: grid_lat_max(:,:),grid_lat_min(:,:) ! Useful limits
    INTEGER :: tot_obs_gp_lcl ! Total number of observations on local PE
    INTEGER, ALLOCATABLE :: temp(:,:,:) ! Gridbox info on observations on local PE
    INTEGER, ALLOCATABLE :: debug_arr(:), debug_arr2(:) ! Arrays useful for debug


! ********************************************
! Compute array of all (lon, lat) points on PE
! ********************************************

    ! First we compute a modified dimension of the local PE, depending on how the
    ! MPI decomposition is performed. This is necessary as due to the E/W and NP
    ! periodicities, certain gridboxes will occur on multiple PEs. Hence, we
    ! introduce a modified dimension so (almost) all gridboxes will only occur
    ! on a single PE. This is necessary so that we don't count the same observation
    ! multiple times when we assign it a gridbox.
    !
    col=MOD(mype_ens,jpni_par)
    row=INT(mype_ens/jpni_par)
    mid_col=INT(jpni_par/2)
    IF( (col == jpni_par-1) .AND. (row == jpnj_par-1) ) THEN
       ilim = idim_pe - 1
       jlim = jdim_pe - 1
    ELSE IF(col == jpni_par-1)  THEN
       ilim = idim_pe - 1
       jlim = jdim_pe
    ELSE IF( (row == jpnj_par - 1) .AND. (col < mid_col) )THEN
       ilim = idim_pe
       jlim = jdim_pe - 2
    ELSE IF( row == jpnj_par - 1 ) THEN
       ilim = idim_pe
       jlim = jdim_pe - 1
    ELSE
       ilim = idim_pe
       jlim = jdim_pe
    END IF

    ! Consider a gridbox (i,j), (i+1,j), (i,j+1), (i+1,j+1).
    ! We index this gridbox by its bottom left coordinate (i,j). There
    ! are (ilim)*(jlim) such coordinates. Each gridbox contains
    ! four vertices. Hence there are 4*(ilim)*(jlim) entries. We
    ! also add 1 extra entry to store original gridbox coordinates for
    ! wraparound case. Hence we have 5*(ilim)*(jlim) entries.
    ALLOCATE( grid_lon(5,ilim,jlim) )
    ALLOCATE( grid_lat(5,ilim,jlim) )

    ! *******************
    ! Non-wraparound case
    ! *******************
    !
    ! Setup array of global coordinates on PE for non-wraparound case.
    DO j = 1, jlim
       DO i = 1, ilim
          ! i,j coordinates for *bottom left* of grid box.
          i_pe = i
          j_pe = j

          ! Compute halo offset.
          i0 = nldi_par - 1
          j0 = nldj_par - 1

          ! bottom left
          grid_lon(1,i,j) = glamt_par( i_pe+i0, j_pe+j0 )
          ! Convert to 0-360 degree format
          IF(grid_lon(1,i,j) < 0.0_pwp) &
               grid_lon(1,i,j) = grid_lon(1,i,j) + 360.0_pwp
          grid_lat(1,i,j) = gphit_par( i_pe+i0, j_pe+j0 )

          ! bottom right
          grid_lon(2,i,j) = glamt_par( i_pe+i0+1, j_pe+j0 )
          ! Convert to 0-360 degree format
          IF(grid_lon(2,i,j) < 0.0_pwp) &
               grid_lon(2,i,j) = grid_lon(2,i,j) + 360.0_pwp
          grid_lat(2,i,j) = gphit_par( i_pe+i0+1, j_pe+j0 )

          ! top left
          grid_lon(3,i,j) = glamt_par( i_pe+i0, j_pe+j0+1 )
          ! Convert to 0-360 degree format
          IF(grid_lon(3,i,j) < 0.0_pwp) &
               grid_lon(3,i,j) = grid_lon(3,i,j) + 360.0_pwp
          grid_lat(3,i,j) = gphit_par( i_pe+i0, j_pe+j0+1 )

          ! top right
          grid_lon(4,i,j) = glamt_par( i_pe+i0+1, j_pe+j0+1 )
          ! Convert to 0-360 degree format
          IF(grid_lon(4,i,j) < 0.0_pwp) &
               grid_lon(4,i,j) = grid_lon(4,i,j) + 360.0_pwp
          grid_lat(4,i,j) = gphit_par( i_pe+i0+1, j_pe+j0+1 )

          ! Store bottom left vertex. Redundant here, but needed later
          ! for wraparound case.
          grid_lon(5,i,j) = i
          grid_lat(5,i,j) = j
       END DO
    END DO

    ! ***************
    ! Wraparound case
    ! ***************
    !
    ! We are required to consider case where part of a gridbox crosses the
    ! prime meridian (Greenwich). For example, consider a regular gridbox
    ! with longitudes 359 and 5. In this case we use two regular gridboxes.
    ! One gridbox has longitudes 359 and 365, the second has -1 and 5.
    !
    ! Begin by creating 'second' gridbox values as described above.
    cnt=0
    DO j = 1, jlim
       DO i = 1, ilim
          ! Crude hack to avoid messy calculations near poles
          IF (grid_lat(1,i,j) < 89.0 .AND. grid_lon(1,i,j) > 310.0) THEN
             lonmin = MINVAL( grid_lon(1:4,i,j) )
             IF ( MAXVAL(grid_lon(1:4,i,j)) > 180.0_pwp+lonmin ) cnt=cnt+1
          END IF
       END DO
    END DO

    IF (cnt > 0) THEN
       ALLOCATE(grid_lon_wrap(5,cnt))
       ALLOCATE(grid_lat_wrap(5,cnt))
       cnt=0
       DO j = 1, jlim
          DO i = 1, ilim
             ! Crude hack to avoid messy calculations near poles
             IF (grid_lat(1,i,j) < 89.0 .AND. grid_lon(1,i,j) > 310.0) THEN
                lonmin = MINVAL( grid_lon(1:4,i,j) )
                IF ( MAXVAL(grid_lon(1:4,i,j)) > 180.0_pwp+lonmin ) THEN
                   cnt=cnt+1
                   ! latitudes stay the same
                   grid_lat_wrap(1,cnt)=grid_lat(1,i,j)
                   grid_lat_wrap(2,cnt)=grid_lat(2,i,j)
                   grid_lat_wrap(3,cnt)=grid_lat(3,i,j)
                   grid_lat_wrap(4,cnt)=grid_lat(4,i,j)
                   ! change the longitudes
                   grid_lon_wrap(1,cnt)=grid_lon(1,i,j)-360.0_pwp
                   grid_lon_wrap(2,cnt)=grid_lon(2,i,j)-360.0_pwp
                   grid_lon_wrap(3,cnt)=grid_lon(3,i,j)-360.0_pwp
                   grid_lon_wrap(4,cnt)=grid_lon(4,i,j)-360.0_pwp
                   ! Don't change longitudes east of prime meridian
                   WHERE (grid_lon_wrap(1:4,cnt) < -180.0_pwp) &
                        grid_lon_wrap(1:4,cnt)=grid_lon_wrap(1:4,cnt)+360.0_pwp
                   ! Retain information about original gridbox location
                   grid_lon_wrap(5,cnt)=grid_lon(5,i,j)
                   grid_lat_wrap(5,cnt)=grid_lat(5,i,j)
                END IF
             END IF
          END DO
       END DO
    END IF

    ! Store number of wraparound gridboxes for later
    wrap_cnt = cnt

    ! Modify original gridbox value for first wraparound case - must be
    ! done *after* second wraparound case above.
    DO j = 1, jlim
       DO i = 1, ilim
          ! Crude hack to avoid messy calculations near poles
          IF (grid_lat(1,i,j) < 89.0 .AND. grid_lon(1,i,j) > 310.0) THEN
             lonmax = MAXVAL( grid_lon(1:4,i,j) )
             WHERE (lonmax-grid_lon(1:4,i,j) > 180.0_pwp) &
                  grid_lon(1:4,i,j)=grid_lon(1:4,i,j)+360.0_pwp
          END IF
       END DO
    END DO


! ************************************
! Determine whether observations on PE
! ************************************

    ! *********************************
    ! Compute some useful bounds/limits
    ! *********************************
    !
    ! Compute coordinates for *regular* gridbox that contains all of gridbox
    ! from NEMO. Obs *must* be in regular gridbox if in gridbox from NEMO.
    ! Use regular gridbox coordinates as upper/lower bounds later.
    ALLOCATE( grid_lat_max(ilim,jlim) )
    ALLOCATE( grid_lat_min(ilim,jlim) )
    ALLOCATE( grid_lon_max(ilim,jlim) )
    ALLOCATE( grid_lon_min(ilim,jlim) )

    DO j=1,jlim
       DO i=1,ilim
          grid_lat_max(i,j) = MAXVAL( grid_lat(1:4,i,j) )
          grid_lat_min(i,j) = MINVAL( grid_lat(1:4,i,j) )
          grid_lon_max(i,j) = MAXVAL( grid_lon(1:4,i,j) )
          grid_lon_min(i,j) = MINVAL( grid_lon(1:4,i,j) )
       END DO
    END DO

    ! *************************************
    ! Enter main loop - non-wraparound case
    ! *************************************
    !
    tot_obs_gp_lcl=0

    ! Temporary array to store if obs in gridbox and if so, to also store
    ! gridbox coordinates.
    ALLOCATE( temp(SIZE(longitude),SIZE(latitude),3) )
    temp(:,:,:) = 0.0_pwp
    DO jj = 1, SIZE(latitude)
       DO ii = 1, SIZE(longitude)
          lon_obs = longitude(ii)
          lat_obs = latitude(jj)

          ! First consider case where obs located on vertex.
          ! We consider the bottom left vertex of each gridbox.
          vertexloop:DO j = 1, jlim
             DO i = 1, ilim
                IF ( ABS(grid_lat(1,i,j)-lat_obs) .LT. 1e-6_pwp ) THEN
                   IF ( ABS(grid_lon(1,i,j)-lon_obs) .LT. 1e-6_pwp ) THEN
                      temp(ii,jj,1)=temp(ii,jj,1)+1.0_pwp ! Indicates obs on PE
                      temp(ii,jj,2)=grid_lon(5,i,j)
                      temp(ii,jj,3)=grid_lat(5,i,j)
                      tot_obs_gp_lcl=tot_obs_gp_lcl+1
                      EXIT vertexloop
                   END IF
                END IF
             END DO
          END DO vertexloop

          ! Next consider case where obs located inside gridbox
          ! *********************************************************
          ! To determine whether obs contained in gridbox, we compute
          ! cross products of co-planar vectors. If all values
          ! are negative, then obs is contained in gridbox. If a value
          ! is zero then it is contained on a gridbox boundary. All
          ! obs that are contained on the left or bottom boundaries
          ! of a gridbox are considered to be contained in a gridbox.
          ! Based on notes from 'NEMO documentation, section:
          ! theoretical details of obs operator'.
          ! *********************************************************

          notvertex:IF (temp(ii,jj,1) .EQ. 0.0_pwp) THEN 
             gridloop:DO j = 1, jlim
                DO i = 1, ilim

                   ! Conditions when obs definitely not in gridbox
                   IF ( ABS(lat_obs) .LT. 85.0_pwp ) THEN ! Included from NEMO
                      ! obs operator source code...
                      IF ( (lat_obs .GT. grid_lat_max(i,j) ) .OR. &
                           & (lat_obs .LT. grid_lat_min(i,j)) ) CYCLE
                   END IF
                   IF ( (lon_obs .GT. grid_lon_max(i,j) ) .OR. &
                        & (lon_obs .LT. grid_lon_min(i,j) ) ) CYCLE

                   ! Cross product anti-commutative - careful with order
                   CALL cross_prod_2d(grid_lon(1,i,j),grid_lat(1,i,j),&
                        grid_lon(3,i,j),grid_lat(3,i,j),lon_obs,lat_obs,&
                        c_prod_13)
                   CALL cross_prod_2d(grid_lon(2,i,j),grid_lat(2,i,j),&
                        grid_lon(1,i,j),grid_lat(1,i,j),lon_obs,lat_obs,&
                        c_prod_21)
                   CALL cross_prod_2d(grid_lon(3,i,j),grid_lat(3,i,j),&
                        grid_lon(4,i,j),grid_lat(4,i,j),lon_obs,lat_obs,&
                        c_prod_34)
                   CALL cross_prod_2d(grid_lon(4,i,j),grid_lat(4,i,j),&
                        grid_lon(2,i,j),grid_lat(2,i,j),lon_obs,lat_obs,&
                        c_prod_42)

                   ! Allow cross products for left hand and bottom
                   ! boundaries to be zero
                   IF (c_prod_13 .LE. 0.0_pwp .AND. c_prod_21 .LE. 0.0_pwp &
                        .AND. c_prod_34 .LT. 0.0_pwp .AND. &
                        c_prod_42 .LT. 0.0_pwp) THEN
                      temp(ii,jj,1)=temp(ii,jj,1)+1.0_pwp ! Indicates obs on PE
                      temp(ii,jj,2)=grid_lon(5,i,j)
                      temp(ii,jj,3)=grid_lat(5,i,j)
                      tot_obs_gp_lcl=tot_obs_gp_lcl+1
                      EXIT gridloop
                   END IF
                END DO
             END DO gridloop
          END IF notvertex
       END DO
    END DO

    ! Clean-up
    DEALLOCATE(grid_lon,grid_lat,grid_lat_max,grid_lat_min,grid_lon_max,&
         grid_lon_min)

    wrap_true:IF (wrap_cnt > 0) THEN
       ! *********************************
       ! Enter main loop - wraparound case
       ! *********************************
       !
       ! Temporary array to store if obs in gridbox and if so, to also store
       ! gridbox coordinates.
       DO jj = 1, SIZE(latitude)
          DO ii = 1, SIZE(longitude)
             lat_obs = latitude(jj)
             lon_obs = longitude(ii)

             ! First consider case where obs located on vertex.
             ! We consider the bottom left vertex of each gridbox.
             wrap_vertexloop:DO i = 1, wrap_cnt
                IF ( ABS(grid_lat_wrap(1,i)-lat_obs) .LT. 1e-6_pwp ) THEN
                   IF ( ABS(grid_lon_wrap(1,j)-lon_obs) .LT. 1e-6_pwp ) THEN
                      temp(ii,jj,1)=temp(ii,jj,1)+1.0_pwp ! Indicates obs on PE
                      temp(ii,jj,2)=grid_lon_wrap(5,i)
                      temp(ii,jj,3)=grid_lat_wrap(5,i)
                      tot_obs_gp_lcl=tot_obs_gp_lcl+1
                      EXIT wrap_vertexloop
                   END IF
                END IF
             END DO wrap_vertexloop

             ! Next consider case where obs located inside gridbox
             ! *********************************************************
             ! To determine whether obs contained in gridbox, we compute
             ! cross products of co-planar vectors. If all values
             ! are negative, then obs is contained in gridbox. If a value
             ! is zero then it is contained on a gridbox boundary. All
             ! obs that are contained on the left or bottom boundaries
             ! of a gridbox are considered to be contained in a gridbox.
             ! Based on notes from 'NEMO documentation, section:
             ! theoretical details of obs operator'.
             ! *********************************************************

             wrap_notvertex:IF (temp(ii,jj,1) .EQ. 0.0_pwp) THEN
                wrap_gridloop:DO i = 1, wrap_cnt
                   ! Cross product anti-commutative - careful with order
                   CALL cross_prod_2d(grid_lon_wrap(1,i), grid_lat_wrap(1,i), &
                        grid_lon_wrap(3,i), grid_lat_wrap(3,i), lon_obs, &
                        lat_obs, c_prod_13)
                   CALL cross_prod_2d(grid_lon_wrap(2,i), grid_lat_wrap(2,i), &
                        grid_lon_wrap(1,i),grid_lat_wrap(1,i), lon_obs, &
                        lat_obs, c_prod_21)
                   CALL cross_prod_2d(grid_lon_wrap(3,i), grid_lat_wrap(3,i), &
                        grid_lon_wrap(4,i), grid_lat_wrap(4,i), lon_obs, &
                        lat_obs, c_prod_34)
                   CALL cross_prod_2d(grid_lon_wrap(4,i), grid_lat_wrap(4,i), &
                        grid_lon_wrap(2,i), grid_lat_wrap(2,i), lon_obs, &
                        lat_obs, c_prod_42)

                   ! Allow cross products for left hand and bottom
                   ! boundaries to be zero
                   IF (c_prod_13 .LE. 0.0_pwp .AND. c_prod_21 .LE. 0.0_pwp &
                        .AND. c_prod_34 .LT. 0.0_pwp .AND. &
                        c_prod_42 .LT. 0.0_pwp) THEN
                      temp(ii,jj,1)=temp(ii,jj,1)+1.0_pwp ! Indicates obs on PE
                      temp(ii,jj,2)=grid_lon_wrap(5,i)
                      temp(ii,jj,3)=grid_lat_wrap(5,i)
                      tot_obs_gp_lcl=tot_obs_gp_lcl+1
                      EXIT wrap_gridloop
                   END IF
                END DO wrap_gridloop
             END IF wrap_notvertex
          END DO
       END DO

       ! Clean-up
       DEALLOCATE(grid_lon_wrap,grid_lat_wrap)

    END IF wrap_true


! *************************
! Format output for PDAFomi
! *************************

    IF(tot_obs_gp_lcl .EQ. 0) THEN
       WRITE(*,*)'ERROR: No obs on local PE. May need to modify algorithm.'
       CALL abort_parallel()
    END IF

    ! Reduce 3d array of obs information to 2d array
    ALLOCATE( obs_gp_lcl(tot_obs_gp_lcl,4) )
    s=0 ! temporary counter
    DO j=1,SIZE(latitude)
       DO i=1,SIZE(longitude)
          IF (temp(i,j,1) .EQ. 1.0_pwp) THEN
             s=s+1
             obs_gp_lcl(s,1)=i
             obs_gp_lcl(s,2)=j
             obs_gp_lcl(s,3)=temp(i,j,2) ! gridbox coordinate
             obs_gp_lcl(s,4)=temp(i,j,3) ! gridbox coordinate
          END IF
       END DO
    END DO

    !! Useful code for debugging (in association with plot.py).
    !! Determines how many times an observation is counted across all PEs.
    !s=0
    !ALLOCATE(debug_arr(SIZE(latitude)*SIZE(longitude)))
    !debug_arr=0
    !DO jj = 1, SIZE(latitude)
    !   DO ii = 1, SIZE(longitude)
    !      s=s+1
    !      debug_arr(s)=temp(ii,jj,1)
    !   END DO
    !END DO
    !ALLOCATE(debug_arr2(SIZE(latitude)*SIZE(longitude)))
    !debug_arr2=0
    !CALL MPI_AllReduce(debug_arr, debug_arr2, s, MPI_INT, MPI_Sum, COMM_filter, i )
    !      !CALL MPI_Barrier(COMM_filter, i)
    !IF (mype_filter == 0) THEN
    !   OPEN(31, file = 'debug_arr2.txt', status = 'new')
    !   DO i =1,s
    !      WRITE (31, *) debug_arr2(i)
    !   END DO
    !   CLOSE(31)
    !   WRITE(*,*) 's test',s,size(debug_arr2),debug_arr2
    !END IF
    !! Hack to make sure that output is correctly written before potential
    !! crash elsehwere!
    !CALL MPI_Barrier(COMM_filter, i)
    !DEALLOCATE(debug_arr,debug_arr2)

    ! Final clean-up
    DEALLOCATE(temp)

  END SUBROUTINE obs_in_grid_lcl


  SUBROUTINE cross_prod_2d(a1,a2,b1,b2,p1,p2,c_prod)
    ! *******************************************************************
    !                     ROUTINE cross_prod_2d
    !
    ! Purpose : Compute cross product of coplanar vectors.
    !
    ! Details : Consider points a=(a1,a2,c), b=(b1,b2,c) and p=(p1,p2,c).
    !           Consider vectors _pa_ and _pb_. We compute cross product:
    !            _pa_ X _pb_ = (0,0,c_prod).
    !
    ! *******************************************************************


    IMPLICIT NONE

    REAL(pwp), INTENT(in)  :: a1,a2,b1,b2
    REAL(pwp), INTENT(in)  :: p1,p2
    REAL(pwp), INTENT(out) :: c_prod

    ! Local variables
    REAL(pwp) :: d1,d2,d3,d4


    d1 = a1-p1
    d2 = b2-p2
    d3 = b1-p1
    d4 = a2-p2
    c_prod = (d1*d2) - (d3*d4)

  END SUBROUTINE cross_prod_2d


  FUNCTION land_vertex_2d(vert_i,vert_j) RESULT(land)
    ! *******************************************************************
    !                     ROUTINE land_vertex_2d
    !
    ! Purpose : Determine whether gridbox contains vertex on land. Input
    !           is **bottom left corner** vertex of gridbox.
    ! *******************************************************************


    USE mod_parallel_pdaf, &
         ONLY: nldi_par, nldj_par
    USE mod_agrif_pdaf, &
         ONLY: tmask_par

    IMPLICIT NONE

    INTEGER, INTENT(in) :: vert_i, vert_j
    LOGICAL :: land

    INTEGER :: i0, j0                        ! Halo offset for local PE


    ! Compute halo offset.
    i0 = nldi_par - 1
    j0 = nldj_par - 1

    land=.FALSE.
    IF(tmask_par(vert_i+i0,vert_j+j0,1) .EQ. 0)     land=.TRUE. ! bottom left vertex
    IF(tmask_par(vert_i+i0,vert_j+j0+1,1) .EQ. 0)   land=.TRUE. ! top left vertex
    IF(tmask_par(vert_i+i0+1,vert_j+j0,1) .EQ. 0)   land=.TRUE. ! bottom right vertex
    IF(tmask_par(vert_i+i0+1,vert_j+j0+1,1) .EQ. 0) land=.TRUE. ! top right vertex

  END FUNCTION land_vertex_2d


  FUNCTION grt_cir_dis(obs_lon,obs_lat,vert_lon,vert_lat) RESULT(arclen)
    ! *******************************************************************
    !                     ROUTINE grt_cir_dis
    !
    ! Purpose : Compute great-circle distance between (obs_lon,obs_lat)
    !           and (vert_lon,vert_lat), where coordinates in *degrees*.
    !
    ! Details : The formula used is a variant of the 'haversine formula'.
    !           Given A=(phi1,lam1) and B=(phi2,lam2) in *radians*, then
    !           the distance between these two points d(AB) is given by:
    !
    !           d(AB)=2*arcsin(sqrt( (1-x)/2 ), where x=
    !           sin(phi1)*sin(phi2) +
    !           cos(phi1)*cos(phi2)*cos(lam1)*cos(lam2) +
    !           cos(phi1)*cos(phi2)*sin(lam1)*sin(lam2)
    !
    ! *******************************************************************


    IMPLICIT NONE

    REAL(pwp), INTENT(in) :: obs_lon, obs_lat   ! Obs coordinates in degrees
    REAL(pwp), INTENT(in) :: vert_lon, vert_lat ! Vertex coordinates in degrees
    REAL(pwp) :: arclen ! distance d(AB)

    ! Local variables
    REAL(pwp) :: obs_lat_rad, obs_lon_rad    ! Radian conversion of obs
    REAL(pwp) :: vert_lat_rad, vert_lon_rad  ! Radian conversion of vertex
    REAL(pwp) :: x     ! Temporary value used for calculating haversine formula


    ! Convert to radians
    obs_lat_rad=obs_lat*rad_conv
    obs_lon_rad=obs_lon*rad_conv
    vert_lat_rad=vert_lat*rad_conv
    vert_lon_rad=vert_lon*rad_conv

    ! Compute variant of haversine formula
    x= ( SIN(vert_lat_rad)*SIN(obs_lat_rad) ) + &
         ( COS(vert_lat_rad)*COS(vert_lon_rad)*COS(obs_lat_rad)*COS(obs_lon_rad) ) + &
         ( COS(vert_lat_rad)*SIN(vert_lon_rad)*COS(obs_lat_rad)*SIN(obs_lon_rad) )

    arclen=2*ASIN(SQRT( (1.0-x)/2 ))

  END FUNCTION grt_cir_dis

!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_obs_ssh_NEMO_pdafomi
