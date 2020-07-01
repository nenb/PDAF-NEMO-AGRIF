MODULE mod_statevector_pdaf
!$AGRIF_DO_NOT_TREAT
! !DESCRIPTION:
! This module provides variables & routines for
! manipulating the state vector.

  ! !USES:
  USE mod_kind_pdaf

  IMPLICIT NONE
  SAVE

! *********************
! Parent grid variables
! *********************
! 2d state vector variables - start index
  INTEGER :: ssh_p_offset_par

! Array holding 2d state variable offsets
  INTEGER :: var2d_p_offset_par(1)

! 2d state vector variables - dimension size
  INTEGER :: ssh_p_dim_state_par

! 2d global state vector variables - dimension size
  INTEGER :: ssh_dim_state_par

! 3d state vector variables - start index
  INTEGER :: t_p_offset_par
  INTEGER :: s_p_offset_par
  INTEGER :: u_p_offset_par
  INTEGER :: v_p_offset_par

! Array holding 3d state variable offsets
  INTEGER :: var3d_p_offset_par(4)

! 3d state vector variables - dimension size
  INTEGER :: t_p_dim_state_par
  INTEGER :: s_p_dim_state_par
  INTEGER :: u_p_dim_state_par
  INTEGER :: v_p_dim_state_par

! 3d global state vector variables - dimension size
  INTEGER :: t_dim_state_par
  INTEGER :: s_dim_state_par
  INTEGER :: u_dim_state_par
  INTEGER :: v_dim_state_par

! Dimensions for MPI subdomain
  INTEGER :: mpi_subd_lat_par
  INTEGER :: mpi_subd_lon_par
  INTEGER :: mpi_subd_vert_par

! ********************
! Child grid variables
! ********************
! 2d state vector variables - start index
  INTEGER :: ssh_p_offset_child

! Array holding 2d state variable offsets
  INTEGER :: var2d_p_offset_child(1)

! 2d state vector variables - dimension size
  INTEGER :: ssh_p_dim_state_child

! 2d global state vector variables - dimension size
  INTEGER :: ssh_dim_state_child

! 3d state vector variables - start index
  INTEGER :: t_p_offset_child
  INTEGER :: s_p_offset_child
  INTEGER :: u_p_offset_child
  INTEGER :: v_p_offset_child

! Array holding 3d state variable offsets
  INTEGER :: var3d_p_offset_child(4)

! 3d state vector variables - dimension size
  INTEGER :: t_p_dim_state_child
  INTEGER :: s_p_dim_state_child
  INTEGER :: u_p_dim_state_child
  INTEGER :: v_p_dim_state_child

! 3d global state vector variables - dimension size
  INTEGER :: t_dim_state_child
  INTEGER :: s_dim_state_child
  INTEGER :: u_dim_state_child
  INTEGER :: v_dim_state_child

! Dimensions for MPI subdomain
  INTEGER :: mpi_subd_lat_child
  INTEGER :: mpi_subd_lon_child
  INTEGER :: mpi_subd_vert_child

  ! Array for 2d/3d state variable .NC ids
  CHARACTER(len=20), DIMENSION(1) :: id2d_list
  CHARACTER(len=20), DIMENSION(4) :: id3d_list

! Fill array of 2d/3d state variable .NC ids
  DATA id2d_list / 'sossheig' /
  DATA id3d_list / 'votemper', 'vosaline', 'vozocrtx', 'vomecrty' /

CONTAINS

  SUBROUTINE calc_mpi_dim()

    ! !DESCRIPTION:
    ! This routine calculates the dimensions of the MPI subdomain
    ! that is used to fill the statevector.

    ! !USES:
    USE mod_parallel_pdaf, ONLY: abort_parallel, jpk_par, jpiglo_par, jpjglo_par, &
         nldi_par, nldj_par, nlei_par, nlej_par, nimpp_par, njmpp_par, jpk_child, &
         jpiglo_child, jpjglo_child, nldi_child, nldj_child, nlei_child, nlej_child, &
         nimpp_child, njmpp_child

    IMPLICIT NONE


    mpi_subd_lon_par = nlei_par-nldi_par+1
    mpi_subd_lat_par = nlej_par-nldj_par+1
    mpi_subd_vert_par = jpk_par

    ! Treat case where MPI subdomain goes outside limits of
    ! global domain. In this case, MPI subdomain is likely
    ! padded with zeros.
    ! (Not sure if this is necessary, but included just in case.)
    IF(nimpp_par + mpi_subd_lon_par - 1 .GT. jpiglo_par) THEN
       WRITE(*,*) 'ERROR: MPI subdomain is padded with zeros along i&
       &dim. Need to modify calc_mpi_dim routine in mod_statevector&
       &to remove this zero padding. Parent grid.'
       CALL abort_parallel()
    END IF
    IF(njmpp_par + mpi_subd_lat_par - 1 .GT. jpjglo_par) THEN
       WRITE(*,*) 'ERROR: MPI subdomain is padded with zeros along j&
       &dim. Need to modify calc_mpi_dim routine in mod_statevector&
       &to remove this zero padding. Parent grid.'
       CALL abort_parallel()
    END IF

#if defined key_agrif
    mpi_subd_lon_child = nlei_child-nldi_child+1
    mpi_subd_lat_child = nlej_child-nldj_child+1
    mpi_subd_vert_child = jpk_child

    ! Treat case where MPI subdomain goes outside limits of
    ! global domain. In this case, MPI subdomain is likely
    ! padded with zeros.
    ! (Not sure if this is necessary, but included just in case.)
    IF(nimpp_child + mpi_subd_lon_child - 1 .GT. jpiglo_child) THEN
       WRITE(*,*) 'ERROR: MPI subdomain is padded with zeros along i&
       &dim. Need to modify calc_mpi_dim routine in mod_statevector&
       &to remove this zero padding. Child grid.'
       CALL abort_parallel()
    END IF
    IF(njmpp_child + mpi_subd_lat_child - 1 .GT. jpjglo_child) THEN
       WRITE(*,*) 'ERROR: MPI subdomain is padded with zeros along j&
       &dim. Need to modify calc_mpi_dim routine in mod_statevector&
       &to remove this zero padding. Child grid.'
       CALL abort_parallel()
    END IF
#endif

  END SUBROUTINE calc_mpi_dim

  SUBROUTINE calc_dim()

    ! !DESCRIPTION:
    ! This routine calculates the dimension of the state variables.

    ! !USES:
    IMPLICIT NONE


    ! Compute MPI subdomain dimensions in case not already done
    CALL calc_mpi_dim()

    ssh_p_dim_state_par = mpi_subd_lat_par*mpi_subd_lon_par
    t_p_dim_state_par = mpi_subd_lat_par*mpi_subd_lon_par*mpi_subd_vert_par
    s_p_dim_state_par = mpi_subd_lat_par*mpi_subd_lon_par*mpi_subd_vert_par
    u_p_dim_state_par = mpi_subd_lat_par*mpi_subd_lon_par*mpi_subd_vert_par
    v_p_dim_state_par = mpi_subd_lat_par*mpi_subd_lon_par*mpi_subd_vert_par
#if defined key_agrif
    ssh_p_dim_state_child = mpi_subd_lat_child*mpi_subd_lon_child
    t_p_dim_state_child = mpi_subd_lat_child*mpi_subd_lon_child*mpi_subd_vert_child
    s_p_dim_state_child = mpi_subd_lat_child*mpi_subd_lon_child*mpi_subd_vert_child
    u_p_dim_state_child = mpi_subd_lat_child*mpi_subd_lon_child*mpi_subd_vert_child
    v_p_dim_state_child = mpi_subd_lat_child*mpi_subd_lon_child*mpi_subd_vert_child
#endif

  END SUBROUTINE calc_dim

  SUBROUTINE calc_offset()

    ! !DESCRIPTION:
    ! This routine calculates the offset values for the state variables.
    ! It then stores the 2d/3d offset values in separate arrays.

    ! !USES:
    IMPLICIT NONE


    ! Compute state vector dimensions in case not already done
    CALL calc_dim()

    ssh_p_offset_par = 0
    t_p_offset_par = ssh_p_offset_par + ssh_p_dim_state_par
    s_p_offset_par = t_p_offset_par + t_p_dim_state_par
    u_p_offset_par = s_p_offset_par + s_p_dim_state_par
    v_p_offset_par = u_p_offset_par + u_p_dim_state_par

    ! Fill array of 2D state variable offsets for local PE
    var2d_p_offset_par(1) = ssh_p_offset_par

    ! Fill array of 3D state variable offsets for local PE
    var3d_p_offset_par(1) = t_p_offset_par
    var3d_p_offset_par(2) = s_p_offset_par
    var3d_p_offset_par(3) = u_p_offset_par
    var3d_p_offset_par(4) = v_p_offset_par

#if defined key_agrif
    ! Continue child variables directly after parent variables
    ssh_p_offset_child = v_p_offset_par + v_p_dim_state_par
    t_p_offset_child = ssh_p_offset_child + ssh_p_dim_state_child
    s_p_offset_child = t_p_offset_child + t_p_dim_state_child
    u_p_offset_child = s_p_offset_child + s_p_dim_state_child
    v_p_offset_child = u_p_offset_child + u_p_dim_state_child

    ! Fill array of 2D state variable offsets for local PE
    var2d_p_offset_child(1) = ssh_p_offset_child

    ! Fill array of 3D state variable offsets for local PE
    var3d_p_offset_child(1) = t_p_offset_child
    var3d_p_offset_child(2) = s_p_offset_child
    var3d_p_offset_child(3) = u_p_offset_child
    var3d_p_offset_child(4) = v_p_offset_child
#endif

  END SUBROUTINE calc_offset

  SUBROUTINE calc_statevector_dim(dim_p, dim_p_par, dim_p_child)

    ! !DESCRIPTION:
    ! This routine calculates the total state vector dimension.

    ! !USES:
    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(inout) :: dim_p       ! Local state vector total dim
    INTEGER, INTENT(inout) :: dim_p_par   ! Local state vector parent grid dim
    INTEGER, INTENT(inout) :: dim_p_child ! Local state vector child grid dim


    ! Calculate state variable dimensions in case not already done
    CALL calc_dim()

#if defined key_agrif
    dim_p = ssh_p_dim_state_par + t_p_dim_state_par + &
         s_p_dim_state_par + u_p_dim_state_par + v_p_dim_state_par + &
         ssh_p_dim_state_child + t_p_dim_state_child +&
         s_p_dim_state_child + u_p_dim_state_child + v_p_dim_state_child
    dim_p_par = ssh_p_dim_state_par + t_p_dim_state_par +&
         s_p_dim_state_par + u_p_dim_state_par + v_p_dim_state_par
    dim_p_child = ssh_p_dim_state_child + t_p_dim_state_child +&
         s_p_dim_state_child + u_p_dim_state_child + v_p_dim_state_child
#else
    dim_p = ssh_p_dim_state_par + t_p_dim_state_par + &
         s_p_dim_state_par + u_p_dim_state_par + v_p_dim_state_par
    dim_p_par = dim_p
    dim_p_child = 0
#endif

  END SUBROUTINE calc_statevector_dim

  SUBROUTINE fill2d_par_ensarray(dim_p, dim_ens, wght, fname, ens_p)

    ! !DESCRIPTION:
    ! Fill local ensemble array with 2d state variables from initial state file.

    ! !USES:
    USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens, nldi_par, nldj_par, &
         nimpp_par, njmpp_par, jpiglo_par,jpjglo_par
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in)   :: dim_p                      ! PE-local state dimension
    INTEGER, INTENT(in)   :: dim_ens                    ! Size of ensemble
    REAL(pwp), INTENT(in)   :: wght(dim_ens)            ! Weights for ensemble
    CHARACTER(len=*), INTENT(in) :: fname               ! Name of NC file
    REAL(pwp), INTENT(inout)   :: ens_p(dim_p, dim_ens) ! PE-local state ensemble

    ! *** local variables ***
    INTEGER :: s, i, j, idx, member     ! Counters
    INTEGER :: ii, jj                   ! Start index for NetCDF file
    INTEGER :: stat(20000)              ! Status flag for NetCDF commands
    INTEGER :: ncid_in                  ! ID for NetCDF file
    INTEGER :: id_2dvar                 ! IDs for fields
    INTEGER :: pos(3),cnt(3)            ! Vectors for 2D reading fields
    REAL(pwp), ALLOCATABLE :: var2d(:,:,:)   ! Array for reading state variables
    REAL(pwp), ALLOCATABLE :: var2d_av(:,:)  ! Ensemble average of var_2d array
    CHARACTER(len=100) :: intervalwrite ! Frequency of snapshots from control


    ! ******************************************
    ! *** Open file containing initial state ***
    ! ******************************************

    s = 1
    stat(s) = NF90_OPEN(fname, NF90_NOWRITE, ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in opening initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

    ! ******************************************
    ! *** Read file containing initial state ***
    ! ******************************************

    ! Calculate offsets in case not already calculated
    CALL calc_offset()

    ALLOCATE (var2d(mpi_subd_lon_par, mpi_subd_lat_par, dim_ens))
    ALLOCATE (var2d_av(mpi_subd_lon_par, mpi_subd_lat_par))

    ! Loop over all state variables using state variable .NC id array
    id2d: DO idx = 1, SIZE(id2d_list)
       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, TRIM(id2d_list(idx)), id_2dvar)

       ! Set the starting position after the halo region
       jj = njmpp_par + nldj_par - 1
       ii = nimpp_par + nldi_par - 1

       ! We increment the time index by member to get a new snapshot for each member.
       pos = (/ ii, jj, 1 /)
       cnt = (/ mpi_subd_lon_par , mpi_subd_lat_par, dim_ens /)

       s = s + 1
       stat(s) = NF90_GET_VAR(ncid_in, id_2dvar, var2d, start=pos, count=cnt)

       ! Output snapshot frequency
       IF (mype_ens == 0) THEN
          s = s + 1
          stat(s) = NF90_INQUIRE_ATTRIBUTE(ncid_in, id_2dvar, "interval_write")
          s = s + 1
          stat(s) = NF90_GET_ATT(ncid_in, id_2dvar, "interval_write", intervalwrite)
          WRITE(*, '(/9x, a, 3x, a)') 'Snapshots currently saved every:',&
               intervalwrite
       END IF

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in reading initial state file, var=', &
                  trim(id2d_list(idx))
             CALL abort_parallel()
          END IF
       END DO

       ! Compute ensemble average
       var2d_av(:,:)=0
       DO member = 1, dim_ens
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                var2d_av(i,j) = var2d_av(i,j) + var2d(i,j,member)
             END DO
          END DO
       END DO
       var2d_av(:,:) = var2d_av(:,:)/dim_ens

       ! Write fields into state vector following method outlined in init_ens_pdaf
       DO member = 1, dim_ens
          DO j = 1, mpi_subd_lat_par
             DO i = 1, mpi_subd_lon_par
                IF (member == 1) THEN  ! Member 1 should not have any perturbation.
                   ens_p(i+(j-1)*mpi_subd_lon_par + var2d_p_offset_par(idx), member) = &
                        var2d(i,j,member)
                ELSE
                   ens_p(i+(j-1)*mpi_subd_lon_par + var2d_p_offset_par(idx), member) = &
                        var2d_av(i,j) + wght(member) * ( var2d(i,j,member) - &
                        var2d_av(i,j) )
                END IF
             END DO
          END DO
       END DO
    END DO id2d

    DEALLOCATE(var2d, var2d_av)

    ! *******************************************
    ! *** Close file containing initial state ***
    ! *******************************************

    s = 1
    stat(s) = NF90_CLOSE(ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in closing initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE fill2d_par_ensarray

  SUBROUTINE fill3d_par_ensarray(dim_p, dim_ens, wght, fname, statevar, ens_p)

    ! !DESCRIPTION:
    ! Fill local ensemble array with 3d state variables from initial state file.

    ! !USES:
    USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens, nldi_par, nldj_par, &
         nimpp_par, njmpp_par
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in)   :: dim_p                     ! PE-local state dimension
    INTEGER, INTENT(in)   :: dim_ens                   ! Size of ensemble
    REAL(pwp), INTENT(in)   :: wght(dim_ens)             ! Weights for ensemble
    CHARACTER(len=*), INTENT(in) :: fname                     ! Name of NC file
    CHARACTER(len=*), INTENT(in) :: statevar           ! Name of state variable
    REAL(pwp), INTENT(inout)   :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble

    ! *** local variables ***
    INTEGER :: s, i, j, k, idx, member   ! Counters
    INTEGER :: lwr_bnd, upr_bnd          ! Index limits for .NC id data
    INTEGER :: ii, jj                    ! Start index for NetCDF file
    INTEGER :: stat(20000)               ! Status flag for NetCDF commands
    INTEGER :: ncid_in                   ! ID for NetCDF file
    INTEGER :: id_3dvar                  ! IDs for fields
    INTEGER :: pos(4),cnt(4)             ! Vectors for 2D reading fields
    REAL(pwp), ALLOCATABLE :: var3d(:,:,:,:)  ! Array for reading state variables
    REAL(pwp), ALLOCATABLE :: var3d_av(:,:,:) ! Ensemble average of var_2d array
    CHARACTER(len=100) :: intervalwrite  ! Frequency of snapshots from control


    ! Determine what .NC ids are needed from id3d_list array
    SELECT CASE ( statevar )
    CASE ('T')
       lwr_bnd = 1
       upr_bnd = 2
    CASE('U')
       lwr_bnd = 3
       upr_bnd = 3
    CASE('V')
       lwr_bnd = 4
       upr_bnd = 4
    END SELECT

    ! ******************************************
    ! *** Open file containing initial state ***
    ! ******************************************

    s = 1
    stat(s) = NF90_OPEN(fname, NF90_NOWRITE, ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in opening initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

    ! ******************************************
    ! *** Read file containing initial state ***
    ! ******************************************

    ! Calculate offsets in case not already calculated
    CALL calc_offset()

    ALLOCATE(var3d(mpi_subd_lon_par, mpi_subd_lat_par, mpi_subd_vert_par, dim_ens))
    ALLOCATE(var3d_av(mpi_subd_lon_par, mpi_subd_lat_par, mpi_subd_vert_par))

    ! Loop over all state variables using state variable .NC id array
    id3d: DO idx = lwr_bnd, upr_bnd
       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, trim(id3d_list(idx)), id_3dvar)

       ! Set the starting position after the halo region
       jj = njmpp_par + nldj_par - 1
       ii = nimpp_par + nldi_par - 1

       ! We increment the time index by member to get a new snapshot for each member.
       pos = (/ ii, jj, 1, 1 /)
       cnt = (/ mpi_subd_lon_par , mpi_subd_lat_par, mpi_subd_vert_par, dim_ens /)
       s = s + 1
       stat(s) = NF90_GET_VAR(ncid_in, id_3dvar, var3d, start=pos, count=cnt)

       ! Output snapshot frequency
       IF (mype_ens == 0) THEN
          s = s + 1
          stat(s) = NF90_INQUIRE_ATTRIBUTE(ncid_in, id_3dvar, "interval_write")
          s = s + 1
          stat(s) = NF90_GET_ATT(ncid_in, id_3dvar, "interval_write", intervalwrite)
          WRITE(*, '(/9x, a, 3x, a)') 'Snapshots currently saved every:',&
               intervalwrite
       END IF

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in reading initial state file, var=', &
                  trim(id3d_list(idx))
             CALL abort_parallel()
          END IF
       END DO

       ! Compute ensemble average
       var3d_av(:,:,:)=0
       DO member = 1, dim_ens
          DO k = 1, mpi_subd_vert_par
             DO j = 1, mpi_subd_lat_par
                DO i = 1, mpi_subd_lon_par
                   var3d_av(i,j,k) = var3d_av(i,j,k) + var3d(i,j,k,member)
                END DO
             END DO
          END DO
       END DO
       var3d_av(:,:,:) = var3d_av(:,:,:)/dim_ens

       ! Write fields into state vector following method outlined in init_ens_pdaf
       DO member = 1, dim_ens
          DO k = 1, mpi_subd_vert_par
             DO j = 1, mpi_subd_lat_par
                DO i = 1, mpi_subd_lon_par
                   IF (member == 1) THEN  ! Member 1 should not have any perturbation.
                      ens_p(i + (j-1)*mpi_subd_lon_par + &
                           (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + &
                           var3d_p_offset_par(idx), member) = var3d(i,j,k,member)
                   ELSE
                      ens_p(i + (j-1)*mpi_subd_lon_par + &
                           (k-1)*mpi_subd_lat_par*mpi_subd_lon_par + &
                           var3d_p_offset_par(idx), member) = var3d_av(i,j,k) + &
                           wght(member) * ( var3d(i,j,k,member) - var3d_av(i,j,k) )
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO id3d

    DEALLOCATE(var3d, var3d_av)

    ! *******************************************
    ! *** Close file containing initial state ***
    ! *******************************************

    s = 1
    stat(s) = NF90_CLOSE(ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in closing initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE fill3d_par_ensarray

  SUBROUTINE fill2d_child_ensarray(dim_p, dim_ens, wght, fname, ens_p)

    ! !DESCRIPTION:
    ! Fill local ensemble array with 2d state variables from initial state file.

    ! !USES:
    USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens, nldi_child, nldj_child, &
         nimpp_child, njmpp_child, jpiglo_child, jpjglo_child
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in)   :: dim_p                      ! PE-local state dimension
    INTEGER, INTENT(in)   :: dim_ens                    ! Size of ensemble
    REAL(pwp), INTENT(in)   :: wght(dim_ens)            ! Weights for ensemble
    CHARACTER(len=*), INTENT(in) :: fname               ! Name of NC file
    REAL(pwp), INTENT(inout)   :: ens_p(dim_p, dim_ens) ! PE-local state ensemble

    ! *** local variables ***
    INTEGER :: s, i, j, idx, member     ! Counters
    INTEGER :: ii, jj                   ! Start index for NetCDF file
    INTEGER :: stat(20000)              ! Status flag for NetCDF commands
    INTEGER :: ncid_in                  ! ID for NetCDF file
    INTEGER :: id_2dvar                 ! IDs for fields
    INTEGER :: pos(3),cnt(3)            ! Vectors for 2D reading fields
    REAL(pwp), ALLOCATABLE :: var2d(:,:,:)   ! Array for reading state variables
    REAL(pwp), ALLOCATABLE :: var2d_av(:,:)  ! Ensemble average of var_2d array
    CHARACTER(len=100) :: intervalwrite ! Frequency of snapshots from control


    ! ******************************************
    ! *** Open file containing initial state ***
    ! ******************************************

    s = 1
    stat(s) = NF90_OPEN(fname, NF90_NOWRITE, ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in opening initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

    ! ******************************************
    ! *** Read file containing initial state ***
    ! ******************************************

    ! Calculate offsets in case not already calculated
    CALL calc_offset()

    ALLOCATE (var2d(mpi_subd_lon_child, mpi_subd_lat_child, dim_ens))
    ALLOCATE (var2d_av(mpi_subd_lon_child, mpi_subd_lat_child))

    ! Loop over all state variables using state variable .NC id array
    id2d: DO idx = 1, SIZE(id2d_list)
       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, TRIM(id2d_list(idx)), id_2dvar)

       ! Set the starting position after the halo region
       jj = njmpp_child + nldj_child - 1
       ii = nimpp_child + nldi_child - 1

       ! We increment the time index by member to get a new snapshot for each member.
       pos = (/ ii, jj, 1 /)
       cnt = (/ mpi_subd_lon_child , mpi_subd_lat_child, dim_ens /)

       s = s + 1
       stat(s) = NF90_GET_VAR(ncid_in, id_2dvar, var2d, start=pos, count=cnt)

       ! Output snapshot frequency
       IF (mype_ens == 0) THEN
          s = s + 1
          stat(s) = NF90_INQUIRE_ATTRIBUTE(ncid_in, id_2dvar, "interval_write")
          s = s + 1
          stat(s) = NF90_GET_ATT(ncid_in, id_2dvar, "interval_write", intervalwrite)
          WRITE(*, '(/9x, a, 3x, a)') 'Snapshots currently saved every:',&
               intervalwrite
       END IF

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in reading initial state file, var=', &
                  trim(id2d_list(idx))
             CALL abort_parallel()
          END IF
       END DO

       ! Compute ensemble average
       var2d_av(:,:)=0
       DO member = 1, dim_ens
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                var2d_av(i,j) = var2d_av(i,j) + var2d(i,j,member)
             END DO
          END DO
       END DO
       var2d_av(:,:) = var2d_av(:,:)/dim_ens

       ! Write fields into state vector following method outlined in init_ens_pdaf
       DO member = 1, dim_ens
          DO j = 1, mpi_subd_lat_child
             DO i = 1, mpi_subd_lon_child
                IF (member == 1) THEN  ! Member 1 should not have any perturbation.
                   ens_p(i+(j-1)*mpi_subd_lon_child + var2d_p_offset_child(idx), member) = &
                        var2d(i,j,member)
                ELSE
                   ens_p(i+(j-1)*mpi_subd_lon_child + var2d_p_offset_child(idx), member) = &
                        var2d_av(i,j) + wght(member) * ( var2d(i,j,member) - &
                        var2d_av(i,j) )
                END IF
             END DO
          END DO
       END DO
    END DO id2d

    DEALLOCATE(var2d, var2d_av)

    ! *******************************************
    ! *** Close file containing initial state ***
    ! *******************************************

    s = 1
    stat(s) = NF90_CLOSE(ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in closing initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE fill2d_child_ensarray

  SUBROUTINE fill3d_child_ensarray(dim_p, dim_ens, wght, fname, statevar, ens_p)

    ! !DESCRIPTION:
    ! Fill local ensemble array with 3d state variables from initial state file.

    ! !USES:
    USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens, nldi_child, nldj_child, &
         nimpp_child, njmpp_child
    USE netcdf

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in)   :: dim_p                     ! PE-local state dimension
    INTEGER, INTENT(in)   :: dim_ens                   ! Size of ensemble
    REAL(pwp), INTENT(in)   :: wght(dim_ens)             ! Weights for ensemble
    CHARACTER(len=*), INTENT(in) :: fname                     ! Name of NC file
    CHARACTER(len=*), INTENT(in) :: statevar           ! Name of state variable
    REAL(pwp), INTENT(inout)   :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble

    ! *** local variables ***
    INTEGER :: s, i, j, k, idx, member   ! Counters
    INTEGER :: lwr_bnd, upr_bnd          ! Index limits for .NC id data
    INTEGER :: ii, jj                    ! Start index for NetCDF file
    INTEGER :: stat(20000)               ! Status flag for NetCDF commands
    INTEGER :: ncid_in                   ! ID for NetCDF file
    INTEGER :: id_3dvar                  ! IDs for fields
    INTEGER :: pos(4),cnt(4)             ! Vectors for 2D reading fields
    REAL(pwp), ALLOCATABLE :: var3d(:,:,:,:)  ! Array for reading state variables
    REAL(pwp), ALLOCATABLE :: var3d_av(:,:,:) ! Ensemble average of var_2d array
    CHARACTER(len=100) :: intervalwrite  ! Frequency of snapshots from control


    ! Determine what .NC ids are needed from id3d_list array
    SELECT CASE ( statevar )
    CASE ('T')
       lwr_bnd = 1
       upr_bnd = 2
    CASE('U')
       lwr_bnd = 3
       upr_bnd = 3
    CASE('V')
       lwr_bnd = 4
       upr_bnd = 4
    END SELECT

    ! ******************************************
    ! *** Open file containing initial state ***
    ! ******************************************

    s = 1
    stat(s) = NF90_OPEN(fname, NF90_NOWRITE, ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in opening initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

    ! ******************************************
    ! *** Read file containing initial state ***
    ! ******************************************

    ! Calculate offsets in case not already calculated
    CALL calc_offset()

    ALLOCATE(var3d(mpi_subd_lon_child, mpi_subd_lat_child, mpi_subd_vert_child, dim_ens))
    ALLOCATE(var3d_av(mpi_subd_lon_child, mpi_subd_lat_child, mpi_subd_vert_child))

    ! Loop over all state variables using state variable .NC id array
    id3d: DO idx = lwr_bnd, upr_bnd
       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, trim(id3d_list(idx)), id_3dvar)

       ! Set the starting position after the halo region
       jj = njmpp_child + nldj_child - 1
       ii = nimpp_child + nldi_child - 1

       ! We increment the time index by member to get a new snapshot for each member.
       pos = (/ ii, jj, 1, 1 /)
       cnt = (/ mpi_subd_lon_child , mpi_subd_lat_child, mpi_subd_vert_child, dim_ens /)
       s = s + 1
       stat(s) = NF90_GET_VAR(ncid_in, id_3dvar, var3d, start=pos, count=cnt)

       ! Output snapshot frequency
       IF (mype_ens == 0) THEN
          s = s + 1
          stat(s) = NF90_INQUIRE_ATTRIBUTE(ncid_in, id_3dvar, "interval_write")
          s = s + 1
          stat(s) = NF90_GET_ATT(ncid_in, id_3dvar, "interval_write", intervalwrite)
          WRITE(*, '(/9x, a, 3x, a)') 'Snapshots currently saved every:',&
               intervalwrite
       END IF

       DO i = 1, s
          IF (stat(i) .NE. NF90_NOERR) THEN
             WRITE(*,'(/9x, a, 3x, a)') &
                  'NetCDF error in reading initial state file, var=', &
                  trim(id3d_list(idx))
             CALL abort_parallel()
          END IF
       END DO

       ! Compute ensemble average
       var3d_av(:,:,:)=0
       DO member = 1, dim_ens
          DO k = 1, mpi_subd_vert_child
             DO j = 1, mpi_subd_lat_child
                DO i = 1, mpi_subd_lon_child
                   var3d_av(i,j,k) = var3d_av(i,j,k) + var3d(i,j,k,member)
                END DO
             END DO
          END DO
       END DO
       var3d_av(:,:,:) = var3d_av(:,:,:)/dim_ens

       ! Write fields into state vector following method outlined in init_ens_pdaf
       DO member = 1, dim_ens
          DO k = 1, mpi_subd_vert_child
             DO j = 1, mpi_subd_lat_child
                DO i = 1, mpi_subd_lon_child
                   IF (member == 1) THEN  ! Member 1 should not have any perturbation.
                      ens_p(i + (j-1)*mpi_subd_lon_child + &
                           (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + &
                           var3d_p_offset_child(idx), member) = var3d(i,j,k,member)
                   ELSE
                      ens_p(i + (j-1)*mpi_subd_lon_child + &
                           (k-1)*mpi_subd_lat_child*mpi_subd_lon_child + &
                           var3d_p_offset_child(idx), member) = var3d_av(i,j,k) + &
                           wght(member) * ( var3d(i,j,k,member) - var3d_av(i,j,k) )
                   END IF
                END DO
             END DO
          END DO
       END DO
    END DO id3d

    DEALLOCATE(var3d, var3d_av)

    ! *******************************************
    ! *** Close file containing initial state ***
    ! *******************************************

    s = 1
    stat(s) = NF90_CLOSE(ncid_in)

    DO i = 1, s
       IF (stat(i) .NE. NF90_NOERR) THEN
          WRITE(*,'(/9x, a, 3x, a)') &
               'NetCDF error in closing initial state file:', fname
          CALL abort_parallel()
       END IF
    END DO

  END SUBROUTINE fill3d_child_ensarray

!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_statevector_pdaf
