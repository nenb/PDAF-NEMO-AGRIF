MODULE mod_statevector

! !DESCRIPTION:
! This module provides variables & routines for
! manipulating the state vector.

! !USES:
  IMPLICIT NONE
  SAVE


! 2d state vector variables - start index
  INTEGER :: ssh_p_offset

! Array holding 2d state variable offsets
  INTEGER :: var2d_p_offset(1)

! 2d state vector variables - dimension size
  INTEGER :: ssh_p_dim_state

! 3d state vector variables - start index
  INTEGER :: t_p_offset
  INTEGER :: s_p_offset
  INTEGER :: u_p_offset
  INTEGER :: v_p_offset
  INTEGER :: w_p_offset

! Array holding 3d state variable offsets
  INTEGER :: var3d_p_offset(5)

! 3d state vector variables - dimension size
  INTEGER :: t_p_dim_state
  INTEGER :: s_p_dim_state
  INTEGER :: u_p_dim_state
  INTEGER :: v_p_dim_state
  INTEGER :: w_p_dim_state

! Dimensions for MPI subdomain
  INTEGER :: mpi_subd_lat
  INTEGER :: mpi_subd_lon
  INTEGER :: mpi_subd_vert

  ! Array for 2d/3d state variable .NC ids
  CHARACTER(len=20), DIMENSION(1) :: id2d_list
  CHARACTER(len=20), DIMENSION(5) :: id3d_list

! Fill array of 2d/3d state variable .NC ids
  DATA id2d_list / 'sossheig' /
  DATA id3d_list / 'votemper', 'vosaline', 'vozocrtx', &
       'vomecrty', 'vovecrtz' /

CONTAINS

  SUBROUTINE calc_mpi_dim()

    ! !DESCRIPTION:
    ! This routine calculates the dimensions of the subdomain of the MPI
    ! domain that is used to fill the statevector. There are two options:
    ! 1) PDAF_optim
    ! Remove  halo regions from the state vector.
    ! 2) *ELSE*
    ! Include halo regions in the state vector.

    ! !USES:
    USE par_oce, ONLY: jpk, jpi, jpj, jpiglo ,jpjglo
    USE dom_oce, ONLY: nldi, nldj, nlei, nlej, nimpp, njmpp

    IMPLICIT NONE


#if defined PDAF_optim
    mpi_subd_lon = nlei-nldi+1
    mpi_subd_lat = nlej-nldj+1
    mpi_subd_vert = jpk
#else
    mpi_subd_lon= jpi
    mpi_subd_lat= jpj
    mpi_subd_vert= jpk

    ! Treat case where MPI subdomain goes outside limits of
    ! global domain. In this case, MPI subdomain is likely
    ! padded with zeros, which we now ignore.
    IF (nimpp + jpi .GT. jpiglo+1) THEN
       mpi_subd_lon = jpiglo+1-nimpp
    END IF
    IF (njmpp + jpj .GT. jpjglo+1) THEN
       mpi_subd_lat = jpjglo+1-njmpp
    END IF
#endif

  END SUBROUTINE calc_mpi_dim

  SUBROUTINE calc_2d_dim()

    ! !DESCRIPTION:
    ! This routine calculates the dimension of the 2d state variables.
    ! There are two options:
    ! 1) PDAF_optim
    ! Land points are excluded from the state vector.
    ! 2) *ELSE*
    ! Land points are included in the state vector.

    ! !USES:
    USE dom_oce, ONLY: ssmask, nldi, nldj

    IMPLICIT NONE

    ! *** local variables ***
    INTEGER :: cnt   ! Counter


    ! Compute MPI subdomain dimensions in case not already done
    CALL calc_mpi_dim()

#if defined PDAF_optim
    cnt=0
    DO j = 1, mpi_subd_lat
       DO i = 1, mpi_subd_lon
          IF (ss_mask(nldi+i-1,nldj+j-1) .EQ. 1) THEN
             cnt=cnt+1
          END IF
       END DO
    END DO
    ssh_p_dim_state=cnt
#else
    ssh_p_dim_state = mpi_subd_lat*mpi_subd_lon
#endif

  END SUBROUTINE calc_2d_dim

  SUBROUTINE calc_3d_dim()

        ! !DESCRIPTION:
    ! This routine calculates the dimension of the 2d state variables.
    ! There are two options:
    ! 1) PDAF_optim
    ! Land points are excluded from the state vector.
    ! 2) *ELSE*
    ! Land points are included in the state vector.

    ! !USES:
    USE dom_oce, ONLY: tmask, umask, vmask, wmask, nldi, nldj

    IMPLICIT NONE

    ! *** local variables ***
    INTEGER :: cnt    ! Counter


    ! Compute MPI subdomain dimensions in case not already done
    CALL calc_mpi_dim()

#if defined PDAF_optim
    cnt=0
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             IF (tmask(nldi+i-1,nldj+j-1,k) .EQ. 1) THEN
                cnt=cnt+1
             END IF
          END DO
       END DO
    END DO
    t_p_dim_state=cnt
    s_p_dim_state=cnt

    cnt=0
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             IF (umask(nldi+i-1,nldj+j-1,k) .EQ. 1) THEN
                cnt=cnt+1
             END IF
          END DO
       END DO
    END DO
    u_p_dim_state=cnt

    cnt=0
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             IF (vmask(nldi+i-1,nldj+j-1,k) .EQ. 1) THEN
                cnt=cnt+1
             END IF
          END DO
       END DO
    END DO
    v_p_dim_state=cnt

    cnt=0
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             IF (wmask(nldi+i-1,nldj+j-1,k) .EQ. 1) THEN
                cnt=cnt+1
             END IF
          END DO
       END DO
    END DO
    w_p_dim_state=cnt
#else
    t_p_dim_state = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
    s_p_dim_state = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
    u_p_dim_state = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
    v_p_dim_state = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
    w_p_dim_state = mpi_subd_lat*mpi_subd_lon*mpi_subd_vert
#endif

  END SUBROUTINE calc_3d_dim

  SUBROUTINE calc_2d_offset()

    ! !DESCRIPTION:
    ! This routine calculates the offset values for the 2d state variables.
    ! It then stores these offset values in a 1D array.

    ! !USES:

    IMPLICIT NONE


    ! Compute state vector dimensions in case not already done
    CALL calc_2d_dim()

    ssh_p_offset = 0

    ! Fill array of state variable offsets for local processor element
    var2d_p_offset(1) = ssh_p_offset

  END SUBROUTINE calc_2d_offset


  SUBROUTINE calc_3d_offset()

    ! !DESCRIPTION:
    ! This routine calculates the offset values for the 3d state variables.
    ! It then stores these offset values in a 1D array.

    ! !USES:

    IMPLICIT NONE


    ! Calculate 2D offsets in case not already done
    CALL calc_2d_offset()

    ! Compute state vector dimensions in case not already done
    CALL calc_3d_dim()

    t_p_offset = ssh_p_offset + ssh_p_dim_state ! Continue from ssh field
    s_p_offset = t_p_offset + t_p_dim_state
    u_p_offset = s_p_offset + s_p_dim_state
    v_p_offset = u_p_offset + s_p_dim_state
    w_p_offset = v_p_offset + v_p_dim_state

    ! Fill array of state variable offsets for local processor element
    var3d_p_offset(1) = t_p_offset
    var3d_p_offset(2) = s_p_offset
    var3d_p_offset(3) = u_p_offset
    var3d_p_offset(4) = v_p_offset
    var3d_p_offset(5) = w_p_offset

  END SUBROUTINE calc_3d_offset


  SUBROUTINE calc_statevector_dim(dim_p)

    ! !DESCRIPTION:
    ! This routine calculates the total state vector dimension.

    ! !USES:
    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(inout) :: dim_p ! Dimension of local state vector

    ! Calculate state variable dimensions in case not already done
    CALL calc_2d_dim()
    CALL calc_3d_dim()

    dim_p = ssh_p_dim_state + t_p_dim_state +&
         s_p_dim_state + u_p_dim_state + v_p_dim_state + w_p_dim_state

  END SUBROUTINE calc_statevector_dim

  SUBROUTINE fill2d_ensarray(dim_p, dim_ens, wght, fname, ens_p)

    ! !DESCRIPTION:
    ! Fill local ensemble array with 2d state variables from initial state file.

    ! !USES:
    USE netcdf
    USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens
    USE dom_oce, ONLY: nldi, nldj, nimpp, njmpp
    USE par_oce, ONLY: jpiglo,jpjglo

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in)   :: dim_p                     ! PE-local state dimension
    INTEGER, INTENT(in)   :: dim_ens                   ! Size of ensemble
    REAL, INTENT(in)   :: wght(dim_ens)             ! Weights for ensemble
    CHARACTER(len=*), INTENT(in) :: fname                     ! Name of NC file
    REAL, INTENT(inout)   :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble

    ! *** local variables ***
    INTEGER :: s, i, j, idx, member     ! Counters
    INTEGER :: ii, jj                   ! Start index for NetCDF file
    INTEGER :: stat(20000)              ! Status flag for NetCDF commands
    INTEGER :: ncid_in                  ! ID for NetCDF file
    INTEGER :: id_2dvar                 ! IDs for fields
    INTEGER :: pos(3),cnt(3)            ! Vectors for 2D reading fields
    REAL, ALLOCATABLE :: var2d(:,:,:)   ! Array for reading state variables
    REAL, ALLOCATABLE :: var2d_av(:,:)  ! Ensemble average of var_2d array
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
    CALL calc_2d_offset()

    ALLOCATE (var2d(mpi_subd_lon, mpi_subd_lat, dim_ens))
    ALLOCATE (var2d_av(mpi_subd_lon, mpi_subd_lat))

    ! Loop over all state variables using state variable .NC id array
    id2d: DO idx = 1, SIZE(id2d_list)
       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, TRIM(id2d_list(idx)), id_2dvar)
#if defined PDAF_optim
       jj = njmpp + nldj - 1
       ii = nimpp + nldi - 1
#else
       jj = njmpp
       ii = nimpp
#endif

       ! We increment the time index by member to get a new snapshot for each member.
       pos = (/ ii, jj, 1 /)
       cnt = (/ mpi_subd_lon , mpi_subd_lat, dim_ens /)

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
          DO j = 1, mpi_subd_lat
             DO i = 1, mpi_subd_lon
                var2d_av(i,j) = var2d_av(i,j) + var2d(i,j,member)
             END DO
          END DO
       END DO
       var2d_av(:,:) = var2d_av(:,:)/dim_ens

       ! Write fields into state vector following method outlined in init_ens_pdaf
       DO member = 1, dim_ens
          DO j = 1, mpi_subd_lat
             DO i = 1, mpi_subd_lon
                IF (member == 1) THEN  ! Member 1 should not have any perturbation.
                   ens_p(i+(j-1)*mpi_subd_lon + var2d_p_offset(idx), member) = &
                        var2d(i,j,member)
                ELSE
                   ens_p(i+(j-1)*mpi_subd_lon + var2d_p_offset(idx), member) = &
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

  END SUBROUTINE fill2d_ensarray

  SUBROUTINE fill2d_statevector(dim_p, state_p)

    ! !DESCRIPTION:
    ! Fill state vector with 2d state variables.

    ! !USES:
    USE dom_oce, ONLY: nldi, nldj
    USE oce, ONLY: sshb

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

    ! *** local variables ***
    INTEGER :: i, j        ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


    ! Calculate offsets in case not already calculated
    CALL calc_2d_offset()

#if defined PDAF_optim
    j0 = nldj - 1
    i0 = nldi - 1
#else
    j0 = 0
    i0 = 0
#endif

    ! Fill state vector with 2d variables
    ! SSH
    DO j = 1, mpi_subd_lat
       DO i = 1, mpi_subd_lon
          state_p(i+(j-1)*mpi_subd_lon + ssh_p_offset) = sshb(i+i0,j+j0)
       END DO
    END DO

  END SUBROUTINE fill2d_statevector

  SUBROUTINE distrib2d_statevector(dim_p, state_p)

    ! !DESCRIPTION:
    ! Distribute state vector 2d state variables.

    ! !USES:
    USE dom_oce, ONLY: nldi, nldj
    USE oce, ONLY: sshb

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

    ! *** local variables ***
    INTEGER :: i, j        ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


    ! Calculate offsets in case not already calculated
    CALL calc_2d_offset()

#if defined PDAF_optim
    j0 = nldj - 1
    i0 = nldi - 1
#else
    j0 = 0
    i0 = 0
#endif

    ! Fill state vector with 2d variables
    ! SSH
    DO j = 1, mpi_subd_lat
       DO i = 1, mpi_subd_lon
          sshb(i+i0,j+j0) = state_p(i+(j-1)*mpi_subd_lon + ssh_p_offset)
       END DO
    END DO

    ! Fill halo regions
#if defined PDAF_optim
    CALL lbc_lnk(sshb, 'T', 1.)
#endif

  END SUBROUTINE distrib2d_statevector

  SUBROUTINE fill3d_ensarray(dim_p, dim_ens, wght, fname, statevar, ens_p)

    ! !DESCRIPTION:
    ! Fill local ensemble array with 3d state variables from initial state file.

    ! !USES:
    USE netcdf
    USE mod_parallel_pdaf, ONLY: abort_parallel, mype_ens
    USE dom_oce, ONLY: nldi, nldj, nimpp, njmpp

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in)   :: dim_p                     ! PE-local state dimension
    INTEGER, INTENT(in)   :: dim_ens                   ! Size of ensemble
    REAL, INTENT(in)   :: wght(dim_ens)             ! Weights for ensemble
    CHARACTER(len=*), INTENT(in) :: fname                     ! Name of NC file
    CHARACTER(len=*), INTENT(in) :: statevar           ! Name of state variable
    REAL, INTENT(inout)   :: ens_p(dim_p, dim_ens)     ! PE-local state ensemble

    ! *** local variables ***
    INTEGER :: s, i, j, k, idx, member   ! Counters
    INTEGER :: lwr_bnd, upr_bnd          ! Index limits for .NC id data
    INTEGER :: ii, jj                    ! Start index for NetCDF file
    INTEGER :: stat(20000)               ! Status flag for NetCDF commands
    INTEGER :: ncid_in                   ! ID for NetCDF file
    INTEGER :: id_3dvar                  ! IDs for fields
    INTEGER :: pos(4),cnt(4)             ! Vectors for 2D reading fields
    REAL, ALLOCATABLE :: var3d(:,:,:,:)  ! Array for reading state variables
    REAL, ALLOCATABLE :: var3d_av(:,:,:) ! Ensemble average of var_2d array
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
    CASE('W')
       lwr_bnd = 5
       upr_bnd = 5
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
    CALL calc_3d_offset()

    ALLOCATE (var3d(mpi_subd_lon, mpi_subd_lat, mpi_subd_vert, dim_ens))
    ALLOCATE (var3d_av(mpi_subd_lon, mpi_subd_lat, mpi_subd_vert))

    ! Loop over all state variables using state variable .NC id array
    id3d: DO idx = lwr_bnd, upr_bnd
       s=1
       stat(s) = NF90_INQ_VARID(ncid_in, trim(id3d_list(idx)), id_3dvar)
#if defined PDAF_optim
       jj = njmpp + nldj - 1
       ii = nimpp + nldi - 1
#else
       jj = njmpp
       ii = nimpp
#endif
       ! We increment the time index by member to get a new snapshot for each member.
       pos = (/ ii, jj, 1, 1 /)
       cnt = (/ mpi_subd_lon , mpi_subd_lat, mpi_subd_vert, dim_ens /)
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
          DO k = 1, mpi_subd_vert
             DO j = 1, mpi_subd_lat
                DO i = 1, mpi_subd_lon
                   var3d_av(i,j,k) = var3d_av(i,j,k) + var3d(i,j,k,member)
                END DO
             END DO
          END DO
       END DO
       var3d_av(:,:,:) = var3d_av(:,:,:)/dim_ens

       ! Write fields into state vector following method outlined in init_ens_pdaf
       DO member = 1, dim_ens
          DO k = 1, mpi_subd_vert
             DO j = 1, mpi_subd_lat
                DO i = 1, mpi_subd_lon
                   IF (member == 1) THEN  ! Member 1 should not have any perturbation.
                      ens_p(i + (j-1)*mpi_subd_lon + (k-1)*mpi_subd_lat*mpi_subd_lon + &
                           var3d_p_offset(idx), member) = var3d(i,j,k,member)
                   ELSE
                      ens_p(i + (j-1)*mpi_subd_lon + (k-1)*mpi_subd_lat*mpi_subd_lon + &
                           var3d_p_offset(idx), member) = var3d_av(i,j,k) + &
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

  END SUBROUTINE fill3d_ensarray

  SUBROUTINE fill3d_statevector(dim_p, state_p)

    ! !DESCRIPTION:
    ! Fill state vector with 2d state variables.

    ! !USES:
    USE dom_oce, ONLY: nldi, nldj
    USE oce, ONLY: tsb, ub, vb, wn
    USE par_oce, ONLY: jp_tem, jp_sal

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

    ! *** local variables ***
    INTEGER :: i, j, k     ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


    ! Calculate offsets in case not already calculated
    CALL calc_3d_offset()

#if defined PDAF_optim
    j0 = nldj - 1
    i0 = nldi - 1
#else
    j0 = 0
    i0 = 0
#endif
    ! Fill state vector with 3d variables
    ! T
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             state_p(i+(j-1)*mpi_subd_lon + (k-1)*mpi_subd_lat*mpi_subd_lon + &
                  t_p_offset) =  tsb(i+i0,j+j0,k,jp_tem)
          END DO
       END DO
    END DO
    ! S
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             state_p(i+(j-1)*mpi_subd_lon + (k-1)*mpi_subd_lat*mpi_subd_lon + &
                  s_p_offset) = tsb(i+i0,j+j0,k,jp_sal)
          END DO
       END DO
    END DO
    ! U
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             state_p(i+(j-1)*mpi_subd_lon + (k-1)*mpi_subd_lat*mpi_subd_lon + &
                  u_p_offset) = ub(i+i0,j+j0,k)
          END DO
       END DO
    END DO
    ! V
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             state_p(i+(j-1)*mpi_subd_lon + (k-1)*mpi_subd_lat*mpi_subd_lon + &
                  v_p_offset) = vb(i+i0,j+j0,k)
          END DO
       END DO
    END DO
    ! W
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             state_p(i+(j-1)*mpi_subd_lon + (k-1)*mpi_subd_lat*mpi_subd_lon + &
                  w_p_offset)  = wn(i+i0,j+j0,k)
          END DO
       END DO
    END DO

  END SUBROUTINE fill3d_statevector

  SUBROUTINE distrib3d_statevector(dim_p, state_p)

    ! !DESCRIPTION:
    ! Distribute state vector 2d state variables.

    ! !USES:
    USE dom_oce, ONLY: nldi, nldj
    USE oce, ONLY: tsb, ub, vb, wn
    USE par_oce, ONLY: jp_tem, jp_sal

    IMPLICIT NONE

    ! !ARGUMENTS
    INTEGER, INTENT(in) :: dim_p                   ! PE-local state dimension
    REAL, INTENT(inout) :: state_p(dim_p)          ! PE-local model state

    ! *** local variables ***
    INTEGER :: i, j, k     ! Counters
    INTEGER :: i0, j0      ! Start index for MPI subdomain


    ! Calculate offsets in case not already calculated
    CALL calc_3d_offset()

#if defined PDAF_optim
    j0 = nldj - 1
    i0 = nldi - 1
#else
    j0 = 0
    i0 = 0
#endif
    ! Fill state vector with 3d variables
    ! T
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             tsb(i+i0,j+j0,k,jp_tem) = &
                  state_p(i+(j-1)*mpi_subd_lon + &
                  (k-1)*mpi_subd_lat*mpi_subd_lon + t_p_offset)
          END DO
       END DO
    END DO
    ! S
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             tsb(i+i0,j+j0,k,jp_sal) = &
                  state_p(i+(j-1)*mpi_subd_lon + &
                  (k-1)*mpi_subd_lat*mpi_subd_lon + s_p_offset)
          END DO
       END DO
    END DO
    ! U
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             ub(i+i0,j+j0,k) = &
                  state_p(i+(j-1)*mpi_subd_lon + &
                  (k-1)*mpi_subd_lat*mpi_subd_lon + u_p_offset)
          END DO
       END DO
    END DO
    ! V
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             vb(i+i0,j+j0,k) = &
                  state_p(i+(j-1)*mpi_subd_lon + &
                  (k-1)*mpi_subd_lat*mpi_subd_lon + v_p_offset)
          END DO
       END DO
    END DO
    ! W
    DO k = 1, mpi_subd_vert
       DO j = 1, mpi_subd_lat
          DO i = 1, mpi_subd_lon
             wn(i+i0,j+j0,k) = &
                  state_p(i+(j-1)*mpi_subd_lon + &
                  (k-1)*mpi_subd_lat*mpi_subd_lon + w_p_offset)
          END DO
       END DO
    END DO

    ! Fill halo regions
#if defined PDAF_optim
    CALL lbc_lnk(tsb, 'T', 1.)
    CALL lbc_lnk(tsb, 'T', 1.)
    CALL lbc_lnk(ub, 'U', -1.)
    CALL lbc_lnk(vb, 'V', -1.)
    CALL lbc_lnk(wn, 'W', 1.)
#endif

  END SUBROUTINE distrib3d_statevector

END MODULE mod_statevector
