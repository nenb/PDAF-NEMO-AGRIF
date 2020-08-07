!$Id: l2g_state_pdaf.F90 332 2019-12-30 09:37:03Z lnerger $
!>  Initialize full state from local analysis
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during the loop over all
!! local analysis domains in PDAF_X_update 
!! after the analysis and ensemble transformation 
!! on a single local analysis domain. It has to 
!! initialize elements of the PE-local full state 
!! vector from the provided analysis state vector 
!! on the local analysis domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! __Revision history:__
!! * 2013-02 - Lars Nerger - Initial code
!! * Later revisions - see repository log
!!
SUBROUTINE l2g_state_pdaf(step, domain_p, dim_l, state_l, dim_p, state_p)
!$AGRIF_DO_NOT_TREAT

  USE mod_kind_pdaf
  USE mod_assimilation_pdaf, &
       ONLY: coords_l, indx_dom_l_par, indx_dom_l_child
  USE mod_statevector_pdaf, &
       ONLY: mpi_subd_vert_child, mpi_subd_vert_par, mpi_subd_lon_child, &
       mpi_subd_lon_par, var2d_p_offset_par, var2d_p_offset_child, &
       var3d_p_offset_par, var3d_p_offset_child, mpi_subd_lat_par, &
       mpi_subd_lat_child
  USE mod_parallel_pdaf, &
       ONLY: nldj_child, nldj_par, nldi_child, nldi_par
  USE mod_agrif_pdaf, &
       ONLY: tmask_par, tmask_child

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in) :: step           !< Current time step
  INTEGER, INTENT(in) :: domain_p       !< Current local analysis domain
  INTEGER, INTENT(in) :: dim_l          !< Local state dimension
  INTEGER, INTENT(in) :: dim_p          !< PE-local full state dimension
  REAL(pwp), INTENT(in)    :: state_l(dim_l) !< State vector on local analysis domain
  REAL(pwp), INTENT(inout) :: state_p(dim_p) !< PE-local full state vector 

    ! *** local variables ***
  INTEGER :: domain_p_child        ! Counter for local analysis domain on child grid
  INTEGER :: i, cnt                ! Counter
  INTEGER :: a, b, c               ! Variables for 3D state variable index
  INTEGER :: tot_dom_l_par         ! Number of local domains on parent grid
  INTEGER :: i_par, j_par          ! Grid coordinates for local analysis domain
  INTEGER :: i_child, j_child      ! Grid coordinates for local analysis domain
  INTEGER :: i0, j0                ! Halo offset for local PE
  INTEGER :: loc_2d             ! 2D coordinate in statevector
  INTEGER :: dim_vert              ! Vertical dimension of local state vector
  REAL(pwp), ALLOCATABLE :: v_coord(:) ! Array for converting vertical coordinate in
                                       ! local state vector to parent grid.


! *************************************
! *** Initialize local state vector ***
! *************************************

#ifndef key_agrif

  ! *****************************************************************
  ! Statevector counts ocean *and* land points, whereas local domains
  ! only count ocean points. Hence we need to introduce a conversion
  ! between statevector and local domains.
  ! *****************************************************************

  ! Compute i,j indices for parent grid.
  i_par = indx_dom_l_par(1, domain_p)
  j_par = indx_dom_l_par(2, domain_p)

  ! Compute horizontal in statevector.
  loc_2d = (j_par - 1)*(mpi_subd_lon_par) + i_par

  ! Compute vertical dimension of local state vector.
  dim_vert = ( dim_l - SIZE(var2d_p_offset_par) )/( SIZE(var3d_p_offset_par) )

  ! Compute halo offset for parent grid.
  i0 = nldi_par - 1
  j0 = nldj_par - 1

  ! Compute array for converting vertical coordinate in local state
  ! vector to vertical coordinate on parent grid. 
  ALLOCATE(v_coord(dim_vert))
  cnt=0
  DO i = 1, mpi_subd_vert_par
     IF(tmask_par(i_par+i0, j_par+j0, i) == 1) THEN
        cnt = cnt + 1
        v_coord(cnt) = i
     END IF
  END DO

  ! The local state vector elements are the 2D state variables as well as
  ! the values at different depths for the 3D state variables.
  DO i = 1, dim_l
     IF( i <= SIZE(var2d_p_offset_par) ) THEN
        state_p( var2d_p_offset_par(i) + loc_2d ) = state_l(i)
     ELSE
        ! Now, fill the local state vector with 3D state variables.
        ! Each 3D state variable has dimension nx*ny*nz. We first break
        ! this into chunks of size nx*ny. Then we identify each chunk as a
        ! different element of the local state vector.
        !
        ! Counter for 3D state variables in local state vector.
        a = i - SIZE(var2d_p_offset_par)
        ! Counter for 3D state variable.
        b = INT( CEILING( REAL(a)/REAL(dim_vert) ) )
        ! Counter for vertical coordinate.
        c = MOD(a-1,dim_vert)
        state_p( var3d_p_offset_par(b) + &
             (mpi_subd_lon_par*mpi_subd_lat_par*(v_coord(c)-1)) + loc_2d ) = &
             state_l(i)
     END IF
  END DO

  !Clean-up
  DEALLOCATE(v_coord)

#else

  WRITE(*,*) 'TEST g2l ERROR, need to remove NEMO grid here'
  CALL abort_parallel()

  ! Child grid
  child_grd: IF(domain_p > tot_dom_l_par) THEN
     ! Convert to local domain value for child grid.
     domain_p_child = domain_p - tot_dom_l_par

     ! Compute i,j indices for child grid.
     i_child = indx_dom_l_child(1, domain_p_child)
     j_child = indx_dom_l_child(2, domain_p_child)

     ! Compute horizontal location of i,j indices in statevector.
     loc_2d = (j_child - 1)*(mpi_subd_lon_child) + i_child

     ! Compute vertical dimension of local state vector.
     dim_vert = ( dim_l - SIZE(var2d_p_offset_child) )/( SIZE(var3d_p_offset_child) )

     ! Compute halo offset for child grid.
     i0 = nldi_child - 1
     j0 = nldj_child - 1

     ! Compute array for converting vertical coordinate in local state
     ! vector to vertical coordinate on parent grid.
     ALLOCATE(v_coord(dim_vert))
     cnt=0
     DO i = 1, mpi_subd_vert_child
        IF(tmask_child(i_child+i0, j_child+j0, i) == 1) THEN
           cnt = cnt + 1
           v_coord(cnt) = i
        END IF
     END DO

     ! The local state vector elements are the 2D state variables as well as
     ! the values at different depths for the 3D state variables.
     DO i = 1, dim_l
        IF( i <= SIZE(var2d_p_offset_child) ) THEN
           state_p( var2d_p_offset_child(i) + loc_2d ) = state_l(i)
        ELSE
           ! Now, fill the local state vector with 3D state variables.
           ! Each 3D state variable has dimension nx*ny*nz. We first break
           ! this into chunks of size nx*ny. Then we identify each chunk as a
           ! different element of the local state vector.
           !
           ! Counter for 3D state variables in local state vector.
           a = i - SIZE(var2d_p_offset_child)
           ! Counter for 3D state variable.
           b = INT( CEILING( REAL(a)/REAL(dim_vert) ) )
           ! Counter for vertical coordinate.
           c = MOD(a-1,dim_vert)
           state_p( var3d_p_offset_child(b) + &
                (mpi_subd_lon_child*mpi_subd_lat_child*(v_coord(c)-1)) + &
                loc_2d ) =  state_l(i)
        END IF
     END DO

     !Clean-up
     DEALLOCATE(v_coord)
  END IF child_grd
#endif

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE l2g_state_pdaf
