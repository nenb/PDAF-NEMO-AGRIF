!$Id: init_dim_l_pdaf.F90 344 2020-01-21 14:47:57Z lnerger $
!>  Set dimension of local model state
!!
!! User-supplied call-back routine for PDAF.
!!
!! Used in the filters: LSEIK/LETKF/LESTKF/LNETF
!!
!! The routine is called during analysis step
!! in PDAF_X_update in the loop over all local
!! analysis domains. It has to set the dimension
!! of the local model  state on the current analysis
!! domain.
!!
!! Implementation for the 2D online example
!! with or without parallelization.
!!
!! 2013-02 - Lars Nerger - Initial code
!! Later revisions - see repository log
!!
SUBROUTINE init_dim_l_pdaf(step, domain_p, dim_l)
!$AGRIF_DO_NOT_TREAT
  USE mod_assimilation_pdaf, &      ! Coordinates of local analysis domain
       ONLY: coords_l
  USE mod_statevector_pdaf, &
       ONLY: mpi_subd_vert_child, mpi_subd_vert_par

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  INTEGER :: i                       ! Counters
  INTEGER :: off_p                   ! Process-local offset in global state vector


! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! Total dimension is sum of number of 2D state variables plus (number
  ! of 3D variables * number of vertical points) on both parent grid
  ! and child grid (if AGRIF used).
  dim_l = 1 + (4*mpi_subd_vert_par)
#if defined key_agrif
  dim_l = dim_l + 1 + (4*mpi_subd_vert_child)
#endif

! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

!!$  ! Coordinates are defined using T longitude/latitude grid values
!!$
!!$  ! First, compute (i,j) grid coordinates of local domain
!!$  j = INT( CEILING(REAL(domain_p)/REAL(nx_global)) )
!!$  i = domain_p - (j-1)*(nx_global)
!!$
!!$  ! Now, convert to T longitude/latitude grid values.
!!$  ! NOTE: tlon and tlat are in radians, and there are ghost cells.
!!$  coords_l(1)=tlon(i+1,j+1,1)
!!$  coords_l(2)=tlat(i+1,j+1,1)

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_dim_l_pdaf
