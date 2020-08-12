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

  USE mod_kind_pdaf
  USE mod_assimilation_pdaf, &
       ONLY: coords_l, indx_dom_l_par, indx_dom_l_child, num_domains_par
  USE mod_statevector_pdaf, &
       ONLY: mpi_subd_vert_child, mpi_subd_vert_par, var2d_p_offset_par, &
       var2d_p_offset_child, var3d_p_offset_par, var3d_p_offset_child
  USE mod_parallel_pdaf, &
       ONLY: nldj_child, nldj_par, nldi_child, nldi_par
  USE mod_agrif_pdaf, &
       ONLY: glamt_child, glamt_par, gphit_child, gphit_par, tmask_par, &
       tmask_child

  IMPLICIT NONE

! *** Arguments ***
  INTEGER, INTENT(in)  :: step     !< Current time step
  INTEGER, INTENT(in)  :: domain_p !< Current local analysis domain
  INTEGER, INTENT(out) :: dim_l    !< Local state dimension

! *** local variables ***
  INTEGER :: cnt, indx        ! Counters
  INTEGER :: i_par, j_par     ! Grid coordinates for local analysis domain
  INTEGER :: i_child, j_child ! Grid coordinates for local analysis domain
  INTEGER :: i0, j0           ! Halo offset for local PE
  INTEGER :: domain_p_child   ! Counter for local analysis domain on child grid
  REAL(pwp) :: lon, lat       ! Longitude, latitude for local analysis domain
  REAL(pwp) :: rad_conv = 3.141592653589793/180.0   ! Degree to radian conversion


! **********************************************
! *** Initialize coordinates of local domain ***
! **********************************************

  ! Determine whether local domain belongs to parent or child grid.
  ! Index for domain_p starts at 0.
  IF(domain_p <= num_domains_par) THEN
     ! Compute halo offset for parent grid.
     i0 = nldi_par - 1
     j0 = nldj_par - 1
     ! Compute coordinates for parent grid.
     i_par = indx_dom_l_par(1, domain_p)
     j_par = indx_dom_l_par(2, domain_p)
     ! Use T-values to get local coordinates.
     lat = gphit_par(i_par+i0,j_par+j0)
     lon = glamt_par(i_par+i0,j_par+j0)
     ! Convert local domain coordinates to radians (as required by PDAF-OMI)
     coords_l(1)=lon*rad_conv
     coords_l(2)=lat*rad_conv
  END IF

#if defined key_agrif
  ! Determine whether local domain belongs to parent or child grid.
  IF(domain_p > num_domains_par) THEN
     ! Convert to local domain value for child grid.
     domain_p_child = domain_p - num_domains_par
     ! Compute halo offset on child grid.
     i0 = nldi_child - 1
     j0 = nldj_child - 1
     ! Compute coordinates for child grid.
     i_child = indx_dom_l_child(1, domain_p_child)
     j_child = indx_dom_l_child(2, domain_p_child)
     ! Use T-values to get local coordinates.
     lat = gphit_child(i_child+i0,j_child+j0)
     lon = glamt_child(i_child+i0,j_child+j0)
     ! Convert local domain coordinates to radians (as required by PDAF-OMI)
     coords_l(1)=lon*rad_conv
     coords_l(2)=lat*rad_conv
  END IF
#endif

! ****************************************
! *** Initialize local state dimension ***
! ****************************************

  ! **********************************************************************
  ! Total dimension is sum of number of 2D state variables plus (number
  ! of 3D variables * number of ocean vertical points) on both parent grid
  ! and child grid (if AGRIF used). We need to perform the calculation for
  ! the number of ocean vertical points here ie we need to determine how
  ! many points in the vertical are ocean and how may are land (we do not
  ! include land points in our local state vector).
  ! **********************************************************************

  ! Initialization.
  cnt=0

  ! Parent grid
  IF(domain_p <= num_domains_par) THEN
     ! Compute halo offset for parent grid.
     i0 = nldi_par - 1
     j0 = nldj_par - 1
     ! Compute coordinates for parent grid.
     i_par = indx_dom_l_par(1, domain_p)
     j_par = indx_dom_l_par(2, domain_p)
     ! Count number of ocean vertical points
     DO indx = 1, mpi_subd_vert_par
        IF(tmask_par(i_par+i0, j_par+j0, indx) == 1) cnt = cnt + 1
     END DO
     dim_l = SIZE(var2d_p_offset_par) + (SIZE(var3d_p_offset_par)*cnt)
  END IF

#if defined key_agrif
  ! Child grid.
  IF(domain_p > num_domains_par) THEN
     ! Convert to local domain value for child grid.
     domain_p_child = domain_p - num_domains_par
     ! Compute halo offset on child grid.
     i0 = nldi_child - 1
     j0 = nldj_child - 1
     ! Compute coordinates for child grid.
     i_child = indx_dom_l_child(1, domain_p_child)
     j_child = indx_dom_l_child(2, domain_p_child)
     ! Count number of ocean vertical points
     DO indx = 1, mpi_subd_vert_child
        IF(tmask_child(i_child+i0, j_child+j0, indx) == 1) cnt = cnt + 1
     END DO
     dim_l = SIZE(var2d_p_offset_child) + (SIZE(var3d_p_offset_child)*cnt)
  END IF
#endif

!$AGRIF_END_DO_NOT_TREAT
END SUBROUTINE init_dim_l_pdaf
