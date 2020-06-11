MODULE kind_pdaf

!  !DESCRIPTION:
!  This module defines the kind of real for PDAF interface/call-back
!  routines.

!  !USES:
   IMPLICIT NONE


   INTEGER, PUBLIC, PARAMETER :: pdp = SELECTED_REAL_KIND(12,307) !: double precision (real 8)
   INTEGER, PUBLIC, PARAMETER :: pwp = pdp                        !: working precision

END MODULE kind_pdaf
