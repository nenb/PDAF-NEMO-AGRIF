MODULE mod_agrif_pdaf
!$AGRIF_DO_NOT_TREAT
!  !DESCRIPTION:
!  This module contains auxiliary functions for exchanging data between
!  NEMO/AGRIF and PDAF.

!  !USES:
   IMPLICIT NONE


   INTEGER :: nitend = 250
   INTEGER :: jpiglo = 722
   INTEGER :: jpjglo = 511
   INTEGER :: jpk = 46
   INTEGER :: nit000 = 0
   INTEGER :: jpnij = 96
   INTEGER :: nldi = 1
   INTEGER :: nlei = 30
   INTEGER :: nldj = 1
   INTEGER :: nlej = 20
   INTEGER :: nimpp = 1
   INTEGER :: njmpp = 1
   INTEGER :: jp_tem = 1
   INTEGER :: jp_sal = 2
   REAL :: rdt = 2160
   REAL :: ssmask(511,722)
   REAL :: tmask(511,722,46)
   REAL :: umask(511,722,46)
   REAL :: vmask(511,722,46)
   REAL :: sshb(511,722)
   REAL :: ub(511,722,46)
   REAL :: vb(511,722,46)
   REAL :: tsb(511,722,46,2)

!$AGRIF_END_DO_NOT_TREAT
END MODULE mod_agrif_pdaf
