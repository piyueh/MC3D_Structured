!**********************************************************************
!   HEADING: MC3D - Structured Grids: VARIABLES MODULE for Initial
!   AUTHOR: MJ Huang, PY Chuang
!   PURPOSE: All the global constants and variables are defined.
!   DATE : 07/10/2009
!**********************************************************************
MODULE mod_VARIABLES
IMPLICIT NONE

!fixed PARAMETERs
!======================================================================
!----------------------------------------------------------------------
REAL(KIND=8), PARAMETER:: M_PI = 3.14159265358979323846D0
REAL(KIND=8), PARAMETER:: M_PI_2 = 6.28318530717958647692D0
REAL(KIND=8), PARAMETER:: M_PI_half = 1.57079632679489661923D0
REAL(KIND=8), PARAMETER:: M_PI_3 = 9.4247779607693797153879301498385D0
REAL(KIND=8), PARAMETER:: kB = 8.617d-2, barh = 6.5822d-1 ! meV/K, meV.ps
REAL(KIND=8), PARAMETER:: zero_tol = 1D-12
INTEGER(KIND=4), PARAMETER:: LRSi = 110, LRGe = 120, LR = 130
INTEGER(KIND=4), PARAMETER:: LW1 = 140, LW2 = 150, LW3 = 160


!material information
!======================================================================
!rho_Ge=5.323d3, rho_Si=2.329d3 [kg/m^3]
!----------------------------------------------------------------------
INTEGER(KIND=4):: iN_Ge1, iN_Ge2
INTEGER(KIND=4):: iN_Si1, iN_Si2
REAL(KIND=8):: Ge_start, Si_start, dU_Ge, dU_Si
REAL(KIND=8):: rho(2) = (/5.323d3, 2.329d3/)
REAL(KIND=8), ALLOCATABLE:: Ge_table(:, :), Si_table(:, :)
INTEGER(KIND=4), ALLOCATABLE:: iCmat(:, :, :)

!domain information
!======================================================================
! option=1: partially specularly and partially diffusely reflected
! option=2: periodic boundary condition
! option=3: leaving and being saved for heat control
!----------------------------------------------------------------------
INTEGER(KIND=4):: option(3)
REAL(KIND=8):: dLdomain(3), dPP(3), dPPB(3)


!cell information
!======================================================================
!----------------------------------------------------------------------
INTEGER(KIND=4):: iNcell(3)
REAL(KIND=8):: dLclen(3), dArea, dVolume  ! nm, nm^3
INTEGER(KIND=4), ALLOCATABLE:: iNnumcell(:, :, :), iNbgcell(:, :, :)
REAL(KIND=8), ALLOCATABLE:: dEcell(:, :, :), dEunit(:, :, :)
REAL(KIND=8), ALLOCATABLE:: dVunit(:, :, :), MFP(:, :, :)
REAL(KIND=8), ALLOCATABLE:: dTemp(:, :, :), dEdiff(:, :, :)


!heatcontrol information
!======================================================================
!----------------------------------------------------------------------
INTEGER(KIND=4):: WAY_HEAT, WAY_DIR
REAL(KIND=8):: TEMP_HEAT(2)
INTEGER(KIND=4), ALLOCATABLE:: mlost(:, :, :), iNemit(:, :, :)
REAL(KIND=8), ALLOCATABLE:: dPpool(:, :, :, :, :)
REAL(KIND=8), ALLOCATABLE:: dtheat(:), dEinject(:, :, :, :)
REAL(KIND=8), ALLOCATABLE:: dVinject(:, :, :, :), dElost(:, :, :)
REAL(KIND=8), ALLOCATABLE:: dEheatflux(:, :, :)


! phonon information
!======================================================================
! 3D simulation:: iNprop=8
! 1: x 2: y 3: z  4: theta 5: phi 6: energy
! 7: velocity 8: energy_material
!----------------------------------------------------------------------
INTEGER(KIND=4):: iNph, newiNph, iNprop, iNmakeup !iNmakeup is pool size
REAL(KIND=8), ALLOCATABLE:: phn(:, :), newphn(:, :)


!simulation parameters
!======================================================================
!----------------------------------------------------------------------
REAL(KIND=8):: dt, bundle(2), time0, time
INTEGER(KIND=4):: iter0


!output data
!======================================================================
!----------------------------------------------------------------------
REAL(KIND=8), ALLOCATABLE:: Tz(:, :, :), qflow(:, :), qctrl(:, :)
CHARACTER(LEN=72):: casename, inputfilename, outputfilename, restartfilename
CHARACTER(LEN=72):: logfilename

END MODULE mod_VARIABLES
