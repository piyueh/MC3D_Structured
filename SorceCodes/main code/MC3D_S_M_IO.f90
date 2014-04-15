!****************************************************************************
!   HEADING: MC3D I/O MODULE
!   AUTHOR: MJ Huang
!   PURPOSE: This module handles the input/output processes.
!   DATE: 2009.7.10
!****************************************************************************
MODULE mod_IO
USE mod_VARIABLES
USE mod_ROUTINES
IMPLICIT NONE
!############################################################################
CONTAINS
!============================================================================
SUBROUTINE readtable
IMPLICIT NONE

OPEN(unit=LRGe,file="Ge_real_table.txt")
OPEN(unit=LRSi,file="Si_real_table.txt")

READ(LRGe,*) iN_Ge1,iN_Ge2
READ(LRSi,*) iN_Si1,iN_Si2

ALLOCATE( Ge_table(iN_Ge1,iN_Ge2) )
ALLOCATE( Si_table(iN_Si1,iN_Si2) )
READ(LRGe,*) Ge_table
READ(LRSi,*) Si_table

Ge_start=Ge_table(3,1)
dU_Ge=Ge_table(3,2)-Ge_table(3,1)
Si_start=Si_table(3,1)
dU_Si=Si_table(3,2)-Si_table(3,1)

PRINT*,'# of Ge data = ',iN_Ge2
PRINT*,'# of Si data = ',iN_Si2
PRINT*,'Ge_start (meV) = ',Ge_start
PRINT*,'dU_Ge    (meV) = ',dU_Ge
PRINT*,'Si_start (meV) = ',Si_start
PRINT*,'Si_Ge    (meV) = ',dU_Si

CLOSE(LRGe)
CLOSE(LRSi)
END SUBROUTINE readtable
!============================================================================
SUBROUTINE initialize
IMPLICIT NONE
INTEGER*4:: idx
real*8:: falseheatflux
! DPP=specular fraction of internal interfaces
! DPPB=specular fraction of computational boundaries
READ(LR,*) bundle,dt,time0,iNcell 
READ(LR,*) dLdomain,dLclen,option,DPP,DPPB
!-------------------------------
 if (cba.eq.1) then
    read(LR,*) dEheatflux0
 else if (cba.eq.2) then
    read(LR,*) falseheatflux
 endif
!------------------------------

READ(LR,*) iNprop,iNph
PRINT*,'domain = ',dLdomain
PRINT*,'Ncells = ',iNcell
PRINT*,'bundle = ',bundle
PRINT*,'option = ',option
PRINT*,'DPP = ',DPP
PRINT*,'DPPB = ',DPPB
PRINT*,'dEheatflux [meV/(ps.nm^2)] = ',dEheatflux0
PRINT*,'dt (ps) = ',dt
PRINT*,'time (ps) = ',time0
!----------------------------------------------------------------------
dArea=dLclen(2)*dLclen(3)
dVolume=dLclen(1)*dArea
!------------------------
if (cba.ne.2) then
    dEheatflux0=dEheatflux0*dArea*dt*DBLE(iNcell(2)*iNcell(3))   ! meV/(ps.nm^2)*(nm^2)*(ps)=meV
endif
!--------------------------
PRINT*,'dEheatflux*dArea*dt (meV/cell) = ',dEheatflux0/DBLE(iNcell(2)*iNcell(3))

ALLOCATE(  iCmat(iNcell(1),iNcell(2),iNcell(3)),iNnumcell(iNcell(1),iNcell(2),iNcell(3)),iNbgcell(iNcell(1),iNcell(2),iNcell(3)) )
ALLOCATE( dEcell(iNcell(1),iNcell(2),iNcell(3)),    dTemp(iNcell(1),iNcell(2),iNcell(3)),  dEunit(iNcell(1),iNcell(2),iNcell(3)) )
ALLOCATE( dEdiff(iNcell(1),iNcell(2),iNcell(3)),   dVunit(iNcell(1),iNcell(2),iNcell(3)),     MFP(iNcell(1),iNcell(2),iNcell(3)) )

IF (option(1).eq.3) THEN
    ALLOCATE( dElost(iNcell(2),iNcell(3),2),dEinject(iNcell(2),iNcell(3),2,2),dVinject(iNcell(2),iNcell(3),2,2),iNemit(iNcell(2),iNcell(3),2) )
    ALLOCATE( qflow(iNcell(2),iNcell(3)),qctrl0(iNcell(2),iNcell(3)),dEheatflux(iNcell(2),iNcell(3),2) )
    if (cba.ne.2) then
        allocate(qctrl(iNcell(2),iNcell(3)))
    endif
ENDIF

ALLOCATE( phn(iNprop,iNph) )

dEdiff=0

READ(LR,*) iCmat
READ(LR,*) phn
READ(LR,*) dEdiff

CALL proc_reorder !聲子排序重整，程式剛運行時此指令無作用??? 只有開始跑程式後才有作用???

CALL cellinfo !得到網格的性質：單位體積能量(聲子束數量整數化後)、該能量對應之群速、MFP、聲子束能量

PRINT*,'# of phonons:: ', MAXVAL(iNnumcell),MINVAL(iNnumcell)
PRINT*,MAXVAL(dTemp)

END SUBROUTINE initialize
!============================================================================
SUBROUTINE restart
IMPLICIT NONE

WRITE(LW3,*) bundle,dt,time,iNcell 
WRITE(LW3,*) dLdomain,dLclen,option,DPP,DPPB
WRITE(LW3,*) dEheatflux0/(dArea*dt*DBLE(iNcell(2)*iNcell(3)))
WRITE(LW3,*) iNprop,iNph
WRITE(LW3,*) iCmat
WRITE(LW3,*) phn
WRITE(LW3,*) dEdiff
IF (WAY_DIR.eq.1) THEN
   WRITE(LW3,*) iNmakeup
   WRITE(LW3,*) dPpool
ENDIF
IF (WAY_HEAT.eq.1) THEN
   WRITE(LW3,*) sumQ
   WRITE(LW3,*) qctrl0
   WRITE(LW3,*) abc
ENDIF

END SUBROUTINE restart
!============================================================================
SUBROUTINE output_pure(iter,n1,n2,noutput)
    IMPLICIT NONE
    INTEGER*4::i,j,k,iter,N,n1,n2,noutput
    REAL*8:: Tmp,tot
    REAL*8,ALLOCATABLE :: Tzoft(:,:,:)   

    IF (n1.eq.1) THEN
        ALLOCATE( Tzoft(iNcell(1),iNcell(2),iNcell(3)) )        
        Tzoft = 0                       
    ENDIF

    DO k=1,iNcell(3)
        DO j=1,iNcell(2)
            DO i=1,iNcell(1)
                tot=tot+dEcell(i,j,k)
                CALL Etable(iCmat(i,j,k),1,dEcell(i,j,k),Tmp)
                IF (n1.eq.1) Tzoft(i,j,k)=Tmp                           
                IF (n2.eq.1) Tz(i,j,k)=Tz(i,j,k)+Tmp
            ENDDO
        ENDDO
    ENDDO

    IF (n1.eq.1) THEN
        WRITE(LW1,*) time,iNcell,dLclen
        WRITE(LW1,*) Tzoft 
        DEALLOCATE( Tzoft )  
        n1=0                     
    ENDIF
                                       
END SUBROUTINE output_pure
!============================================================================
END MODULE mod_IO
