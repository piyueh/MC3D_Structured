!****************************************************************************
!   HEADING: MC3D - Structured Solver: I/O MODULE
!   AUTHOR: MJ Huang, PY Chuang
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

        OPEN(UNIT = LRGe, FILE = "Ge_real_table.txt")
        OPEN(UNIT = LRSi, FILE = "Si_real_table.txt")

        READ(LRGe, *) iN_Ge1, iN_Ge2
        READ(LRSi, *) iN_Si1, iN_Si2

        ALLOCATE(Ge_table(iN_Ge1, iN_Ge2))
        ALLOCATE(Si_table(iN_Si1, iN_Si2))
        READ(LRGe,*) Ge_table
        READ(LRSi,*) Si_table

        Ge_start = Ge_table(3, 1)
        dU_Ge = Ge_table(3, 2) - Ge_table(3, 1)
        Si_start = Si_table(3, 1)
        dU_Si = Si_table(3, 2) - Si_table(3, 1)

        WRITE(*, *) '# of Ge data = ',iN_Ge2
        WRITE(*, *) '# of Si data = ',iN_Si2
        WRITE(*, *) 'Ge_start (meV) = ',Ge_start
        WRITE(*, *) 'dU_Ge    (meV) = ',dU_Ge
        WRITE(*, *) 'Si_start (meV) = ',Si_start
        WRITE(*, *) 'Si_Ge    (meV) = ',dU_Si

        CLOSE(LRGe)
        CLOSE(LRSi)

    END SUBROUTINE readtable
!============================================================================
!============================================================================
    SUBROUTINE initialize
    IMPLICIT NONE
    INTEGER*4:: idx
    ! DPP=specular fraction of internal interfaces
    ! DPPB=specular fraction of computational boundaries

        OPEN(LR, FILE = inputfilename)

        READ(LR, *) bundle,dt,time0,iNcell
        READ(LR, *) dLdomain,dLclen,option,DPP,DPP
        READ(LR, *) iNprop,iNph

        dArea = dLclen(2) * dLclen(3)
        dVolume = dLclen(1) * dArea
        dEheatflux0 = dEheatflux0 * dArea * dt * DBLE(iNcell(2) * iNcell(3))
        ! meV/(ps.nm^2)*(nm^2)*(ps)=meV
!--------------------------


        ALLOCATE( iCmat(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( iNnumcell(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( iNbgcell(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dEcell(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dTemp(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dEunit(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dEdiff(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dVunit(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( MFP(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( phn(iNprop, iNph) )

        READ(LR,*) iCmat
        READ(LR,*) phn
        READ(LR,*) dEdiff

        IF (option(1).eq.3) THEN
            ALLOCATE( dElost(iNcell(2), iNcell(3), 2) )
            ALLOCATE( dEinject(iNcell(2), iNcell(3), 2, 2) )
            ALLOCATE( dVinject(iNcell(2), iNcell(3), 2, 2) )
            ALLOCATE( iNemit(iNcell(2), iNcell(3), 2) )
            ALLOCATE( qflow(iNcell(2), iNcell(3)) )
            ALLOCATE( qctrl0(iNcell(2), iNcell(3)) )
            ALLOCATE( dEheatflux(iNcell(2), iNcell(3), 2) )
            ALLOCATE( qctrl(iNcell(2), iNcell(3)) )

            IF (WAY_DIR.eq.1) THEN
                READ(LR,*) iNmakeup
                ALLOCATE( mlost(iNcell(2), iNcell(3), 2) )
                ALLOCATE( dPpool(6, iNmakeup, iNcell(2), iNcell(3), 2) )
                mlost=0
                dPpool=0
                READ(LR,*) dPpool
            ENDIF

            IF (WAY_HEAT.eq.1) THEN !1: phonons injected at specified heat flux
                sumQ=0
                qctrl0=0
                READ(LR,*) sumQ
                READ(LR,*) qctrl0
                read(LR,*) abc
                qctrl=qctrl0/sumQ
            ENDIF
        ENDIF

        dEdiff=0



        CALL proc_reorder
        CALL cellinfo

        CLOSE(LR)

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
