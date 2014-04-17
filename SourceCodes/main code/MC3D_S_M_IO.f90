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
    SUBROUTINE initialize(NorR)
    IMPLICIT NONE
    INTEGER*4:: idx, NorR
    ! DPP=specular fraction of internal interfaces
    ! DPPB=specular fraction of computational boundaries

        SELECTCASE(NorR)
        CASE(1)
            WRITE(*, *) "Start to Read Initialization Data..."
            OPEN(LR, FILE = inputfilename)
        CASE(2)
            WRITE(*, *) "Start to Read Restart Data..."
            OPEN(LR, FILE = restartfilename)
        ENDSELECT

        WRITE(*, *) "   Reading dt..."
        READ(LR, *) dt
        WRITE(*, *) "   Reading iter0..."
        READ(LR, *) iter0
        WRITE(*, *) "   Reading dLdomain..."
        READ(LR, *) dLdomain
        WRITE(*, *) "   Reading iNcell..."
        READ(LR, *) iNcell
        WRITE(*, *) "   Reading dLclen..."
        READ(LR, *) dLclen
        WRITE(*, *) "   Reading bundle..."
        READ(LR, *) bundle
        WRITE(*, *) "   Reading iNprop & iNph..."
        READ(LR, *) iNprop,iNph
        WRITE(*, *) "   Reading iNmakeup..."
        READ(LR, *) iNmakeup

        ALLOCATE( iCmat(iNcell(1), iNcell(2), iNcell(3)) )
        WRITE(*, *) "   Reading iCmat..."
        READ(LR, *) iCmat

        ALLOCATE( phn(iNprop, iNph) )
        WRITE(*, *) "   Reading phn..."
        READ(LR, *) phn

        ALLOCATE( dEdiff(iNcell(1), iNcell(2), iNcell(3)) )
        WRITE(*, *) "   Reading dEdiff..."
        READ(LR, *) dEdiff

        ALLOCATE( dPpool(6, iNmakeup, iNcell(2), iNcell(3), 2) )
        WRITE(*, *) "   Reading dPpool..."
        READ(LR, *) dPpool

        ALLOCATE( mlost(iNcell(2), iNcell(3), 2) )
        WRITE(*, *) "   Reading mlost..."
        READ(LR, *) mlost

        CLOSE(LR)

        ALLOCATE( iNnumcell(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( iNbgcell(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dEcell(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dTemp(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dEunit(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( dVunit(iNcell(1), iNcell(2), iNcell(3)) )
        ALLOCATE( MFP(iNcell(1), iNcell(2), iNcell(3)) )


        dArea = dLclen(2) * dLclen(3)
        dVolume = dLclen(1) * dArea

        IF (option(1).eq.3) THEN
            ALLOCATE( dElost(iNcell(2), iNcell(3), 2) )
            ALLOCATE( dEinject(iNcell(2), iNcell(3), 2, 2) )
            ALLOCATE( dVinject(iNcell(2), iNcell(3), 2, 2) )
            ALLOCATE( iNemit(iNcell(2), iNcell(3), 2) )
            ALLOCATE( qflow(iNcell(2), iNcell(3)) )
            ALLOCATE( qctrl0(iNcell(2), iNcell(3)) )
            ALLOCATE( dEheatflux(iNcell(2), iNcell(3), 2) )
            ALLOCATE( qctrl(iNcell(2), iNcell(3)) )
            
            dElost = 0d0
            
            SELECTCASE(WAY_HEAT)
            CASE(1)
            ! phonons injected at specified heat flux
            !==========================================================
                sumQ=0
                qctrl0=0
                READ(LR,*) sumQ
                READ(LR,*) qctrl0
                read(LR,*) abc
                qctrl=qctrl0/sumQ
            CASE(2)
            ! phonons injected at specified temperatures
            !==========================================================
                dEheatflux(:,:,1)=TEMP_HEAT(1)
                dEheatflux(:,:,2)=TEMP_HEAT(2)

                Call proc_BC( dEheatflux(1:iNcell(2), 1:iNcell(3), 1), &
                                        dEheatflux(1:iNcell(2), 1:iNcell(3), 2) )

                DO k = 1, iNcell(3)
                    DO j = 1, iNcell(2)
                        CALL proc_energy( iCmat(1,j,k), TEMP_HEAT(1), tmp)
                        dEheatflux(j, k, 1) = tmp * dVinject(j, k, iCmat(1,j,k), 1)
                        CALL proc_energy( iCmat(iNcell(1), j, k), TEMP_HEAT(2), tmp)
                        dEheatflux(j, k, 2) = tmp * dVinject(j, k, iCmat(iNcell(1), j, k), 2)
                    ENDDO
                ENDDO
                dEheatflux = dEheatflux / 4d0 * dArea * dt
                ! dEheatflux:
                !       The energy that should be injected at each time step
                !       and each boundary element.
            ENDSELECT
            !==========================================================
        ENDIF

        CALL proc_reorder
        CALL cellinfo

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
END MODULE mod_IO
