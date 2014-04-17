!****************************************************************************
!   HEADING: MC3D - Structured Solver: Main PROGRAM
!   AUTHOR: MJ Huang, PY Chuang
!   PURPOSE: To simulate thermal conductivity of nanostructured material
!   DATE: 2009.7.10
!****************************************************************************
PROGRAM DSMC3D
USE mod_VARIABLES
USE mod_heatcontrol
USE mod_ROUTINES
USE mod_ADVANCE
USE mod_IO
IMPLICIT NONE
INTEGER*4:: iter, NorR

    CALL RANDOM_SEED()

    WRITE(*,*) 'INPUT THE CASE NAME(NO .''TXT'', NO ''_GRIDFILE''): '
    READ(*,"(A72)") casename
    inputfilename = casename(1:LEN_TRIM(casename))//'_initial.txt'
    restartfilename = casename(1:LEN_TRIM(casename))//'_restart.txt'
    logfilename = casename(1:LEN_TRIM(casename))//'_log.txt'

    iterations=5000000
    noutput=2500

    WRITE(*, *) "Enter 1 for a New Simulation or 2 to Restart a Existing Simulation:"
    READ(*, *) NorR

    SELECTCASE(NorR)
    CASE(1)
        OPEN(UNIT=LW1, FILE=logfilename)

        WRITE(*, *) "Enter the type of BCx, BCy, BCz in Each Direction:"
        WRITE(*, *) "   1: partially specularly and partially diffusely reflected"
        WRITE(*, *) "   2: periodic boundary condition"
        WRITE(*, *) "   3: leaving and being saved for heat control"
        READ(*, *) option

        WRITE(LW1, *) option

        WRITE(*, *) "Enter the Specularity of Domain Boundaries: "
        READ(*, *) dPPB
        WRITE(*, *) "Enter the Specularity of Interface: "
        READ(*, *) dPP

        WRITE(LW1, *) dPPB, dPP

        IF ((option(2).eq.3).OR.(option(3).eq.3)) THEN
            WRITE(*, *) "BCy = 3 or BCz = 3 is NOT Supported in This Version."
            WRITE(*, *) "The Program is Going to Shut Down in 5 Seconds."
            CALL SLEEP(5)
            STOP
        ENDIF

        IF (option(1).eq.3) THEN

            WRITE(*, *) 'Way to inject heat phonons:'
            WRITE(*, *) '    1: phonons injected at specified heat flux (Discarded in This Version)'
            WRITE(*, *) '    2: phonons injected at specified temperatures'
            WRITE(*, *) 'Enter the # of option and press <ENTER>: '
            READ(*, *) WAY_HEAT

            WRITE(LW1, *) WAY_HEAT

            IF (WAY_HEAT.eq.1) THEN
                !WRITE(*, *) "Then, the specific heat flux (meV / nm^2-ps) = "
                !READ(*, *) dEheatflux0
                !WRITE(LW1, *) dEheatflux0
                !dEheatflux0 = dEheatflux0 * dArea * dt * DBLE(iNcell(2) * iNcell(3))
                WRITE(*, *) "Phonons Injected at Specified Heat Flux is Discarded in This Version."
                WRITE(*, *) "The Program is Going to Shut Down in 5 Seconds."
                CALL SLEEP(5)
                STOP
            ELSEIF (WAY_HEAT.eq.2) THEN
                WRITE(*, *) 'Then, the specified temperatures (K): TL,TR = '
                READ(*, *) TEMP_HEAT
                WRITE(LW1, *) TEMP_HEAT
            ELSE
                WRITE(*, *) 'Wrong Option!'
                WRITE(*, *) "The Program is Going to Shut Down in 5 Seconds."
                CALL SLEEP(5)
                STOP
            ENDIF

            WRITE(*, *) 'Way to assign the directions of incident phonons'
            WRITE(*, *) '    1: periodically assigned '
            WRITE(*, *) '    2: randomly assigned'
            WRITE(*, *) 'Enter the # of option and press <ENTER>: '
            READ(*, *) WAY_DIR
            WRITE(LW1, *) WAY_DIR
        ENDIF

        CLOSE(LW1)

    CASE(2)

        OPEN(UNIT=LR, FILE=logfilename)
        READ(LR, *) option
        READ(LR, *) dPPB, dPP
        IF (option(1).eq.3) THEN
            READ(LR, *) WAY_HEAT
            READ(LR, *) TEMP_HEAT
            READ(LR, *) WAY_DIR
        ENDIF
        CLOSE(LR)

    CASE DEFAULT

        WRITE(*, *) 'Wrong Option!'
        WRITE(*, *) "The Program is Going to Shut Down in 5 Seconds."
        CALL SLEEP(5)
        STOP

    ENDSELECT

    ! Read Material Property Table
    CALL readtable

    ! Initialization
    CALL initialize(NorR)

    WRITE(*, *) 'PREPROCESSING SESSION HAS FINISHED!!'

    DO iter = iter0,iterations

        IF (n2.ne.0.and.mod(iter-iter0,nsteady).eq.0) THEN !3

            WRITE(LW2,*) iNcell,dLclen,dt,dLdomain,nsteady
            WRITE(LW2,*) qflow/DBLE(nsteady)
            WRITE(LW2,*) Tz/DBLE(nsteady)

            IF (WAY_HEAT.eq.1) THEN !3-1
                sumQ=sumQ+SUM( qflow )
                qctrl0=qctrl0+qflow
                qctrl=qctrl0/sumQ
            !-----------------------------
                if (abc.eq.1.and.iter.eq.10000) then !3-1-1
                    open(LR,file='initial.txt')
                    cba=2
                    DEALLOCATE( iCmat,iNnumcell,iNbgcell )
                    DEALLOCATE( dEcell, dTemp, dEunit )
                    DEALLOCATE( dEdiff,  dVunit, MFP )
                    DEALLOCATE( dElost,dEinject,dVinject,iNemit )
                    DEALLOCATE( qflow,qctrl0,dEheatflux,phn )
                    call initialize
                    close(LR)

                    dEheatflux0=dEheatflux0*50d0/(sum(Tz(1,:,:))-sum(Tz(iNcell(1),:,:)))*dble(iNcell(2)*iNcell(3)*nsteady)
                    sumQ=0
                    qctrl0=0
                    write(1,*) dEheatflux0/dArea/dt/DBLE(iNcell(2)*iNcell(3))
                    abc=2
                    close(1)
                endif !3-1-1
	            !-------------------------------
            ENDIF !3-1

            qflow=0
            Tz=0
        ENDIF !3

        IF (outputfilename2.ne.'-'.and.iter.eq.iter0) THEN !4
            IF (WAY_HEAT.eq.1) THEN !4-1
                sumQ=sumQ+SUM( qflow )
                qctrl0=qctrl0+qflow
                qctrl=qctrl0/sumQ
            ENDIF !4-1
            n2=1
            qflow=0
        ENDIF !4

        CALL advance(iter,iter0)
        time=time+dt

        IF (restartfilename.ne.'-') THEN !5
            IF (iter.eq.iterations.or.MOD(iter,nsteady).eq.0) THEN !5-1
                OPEN(LW3,file=restartfilename)
                CALL restart
                CLOSE(LW3)
            ENDIF !5-1
        ENDIF !5

        WRITE(*, *)  iter,iNmakeup,MINVAL(mlost)
ENDDO

IF (n2.ne.0) THEN
   CLOSE(LW2)
   DEALLOCATE( Tz )
ENDIF

DEALLOCATE( Ge_table,Si_table )
DEALLOCATE( iCmat,iNnumcell,iNbgcell,dEcell,dEdiff,dEunit,dVunit,MFP,dTemp,phn,qflow )
IF (option(1).eq.3) DEALLOCATE( dElost,dEinject,dVinject,dEheatflux,iNemit,qctrl,qctrl0,qflow )
IF (WAY_DIR.eq.1) DEALLOCATE( dPpool,mlost )

END PROGRAM DSMC3D
!============================================================================

!============================================================================

!==================================================================================
