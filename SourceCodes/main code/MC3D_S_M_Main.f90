!**********************************************************************
!   HEADING: MC3D - Structured Solver: Main Program
!   AUTHOR: MJ Huang, PY Chuang
!   PURPOSE: To simulate thermal conductivity of nanostructured
!            material
!   DATE: 2009.7.10
!**********************************************************************
PROGRAM DSMC3D
USE mod_VARIABLES
USE mod_heatcontrol
USE mod_ROUTINES
USE mod_ADVANCE
USE mod_IO
IMPLICIT NONE
INTEGER(KIND=4):: iter, NorR

    CALL RANDOM_SEED()

    WRITE(*,*) 'INPUT THE CASE NAME(NO .''TXT''): '
    READ(*,"(A72)") casename
    inputfilename = casename(1:LEN_TRIM(casename))//'_initial.txt'
    restartfilename = casename(1:LEN_TRIM(casename))//'_restart.txt'
    logfilename = casename(1:LEN_TRIM(casename))//'_log.txt'

    iterations = 5000000
    noutput = 2500
    backupstep = 5000

    WRITE(*, *) "Enter 1 for a New Simulation or 2 to Restart a " // &
                "Existing Simulation:"
    READ(*, *) NorR

    SELECTCASE(NorR)
    CASE(1)
        OPEN( UNIT = LW1, FILE = logfilename)

        WRITE(*, *) "Enter the type of BCx, BCy, BCz:"
        WRITE(*, *) "   1: partially specularly and partially " // &
                    "diffusely reflected"
        WRITE(*, *) "   2: periodic boundary condition"
        WRITE(*, *) "   3: leaving and being saved for heat control"
        READ(*, *) option

        WRITE(LW1, *) option

        WRITE(*, *) "Enter the Specularity of Domain Boundaries: "
        READ(*, *) dPPB
        WRITE(*, *) "Enter the Specularity of Interface: "
        READ(*, *) dPP

        WRITE(LW1, *) dPPB, dPP

        IF ( (option(2).eq.3) .OR. (option(3).eq.3) ) THEN
            WRITE(*, *) "BCy = 3 or BCz = 3 is NOT Supported in " // &
                        "This Version."
            WRITE(*, *) "The Program Will Be Shut Down in 5 Seconds."
            CALL SLEEP(5)
            STOP
        ENDIF

        IF (option(1).eq.3) THEN

            WRITE(*, *) 'Way to inject heat phonons:'
            WRITE(*, *) '    1: phonons injected at specified ' // &
                        'heat flux (Discarded in This Version)'
            WRITE(*, *) '    2: phonons injected at specified ' // &
                        'temperatures'
            WRITE(*, *) 'Enter the # of option and press <ENTER>: '
            READ(*, *) WAY_HEAT

            WRITE(LW1, *) WAY_HEAT

            IF (WAY_HEAT.eq.1) THEN
                !WRITE(*, *) "Then, the specific heat flux " // &
                !            "(meV / nm^2-ps) = "
                !READ(*, *) dEheatflux0
                !WRITE(LW1, *) dEheatflux0
                !dEheatflux0 = dEheatflux0 * dArea * dt * &
                !              DBLE(iNcell(2) * iNcell(3))
                WRITE(*, *) "Phonons Injected at Specified Heat" // &
                            " Flux is Discarded in This Version."
                WRITE(*, *) "The Program is Going to Shut Down in" // &
                            " 5 Seconds."
                CALL SLEEP(5)
                STOP
            ELSEIF (WAY_HEAT.eq.2) THEN
                WRITE(*, *) "Then, the specified temperatures (K):" //&
                            " TL,TR = "
                READ(*, *) TEMP_HEAT
                WRITE(LW1, *) TEMP_HEAT
            ELSE
                WRITE(*, *) 'Wrong Option!'
                WRITE(*, *) "The Program is Going to Shut Down in" // &
                            " 5 Seconds."
                CALL SLEEP(5)
                STOP
            ENDIF

            WRITE(*, *) 'Way to assign the directions of incident' // &
                        ' phonons'
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

    DO iter = iter0 + 1, iterations

        CALL advance
        time = time + dt
        Tz = Tz + dTemp
        Ez = Ez + dEcell
        ct = ct + 1

        IF ( MOD( iter, noutput ).eq.0 ) THEN
            WRITE(outputfilename, "(I8.8, '.txt')") iter
            OPEN( UNIT = LW1, FILE = outputfilename )
            WRITE(LW1, *) iter, time, ct
            WRITE(LW1, *) SUM( qflowL ) / DBLE( ct )
            WRITE(LW1, *) SUM( qflow ) / DBLE( ct )
            WRITE(LW1, *) SUM( qflowR ) / DBLE( ct )
            WRITE(LW1, *) qflowL
            WRITE(LW1, *) qflow
            WRITE(LW1, *) qflowR
            WRITE(LW1, *) Tz / DBLE( ct )
            WRITE(LW1, *) Ez / DBLE( ct )
            CLOSE( LW1 )
            qflow = 0d0
            qflowL = 0d0
            qflowR = 0d0
            Tz = 0d0
            Ez = 0d0
            ct = 0d0
        ENDIF

        IF ( MOD( iter, backupstep) ) THEN
            CALL restart
        ENDIF

        WRITE(*, *)  iter, time, ct
    ENDDO


END PROGRAM DSMC3D
!======================================================================
!======================================================================
