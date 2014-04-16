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
INTEGER*4:: iter,iter0,iterations,n1,n2,noutput,nsteady,j,k,m,s
INTEGER*4:: NorR
LOGICAL:: TF

    CALL RANDOM_SEED()

    WRITE(*,*) 'INPUT THE CASE NAME(NO .''TXT'', NO ''_GRIDFILE''): '
    READ(*,"(A72)") casename
    inputfilename = casename(1:LEN_TRIM(casename))//'_initial.txt'
    restartfilename = casename(1:LEN_TRIM(casename))//'_restart.txt'
    logfilename = casename(1:LEN_TRIM(casename))//'_log.txt'

    iterations=400000
    noutput=10000000
    nsteady=5000
    iter0=5000

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
    CALL initialize

    WRITE(*, *) 'PREPROCESSING HAS FINISHED!!'
!----------------------------------------------------------------------------------------------------
n1=0
IF (outputfilename1.ne.'-') THEN
   OPEN(LW1,file=outputfilename1)
   WRITE(LW1,*) iterations/noutput
ENDIF
n2=0
IF (outputfilename2.ne.'-') THEN
   OPEN(LW2,file=outputfilename2)
   ALLOCATE( Tz(iNcell(1),iNcell(2),iNcell(3)) )
   Tz=0
   qflow=0
ENDIF

time = time0

DO iter=0,iterations

    IF (mod(iter+1,noutput).eq.0) n1=1  !1(第一個IF)) !noutput：每幾步要output資料一次

    IF (n1+n2.ne.0) CALL output_pure(iter,n1,n2,noutput) !2

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
SUBROUTINE advance(iter,iter0) !所有聲子實際運動過程
    USE mod_VARIABLES
    USE mod_ADVANCE
    USE mod_ROUTINES
    USE mod_heatcontrol
    IMPLICIT NONE
    INTEGER*4::iter,iter0,i,j,k,bg,ed
    REAL*8,ALLOCATABLE::cellbdy(:,:)
    !------------------------------------------------
    ALLOCATE( cellbdy(2,3) )
    DO k=1,iNcell(3)
        DO j=1,iNcell(2)
            DO i=1,iNcell(1)
                cellbdy(1,1)=dLclen(1)*DBLE(i-1)
                cellbdy(2,1)=dLclen(1)*DBLE(i)
                cellbdy(1,2)=dLclen(2)*DBLE(j-1)
                cellbdy(2,2)=dLclen(2)*DBLE(j)
                cellbdy(1,3)=dLclen(3)*DBLE(k-1)
                cellbdy(2,3)=dLclen(3)*DBLE(k)
                bg=iNbgcell(i,j,k)+1
                ed=iNbgcell(i,j,k)+iNnumcell(i,j,k) !所以第(i,j,k)網格裡的聲子編號是第bg個~第ed個
                CALL proc_advection(i,j,k,cellbdy,iNnumcell(i,j,k),phn(1:iNprop,bg:ed),1)
            ENDDO
        ENDDO
    ENDDO

    ! So far, the movement, boundary scattering, and interfacial scattering procedures have been finished.
    IF (option(1).eq.3) CALL proc_heatcontrol(iter,iter0)
    !-------build cell information for future use

    CALL proc_reorder
    CALL proc_createdelete
    CALL cellinfo
    DEALLOCATE( cellbdy )
END SUBROUTINE advance
!============================================================================
SUBROUTINE cellinfo !所有網格內所有聲子的能量總和(單位體積)、該能量對應之溫度、每顆聲子(束)的能量、平均聲子群速、MFP
USE mod_VARIABLES
USE mod_ROUTINES
IMPLICIT NONE
INTEGER*4::i,j,k,bg,ed

DO k=1,iNcell(3)
    DO j=1,iNcell(2)
        DO i=1,iNcell(1)
            bg=iNbgcell(i,j,k)+1
            ed=iNbgcell(i,j,k)+iNnumcell(i,j,k)
            dEcell(i,j,k) = SUM( phn(6,bg:ed) )/dVolume  ! total energy in this cell(是聲子束整數化後的能量))
            CALL Etable(iCmat(i,j,k),2,dEcell(i,j,k),dEunit(i,j,k)) ! compute the equilibrium number of phonons
            CALL Etable(iCmat(i,j,k),4,dEcell(i,j,k),dVunit(i,j,k)) ! compute the average phonon group velocity
            CALL Etable(iCmat(i,j,k),5,dEcell(i,j,k),MFP(i,j,k)) !calculate MFP in this cell
            CALL Etable(iCmat(i,j,k),1,dEcell(i,j,k),dTemp(i,j,k)) !該單位體積能量對應之溫度
            dEunit(i,j,k)=dEcell(i,j,k)/dEunit(i,j,k)*bundle(iCmat(i,j,k))  ! average energy per phonon bundle
        ENDDO
    ENDDO
ENDDO

END SUBROUTINE cellinfo
!==================================================================================
