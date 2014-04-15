!****************************************************************************
!   HEADING: MC3D MAIN PROGRAM
!   AUTHOR: MJ Huang
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
INTEGER*4::iter,iter0,iterations,n1,n2,noutput,nsteady,j,k,m,s
CHARACTER*40::inputfilename, outputfilename1, outputfilename2, restartfilename
REAL*8:: tmp
!PRINT*,'input filename = '
!READ*, inputfilename
!PRINT*,'output filename (transient) = '
!READ*, outputfilename1
!PRINT*,'output filename (stationary) = '
!READ*, outputfilename2
!PRINT*,'restart filename = '
!READ*, restartfilename
!PRINT*,'noutput = '
!READ*,noutput
!PRINT*, 'number of iterations to execute = '
!READ*, iterations
inputfilename='initial.txt'
outputfilename1='out1.txt'
outputfilename2='out2.txt'
restartfilename='restart.txt'
iterations=400000
noutput=10000000
nsteady=5000
iter0=5000
cba=1
OPEN(LR,file=inputfilename) !inputfilename='initial.txt'
open(1,file='heatflux.txt')
!open(999,file='test.txt')
!read Ge and Si table from files
CALL readtable
!initialization	
CALL initialize
!-------------------------------------------------------------
CALL random_seed() !要先呼叫random_seed()這函數，呼叫之後就可以使用random_number(欲儲存亂數的變數)這含個函數來得到變數

WAY_HEAT=0 !1為定熱通量，2為固定邊界溫度
WAY_DIR=0 !1為週期入射法，2為亂數入設法

IF (option(1).eq.3) THEN !option=3的話是熱流控制邊界條件
    WAY_HEAT=1
        !PRINT*,'Way to inject heat phonons:'
        !PRINT*,'1: phonons injected at specified heat flux'
        !PRINT*,'2: phonons injected at specified temperatures'
        !READ*,WAY_HEAT
    WAY_DIR=1
        !PRINT*,'Way to assign the directions of incident phonons'
        !PRINT*,'1: periodically assigned '
        !PRINT*,'2: randomly assigned'
        !READ*,WAY_DIR
    
    !----------------------------------------------------------------------------------------------------      
    IF (WAY_HEAT.eq.2) THEN !2: phonons injected at specified temperatures
        PRINT*,'the specified temperatures (K): TL,TR = '
	    READ*, TEMP_HEAT !要自行從螢幕輸入
	    
	    dEheatflux(:,:,1)=TEMP_HEAT(1)
	    dEheatflux(:,:,2)=TEMP_HEAT(2)
	    
	    Call proc_BC(dEheatflux(1:iNcell(2),1:iNcell(3),1),dEheatflux(1:iNcell(2),1:iNcell(3),2))!得到邊界入射聲子性質
	  
	    DO k=1,iNcell(3)
	        DO j=1,iNcell(2)
	     
	            CALL proc_energy(iCmat(1,j,k),TEMP_HEAT(1),tmp) !此行tmp:TEMP_HEAT(1)對應的單位體積能量
	            
	            dEheatflux(j,k,1)=tmp*dVinject(j,k,iCmat(1,j,k),1)
		        
		        CALL proc_energy(iCmat(iNcell(1),j,k),TEMP_HEAT(2),tmp) !此行tmp:TEMP_HEAT(2)對應的單位體積能量
	        
	            dEheatflux(j,k,2)=tmp*dVinject(j,k,iCmat(iNcell(1),j,k),2)
	        ENDDO
	    ENDDO
	    dEheatflux=dEheatflux/4d0*dArea*dt  !每個網格每個步進時間最後應該要入射或射出的能量
    ENDIF

    IF (WAY_DIR.eq.1) THEN
        READ(LR,*) iNmakeup !iNmakeup為pool裡每個邊界網格可以儲存的離開聲子束數量
        ALLOCATE( mlost(iNcell(2),iNcell(3),2),dPpool(6,iNmakeup,iNcell(2),iNcell(3),2) )
        !dPpool(x,：,：,：,：)：
        !x=1，為剩餘時間、2及3為非熱流方向的其餘兩方向位置、4及5為移動方向、6為能量材料屬性
      
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
    
    dElost=0

ENDIF !(option(1).eq.3)

PRINT*,'iNmakeup = ',iNmakeup

CLOSE(LR)

PRINT*,'PREPROCESSING DONE'
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
    
    CALL random_seed()

    CALL advance(iter,iter0)   
    time=time+dt

    IF (restartfilename.ne.'-') THEN !5
        IF (iter.eq.iterations.or.MOD(iter,nsteady).eq.0) THEN !5-1
            OPEN(LW3,file=restartfilename)
	        CALL restart
	        CLOSE(LW3)
        ENDIF !5-1
    ENDIF !5
            
    PRINT*, iter,iNmakeup,MINVAL(mlost)
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
