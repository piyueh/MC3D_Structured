!****************************************************************************
!   HEADING: MC3D ADVANCE MODULE
!   AUTHOR: MJ Huang
!   PURPOSE: This module time marches the simulation.
!   DATE : 2009.7.10
!****************************************************************************
MODULE mod_ADVANCE
USE mod_VARIABLES
USE mod_ROUTINES
IMPLICIT NONE
!****************************************************************************
CONTAINS
!============================================================================
SUBROUTINE proc_advection(i0,j0,k0,cellbdy0,N,ph,nc) 
    IMPLICIT NONE
    INTEGER*4::i0,j0,k0,hit,true,nc
    INTEGER*4::N,m,i,phcell(3),face(3),idxt(1)
    REAL*8:: cellbdy0(2,3), ph(iNprop,N)
    REAL*8::dtremain,dtused
    REAL*8,ALLOCATABLE:: cellbdy(:,:),vel(:),ds(:)
    !phcell：聲子所在網格
    !-----------------------------------
    ALLOCATE( cellbdy(2,3),vel(3),ds(3) )
    !-----------------------------------
    DO m=1,N 
        cellbdy=cellbdy0 !SO此處的cellbdy會等於外界的cellbdy (in main code裡的副程式advance)
        phcell(1)=i0
        phcell(2)=j0
        phcell(3)=k0
        vel(3)=ph(7,m)*DSQRT(1d0-ph(4,m)**2)  ! Vg*sin(theta)
        vel(1)=ph(7,m)*ph(4,m) ! Vg*cos(theta)
        vel(2)=vel(3)*DCOS(ph(5,m))  ! Vg*sin(theta)*cos(phi)
        vel(3)=vel(3)*DSIN(ph(5,m))  ! Vg*sin(theta)*sin(phi)

        IF (nc.eq.1) THEN !so nc=1表示是聲子運動步驟，nc=-1表示是熱流控制步驟??
            dtremain = dt
        ELSE 
            dtremain = dtheat(m)
        ENDIF
        
        !-------------------------------------------
        DO WHILE (dtremain.gt.0d0)	
	        DO i=1,3
	            IF (vel(i).gt.0d0) THEN
	                ds(i)=cellbdy(2,i)
			        face(i)=1
                ELSE
		            ds(i)=cellbdy(1,i)
			        face(i)=-1
                ENDIF
	        ENDDO
	        ds=DABS((ds-ph(1:3,m))/vel) !現在ds變成移動到碰上網格邊界所需的時間
	        idxt=MINLOC(ds) !傳回三個方向中，最先碰上哪個方向的網格邊界，idxt=1 or 2 or 3
	        hit=idxt(1)
            dtused=MIN(ds(hit),dtremain) !比較剩餘時間是否足夠聲子碰上網格邊界
            !-------------------------------------------
            true = 0
	        ph(1:3,m)=ph(1:3,m)+dtused*vel ! movement

            CALL proc_intrinsicscattering(phcell,ph(1:iNprop,m),dtused,true) !判斷會不會發生本質散射，若有散射，會改變聲子性質
	        IF (true.eq.1) THEN !true=1表示有本質散射
	            vel(3)=ph(7,m)*DSQRT(1d0-ph(4,m)**2)
	            vel(1)=ph(7,m)*ph(4,m)
                vel(2)=vel(3)*DCOS(ph(5,m))
		        vel(3)=vel(3)*DSIN(ph(5,m))
            ENDIF
            !-------------------------------------------
            IF (dtused.ge.dtremain) THEN  
	            ! IF (ds(hit).ge.dtremain) THEN  
	            ! if (ds(hit).eq.dtremain) then
	            ! write(999,*) '1'
	            ! endif
          
		        dtremain=0
	        ELSE
	            dtremain = dtremain-dtused
		        true=0
		        
		        IF (face(hit)*vel(hit).gt.0) THEN !face(hit)*vel(hit)一定大於零阿???不一定!!因為發生過散射了，可能vel(hit)改變了!所以此行是在判斷經散射後是否還會穿透網格邊界
		            CALL proc_outdomain(phcell,cellbdy,hit,face(hit),ph(1:iNprop,m),dtremain,true) 
		            !-----check and handle if this phonon hits the computational domain
		            IF (true.eq.0) CALL proc_transmissivity(phcell,cellbdy,hit,face(hit),ph(1:iNprop,m),vel,true)
	            ENDIF
	            
		        IF (true.eq.1) THEN !若有發生邊界反射，或介面鏡/亂穿透or反射則true為1，且在proc_outdomain與proc_transmissivity只決定方向和群速，而沒決定速度分量
	                vel(3)=ph(7,m)*DSQRT(1d0-ph(4,m)**2)
	                vel(1)=ph(7,m)*ph(4,m)
                    vel(2)=vel(3)*DCOS(ph(5,m))
	                vel(3)=vel(3)*DSIN(ph(5,m))
		        ENDIF
	        ENDIF
        ENDDO
    ENDDO
 
    DEALLOCATE( cellbdy,vel,ds )

END SUBROUTINE proc_advection
!============================================================================
SUBROUTINE proc_intrinsicscattering(phcell,phm,dt1,true) !判斷會不會發生本質散射，若有散射，會改變聲子性質
IMPLICIT NONE
INTEGER*4::phcell(3),true !phcell：該聲子所在網格，phm：該聲子的所有性質，dt1：此段運動所需時間，true：將結果存回true，1表示有本質散射
REAL*8::phm(iNprop),dt1,prob
REAL*8,ALLOCATABLE:: rannum1(:)

ALLOCATE( rannum1(3) )
CALL random_number(rannum1)

prob=1d0-DEXP(-dt1*phm(7)/MFP(phcell(1),phcell(2),phcell(3))) !在dt1時間內散射的機率

IF (rannum1(1).le.prob) THEN
   phm(4)=2D0*rannum1(2)-1D0
   phm(5)=M_PI_2*rannum1(3)
   dEdiff(phcell(1),phcell(2),phcell(3))=dEdiff(phcell(1),phcell(2),phcell(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
   !散射後聲子束能量為當前網格溫度下的能量，所以要計算原本聲子能量與散射後生子能量的差值，未來做能量守恆用途
   phm(6)=dEunit(phcell(1),phcell(2),phcell(3)) !本質散射後，聲子能量為當前網格平均聲子能量
   phm(7)=dVunit(phcell(1),phcell(2),phcell(3)) !本質散射後，群速為當前網格平均聲子群速
   phm(8)=iCmat(phcell(1),phcell(2),phcell(3)) !所以散射後，聲子能量材料屬性也變成當前網格的材料
   true=1
ENDIF

DEALLOCATE( rannum1 )
END SUBROUTINE proc_intrinsicscattering
!============================================================================
SUBROUTINE proc_transmissivity(phcell,cellbdy,hit,face0,phm,vel,true) !判斷為鏡/亂穿射or反射，並得到穿/反射後的聲子性質
IMPLICIT NONE
INTEGER*4::phcell(3),neighbor(3),hit,face0,true,i
REAL*8::phm(iNprop),cellbdy(2,3),vel(3) 
REAL*8::tau12,ratio,dcosth2,dsinth2,rho1,rho2,neighborE,neighborV
REAL*8::tau21
REAL*8,ALLOCATABLE:: rannum1(:)
    !phcell：該聲子目前所在網格，cellbdy：所在網格的6個截面在模擬區域的位置，hit：聲子要碰撞的面的方向(1/2/3)，
    !face0：聲子要碰撞的面(1代表正hit方向，-1代表負hit方向)，phm：該聲子的所有性質，vel：速度分量，
    !true：最後傳回true(1表此聲子遇上絕熱邊界反射了，-1表遇上週期性邊界，被移到另一邊邊界，或是穿透到domain外)

    !------------------
    neighbor=phcell
    neighbor(hit)=neighbor(hit)+face0 !所以neighbor變成聲子即將移動過去的網格(OR SAY 新網格)
    
    IF (iCmat(neighbor(1),neighbor(2),neighbor(3)).eq.iCmat(phcell(1),phcell(2),phcell(3))) THEN !兩個網格材料相同，無折射現象
        IF (hit.eq.1) CALL Compute_qflux(phcell,phm,vel) !若聲子是穿過熱流方向的截面，則要先判斷是不是有穿過domain中間截面，因為要記錄中間截面處的熱通量
	    tau12 = 1d0 !穿透率100%
	    phcell = neighbor !聲子所在網格改變為新網格
	    cellbdy(1:2,hit)=cellbdy(1:2,hit)+dLclen(hit)*DBLE(face0) !聲子所在網格的兩端截面(聲子穿透方向)在domain位置相應改變
    
    ELSE ! hit the interface
        ALLOCATE( rannum1(2) )
	    CALL random_number(rannum1)
	    CALL proc_Energy(iCmat(neighbor(1),neighbor(2),neighbor(3)),dTemp(phcell(1),phcell(2),phcell(3)),neighborE) !用當前網格溫度，計算移動目的網格的單位體積能量(U2)
        CALL Etable(iCmat(neighbor(1),neighbor(2),neighbor(3)),4,neighborE,neighborV) !用上一步驟計算的能量求得，當前網格溫度但材料為目的網格材料時的群速(v2)
        !----------------------------------------------------------------------------
	    IF (rannum1(1).le.DPP(hit)) THEN ! 鏡穿/反射
	        
	        ratio=(dEcell(phcell(1),phcell(2),phcell(3))*dVunit(phcell(1),phcell(2),phcell(3)))/(neighborE*neighborV) !(U1*v1)/(U2*v2)
            dcosth2=vel(hit)/phm(7) !穿透方向的速度分量/速度，即入射角的cos值
	        dsinth2=DSQRT((1d0-dcosth2**2)*ratio) !折射角的sin值
	        tau12=0 !若全反射，tau12不會被更改，因此後續會被判斷為鏡反射
	        
	        IF (dsinth2.lt.1d0) THEN !若dsinth2 > 1，數學上無意義，物理上則為全反射
	        
	            rho1=rho(iCmat(phcell(1),phcell(2),phcell(3))) !聲子所在網格的材料密度
		        rho2=rho(iCmat(neighbor(1),neighbor(2),neighbor(3))) !目的網格的材料密度
		        dcosth2=DSQRT(1d0-dsinth2**2) !此時dcosth2變成折射角的cos值
		        tau12=(rho2*neighborV*dcosth2)/DABS(rho1*vel(hit))
		        tau12=1d0-((1d0-tau12)/(1d0+tau12))**2  !鏡穿透率
	        
	        ENDIF
	        
	        IF (rannum1(2).lt.tau12) THEN !! 鏡穿透 (IAMM)
	            IF (hit.eq.1) CALL Compute_qflux(phcell,phm,vel)
	            CALL Snells( ratio,phm,vel,hit ) !決定折射(穿透)後的移動方向
	            cellbdy(1:2,hit)=cellbdy(1:2,hit)+dLclen(hit)*DBLE(face0) 
		        !dEdiff(neighbor(1),neighbor(2),neighbor(3))=dEdiff(neighbor(1),neighbor(2),neighbor(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
                !phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
                phm(7)=dVunit(neighbor(1),neighbor(2),neighbor(3)) !聲子只會改變群速!
		        phcell = neighbor !聲子所屬網格改變
	        ELSE !! 鏡反射
	            IF (hit.eq.1) THEN
		            phm(4)=-phm(4)
		        ELSE IF (hit.eq.2) THEN
		            phm(5)=M_PI-phm(5)
			        IF (phm(5).lt.0) phm(5)=phm(5)+M_PI_2
		        ELSE ! hit.eq.3
		  	        phm(5)=M_PI_2-phm(5)
		        ENDIF
		        !鏡反射除了移動方向改變，其他全部都不會變
	        ENDIF
	    !----------------------------------------------------------------------------
	    ELSE ! 亂穿/反射
            tau12=(neighborE*neighborV)/(dEcell(phcell(1),phcell(2),phcell(3))*dVunit(phcell(1),phcell(2),phcell(3))+neighborE*neighborV) !亂射時的穿透率
	        IF (rannum1(2).le.tau12) THEN !! 亂穿透 (DAMM)		            
	            IF (hit.eq.1) CALL Compute_qflux(phcell,phm,vel)
	            cellbdy(1:2,hit)=cellbdy(1:2,hit)+dLclen(hit)*DBLE(face0) !聲子所屬網格改變
		        !dEdiff(neighbor(1),neighbor(2),neighbor(3))=dEdiff(neighbor(1),neighbor(2),neighbor(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
		        CALL diffuseB( phm,hit,face0,1 ) !決定亂反/穿射後的方向
                phm(7)=dVunit(neighbor(1),neighbor(2),neighbor(3)) !只會改變群速!!!!!
		        phcell = neighbor !聲子所屬網格改變
            !----------------------------------------------------------------------------		  
	        ELSE !! 亂反射 (DAMM)
	            !dEdiff(phcell(1),phcell(2),phcell(3))=dEdiff(phcell(1),phcell(2),phcell(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
		        CALL diffuseB( phm,hit,face0,-1 ) !決定亂反/穿射後的方向
                phm(7)=dVunit(phcell(1),phcell(2),phcell(3)) !只會改變群速!!!!!
	        ENDIF
	    ENDIF
	    
	    !無論是鏡亂/反射穿射，聲子能量通通都不會變，而除了鏡反射只改變方向外，其他三種情況會改變方向與群速
!----------------------------------------------------------------------------
        true = 1
	    DEALLOCATE( rannum1 )
    ENDIF

END SUBROUTINE proc_transmissivity
!============================================================================
SUBROUTINE proc_outdomain(phcell,cellbdy,hit,face0,phm,dtremain,true) !判斷聲子是否會離開domain，若離開的話，依照bc不同，改變性質
    IMPLICIT NONE
    INTEGER*4::phcell(3),hit,face0,true,judge,j,k
    REAL*8::phm(iNprop),cellbdy(2,3),dtremain,tmp,rannum
    !phcell：該聲子目前所在網格，cellbdy：所在網格的6個截面在模擬區域的位置，hit：聲子要碰撞的面的方向(1/2/3)，
    !face0：生子要碰撞的面(1代表正的hit方向，-1代表負的hit方向)，phm：該聲子的所有性質，dtremail：聲子剩餘運動時間，
    !true：最後傳回true(1表此聲子遇上絕熱邊界反射了，-1表遇上週期性邊界，被移到另一邊邊界，或是穿透到domain外)

    CALL random_number(rannum)
    judge=0
    IF (phcell(hit).eq.iNcell(hit).and.face0.gt.0) judge=1 !此網格在模擬區域(若是熱流方向，就是出口邊界)邊界上，且聲子正往邊界移動
    IF (phcell(hit).eq.1.and.face0.lt.0) judge = -1 !此網格在模擬區域邊界上(若是熱流方向，就是入口邊界)邊界上，且聲子正往邊界移動

    IF (judge.ne.0) THEN !若等於0就表示不會穿過介面
        ! option=1: partially specularly and partially diffusely reflected
        ! option=2: periodic boundary condition
        ! option=3: leaving and being saved for heat control
        !-----------------------------------------------------------------------------------------------------------------
        IF (option(hit).eq.1) THEN !若是反射(絕熱)
            IF (rannum.le.DPPB(hit)) THEN ! 鏡反射 
	            IF (hit.eq.1) THEN 
		            phm(4)=-phm(4)
	            ELSE IF (hit.eq.2) THEN
		            phm(5)=M_PI-phm(5)
		            IF (phm(5).lt.0) phm(5)=phm(5)+M_PI_2
	            ELSE ! hit.eq.3
		            phm(5)=M_PI_2-phm(5)
		        ENDIF 
		        !邊界鏡反射只改變方向，其他不改變???
	        ELSE !亂反射
                CALL diffuseB( phm,hit,judge,-1 ) !!決定亂反/鏡射後的方向
		        dEdiff(phcell(1),phcell(2),phcell(3))=dEdiff(phcell(1),phcell(2),phcell(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3)) !此時phm(6)是還未散射前的聲子能量，dEunit(phcell(1),phcell(2),phcell(3))是散射後的聲子能量(為該網格平均聲子能量)
                phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
                phm(7)=dVunit(phcell(1),phcell(2),phcell(3))
		        phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
		        !邊界亂反射除了改變方向，聲子能量、群速、材料屬性也會改變
	        ENDIF 
	        
	        true=1
        !-----------------------------------------------------------------------------------------------------------------
        ELSE IF (option(hit).eq.2) THEN
        !!! The first and the last cells must belong to the same material.
            IF (judge.eq.1) THEN !1-1-2
	            phm(hit)=0d0
		        phcell(hit)=1
		        cellbdy(1,hit)=0
		        cellbdy(2,hit)=dLclen(hit)
	        ELSE IF (judge.eq.-1) THEN
	            phm(hit)=dLdomain(hit)
	            phcell(hit)=iNcell(hit)
		        cellbdy(2,hit)=dLdomain(hit)
		        cellbdy(1,hit)=dLdomain(hit)-dLclen(hit)
	        ENDIF !1-1-2
	        
	        IF (true.eq.0) true=-1
!-----------------------------------------------------------------------------------------------------------------
        ELSE IF (option(hit).eq.3) THEN ! possible only if hit=1
	        j=phcell(2)
	        k=phcell(3)
	        IF (WAY_DIR.eq.1) THEN !WAY_DIR=1為週期入射法，2為亂數入設法
	            IF (judge.eq.1) THEN !judge=1表示從正邊界離開(2面)
	                mlost(j,k,1)=mlost(j,k,1)+1 !從2面離開的聲子會被1面的相同位置網格使用，所以直接記成1，代表這是給1面的邊界網格使用的
			        IF (mlost(j,k,1).gt.iNmakeup) mlost(j,k,1)=1
		            dPpool(1,mlost(j,k,1),j,k,1)=dtremain
		            dPpool(2:5,mlost(j,k,1),j,k,1)=phm(2:5)
			        dPpool(6,mlost(j,k,1),j,k,1)=phm(8)
	            ELSE 
	                mlost(j,k,2)=mlost(j,k,2)+1
			        IF (mlost(j,k,2).gt.iNmakeup) mlost(j,k,2)=1
		            dPpool(1,mlost(j,k,2),j,k,2)=dtremain
		            dPpool(2:5,mlost(j,k,2),j,k,2)=phm(2:5)
			        dPpool(6,mlost(j,k,2),j,k,2)=phm(8)
	            ENDIF
	        ENDIF !1-1-3
	        
	        IF (WAY_HEAT.eq.1) THEN !WAY_HEAT=1為定熱通量，2為固定邊界溫度
	            IF (judge.eq.+1) dElost(j,k,2)=dElost(j,k,2)+phm(6)
	            IF (judge.eq.-1) dElost(j,k,1)=dElost(j,k,1)+phm(6)
	        ENDIF !1-1-4
	        
	        phm(6)=0
	        dtremain=0
	        true=-1
        ENDIF !1-1
    ENDIF !1

END SUBROUTINE proc_outdomain
!============================================================================
SUBROUTINE proc_createdelete !用來進行網格能量守恆時在網格內增或減聲子(內已包含重新整理)
    IMPLICIT NONE
    INTEGER*4::i,j,k,m,s,bg,ed,true
    REAL*8::random1,tmp
    INTEGER*4,ALLOCATABLE::nadd(:,:,:)
    REAL*8,ALLOCATABLE::rannum(:,:),newphn(:,:)
    ALLOCATE( nadd(iNcell(1),iNcell(2),iNcell(3)) )
    nadd=0
    true=0
    DO k=1,iNcell(3)
        DO j=1,iNcell(2)
            DO i=1,iNcell(1)
                tmp=0.5d0*dEunit(i,j,k)
                IF (dEdiff(i,j,k).ge.tmp) THEN
                    nadd(i,j,k) = INT(dEdiff(i,j,k)/dEunit(i,j,k)+0.5d0)
	                dEdiff(i,j,k)=dEdiff(i,j,k)-nadd(i,j,k)*dEunit(i,j,k)
	                true=1
                ELSE IF (dEdiff(i,j,k).lt.-tmp) THEN
                    DO WHILE (dEdiff(i,j,k).lt.-tmp)
	                    CALL random_number(random1)
		                m=MIN(INT(iNbgcell(i,j,k)+random1*iNnumcell(i,j,k)+1),iNbgcell(i,j,k)+iNnumcell(i,j,k))
		                IF (phn(6,m).gt.0) THEN
		                    dEdiff(i,j,k)=dEdiff(i,j,k)+phn(6,m)
		                    phn(6,m)=0
		                    nadd(i,j,k)=nadd(i,j,k)-1
	                    ENDIF
	                ENDDO
	                true=1
                ENDIF
            ENDDO
        ENDDO
    ENDDO
    !-----------------------------
    IF (true.eq.1) THEN
        newiNph=iNph+SUM(nadd)
        ALLOCATE( newphn(iNprop,newiNph) )
        s=0
        DO k=1,iNcell(3)
            DO j=1,iNcell(2)
                DO i=1,iNcell(1)
                    bg=iNbgcell(i,j,k)+1
                    ed=iNbgcell(i,j,k)+iNnumcell(i,j,k)
                    IF (nadd(i,j,k).eq.0) THEN    
	                    newphn(:,s+1:s+iNnumcell(i,j,k))=phn(:,bg:ed)
	                    s=s+iNnumcell(i,j,k)	
	                ELSE IF (nadd(i,j,k).gt.0) THEN
	                    newphn(:,s+1:s+iNnumcell(i,j,k))=phn(:,bg:ed)
	                    s=s+iNnumcell(i,j,k)
                        ALLOCATE( rannum(nadd(i,j,k),5) )
	                    CALL random_number(rannum)
	                    newphn(1,s+1:s+nadd(i,j,k))=(DBLE(i-1)+rannum(:,1))*dLclen(1)
		                newphn(2,s+1:s+nadd(i,j,k))=(DBLE(j-1)+rannum(:,2))*dLclen(2)
		                newphn(3,s+1:s+nadd(i,j,k))=(DBLE(k-1)+rannum(:,3))*dLclen(3)
	                    newphn(4,s+1:s+nadd(i,j,k))=2d0*rannum(:,4)-1d0
		                newphn(5,s+1:s+nadd(i,j,k))=M_PI_2*rannum(:,5)
	                    newphn(6,s+1:s+nadd(i,j,k))=dEunit(i,j,k)
                        newphn(7,s+1:s+nadd(i,j,k))=dVunit(i,j,k)
		                newphn(8,s+1:s+nadd(i,j,k))=iCmat(i,j,k)
	                    s=s+nadd(i,j,k)
		                DEALLOCATE( rannum )
                    ELSE IF (nadd(i,j,k).lt.0) THEN
                        DO m=bg,ed
	                        IF (phn(6,m).gt.0) THEN
		                        s=s+1
	                            newphn(:,s)=phn(:,m)
	                        ENDIF
	                    ENDDO
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        DEALLOCATE( phn )
        ALLOCATE( phn(iNprop,newiNph) )
        iNph=newiNph
        phn = newphn
        DEALLOCATE( newphn )
!----------------------------
        iNbgcell=0
        m=0
        DO k=1,iNcell(3)
            DO j=1,iNcell(2)
                DO i=1,iNcell(1)
                    iNbgcell(i,j,k)=m
                    iNnumcell(i,j,k)=iNnumcell(i,j,k)+nadd(i,j,k)
                    m=m+iNnumcell(i,j,k)
                ENDDO
            ENDDO
        ENDDO
    ENDIF
 
    DEALLOCATE( nadd )
END SUBROUTINE proc_createdelete
!============================================================================
SUBROUTINE Compute_qflux(phcell,phm,vel) !計算中間截面處的熱通量
    IMPLICIT NONE
    INTEGER*4:: phcell(3)
    REAL*8:: phm(7),vel(3)

    IF (phcell(1).eq.iNcell(1)/2.and.vel(1).gt.0d0) THEN !聲子往正方向移動通過中間截面
        qflow(phcell(2),phcell(3))=qflow(phcell(2),phcell(3))+phm(6)
    ELSE IF (phcell(1).eq.iNcell(1)/2+1.and.vel(1).lt.0d0) THEN !聲子往負方向移動穿過中間截面
        qflow(phcell(2),phcell(3))=qflow(phcell(2),phcell(3))-phm(6)
    ENDIF

END SUBROUTINE Compute_qflux
!============================================================================
END MODULE mod_ADVANCE