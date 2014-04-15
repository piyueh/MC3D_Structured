!****************************************************************************
!   HEADING: MC3D HEAT CONTROL MODULE
!   AUTHOR: MJ Huang
!   PURPOSE: This module maintain a constant heat flow rate along one direction.
!   DATE : 2009.7.10
!****************************************************************************
MODULE mod_heatcontrol
USE mod_VARIABLES
USE mod_ROUTINES
USE mod_ADVANCE
IMPLICIT NONE
!****************************************************************************
CONTAINS
!============================================================================
SUBROUTINE proc_heatcontrol(iter,iter0)
IMPLICIT NONE
INTEGER*4::iter,iter0,Nmax,s,m,k,ip,j,i,npool
REAL*8,ALLOCATABLE:: phnheat(:,:),cellheat(:,:),rannum(:,:) 
INTEGER*4,ALLOCATABLE:: im(:,:,:,:)
REAL*8:: tmp

IF (WAY_HEAT.eq.1) THEN !定熱通量

    CALL proc_BC( dTemp(1,1:iNcell(2),1:iNcell(3)),dTemp(iNcell(1),1:iNcell(2),1:iNcell(3)) ) !得到邊界處各網格當下的聲子群速、每顆聲子能量
   
    IF (sumQ.eq.0.and.abc.eq.1) THEN
        dEheatflux(:,:,1)=(dElost(:,:,1)/SUM(dElost(:,:,1))+dElost(:,:,2)/SUM(dElost(:,:,2)))/2d0
        !dElost(:,:,1)，1邊界面各網格離開的能量
        !dElost(:,:,2)，2邊界面各網格離開的能量
        dEheatflux(:,:,1)=dEheatflux0*dEheatflux(:,:,1)
        !dEheatflux0：固定"能量"(目前不是熱通量)
        dEheatflux(:,:,2)=dEheatflux(:,:,1)
    ELSE
        dEheatflux(:,:,1)=dEheatflux0*qctrl(:,:)
        dEheatflux(:,:,2)=dEheatflux(:,:,1)
    ENDIF
    dElost(:,:,1)=dElost(:,:,1)+dEheatflux(:,:,1) !現在dElost變成每個網格要入射的"能量"
    dElost(:,:,2)=dElost(:,:,2)-dEheatflux(:,:,2)
ELSE
    dElost=dElost+dEheatflux
ENDIF
!---
tmp=MAX( 0d0,MAXVAL(dElost)/MINVAL(dEinject) ) !dEinject：該網格入射聲子的能量，MAXVAL(dElost)/MINVAL(dEinject)為最大要入射的聲子數量
i=INT(tmp)+1
ALLOCATE( cellheat(2,3) )
ALLOCATE( im(i,iNcell(2),iNcell(3),2) )
im=0

iNemit=0
DO s=1,2  ! left or right boundary 
   
    IF (s.eq.1) THEN
      ip=1
    ELSE
      ip=iNcell(1)
    ENDIF
   
    DO k=1,iNcell(3)
        DO j=1,iNcell(2)
      
            iNemit(j,k,s)=0 !iNemit紀錄在s面第j,k網格被射入的聲子數量
      
            IF (WAY_DIR.eq.1) THEN
		        npool=mlost(j,k,s)
	            IF (dPpool(6,iNmakeup,j,k,s).gt.0) npool=iNmakeup !所以若第s面的第j,k網格的pool槽的最大編號聲子存在，則npool等於最大聲子編號，否則為目前pool槽裡最新的聲編號mlost(j,k,s)
	            
	            !npool應該要用另一面相對應網格的pool吧???若目前網格在s=1這面，則應該要用s=2那面的網格的pool啊!!!?????????
	            !Ans：在proc_outdomain程序裡已經自動把從第2面第j,k網格離開的聲子存入第1面的第j,k網格的pool，所以第1面網格要入射聲子時，是從第1面第j,k網格的pool選取，而不是從第二面
	            
	        ENDIF
	  
            IF (WAY_DIR.eq.1.and.npool.gt.0) THEN !週期入射法，決定入射聲子是另一面相對應網格的pool槽的哪顆聲子，及此網格總入射聲子數
	            DO WHILE(dElost(j,k,s).ge.0.5d0*dEinject(j,k,iCmat(ip,j,k),s)) !一直執行到此網格"需要入射的能量"小於"一半聲子能量"
		            CALL random_number(tmp) 
		            m=MIN(INT(tmp*npool)+1,npool) !從pool槽裡隨機選一個編號
		            iNemit(j,k,s)=iNemit(j,k,s)+1 !最後的iNmet就是這個網格入射的總聲子數
		            im(iNemit(j,k,s),j,k,s)=m     !紀錄入射聲子在pool槽裡原本的編號
		            dElost(j,k,s)=dElost(j,k,s)-dEinject(j,k,dPpool(6,m,j,k,s),s) !dElost為"還需要射入此網格的能量"
		        ENDDO
            ELSE 
		        iNemit(j,k,s)=MAX(0,INT( (dElost(j,k,s)/dEinject(j,k,iCmat(ip,j,k),s))+0.5d0 ))
		        dElost(j,k,s)=dElost(j,k,s)-dEinject(j,k,iCmat(ip,j,k),s)*DBLE(iNemit(j,k,s))
            ENDIF
        ENDDO
    ENDDO
    
ENDDO

newiNph = SUM( iNemit ) !在domain裡總共新增的聲子數
ALLOCATE( phnheat(iNprop,iNph+newiNph) )
phnheat(1:iNprop,1:iNph)=phn(1:iNprop,1:iNph)
DEALLOCATE( phn )
 
m=iNph

DO s=1,2
    DO k=1,iNcell(3)
        DO j=1,iNcell(2)
            IF (iNemit(j,k,s).gt.0) THEN !若此網格有入射聲子
                ALLOCATE( dtheat(iNemit(j,k,s)) )
                cellheat(1,2)=dLclen(2)*DBLE(j-1)
                cellheat(2,2)=dLclen(2)*DBLE(j)
                cellheat(1,3)=dLclen(3)*DBLE(k-1)
                cellheat(2,3)=dLclen(3)*DBLE(k)
                IF (s.eq.1) THEN
                    cellheat(1,1)=0d0
	                cellheat(2,1)=dLclen(1)
	                ip=1
	                phnheat(1,m+1:m+iNemit(j,k,s))=0d0
                ELSE IF (s.eq.2) THEN
                    cellheat(1,1)=dLdomain(1)-dLclen(1)
                    cellheat(2,1)=dLdomain(1)
	                ip=iNcell(1)
	                phnheat(1,m+1:m+iNemit(j,k,s))=dLdomain(1)
                ENDIF
                !-------------------------------------------------------------------------------------
                ! One way to determine the properties of emitted phonons
                IF (WAY_DIR.eq.1) THEN
                    npool=mlost(j,k,s)
	                IF (dPpool(6,iNmakeup,j,k,s).gt.0) npool=iNmakeup
                ENDIF
                
                IF (WAY_DIR.eq.1.and.mlost(j,k,s).gt.0) THEN !賦予入射性質
	                phnheat(2:5,m+1:m+iNemit(j,k,s))=dPpool(2:5,im(1:iNemit(j,k,s),j,k,s),j,k,s) !把pool裡被選中的聲子資訊複製到phnheat裡
	                phnheat(8,m+1:m+iNemit(j,k,s))=dPpool(6,im(1:iNemit(j,k,s),j,k,s),j,k,s)
	                phnheat(6,m+1:m+iNemit(j,k,s))=dEinject(j,k,phnheat(8,m+1:m+iNemit(j,k,s)),s)	
	                phnheat(7,m+1:m+iNemit(j,k,s))=dVinject(j,k,iCmat(ip,j,k),s)
                    dtheat=dPpool(1,im(1:iNemit(j,k,s),j,k,s),j,k,s) !dtheat：熱流控制週期性入射聲子的殘餘運動時間
                    m=m+iNemit(j,k,s)
                ELSE
                !-------------------------------------------------------------------------------------
                ! Second way to determine the properties of emitted phonons
                    ALLOCATE( rannum(iNemit(j,k,s),5) )
                    CALL random_number( rannum )
	                phnheat(2,m+1:m+iNemit(j,k,s))=cellheat(1,2)+rannum(:,1)*dLclen(2)
                    phnheat(3,m+1:m+iNemit(j,k,s))=cellheat(1,3)+rannum(:,2)*dLclen(3)
                    phnheat(4,m+1:m+iNemit(j,k,s))=DSIGN( DSQRT(rannum(:,3)),DBLE(2-ip)) ! propability proptional to the cosine value when the heat flux is nonzero.
	                phnheat(5,m+1:m+iNemit(j,k,s))=M_PI_2*rannum(:,4) 
	                phnheat(6,m+1:m+iNemit(j,k,s))=dEinject(j,k,iCmat(ip,j,k),s)
	                phnheat(7,m+1:m+iNemit(j,k,s))=dVinject(j,k,iCmat(ip,j,k),s) 
	                phnheat(8,m+1:m+iNemit(j,k,s))=iCmat(ip,j,k)
	                dtheat(:)=dt*rannum(:,5)
	                m=m+iNemit(j,k,s)
	                DEALLOCATE( rannum )
	            ENDIF
!-------------------------------------------------------------------------------------   
                CALL proc_advection(ip,j,k,cellheat,iNemit(j,k,s),phnheat(1:iNprop,m-iNemit(j,k,s)+1:m),-1) 
   
                DEALLOCATE( dtheat )
            ENDIF
        ENDDO
    ENDDO
ENDDO

iNph=iNph+newiNph
ALLOCATE( phn(iNprop,iNph) )
phn=phnheat

DEALLOCATE( phnheat,cellheat,im )
 
END SUBROUTINE proc_heatcontrol
!============================================================================
SUBROUTINE proc_BC(TL,TR) !計算出在固定邊界溫度為TL及TR時，邊界處各網格在不同材料時的入射聲子數群速、能量
    IMPLICIT NONE
    INTEGER*4:: j,k,m
    REAL*8:: TL(iNcell(2),iNcell(3)),TR(iNcell(2),iNcell(3)),dNum,dUL,dUR

    DO m=1,2
        DO k=1,iNcell(3)
            DO j=1,iNcell(2)
                CALL proc_energy(m,TL(j,k),dUL)
                CALL Etable(m,4,dUL,dVinject(j,k,m,1)) !j,k代表網格，m為材料，1為出口或入口
                CALL Etable(m,2,dUL,dNum)
                dEinject(j,k,m,1)=dUL/dNum*bundle(m)
   
                CALL proc_energy(m,TR(j,k),dUR)
                CALL Etable(m,4,dUR,dVinject(j,k,m,2))
                CALL Etable(m,2,dUR,dNum)
                dEinject(j,k,m,2)=dUR/dNum*bundle(m)
            ENDDO
        ENDDO
    ENDDO

END SUBROUTINE proc_BC
!============================================================================
END MODULE mod_heatcontrol