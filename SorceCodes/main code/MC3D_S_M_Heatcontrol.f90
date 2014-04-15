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

IF (WAY_HEAT.eq.1) THEN !�w���q�q

    CALL proc_BC( dTemp(1,1:iNcell(2),1:iNcell(3)),dTemp(iNcell(1),1:iNcell(2),1:iNcell(3)) ) !�o����ɳB�U�����U���n�l�s�t�B�C���n�l��q
   
    IF (sumQ.eq.0.and.abc.eq.1) THEN
        dEheatflux(:,:,1)=(dElost(:,:,1)/SUM(dElost(:,:,1))+dElost(:,:,2)/SUM(dElost(:,:,2)))/2d0
        !dElost(:,:,1)�A1��ɭ��U�������}����q
        !dElost(:,:,2)�A2��ɭ��U�������}����q
        dEheatflux(:,:,1)=dEheatflux0*dEheatflux(:,:,1)
        !dEheatflux0�G�T�w"��q"(�ثe���O���q�q)
        dEheatflux(:,:,2)=dEheatflux(:,:,1)
    ELSE
        dEheatflux(:,:,1)=dEheatflux0*qctrl(:,:)
        dEheatflux(:,:,2)=dEheatflux(:,:,1)
    ENDIF
    dElost(:,:,1)=dElost(:,:,1)+dEheatflux(:,:,1) !�{�bdElost�ܦ��C�Ӻ���n�J�g��"��q"
    dElost(:,:,2)=dElost(:,:,2)-dEheatflux(:,:,2)
ELSE
    dElost=dElost+dEheatflux
ENDIF
!---
tmp=MAX( 0d0,MAXVAL(dElost)/MINVAL(dEinject) ) !dEinject�G�Ӻ���J�g�n�l����q�AMAXVAL(dElost)/MINVAL(dEinject)���̤j�n�J�g���n�l�ƶq
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
      
            iNemit(j,k,s)=0 !iNemit�����bs����j,k����Q�g�J���n�l�ƶq
      
            IF (WAY_DIR.eq.1) THEN
		        npool=mlost(j,k,s)
	            IF (dPpool(6,iNmakeup,j,k,s).gt.0) npool=iNmakeup !�ҥH�Y��s������j,k���檺pool�Ѫ��̤j�s���n�l�s�b�A�hnpool����̤j�n�l�s���A�_�h���ثepool�Ѹ̷̳s���n�s��mlost(j,k,s)
	            
	            !npool���ӭn�Υt�@���۹������檺pool�a???�Y�ثe����bs=1�o���A�h���ӭn��s=2���������檺pool��!!!?????????
	            !Ans�G�bproc_outdomain�{�Ǹ̤w�g�۰ʧ�q��2����j,k�������}���n�l�s�J��1������j,k���檺pool�A�ҥH��1������n�J�g�n�l�ɡA�O�q��1����j,k���檺pool����A�Ӥ��O�q�ĤG��
	            
	        ENDIF
	  
            IF (WAY_DIR.eq.1.and.npool.gt.0) THEN !�g���J�g�k�A�M�w�J�g�n�l�O�t�@���۹������檺pool�Ѫ������n�l�A�Φ������`�J�g�n�l��
	            DO WHILE(dElost(j,k,s).ge.0.5d0*dEinject(j,k,iCmat(ip,j,k),s)) !�@������즹����"�ݭn�J�g����q"�p��"�@�b�n�l��q"
		            CALL random_number(tmp) 
		            m=MIN(INT(tmp*npool)+1,npool) !�qpool�Ѹ��H����@�ӽs��
		            iNemit(j,k,s)=iNemit(j,k,s)+1 !�̫᪺iNmet�N�O�o�Ӻ���J�g���`�n�l��
		            im(iNemit(j,k,s),j,k,s)=m     !�����J�g�n�l�bpool�Ѹ̭쥻���s��
		            dElost(j,k,s)=dElost(j,k,s)-dEinject(j,k,dPpool(6,m,j,k,s),s) !dElost��"�ٻݭn�g�J�����檺��q"
		        ENDDO
            ELSE 
		        iNemit(j,k,s)=MAX(0,INT( (dElost(j,k,s)/dEinject(j,k,iCmat(ip,j,k),s))+0.5d0 ))
		        dElost(j,k,s)=dElost(j,k,s)-dEinject(j,k,iCmat(ip,j,k),s)*DBLE(iNemit(j,k,s))
            ENDIF
        ENDDO
    ENDDO
    
ENDDO

newiNph = SUM( iNemit ) !�bdomain���`�@�s�W���n�l��
ALLOCATE( phnheat(iNprop,iNph+newiNph) )
phnheat(1:iNprop,1:iNph)=phn(1:iNprop,1:iNph)
DEALLOCATE( phn )
 
m=iNph

DO s=1,2
    DO k=1,iNcell(3)
        DO j=1,iNcell(2)
            IF (iNemit(j,k,s).gt.0) THEN !�Y�����榳�J�g�n�l
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
                
                IF (WAY_DIR.eq.1.and.mlost(j,k,s).gt.0) THEN !�ᤩ�J�g�ʽ�
	                phnheat(2:5,m+1:m+iNemit(j,k,s))=dPpool(2:5,im(1:iNemit(j,k,s),j,k,s),j,k,s) !��pool�̳Q�襤���n�l��T�ƻs��phnheat��
	                phnheat(8,m+1:m+iNemit(j,k,s))=dPpool(6,im(1:iNemit(j,k,s),j,k,s),j,k,s)
	                phnheat(6,m+1:m+iNemit(j,k,s))=dEinject(j,k,phnheat(8,m+1:m+iNemit(j,k,s)),s)	
	                phnheat(7,m+1:m+iNemit(j,k,s))=dVinject(j,k,iCmat(ip,j,k),s)
                    dtheat=dPpool(1,im(1:iNemit(j,k,s),j,k,s),j,k,s) !dtheat�G���y����g���ʤJ�g�n�l���ݾl�B�ʮɶ�
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
SUBROUTINE proc_BC(TL,TR) !�p��X�b�T�w��ɷū׬�TL��TR�ɡA��ɳB�U����b���P���Ʈɪ��J�g�n�l�Ƹs�t�B��q
    IMPLICIT NONE
    INTEGER*4:: j,k,m
    REAL*8:: TL(iNcell(2),iNcell(3)),TR(iNcell(2),iNcell(3)),dNum,dUL,dUR

    DO m=1,2
        DO k=1,iNcell(3)
            DO j=1,iNcell(2)
                CALL proc_energy(m,TL(j,k),dUL)
                CALL Etable(m,4,dUL,dVinject(j,k,m,1)) !j,k�N�����Am�����ơA1���X�f�ΤJ�f
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