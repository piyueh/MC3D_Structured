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
    !phcell�G�n�l�Ҧb����
    !-----------------------------------
    ALLOCATE( cellbdy(2,3),vel(3),ds(3) )
    !-----------------------------------
    DO m=1,N 
        cellbdy=cellbdy0 !SO���B��cellbdy�|����~�ɪ�cellbdy (in main code�̪��Ƶ{��advance)
        phcell(1)=i0
        phcell(2)=j0
        phcell(3)=k0
        vel(3)=ph(7,m)*DSQRT(1d0-ph(4,m)**2)  ! Vg*sin(theta)
        vel(1)=ph(7,m)*ph(4,m) ! Vg*cos(theta)
        vel(2)=vel(3)*DCOS(ph(5,m))  ! Vg*sin(theta)*cos(phi)
        vel(3)=vel(3)*DSIN(ph(5,m))  ! Vg*sin(theta)*sin(phi)

        IF (nc.eq.1) THEN !so nc=1��ܬO�n�l�B�ʨB�J�Anc=-1��ܬO���y����B�J??
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
	        ds=DABS((ds-ph(1:3,m))/vel) !�{�bds�ܦ����ʨ�I�W������ɩһݪ��ɶ�
	        idxt=MINLOC(ds) !�Ǧ^�T�Ӥ�V���A�̥��I�W���Ӥ�V��������ɡAidxt=1 or 2 or 3
	        hit=idxt(1)
            dtused=MIN(ds(hit),dtremain) !����Ѿl�ɶ��O�_�����n�l�I�W�������
            !-------------------------------------------
            true = 0
	        ph(1:3,m)=ph(1:3,m)+dtused*vel ! movement

            CALL proc_intrinsicscattering(phcell,ph(1:iNprop,m),dtused,true) !�P�_�|���|�o�ͥ��贲�g�A�Y�����g�A�|�����n�l�ʽ�
	        IF (true.eq.1) THEN !true=1��ܦ����贲�g
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
		        
		        IF (face(hit)*vel(hit).gt.0) THEN !face(hit)*vel(hit)�@�w�j��s��???���@�w!!�]���o�͹L���g�F�A�i��vel(hit)���ܤF!�ҥH����O�b�P�_�g���g��O�_�ٷ|��z�������
		            CALL proc_outdomain(phcell,cellbdy,hit,face(hit),ph(1:iNprop,m),dtremain,true) 
		            !-----check and handle if this phonon hits the computational domain
		            IF (true.eq.0) CALL proc_transmissivity(phcell,cellbdy,hit,face(hit),ph(1:iNprop,m),vel,true)
	            ENDIF
	            
		        IF (true.eq.1) THEN !�Y���o����ɤϮg�A�Τ�����/�ì�zor�Ϯg�htrue��1�A�B�bproc_outdomain�Pproc_transmissivity�u�M�w��V�M�s�t�A�ӨS�M�w�t�פ��q
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
SUBROUTINE proc_intrinsicscattering(phcell,phm,dt1,true) !�P�_�|���|�o�ͥ��贲�g�A�Y�����g�A�|�����n�l�ʽ�
IMPLICIT NONE
INTEGER*4::phcell(3),true !phcell�G���n�l�Ҧb����Aphm�G���n�l���Ҧ��ʽ�Adt1�G���q�B�ʩһݮɶ��Atrue�G�N���G�s�^true�A1��ܦ����贲�g
REAL*8::phm(iNprop),dt1,prob
REAL*8,ALLOCATABLE:: rannum1(:)

ALLOCATE( rannum1(3) )
CALL random_number(rannum1)

prob=1d0-DEXP(-dt1*phm(7)/MFP(phcell(1),phcell(2),phcell(3))) !�bdt1�ɶ������g�����v

IF (rannum1(1).le.prob) THEN
   phm(4)=2D0*rannum1(2)-1D0
   phm(5)=M_PI_2*rannum1(3)
   dEdiff(phcell(1),phcell(2),phcell(3))=dEdiff(phcell(1),phcell(2),phcell(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
   !���g���n�l����q����e����ūפU����q�A�ҥH�n�p��쥻�n�l��q�P���g��ͤl��q���t�ȡA���Ӱ���q�u��γ~
   phm(6)=dEunit(phcell(1),phcell(2),phcell(3)) !���贲�g��A�n�l��q����e���業���n�l��q
   phm(7)=dVunit(phcell(1),phcell(2),phcell(3)) !���贲�g��A�s�t����e���業���n�l�s�t
   phm(8)=iCmat(phcell(1),phcell(2),phcell(3)) !�ҥH���g��A�n�l��q�����ݩʤ]�ܦ���e���檺����
   true=1
ENDIF

DEALLOCATE( rannum1 )
END SUBROUTINE proc_intrinsicscattering
!============================================================================
SUBROUTINE proc_transmissivity(phcell,cellbdy,hit,face0,phm,vel,true) !�P�_����/�ì�gor�Ϯg�A�ño���/�Ϯg�᪺�n�l�ʽ�
IMPLICIT NONE
INTEGER*4::phcell(3),neighbor(3),hit,face0,true,i
REAL*8::phm(iNprop),cellbdy(2,3),vel(3) 
REAL*8::tau12,ratio,dcosth2,dsinth2,rho1,rho2,neighborE,neighborV
REAL*8::tau21
REAL*8,ALLOCATABLE:: rannum1(:)
    !phcell�G���n�l�ثe�Ҧb����Acellbdy�G�Ҧb���檺6�ӺI���b�����ϰ쪺��m�Ahit�G�n�l�n�I����������V(1/2/3)�A
    !face0�G�n�l�n�I������(1�N��hit��V�A-1�N��thit��V)�Aphm�G���n�l���Ҧ��ʽ�Avel�G�t�פ��q�A
    !true�G�̫�Ǧ^true(1���n�l�J�W������ɤϮg�F�A-1��J�W�g������ɡA�Q����t�@����ɡA�άO��z��domain�~)

    !------------------
    neighbor=phcell
    neighbor(hit)=neighbor(hit)+face0 !�ҥHneighbor�ܦ��n�l�Y�N���ʹL�h������(OR SAY �s����)
    
    IF (iCmat(neighbor(1),neighbor(2),neighbor(3)).eq.iCmat(phcell(1),phcell(2),phcell(3))) THEN !��Ӻ�����ƬۦP�A�L��g�{�H
        IF (hit.eq.1) CALL Compute_qflux(phcell,phm,vel) !�Y�n�l�O��L���y��V���I���A�h�n���P�_�O���O����Ldomain�����I���A�]���n�O�������I���B�����q�q
	    tau12 = 1d0 !��z�v100%
	    phcell = neighbor !�n�l�Ҧb������ܬ��s����
	    cellbdy(1:2,hit)=cellbdy(1:2,hit)+dLclen(hit)*DBLE(face0) !�n�l�Ҧb���檺��ݺI��(�n�l��z��V)�bdomain��m��������
    
    ELSE ! hit the interface
        ALLOCATE( rannum1(2) )
	    CALL random_number(rannum1)
	    CALL proc_Energy(iCmat(neighbor(1),neighbor(2),neighbor(3)),dTemp(phcell(1),phcell(2),phcell(3)),neighborE) !�η�e����ūסA�p�Ⲿ�ʥت����檺�����n��q(U2)
        CALL Etable(iCmat(neighbor(1),neighbor(2),neighbor(3)),4,neighborE,neighborV) !�ΤW�@�B�J�p�⪺��q�D�o�A��e����ūצ����Ƭ��ت�������Ʈɪ��s�t(v2)
        !----------------------------------------------------------------------------
	    IF (rannum1(1).le.DPP(hit)) THEN ! ���/�Ϯg
	        
	        ratio=(dEcell(phcell(1),phcell(2),phcell(3))*dVunit(phcell(1),phcell(2),phcell(3)))/(neighborE*neighborV) !(U1*v1)/(U2*v2)
            dcosth2=vel(hit)/phm(7) !��z��V���t�פ��q/�t�סA�Y�J�g����cos��
	        dsinth2=DSQRT((1d0-dcosth2**2)*ratio) !��g����sin��
	        tau12=0 !�Y���Ϯg�Atau12���|�Q���A�]������|�Q�P�_����Ϯg
	        
	        IF (dsinth2.lt.1d0) THEN !�Ydsinth2 > 1�A�ƾǤW�L�N�q�A���z�W�h�����Ϯg
	        
	            rho1=rho(iCmat(phcell(1),phcell(2),phcell(3))) !�n�l�Ҧb���檺���ƱK��
		        rho2=rho(iCmat(neighbor(1),neighbor(2),neighbor(3))) !�ت����檺���ƱK��
		        dcosth2=DSQRT(1d0-dsinth2**2) !����dcosth2�ܦ���g����cos��
		        tau12=(rho2*neighborV*dcosth2)/DABS(rho1*vel(hit))
		        tau12=1d0-((1d0-tau12)/(1d0+tau12))**2  !���z�v
	        
	        ENDIF
	        
	        IF (rannum1(2).lt.tau12) THEN !! ���z (IAMM)
	            IF (hit.eq.1) CALL Compute_qflux(phcell,phm,vel)
	            CALL Snells( ratio,phm,vel,hit ) !�M�w��g(��z)�᪺���ʤ�V
	            cellbdy(1:2,hit)=cellbdy(1:2,hit)+dLclen(hit)*DBLE(face0) 
		        !dEdiff(neighbor(1),neighbor(2),neighbor(3))=dEdiff(neighbor(1),neighbor(2),neighbor(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
                !phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
                phm(7)=dVunit(neighbor(1),neighbor(2),neighbor(3)) !�n�l�u�|���ܸs�t!
		        phcell = neighbor !�n�l���ݺ������
	        ELSE !! ��Ϯg
	            IF (hit.eq.1) THEN
		            phm(4)=-phm(4)
		        ELSE IF (hit.eq.2) THEN
		            phm(5)=M_PI-phm(5)
			        IF (phm(5).lt.0) phm(5)=phm(5)+M_PI_2
		        ELSE ! hit.eq.3
		  	        phm(5)=M_PI_2-phm(5)
		        ENDIF
		        !��Ϯg���F���ʤ�V���ܡA��L���������|��
	        ENDIF
	    !----------------------------------------------------------------------------
	    ELSE ! �ì�/�Ϯg
            tau12=(neighborE*neighborV)/(dEcell(phcell(1),phcell(2),phcell(3))*dVunit(phcell(1),phcell(2),phcell(3))+neighborE*neighborV) !�îg�ɪ���z�v
	        IF (rannum1(2).le.tau12) THEN !! �ì�z (DAMM)		            
	            IF (hit.eq.1) CALL Compute_qflux(phcell,phm,vel)
	            cellbdy(1:2,hit)=cellbdy(1:2,hit)+dLclen(hit)*DBLE(face0) !�n�l���ݺ������
		        !dEdiff(neighbor(1),neighbor(2),neighbor(3))=dEdiff(neighbor(1),neighbor(2),neighbor(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
		        CALL diffuseB( phm,hit,face0,1 ) !�M�w�ä�/��g�᪺��V
                phm(7)=dVunit(neighbor(1),neighbor(2),neighbor(3)) !�u�|���ܸs�t!!!!!
		        phcell = neighbor !�n�l���ݺ������
            !----------------------------------------------------------------------------		  
	        ELSE !! �äϮg (DAMM)
	            !dEdiff(phcell(1),phcell(2),phcell(3))=dEdiff(phcell(1),phcell(2),phcell(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
		        !phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
		        CALL diffuseB( phm,hit,face0,-1 ) !�M�w�ä�/��g�᪺��V
                phm(7)=dVunit(phcell(1),phcell(2),phcell(3)) !�u�|���ܸs�t!!!!!
	        ENDIF
	    ENDIF
	    
	    !�L�׬O���/�Ϯg��g�A�n�l��q�q�q�����|�ܡA�Ӱ��F��Ϯg�u���ܤ�V�~�A��L�T�ر��p�|���ܤ�V�P�s�t
!----------------------------------------------------------------------------
        true = 1
	    DEALLOCATE( rannum1 )
    ENDIF

END SUBROUTINE proc_transmissivity
!============================================================================
SUBROUTINE proc_outdomain(phcell,cellbdy,hit,face0,phm,dtremain,true) !�P�_�n�l�O�_�|���}domain�A�Y���}���ܡA�̷�bc���P�A���ܩʽ�
    IMPLICIT NONE
    INTEGER*4::phcell(3),hit,face0,true,judge,j,k
    REAL*8::phm(iNprop),cellbdy(2,3),dtremain,tmp,rannum
    !phcell�G���n�l�ثe�Ҧb����Acellbdy�G�Ҧb���檺6�ӺI���b�����ϰ쪺��m�Ahit�G�n�l�n�I����������V(1/2/3)�A
    !face0�G�ͤl�n�I������(1�N����hit��V�A-1�N��t��hit��V)�Aphm�G���n�l���Ҧ��ʽ�Adtremail�G�n�l�Ѿl�B�ʮɶ��A
    !true�G�̫�Ǧ^true(1���n�l�J�W������ɤϮg�F�A-1��J�W�g������ɡA�Q����t�@����ɡA�άO��z��domain�~)

    CALL random_number(rannum)
    judge=0
    IF (phcell(hit).eq.iNcell(hit).and.face0.gt.0) judge=1 !������b�����ϰ�(�Y�O���y��V�A�N�O�X�f���)��ɤW�A�B�n�l������ɲ���
    IF (phcell(hit).eq.1.and.face0.lt.0) judge = -1 !������b�����ϰ���ɤW(�Y�O���y��V�A�N�O�J�f���)��ɤW�A�B�n�l������ɲ���

    IF (judge.ne.0) THEN !�Y����0�N��ܤ��|��L����
        ! option=1: partially specularly and partially diffusely reflected
        ! option=2: periodic boundary condition
        ! option=3: leaving and being saved for heat control
        !-----------------------------------------------------------------------------------------------------------------
        IF (option(hit).eq.1) THEN !�Y�O�Ϯg(����)
            IF (rannum.le.DPPB(hit)) THEN ! ��Ϯg 
	            IF (hit.eq.1) THEN 
		            phm(4)=-phm(4)
	            ELSE IF (hit.eq.2) THEN
		            phm(5)=M_PI-phm(5)
		            IF (phm(5).lt.0) phm(5)=phm(5)+M_PI_2
	            ELSE ! hit.eq.3
		            phm(5)=M_PI_2-phm(5)
		        ENDIF 
		        !�����Ϯg�u���ܤ�V�A��L������???
	        ELSE !�äϮg
                CALL diffuseB( phm,hit,judge,-1 ) !!�M�w�ä�/��g�᪺��V
		        dEdiff(phcell(1),phcell(2),phcell(3))=dEdiff(phcell(1),phcell(2),phcell(3))+phm(6)-dEunit(phcell(1),phcell(2),phcell(3)) !����phm(6)�O�٥����g�e���n�l��q�AdEunit(phcell(1),phcell(2),phcell(3))�O���g�᪺�n�l��q(���Ӻ��業���n�l��q)
                phm(6)=dEunit(phcell(1),phcell(2),phcell(3))
                phm(7)=dVunit(phcell(1),phcell(2),phcell(3))
		        phm(8)=iCmat(phcell(1),phcell(2),phcell(3))
		        !��ɶäϮg���F���ܤ�V�A�n�l��q�B�s�t�B�����ݩʤ]�|����
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
	        IF (WAY_DIR.eq.1) THEN !WAY_DIR=1���g���J�g�k�A2���üƤJ�]�k
	            IF (judge.eq.1) THEN !judge=1��ܱq��������}(2��)
	                mlost(j,k,1)=mlost(j,k,1)+1 !�q2�����}���n�l�|�Q1�����ۦP��m����ϥΡA�ҥH�����O��1�A�N��o�O��1������ɺ���ϥΪ�
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
	        
	        IF (WAY_HEAT.eq.1) THEN !WAY_HEAT=1���w���q�q�A2���T�w��ɷū�
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
SUBROUTINE proc_createdelete !�ΨӶi������q�u��ɦb���椺�W�δ��n�l(���w�]�t���s��z)
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
SUBROUTINE Compute_qflux(phcell,phm,vel) !�p�⤤���I���B�����q�q
    IMPLICIT NONE
    INTEGER*4:: phcell(3)
    REAL*8:: phm(7),vel(3)

    IF (phcell(1).eq.iNcell(1)/2.and.vel(1).gt.0d0) THEN !�n�l������V���ʳq�L�����I��
        qflow(phcell(2),phcell(3))=qflow(phcell(2),phcell(3))+phm(6)
    ELSE IF (phcell(1).eq.iNcell(1)/2+1.and.vel(1).lt.0d0) THEN !�n�l���t��V���ʬ�L�����I��
        qflow(phcell(2),phcell(3))=qflow(phcell(2),phcell(3))-phm(6)
    ENDIF

END SUBROUTINE Compute_qflux
!============================================================================
END MODULE mod_ADVANCE