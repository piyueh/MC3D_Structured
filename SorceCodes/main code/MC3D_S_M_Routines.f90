!****************************************************************************
!   HEADING: MC3D ROUTINES MODULE
!   AUTHOR: MJ Huang
!   PURPOSE: This module contains all the routines for transmission and reflection.
!   DATE: 2009.7.10
!****************************************************************************
MODULE mod_ROUTINES
USE mod_VARIABLES
IMPLICIT NONE
!****************************************************************************
CONTAINS
!============================================================================
SUBROUTINE Snells( sinratio,phm,vel,idx ) !決定折射後的方向
    IMPLICIT NONE
    REAL*8:: phm(iNprop),vel(3),sinratio
    REAL*8:: dsinth1,dcosth1,dsinth2,dcosth2,tmp
    INTEGER*4:: idx
    ! Inelastic Acoustic Mismatch Model (IAMM)
    ! sinratio：(sin 折射角)/(sin 入射角)、phm：聲子性質、vel：三個方向的速度分量
    ! idx：穿透面的方向(1/2/3)
    
    IF (idx.eq.3) THEN
    
        tmp=DCOS(phm(5))                               ! phm(5)為23平面上的投影與2軸的夾角
        dsinth1=1d0-(vel(3)/phm(7))**2                 ! sin(th1)的平方 、 th1：入射角
        dsinth2=sinratio*dsinth1                       ! sin(th2)的平方 、 th2：折射角
        dcosth2=DSIGN(DSQRT(1d0-dsinth2),vel(3))       ! cos(th2) 、 
        phm(4)=phm(4)*DSQRT(sinratio)                  ! 原本的phm(4)是與1軸夾角的cos值，新phm(4)=vel2(1)/v2，折射後移動方向與1軸的夾角COS值
        dsinth2=DSQRT(1d0-phm(4)*phm(4))               ! dsinth2為折射後的方向向量在23平面上的投影
        IF (DABS(dcosth2/dsinth2).gt.1) dcosth2=DSIGN(dsinth2,dcosth2) ! in case tmp=0，phm(5)=90/180度，則速度向量在13平面上，無2方向分量
        ! 上面那行???????????????????????????????????????????????????????
        phm(5)=DASIN(dcosth2/dsinth2)                  !因為(sin phm(5))=(方向向量在3軸上的分量)/(方向向量在23平面上的投影)
        IF (tmp.lt.0) phm(5)=M_PI-phm(5)
        IF (phm(5).lt.0) phm(5)=phm(5)+M_PI_2  
    
    ELSE IF (idx.eq.2) THEN
        tmp=DSIN(phm(5))
        dsinth1=1d0-(vel(2)/phm(7))**2                 ! sin(th1)**2
        dsinth2=sinratio*dsinth1                       ! sin(th2)**2
        dcosth2=DSIGN(DSQRT(1d0-dsinth2),vel(2))       ! cos(th2)
        phm(4)=phm(4)*DSQRT(sinratio)
        dsinth2=DSQRT(1d0-phm(4)*phm(4))
        IF (DABS(dcosth2/dsinth2).gt.1) dcosth2=DSIGN(dsinth2,dcosth2) ! in case tmp=0
        phm(5)=DACOS(dcosth2/dsinth2)
        IF (tmp.lt.0) phm(5)=M_PI_2-phm(5)  
    ELSE IF (idx.eq.1) THEN
        dsinth1=1d0-phm(4)**2                          ! sin(th1)**2
        dsinth2=sinratio*dsinth1                       ! sin(th2)**2
        dcosth2=DSIGN(DSQRT(1d0-dsinth2),vel(1))       ! cos(th2)
        phm(4)=dcosth2                                 ! phm(4) not changed 
    ENDIF

END SUBROUTINE Snells
!============================================================================
SUBROUTINE diffuseB( phm,idx,face0,nc )  !決定亂反/穿射後的方向
! nc=1 : transmission
! nc=-1 : reflection
! face > 0 : the right face of the cell
! face < 0 : the left face of the cell
    IMPLICIT NONE
    INTEGER*4:: face0,nc,idx,i
    REAL*8:: phm(iNprop),tmp,tmp2,tmp3
    REAL*8,ALLOCATABLE:: rannum(:)
    ALLOCATE( rannum(2) )
    CALL random_number( rannum )

    phm(4)=DSQRT(rannum(1))*DBLE(nc*face0)
    phm(5)=M_PI_2*rannum(2)

    IF (idx.eq.2) THEN ! U2=phm(4); U1=SQRT(1-U2*U2)*COS(phm(5)); U3=SQRT(1-U2*U2)*SIN(phm(5))
        tmp3=DSIN(phm(5)) ! SIGN(tmp3)=SIGN(U3)
        tmp=DSQRT(1d0-phm(4)**2)*DCOS(phm(5))   ! = V1/V  
        tmp2=DSQRT(1d0-tmp*tmp)
        IF (DABS(phm(4)/tmp2).gt.1d0) phm(4)=DSIGN(tmp2,phm(4)) ! in case tmp3=0
        phm(5)=DACOS( phm(4)/tmp2 )  
        IF (tmp3.lt.0d0) phm(5)=M_PI_2-phm(5)
        phm(4)=tmp
    ELSE IF (idx.eq.3) THEN  ! U3=phm(4); U1=SQRT(1-U3*U3)*COS(phm(5)); U2=SQRT(1-U3*U3)*SIN(phm(5))
        tmp3=DSIN(phm(5))     ! SIGN(tmp3)=SIGN(U2)
        tmp=DSQRT(1d0-phm(4)**2)*DCOS(phm(5))   ! = V1/V 
        tmp2=DSQRT(1d0-tmp*tmp)
        IF (DABS(phm(4)/tmp2).gt.1d0) phm(4)=DSIGN(tmp2,phm(4)) ! in case tmp3=0 
        phm(5)=DASIN( phm(4)/tmp2 )  
        IF (tmp3.lt.0d0) phm(5)=M_PI-phm(5)
        IF (phm(5).lt.0d0) phm(5)=phm(5)+M_PI_2
        phm(4)=tmp
    ENDIF

    DEALLOCATE( rannum )
END SUBROUTINE diffuseB
!============================================================================
    SUBROUTINE proc_reorder
    IMPLICIT NONE
    INTEGER*4:: i, j, k, tot, m
    REAL*8, ALLOCATABLE:: lcr(:, :), newp(:, :)

        iNnumcell=0
        iNbgcell=0
        ALLOCATE(lcr(iNph,3))

        DO m = 1, iNph
            IF ( phn(6, m).gt.0 ) THEN
                i = INT( phn(1, m) / dLclen(1) ) + 1
                j = INT( phn(2, m) / dLclen(2) ) + 1
                k = INT( phn(3, m) / dLclen(3) ) + 1
                iNnumcell(i, j, k) = iNnumcell(i, j, k) + 1
                lcr(m, 1) = i
                lcr(m, 2) = j
                lcr(m, 3) = k
            ENDIF
        ENDDO

        tot=0

        DO k = 1, iNcell(3)
            DO j = 1, iNcell(2)
                DO i = 1, iNcell(1)
                    iNbgcell(i, j, k) = tot
                    tot = tot + iNnumcell(i,j,k)
                ENDDO
            ENDDO
        ENDDO

        iNnumcell=0

        ALLOCATE( newp(iNprop, tot) )

        DO m = 1, iNph
            IF ( phn(6, m).gt.0 ) THEN
                i = lcr(m, 1)
                j = lcr(m, 2)
                k = lcr(m, 3)
                iNnumcell(i, j, k) = iNnumcell(i, j, k)+1 !最後會得到第(i,j,k)網格的聲子束數量
                newp(:, iNbgcell(i, j, k) + iNnumcell(i, j, k)) = phn(:, m)
            ENDIF
        ENDDO

        DEALLOCATE( phn )
        ALLOCATE( phn(iNprop, tot) )
    
        phn = newp
        iNph = tot

        DEALLOCATE(lcr,newp)
    
    END SUBROUTINE proc_reorder
!============================================================================
    SUBROUTINE cellinfo
    IMPLICIT NONE
    INTEGER*4:: i, j, k, bg, ed

        DO k = 1, iNcell(3)
            DO j = 1, iNcell(2)
                DO i = 1, iNcell(1)
                    bg = iNbgcell(i, j, k) + 1
                    ed = iNbgcell(i, j, k) + iNnumcell(i, j, k)
                    dEcell(i, j, k) = SUM( phn(6, bg:ed) ) / dVolume  ! energy density in this cell
                    CALL Etable( iCmat(i,j,k), 2, dEcell(i,j,k), dEunit(i,j,k) ) ! number density of phonons
                    CALL Etable( iCmat(i,j,k), 4, dEcell(i,j,k), dVunit(i,j,k) ) ! phonon group velocity
                    CALL Etable( iCmat(i,j,k), 5, dEcell(i,j,k), MFP(i,j,k) ) ! MFP in this cell
                    CALL Etable( iCmat(i,j,k), 1, dEcell(i,j,k), dTemp(i,j,k)) ! temperature
                    dEunit(i, j, k) = dEcell(i, j, k) / dEunit(i,j,k) * bundle( iCmat(i, j, k) )
                    ! average energy per phonon bundle
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE cellinfo
!============================================================================
SUBROUTINE Etable(mat,idx,E,out) !mat=1:Ge, mat=2:Si, idx=1:temperature, 2:number density, 4:velocity, 5:MFP, 6:specific heat
IMPLICIT NONE
INTEGER*4::mat,idx,iU
REAL*8::E,out,out_a,out_b

IF (mat.eq.1) THEN
	iU = INT( (E-Ge_start)/dU_Ge )+1
	IF (iU.ge.iN_Ge2.or.iU.lt.0) PAUSE 'out of table!'
	out_a = Ge_table(idx,iU)
	out_b = Ge_table(idx,iU+1)
	out = out_a + (out_b-out_a)*(E-Ge_table(3,iU))/dU_Ge
ELSE IF (mat.eq.2) THEN
	iU = int( (E-Si_start)/dU_Si )+1
	IF (iU.ge.iN_Si2.or.iU.lt.0) PAUSE 'out of table'
	out_a = Si_table(idx,iU)
	out_b = Si_table(idx,iU+1)
	out = out_a + (out_b-out_a)*(E-Si_table(3,iU))/dU_Si
END IF

END SUBROUTINE Etable
!============================================================================
SUBROUTINE proc_energy(nc,T0,Eout)
    USE mod_VARIABLES
    IMPLICIT NONE
    INTEGER*4 :: nc
    REAL*8 :: T0,Eout,Ea,Eb,Ta,Tb,Tout

    IF (nc.eq.1) THEN
        Ea = Ge_table(3,1)
        Eb = Ge_table(3,iN_Ge2)
        Ta = Ge_table(1,1)-T0
        Tb = Ge_table(1,iN_Ge2)-T0
    ELSE IF (nc.eq.2) THEN
        Ea = Si_table(3,1)
        Eb = Si_table(3,iN_Si2)
        Ta = Si_table(1,1)-T0
        Tb = Si_table(1,iN_Si2)-T0
    ENDIF

    Eout=(Ea+Eb)/2d0

    CALL ETable(nc,1,Eout,Tout)
    
    Tout=Tout-T0
    
    DO WHILE(ABS(Tout).gt.zero_tol)
        IF (Tout*Ta.lt.0d0) THEN
            Eb=Eout
	        Tb=Tout
        ELSE
	        Ea=Eout
	        Ta=Tout
        ENDIF
        Eout=(Ea+Eb)/2d0
        CALL ETable(nc,1,Eout,Tout)
        Tout=Tout-T0
    ENDDO
END SUBROUTINE proc_energy
!============================================================================
END MODULE mod_ROUTINES