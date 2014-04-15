!**********************************************************************
!   HEADING: MC3D - Structured Grids: Initialization Program
!   AUTHOR: MJ Huang, PY Chuang
!   PURPOSE: All the global constants and variables are defined.
!   DATE : 07/10/2009
!**********************************************************************
PROGRAM initial3D
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4 :: NperCell

!WRITE(*, *)'initial temperature T(Ge), T(Si) in K = '
!READ*,T1,T2
!WRITE(*, *) 'time increment dt (ps) = '
!READ*,dt
!WRITE(*, *)'dLdomain(1:2) (nm),dwidth (nm),iNcell(1:2) ='
!READ*,dLdomain,dwidth,iNcell

! for test01
    abc=1
    dt=1d0
    dLdomain(1)=440d0
    dLdomain(2)=110d0
    dLdomain(3)=110d0
    iNcell(1)=40
    iNcell(2)=22
    iNcell(3)=22

    dLclen=dLdomain/DBLE(iNcell)

    ALLOCATE( iCmat(iNcell(1),iNcell(2),iNcell(3)) )
    ALLOCATE( dEdiff(iNcell(1),iNcell(2),iNcell(3)) )
    ALLOCATE( dTemp(iNcell(1),iNcell(2),iNcell(3)) )
    ALLOCATE( dVunit(iNcell(1),iNcell(2),iNcell(3)) )
    ALLOCATE( dEunit(iNcell(1),iNcell(2),iNcell(3)) )
    ALLOCATE( qctrl(iNcell(2),iNcell(3)) )

    CALL readtable
    CALL domain(NperCell)
    CALL nanostructure
    CALL properties(NperCell)
    CALL output

    DEALLOCATE( iCmat, dEdiff, dTemp, dVunit )
    DEALLOCATE( dEunit, Ge_table, Si_Table, phn )

END PROGRAM initial3D
!======================================================================
!======================================================================
SUBROUTINE readtable
USE mod_VARIABLES
IMPLICIT NONE

    WRITE(*, *) "Reading Material Tables..."

    OPEN(unit=LRGe,file="Ge_real_table.txt")
    OPEN(unit=LRSi,file="Si_real_table.txt")

    READ(LRGe, *) iN_Ge1, iN_Ge2
    READ(LRSi, *) iN_Si1, iN_Si2
    WRITE(*, *) "iN_Ge1 = ", iN_Ge1, "iN_Ge2 = ", iN_Ge2
    WRITE(*, *) "iN_Si1 = ", iN_Si1, "iN_Si2 = ", iN_Si2

    ALLOCATE( Ge_table(iN_Ge1, iN_Ge2) )
    ALLOCATE( Si_table(iN_Si1, iN_Si2) )

    READ(LRGe,*) Ge_table
    READ(LRSi,*) Si_table

    CLOSE(LRGe)
    CLOSE(LRSi)
    Ge_start = Ge_table(3,1)
    dU_Ge=Ge_table(3,2)-Ge_table(3,1)
    Si_start=Si_table(3,1)
    dU_Si=Si_table(3,2)-Si_table(3,1)
    WRITE(*, *) "Ge_start = ", Ge_start, "dU_Ge = ", dU_Ge
    WRITE(*, *) "Si_start = ", Si_start, "dU_Si = ", dU_Si

END SUBROUTINE readtable
!======================================================================
!======================================================================
SUBROUTINE domain(NperCell)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4 :: NperCell

    option(1)=3
    option(2)=2
    option(3)=2

    DPP=0
    DPPB=0
    dEheatflux0=20d0 ! meV/(ps.nm^2)

    iNprop=8
    NperCell=200

    time0=0d0

END SUBROUTINE domain
!======================================================================
!======================================================================
SUBROUTINE nanostructure
USE mod_VARIABLES
IMPLICIT NONE

    iCmat=1
    iCmat(:,7:16,7:16)=2

    dTemp=330d0

END SUBROUTINE nanostructure
!======================================================================
!======================================================================
SUBROUTINE properties(NperCell)
USE mod_VARIABLES
IMPLICIT NONE
REAL*8 :: Eng0(2), N0(2)
INTEGER*4 :: NperCell, i, j, k, m, m_min(2), s
REAL*8,ALLOCATABLE::rannum(:,:)

    WRITE(*, *) "Assigning Properties to Each Mesh..."

    dArea = dLclen(2) * dLclen(3)
    dVolume = dLclen(1) * dArea

    m_min = 100000000

    DO k=1,iNcell(3)
        DO j=1,iNcell(2)
            DO i=1,iNcell(1)
                s=iCmat(i,j,k)
                CALL energy(s, dTemp(i,j,k), dEdiff(i,j,k))     ! U
                CALL ETable(s, 4, dEdiff(i,j,k), dVunit(i,j,k)) ! V
                CALL ETable(s, 2, dEdiff(i,j,k), dEunit(i,j,k)) ! N
                CALL ETable(1, 5, dEdiff(i,j,k), dTemp(i,j,k))  ! MFP

                m=INT( dEunit(i,j,k) * dVolume + 0.5d0)
                IF (m.le.m_min(s)) THEN
                    m_min(s)=m
                    N0(s)=dEunit(i,j,k)
                    Eng0(s)=dEdiff(i,j,k)
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    WRITE(*, *) "   Model Check:"
    WRITE(*, *) "   dt, dLclen = ", dt, dLclen
    WRITE(*, *) "   MFP = ", MINVAL(dTemp)
    WRITE(*, *) "   tau = ", MINVAL(dTemp) / MAXVAL(dVunit)
    WRITE(*, *) "   dz/vel = ", dLclen/MAXVAL(dVunit)
    WRITE(*, *) "   vel = ", MAXVAL(dVunit)
    WRITE(*, *) " "

    !! choose the weighting number under the condition that there are
    !! at least NperCell phonons in a cell if necessary
    IF (m_min(1).le.m_min(2)) THEN
        m=m_min(1)
        iNph=1
    ELSE
        m=m_min(2)
        iNph=2
    ENDIF

    !!-----------------------------------------------------------------
    !! choice 1: same number of phonons per phonon bundle
    !! choice 2: same energy per phonon bundle
    !! choice 3: same number of phonon bundles per cell
    bundle=MAX( 1d0, DBLE(m)/DBLE(Npercell) )
    !! Uncomment the following IF section for choice 2 or 3
    !!IF (iNph.eq.1) THEN
        !!bundle(2)=(bundle(1)*Eng0(1)/N0(1))/(Eng0(2)/N0(2)) !choice 2
        !!bundle(2)=(bundle(1)*N0(2))/N0(1)                   !choice 3
    !!ELSE
        !!bundle(1)=(bundle(2)*Eng0(2)/N0(2))/(Eng0(1)/N0(1)) !choice 2
        !!bundle(1)=(bundle(2)*N0(1))/N0(2)                   !choice 3
    !!ENDIF

    WRITE(*, *)'# of phonons per bundle = ',bundle
    WRITE(*, *)'dEcell = ',Eng0
    WRITE(*, *)'# of phonon bundles per cell = ',INT(N0*dVolume/bundle+0.5d0)
    WRITE(*, *)'energy per phonon bundle = ',Eng0/N0*bundle

    iNmakeup = MAXVAL( INT( N0 * dVolume / bundle + 0.5d0 ) )*2
    ALLOCATE( dPpool(6,iNmakeup,iNcell(2),iNcell(3),2) )
    dPpool=0

    dTemp=dEdiff ! U now
    dEdiff=0
    iNph=0

    DO k = 1, iNcell(3)
        DO j = 1, iNcell(2)
            DO i = 1, iNcell(1)
                s = iCmat(i,j,k)
                m = INT( dEunit(i,j,k) * dVolume / bundle(s) + 0.5d0 )
                IF (m.lt.1) PAUSE
                dEdiff(i,j,k) = &
                    dTemp(i,j,k) * dVolume - &
                    DBLE(m) * bundle(s) * dTemp(i,j,k) / dEunit(i,j,k)
                iNph=iNph+m
            ENDDO
        ENDDO
    ENDDO

    WRITE(*, *)'total # of phonons = ',iNph
    WRITE(*, *)'iNprop = ', iNprop
!------------------------------------
    ALLOCATE( phn(iNprop, iNph) )
    phn=0d0
    iNph = 0
    DO k = 1, iNcell(3)
        DO j = 1, iNcell(2)
            DO i = 1, iNcell(1)
                s=iCmat(i,j,k)
                m=INT(dEunit(i,j,k)*dVolume/bundle(s)+0.5d0)
                phn(6,iNph+1:iNph+m)=dTemp(i,j,k)/dEunit(i,j,k)*bundle(s)
                phn(7,iNph+1:iNph+m)=dVunit(i,j,k)
                phn(8,iNph+1:iNph+m)=iCmat(i,j,k)
                ALLOCATE(rannum(m,5))
                CALL RANDOM_NUMBER(rannum)
                phn(1,iNph+1:iNph+m)=(DBLE(i-1)+rannum(1:m,1))*dLclen(1)
                phn(2,iNph+1:iNph+m)=(DBLE(j-1)+rannum(1:m,2))*dLclen(2)
                phn(3,iNph+1:iNph+m)=(DBLE(k-1)+rannum(1:m,3))*dLclen(3)
                phn(4,iNph+1:iNph+m)=2d0*rannum(1:m,4)-1d0
                phn(5,iNph+1:iNph+m)=M_PI_2*rannum(1:m,5)
                iNph=iNph+m
                DEALLOCATE( rannum )
            ENDDO
        ENDDO
    ENDDO

    WRITE(*, *) "Assigning Properties Has Finished."

END SUBROUTINE properties
!======================================================================
!======================================================================
SUBROUTINE output
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4::i, j

    OPEN(LW1, file='initial.txt')
    WRITE(LW1, *) bundle,dt,time0,iNcell
    WRITE(LW1, *) dLdomain,dLclen,option,DPP,DPPB
    WRITE(LW1, *) dEheatflux0
    WRITE(LW1, *) iNprop,iNph
    WRITE(LW1, *) iCmat
    WRITE(LW1, *) phn(:, :5250000)
    WRITE(LW1, *) dEdiff
    WRITE(LW1, *) iNmakeup
    WRITE(LW1, *) dPpool
    WRITE(LW1, *) 0
    WRITE(LW1, *) qctrl
    WRITE(LW1, *) abc
    CLOSE(LW1)

END SUBROUTINE output
!======================================================================
!======================================================================
SUBROUTINE energy(nc,T0,Eout)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4 :: nc
REAL*8 :: T0, Eout, E1, E2, T1, T2, Tout

    IF (nc.eq.1) THEN
        E1 = Ge_table(3, 1)
        E2 = Ge_table(3, iN_Ge2)
        T1 = Ge_table(1, 1) - T0
        T2 = Ge_table(1, iN_Ge2) - T0
    ELSE IF (nc.eq.2) THEN
        E1 = Si_table(3, 1)
        E2 = Si_table(3, iN_Si2)
        T1 = Si_table(1, 1) - T0
        T2 = Si_table(1, iN_Si2) - T0
    ENDIF

    Eout = (E1 + E2) / 2d0
    CALL ETable(nc, 1, Eout, Tout)
    Tout = Tout - T0

    DO WHILE(ABS(Tout).gt.zero_tol)
        IF (Tout*T1.lt.0d0) THEN
            E2 = Eout
            T2 = Tout
        ELSE
            E1 = Eout
            T1 = Tout
        ENDIF
        Eout = (E1 + E2) / 2d0
        CALL ETable(nc, 1, Eout, Tout)
        Tout = Tout - T0
    ENDDO

END SUBROUTINE energy
!======================================================================
!======================================================================
SUBROUTINE Etable(mat,idx,E,out)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4::mat,idx,iU
REAL*8::E,out,out_a,out_b
!!mat:
!!   1:Ge, mat=2:Si
!!idx:
!!   1:temperature, 2:number density, 4:velocity,
!!   5:MFP, 6:specific heat

    IF (mat.eq.1) THEN
        iU = INT( (E - Ge_start) / dU_Ge ) + 1
        out_a = Ge_table(idx, iU)
        out_b = Ge_table(idx, iU+1)
        out = out_a + (out_b - out_a) * (E - Ge_table(3, iU)) / dU_Ge
    ELSE IF (mat.eq.2) THEN
        iU = int( (E - Si_start) / dU_Si ) + 1
        out_a = Si_table(idx, iU)
        out_b = Si_table(idx, iU+1)
        out = out_a + (out_b - out_a) * (E - Si_table(3, iU)) / dU_Si
    ENDIF

END SUBROUTINE Etable
!======================================================================
!======================================================================
