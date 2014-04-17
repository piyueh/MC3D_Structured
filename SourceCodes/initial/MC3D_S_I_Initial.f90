!**********************************************************************
!   HEADING: MC3D - Structured Grids: Initialization Program
!   AUTHOR: MJ Huang, PY Chuang
!   PURPOSE: All the global constants and variables are defined.
!   DATE : 07/10/2009
!**********************************************************************
PROGRAM initial3D
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4:: NperCell

    WRITE(*,*) 'INPUT THE CASE NAME(NO .''TXT'', NO ''_GRIDFILE''): '
    READ(*,"(A72)") casename

    ! Read material tables
    CALL readtable

    ! Set up the computational domain
    CALL domain(NperCell)

    ! Allocate memory
    CALL alloc_variables

    ! Set up the nanostructure or material distribution
    CALL nanostructure

    ! Set up initial temperature field
    CALL temperature

    ! Calculate properties of phonons and elements
    CALL properties(NperCell)

    ! Output to initial file
    CALL output

    ! Deallocate memory
    CALL Dealloc_variables

END PROGRAM initial3D
!======================================================================
!======================================================================
SUBROUTINE alloc_variables
USE mod_VARIABLES
IMPLICIT NONE

    WRITE(*, *) 'Allocating Memory...'
    ALLOCATE( iCmat(iNcell(1), iNcell(2), iNcell(3)) )
    ALLOCATE( dEcell(iNcell(1), iNcell(2), iNcell(3)) )
    ALLOCATE( dEunit(iNcell(1), iNcell(2), iNcell(3)) )
    ALLOCATE( dVunit(iNcell(1), iNcell(2), iNcell(3)) )
    ALLOCATE( MFP(iNcell(1), iNcell(2), iNcell(3)) )
    ALLOCATE( dTemp(iNcell(1), iNcell(2), iNcell(3)) )
    ALLOCATE( dEdiff(iNcell(1), iNcell(2), iNcell(3)) )
    ALLOCATE( iNnumcell(iNcell(1), iNcell(2), iNcell(3)) )
    WRITE(*, *) 'Allocating Memory Has Finished.'

END SUBROUTINE alloc_variables
!======================================================================
!======================================================================
SUBROUTINE Dealloc_variables
USE mod_VARIABLES
IMPLICIT NONE

    WRITE(*, *) 'Deallocating Memory...'
    DEALLOCATE( iCmat )
    DEALLOCATE( dEcell )
    DEALLOCATE( dEunit )
    DEALLOCATE( dVunit )
    DEALLOCATE( MFP )
    DEALLOCATE( dTemp )
    DEALLOCATE( dEdiff )
    DEALLOCATE( iNnumcell )

    DEALLOCATE( Ge_table, Si_Table, phn, dPpool, mlost )
    WRITE(*, *) 'Deallocating Memory Has Finished.'

END SUBROUTINE Dealloc_variables
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
    WRITE(*, *) "   iN_Ge1 = ", iN_Ge1, "iN_Ge2 = ", iN_Ge2
    WRITE(*, *) "   iN_Si1 = ", iN_Si1, "iN_Si2 = ", iN_Si2

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
    WRITE(*, *) "   Ge_start = ", Ge_start, "dU_Ge = ", dU_Ge
    WRITE(*, *) "   Si_start = ", Si_start, "dU_Si = ", dU_Si

    WRITE(*, *) "Reading Material Tables Has Finished."

END SUBROUTINE readtable
!======================================================================
!======================================================================
SUBROUTINE domain(NperCell)
USE mod_VARIABLES
IMPLICIT NONE
INTEGER*4 :: NperCell

    WRITE(*, *) 'Enter the dimension Lx, Ly, Lz of the model:'
    READ(*, *) dLdomain

    WRITE(*, *) 'Enter the number of elements on x, y, z direction:'
    READ(*, *) iNcell

    dLclen = dLdomain / DBLE(iNcell)
    dArea = dLclen(2) * dLclen(3)
    dVolume = dLclen(1) * dArea

    iNprop = 8

    WRITE(*, *) 'Enter the min. number of phonon group(NperCell) in elements:'
    READ(*, *) NperCell

END SUBROUTINE domain
!======================================================================
!======================================================================
SUBROUTINE nanostructure
USE mod_VARIABLES
IMPLICIT NONE

    WRITE(*, *) 'Setting Up Nanostructure or Material Distribution...'
    iCmat = 1
    WRITE(*, *) 'Finished.'

END SUBROUTINE nanostructure
!======================================================================
!======================================================================
SUBROUTINE temperature
USE mod_VARIABLES
IMPLICIT NONE

    WRITE(*, *) 'Setting Up Initial Temperature Field...'
    dTemp = 330d0
    WRITE(*, *) 'Finished.'

END SUBROUTINE temperature
!======================================================================
!======================================================================
SUBROUTINE properties(NperCell)
USE mod_VARIABLES
IMPLICIT NONE
REAL*8 :: Eng0(2), N0(2)
INTEGER*4 :: NperCell, i, j, k, m, m_min(2), s
REAL*8,ALLOCATABLE::rannum(:,:)

    WRITE(*, *) "Assigning Properties to Each Mesh..."

    m_min = 100000000

    DO k = 1, iNcell(3)
        DO j = 1, iNcell(2)
            DO i = 1, iNcell(1)
                s = iCmat(i, j, k)
                CALL energy( s, dTemp(i, j, k), dEcell(i, j, k) )     ! U
                CALL ETable( s, 4, dEcell(i, j, k), dVunit(i, j, k) ) ! V
                CALL ETable( s, 2, dEcell(i, j, k), dEunit(i,j,k)) ! N
                CALL ETable( s, 5, dEcell(i, j, k), MFP(i,j,k))  ! MFP
                !dEcell is energy density rather than energy

                m = INT( dEunit(i, j, k) * dVolume + 0.5d0)
                !dEunit is number density currently.
                !m is number of phonons in an element

                IF ( m.le.m_min(s) ) THEN
                    m_min(s) = m
                    N0(s) = dEunit(i, j, k) !dEunit is number density
                    Eng0(s) = dEcell(i, j, k)
                ENDIF
            ENDDO
        ENDDO
    ENDDO

    dt = MAXVAL( dLclen ) * 0.5 / MAXVAL( dVunit )

    bundle=MAX( 1d0, DBLE(MINVAL(m_min)) / DBLE(Npercell) )

    iNmakeup = MAXVAL( INT( N0 * dVolume / bundle + 0.5d0 ) ) * 5
    ALLOCATE( dPpool(6, iNmakeup, iNcell(2), iNcell(3), 2) )
    ALLOCATE( mlost(iNcell(2), iNcell(3), 2) )
    dPpool = 0d0
    mlost = 0

    iNph=0
    DO k = 1, iNcell(3)
        DO j = 1, iNcell(2)
            DO i = 1, iNcell(1)
                s = iCmat(i, j, k)
                iNnumcell(i, j, k) = INT( dEunit(i, j, k) * dVolume / bundle(s) + 0.5d0 )
                IF ( iNnumcell(i, j, k).lt.1 ) PAUSE
                dEdiff(i, j, k) = &
                    dEcell(i, j, k) * dVolume - &
                    DBLE(m) * bundle(s) * dEcell(i, j, k) / dEunit(i, j, k)
            ENDDO
        ENDDO
    ENDDO

    iNph = SUM( iNnumcell )
    ALLOCATE( phn(iNprop, iNph) )
    phn=0d0
    iNph = 0
    DO k = 1, iNcell(3)
        DO j = 1, iNcell(2)
            DO i = 1, iNcell(1)
                s = iCmat(i, j, k)
                m = INT( dEunit(i, j, k) * dVolume / bundle(s) + 0.5d0 )
                phn(6, iNph+1:iNph+m) = dEcell(i, j, k) / dEunit(i, j, k) * bundle(s)
                phn(7, iNph+1:iNph+m) = dVunit(i, j, k)
                phn(8, iNph+1:iNph+m) = iCmat(i, j, k)
                ALLOCATE( rannum(m, 5) )
                CALL RANDOM_NUMBER( rannum )
                phn(1, iNph+1:iNph+m) = ( DBLE(i-1) + rannum(:, 1) ) * dLclen(1)
                phn(2, iNph+1:iNph+m) = ( DBLE(j-1) + rannum(:, 2) ) * dLclen(2)
                phn(3, iNph+1:iNph+m) = ( DBLE(k-1) + rannum(:, 3) ) * dLclen(3)
                phn(4, iNph+1:iNph+m) = 2d0 * rannum(:, 4)-1d0
                phn(5, iNph+1:iNph+m) = M_PI_2 * rannum(:, 5)
                iNph=iNph+m
                DEALLOCATE( rannum )
            ENDDO
        ENDDO
    ENDDO

    IF ( iNph.ne.SUM(iNnumcell) ) PAUSE

    WRITE(*, *) "   Model Check:"
    WRITE(*, *) "   dt = ", dt
    WRITE(*, *) "   dLclen = ", dLclen
    WRITE(*, *) "   min. MFP = ", MINVAL(MFP)
    WRITE(*, *) "   MAX. Vg = ",  MAXVAL(dVunit)
    WRITE(*, *) "   # of phonons per bundle = ", bundle
    WRITE(*, *) "   min. # of phonon bundles per cell = ", MINVAL( iNnumcell )
    WRITE(*, *) "   MAX. # of phonon bundles per cell = ", MAXVAL( iNnumcell )
    WRITE(*, *) "   Total # of phonons = ", iNph
    WRITE(*, *) "   iNmakeup = ", iNmakeup
    WRITE(*, *) "Assigning Properties Has Finished."

END SUBROUTINE properties
!======================================================================
!======================================================================
SUBROUTINE output
USE mod_VARIABLES
IMPLICIT NONE

    WRITE(*, *) "Start to Output..."
    OPEN(LW1, FILE=casename(1:LEN_TRIM(casename))//'_initial.txt')
    WRITE(*, *) "   dt..."
    WRITE(LW1, *) dt
    WRITE(*, *) "   iter0..."
    WRITE(LW1, *) 0
    WRITE(*, *) "   dLdomain..."
    WRITE(LW1, *) dLdomain
    WRITE(*, *) "   iNcell..."
    WRITE(LW1, *) iNcell
    WRITE(*, *) "   dLclen..."
    WRITE(LW1, *) dLclen
    WRITE(*, *) "   bundle..."
    WRITE(LW1, *) bundle
    WRITE(*, *) "   iNprop, iNph..."
    WRITE(LW1, *) iNprop, iNph
    WRITE(*, *) "   iNmakeup..."
    WRITE(LW1, *) iNmakeup
    WRITE(*, *) "   iCmat..."
    WRITE(LW1, *) iCmat
    WRITE(*, *) "   phn..."
    WRITE(LW1, *) phn
    WRITE(*, *) "   dEdiff..."
    WRITE(LW1, *) dEdiff
    WRITE(*, *) "   dPpool..."
    WRITE(LW1, *) dPpool
    WRITE(*, *) "   mlost..."
    WRITE(LW1, *) mlost
    CLOSE(LW1)
    WRITE(*, *) "Finished."

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
