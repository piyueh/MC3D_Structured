!**********************************************************************
!   HEADING: MC3D - Structured Grids:
!                                 routines related to phonon movements
!   AUTHOR: MJ Huang, PY Chuang
!   PURPOSE: This module time marches the simulation.
!   DATE : 07/10/2009
!**********************************************************************
MODULE mod_ADVANCE
USE mod_VARIABLES
USE mod_ROUTINES
IMPLICIT NONE
CONTAINS
!======================================================================
    SUBROUTINE advance
    USE mod_heatcontrol
    IMPLICIT NONE
    INTEGER*4:: i, j, k, bg, ed
    REAL*8:: cellbdy(2, 3)

        DO k = 1, iNcell(3)
            DO j = 1, iNcell(2)
                DO i = 1, iNcell(1)
                    cellbdy(1, 1) = dLclen(1) * DBLE(i-1)
                    cellbdy(2, 1) = dLclen(1) * DBLE(i)
                    cellbdy(1, 2) = dLclen(2) * DBLE(j-1)
                    cellbdy(2, 2) = dLclen(2) * DBLE(j)
                    cellbdy(1, 3) = dLclen(3) * DBLE(k-1)
                    cellbdy(2, 3) = dLclen(3) * DBLE(k)
                    bg = iNbgcell(i, j, k) + 1
                    ed = iNbgcell(i, j, k) + iNnumcell(i, j, k)
                    CALL proc_advection(i, j, k, cellbdy, &
                                        iNnumcell(i, j, k), &
                                        phn(1:iNprop, bg:ed), 1)
                ENDDO
            ENDDO
        ENDDO
        ! So far, the movement, boundary scattering, and
        ! interfacial scattering procedures have been finished.

        IF (option(1).eq.3) CALL proc_heatcontrol

        !build cell information for future use
        CALL proc_reorder
        CALL proc_createdelete
        CALL cellinfo

    END SUBROUTINE advance
!======================================================================
!======================================================================
    SUBROUTINE proc_advection(i0, j0, k0, cellbdy0, N, ph, nc)
    IMPLICIT NONE
    INTEGER*4:: i0, j0, k0, hit, true, nc
    INTEGER*4:: N, m, i, phcell(3), face(3), idxt(1)
    REAL*8:: cellbdy0(2, 3), ph(iNprop, N)
    REAL*8:: dtremain, dtused
    REAL*8:: cellbdy(2, 3), vel(3), ds(3)
    ! i0, j0, k0: the initial cell index of the phonons
    ! phcell: the cell index follows with the target phonon
    !------------------------------------------------------------------

        DO m=1,N


            cellbdy = cellbdy0
            phcell = (/i0, j0, k0/)


            !----------------------------------------------------------
            ! Vx = vel(1) = Vg * cos(theta)
            ! Vy = vel(2) = Vg * sin(theta) * cos(phi)
            ! Vz = vel(3) = Vg * sin(theta) * sin(phi)
            !----------------------------------------------------------
            vel(3) = ph(7, m) * DSQRT( 1d0 - ph(4, m)**2 )
            vel(1) = ph(7, m) * ph(4, m)
            vel(2) = vel(3) * DCOS( ph(5, m) )
            vel(3) = vel(3) * DSIN( ph(5, m) )


            IF (nc.eq.1) THEN
                dtremain = dt
            ELSE
                dtremain = dtheat(m)
            ENDIF


            DO WHILE ( DABS( dtremain ).le.zero_tol )

                DO i = 1, 3
                    IF ( vel(i).gt.0d0 ) THEN
                        ds(i) = cellbdy(2, i)
                        face(i) = 1
                    ELSE
                        ds(i) = cellbdy(1, i)
                        face(i) = -1
                    ENDIF
                ENDDO

                !------------------------------------------------------
                ! ds: times to reach element's surface
                ! idxt and hit: surface the phonon will reach actually
                ! dtused: the time used for current movement.
                !         If dtremain is less than dt(hit), the phonon
                !         will not reach the element surface actually.
                ! The phonon will move to the new position using dtused
                ! and be adjusted whether intrinsic scattering occurs.
                !------------------------------------------------------
                ds = DABS( ( ds - ph(1:3,m) ) / vel )
                idxt = MINLOC( ds )
                hit = idxt(1)
                dtused = MIN( ds(hit), dtremain )
                ph(1:3, m) = ph(1:3, m) + dtused * vel
                CALL proc_intrinsicscattering( phcell, &
                                               ph(1:iNprop, m), &
                                               dtused, true )

                !------------------------------------------------------
                ! true = 1 represents that the intrinsic scattering
                ! occurs. The velocity will be changed.
                !------------------------------------------------------
                IF (true.eq.1) THEN
                    vel(3) = ph(7, m) * DSQRT( 1d0 - ph(4, m)**2 )
                    vel(1) = ph(7, m) * ph(4, m)
                    vel(2) = vel(3) * DCOS( ph(5, m) )
                    vel(3) = vel(3) * DSIN( ph(5, m) )
                ENDIF

                IF ( dtused.lt.dtremain ) THEN

                    dtremain = dtremain - dtused
                    true=0

                    !--------------------------------------------------
                    ! Adjust whether the phonon still go through the
                    ! element surface to a neighbor element (true = 0),
                    ! ,go through the boundary of computational domain 
                    ! (true = -1), or will be reflected back to the 
                    ! same element (true = 1).
                    !--------------------------------------------------
                    IF ( ( face(hit) * vel(hit) ).gt.0 ) THEN
                        CALL proc_outdomain( phcell, cellbdy, hit, &
                                             face(hit), &
                                             ph(1:iNprop, m), &
                                             dtremain, true )
                        !----------------------------------------------
                        ! true = 0, the phonon still go through the
                        ! element surface to a neighbor element.
                        ! proc_transmissivity will adjust whether the
                        ! boundary of this element is a material/grain 
                        ! interface
                        !----------------------------------------------
                        IF ( true.eq.0 ) &
                            CALL proc_transmissivity( phcell, cellbdy,&
                                                      hit, face(hit), &
                                                      ph(1:iNprop, m),&
                                                      vel, true )
                    ENDIF

                    !--------------------------------------------------
                    ! If true is 1, it represents the direction is 
                    ! changed.  The velocity vector therefore must be 
                    ! changed, too.
                    !--------------------------------------------------
                    IF ( true.eq.1 ) THEN
                        vel(3) = ph(7, m) * DSQRT( 1d0 - ph(4,m)**2 )
                        vel(1) = ph(7, m) * ph(4, m)
                        vel(2) = vel(3) * DCOS( ph(5, m) )
                        vel(3) = vel(3) * DSIN( ph(5, m) )
                    ENDIF
                    
                ELSE
                    dtremain = 0d0
                ENDIF

            ENDDO
        ENDDO

    END SUBROUTINE proc_advection
!======================================================================
!======================================================================
    SUBROUTINE proc_intrinsicscattering( phcell, phm, dt1, true )
    IMPLICIT NONE
    INTEGER*4:: phcell(3), true
    REAL*8:: phm(iNprop), dt1, prob, rannum1(3)
    !------------------------------------------------------------------
    ! This subroutine is used to adjust whether the intrinsic
    ! scattering occurs.  If occurs, the properties of the phonon will
    ! be changed
    !
    ! phcell: the element index of the phonon
    ! phm: the properties of the phonon
    ! dt1: the time needed for the movement
    ! true: the subroutine will return value 1 to this parameter if the
    !       scattering occurs and value 0 otherwise
    !------------------------------------------------------------------

        CALL RANDOM_NUMBER( rannum1 )

        ! prob represents the probability of scattering occured during
        ! time interval dt1
        prob = 1d0 - &
               DEXP( -dt1 * phm(7) / &
               MFP( phcell(1), phcell(2), phcell(3) ) )

        IF ( rannum1(1).le.prob ) THEN
            phm(4) = 2D0 * rannum1(2) - 1D0
            phm(5) = M_PI_2 * rannum1(3)
            dEdiff(phcell(1), phcell(2), phcell(3)) = &
                            dEdiff(phcell(1), phcell(2), phcell(3)) + &
                            phm(6) - &
                            dEunit(phcell(1), phcell(2), phcell(3))
            phm(6) = dEunit(phcell(1), phcell(2), phcell(3))
            phm(7) = dVunit(phcell(1), phcell(2), phcell(3))
            phm(8) = iCmat(phcell(1), phcell(2), phcell(3))
            true = 1
        ELSEIF
            true = 0
        ENDIF

    END SUBROUTINE proc_intrinsicscattering
!======================================================================
!======================================================================
    SUBROUTINE proc_transmissivity( phcell, cellbdy, hit, &
                                                face0, phm, vel, true)
    IMPLICIT NONE
    INTEGER*4:: phcell(3), neighbor(3), hit, face0, true, i
    REAL*8:: phm(iNprop), cellbdy(2,3), vel(3)
    REAL*8:: ratio, dcosth2, dsinth2, rho1, rho2
    REAL*8:: neighborE, neighborV
    REAL*8:: tau21, tau12
    REAL*8:: rannum1(2)
    !------------------------------------------------------------------
    ! This subroutine is used to determine whether the phonon encounter
    ! a material/grain interface when it in going through a element's
    ! boundary.  And determine whether the response is diffused or
    ! specular.
    !
    ! phcell: the index of the element in which the target phonon is.
    ! cellbdy: the coordinate of the 6 surfaces of the element
    ! hit: the direction which the target phonon will transmit through
    ! face0: the surface which the target phonon will transmit through
    ! phm: the properties of target phonon
    ! vel: the velocity vector of the target phonon
    ! true: the subroutine will return true = 1 if the element's 
    !       boundary is a material/grain interface
    !------------------------------------------------------------------

        neighbor = phcell
        neighbor(hit) = neighbor(hit) + face0
        ! neighbor now is the element which the phonon will go into.

        IF ( iCmat(neighbor(1), neighbor(2), neighbor(3)).eq. &
                          iCmat(phcell(1), phcell(2), phcell(3)) ) THEN
        
            IF ( hit.eq.1 ) CALL Compute_qflux(phcell,phm,vel)
            tau12 = 1d0
            phcell = neighbor
            cellbdy(1:2, hit) = cellbdy(1:2, hit) + &
                                              dLclen(hit) * DBLE(face0)

        ELSE ! hit the interface
    
            CALL RANDOM_NUMBER( rannum1 )
            !----------------------------------------------------------
            ! Note:
            !   Use the temperature of current element and the material
            !   of the neighbor element to calculate the properties of 
            !   the neighbor element which the phonon will go into.
            !----------------------------------------------------------
            CALL proc_Energy( &
                        iCmat(neighbor(1), neighbor(2), neighbor(3)), &
                        dTemp(phcell(1), phcell(2), phcell(3)), &
                        neighborE ) 
            CALL Etable( iCmat(neighbor(1), neighbor(2), neighbor(3)),&
                         4, neighborE, neighborV )
                     
            IF ( rannum1(1).le.DPP(hit) ) THEN ! specular response

                ratio = (dEcell(phcell(1), phcell(2), phcell(3)) * &
                        dVunit(phcell(1), phcell(2), phcell(3))) / &
                        (neighborE * neighborV) !(U1*v1)/(U2*v2)
                dcosth2 = vel(hit) / phm(7) ! cos(theta1)
                dsinth2 = DSQRT( (1d0 - dcosth2**2) * ratio )
                tau12=0d0
            
                !------------------------------------------------------
                ! If dsinth2 > 1, it represents total reflection. Then, 
                ! tau12 keeps 0. Otherwise if dsinth2 < 1, refraction
                ! occurs.
                !------------------------------------------------------
                IF ( dsinth2.lt.1d0 ) THEN
                    rho1 = rho(iCmat(phcell(1), phcell(2), phcell(3)))
                    rho2 = rho(iCmat(neighbor(1), neighbor(2), &
                                                          neighbor(3)))
                    dcosth2 = DSQRT( 1d0 - dsinth2**2 ) ! cos(theta2)
                    tau12 = (rho2 * neighborV * dcosth2) / &
                                                DABS( rho1 * vel(hit) )
                    tau12 = 1d0 - ((1d0 - tau12) / (1d0 + tau12))**2
                ENDIF

                !------------------------------------------------------
                ! If the random number is smaller tau12, the specular 
                ! transmission occurs and IAMM is applied. Otherwise 
                ! specular reflection occurs.  Only the direction and 
                ! group velocity will be changed during specular 
                ! transmission.  And only the direction will be changed
                ! during specular reflaction.
                !------------------------------------------------------
                IF ( rannum1(2).lt.tau12 ) THEN
                    IF ( hit.eq.1 ) &
                                 CALL Compute_qflux( phcell, phm, vel )
                    CALL Snells( ratio, phm, vel, hit )
                    cellbdy(1:2, hit) = cellbdy(1:2, hit) + &
                                              dLclen(hit) * DBLE(face0)
                    phm(7) = dVunit( neighbor(1), neighbor(2), &
                                                          neighbor(3) ) 
                    phcell = neighbor
                ELSE
                    SELECTCASE(hit)
                    CASE(1)
                        phm(4) = -phm(4)
                    CASE(2)
                        phm(5) = M_PI - phm(5)
                        IF ( phm(5).lt.0 ) phm(5) = phm(5) + M_PI_2
                    CASE(3)
                        phm(5) = M_PI_2 - phm(5)
                    ENDIF
                ENDIF
	    
            ELSE ! diffused response
                tau12 = (neighborE * neighborV) / &
                        (dEcell(phcell(1), phcell(2), phcell(3)) * &
                        dVunit(phcell(1), phcell(2), phcell(3)) + &
                        neighborE * neighborV)
                !------------------------------------------------------
                ! If the random number is smaller tau12, the diffused 
                ! transmission occurs and DMM is applied. Otherwise 
                ! diffused reflection occurs.  The direction and 
                ! group velocity will be changed during both diffused 
                ! transmission and reflaction.
                !------------------------------------------------------
                IF ( rannum1(2).le.tau12 ) THEN
                    IF ( hit.eq.1 ) &
                                 CALL Compute_qflux( phcell, phm, vel )
                    cellbdy(1:2, hit) = cellbdy(1:2, hit) + &
                                                dLclen(hit)*DBLE(face0)
                    CALL diffuseB( phm, hit, face0, 1 )
                    phm(7) = dVunit(neighbor(1), neighbor(2), neighbor(3))
                    phcell = neighbor
                ELSE
                    CALL diffuseB( phm, hit, face0, -1 )
                    phm(7) = dVunit(phcell(1), phcell(2), phcell(3))
                ENDIF
            ENDIF
	    
            true = 1
        
        ENDIF

    END SUBROUTINE proc_transmissivity
!======================================================================
!======================================================================
    SUBROUTINE proc_outdomain( phcell, cellbdy, hit, &
                                            face0, phm, dtremain, true)
    IMPLICIT NONE
    INTEGER*4:: phcell(3), hit, face0, true, judge, j, k
    REAL*8:: phm(iNprop), cellbdy(2,3), dtremain, tmp, rannum
    !------------------------------------------------------------------
    ! This subroutine adjust whether the target phonon will leave the
    ! computational domain.  If this occurs, it will modify some
    ! properties of the phonon according to the boundary conditions.
    !
    ! phcell: the index of the element in which the target phonon is.
    ! cellbdy: the coordinate of the 6 surfaces of the element
    ! hit: the direction which the target phonon will transmit through
    ! face0: the surface which the target phonon will transmit through
    !        1 represent positive surface, -1 otherwise
    ! phm: the properties of target phonon
    ! dtremail: the remaining drift time of target phonon
    ! true: the subroutine will return 1, which represents the phonon 
    !       has been reflected into the computational domain again. 
    !       Othereise, it will return -1, which represents the phonon
    !       has transmitted through the domain boundary.
    !
    ! P.S. In current version. The first and the last cells must belong
    !      to the same material.
    !------------------------------------------------------------------
    
        CALL RANDOM_NUMBER( rannum )

        judge=0

        IF ( ( phcell(hit).eq.iNcell(hit) ).and.( face0.gt.0 ) ) &
                                                              judge = 1

        IF ( ( phcell(hit).eq.1 ).and.( face0.lt.0 ) judge = -1
        
        !--------------------------------------------------------------
        ! The phonon will transmit through the domain 
        ! boundary if judge.ne.0
        !--------------------------------------------------------------
        IF ( judge.ne.0 ) THEN

            SELECTCASE( option(hit) )
            !----------------------------------------------------------
            ! option(hit):
            !       1: adiabatic BC (the phonon will be reflected)
            !       2: periodic BC (the phonon will be moved to the 
            !                       other side)
            !       3: thermal control BC
            !----------------------------------------------------------
            CASE(1)
                !------------------------------------------------------
                ! rannum <= dPPB represents specular reflection, and
                ! rannum > dPPPB represents diffused reflection.
                !------------------------------------------------------
                IF ( rannum.le.dPPB(hit) ) THEN
                
                    SELECTCASE(hit)
                    CASE(1)
                        phm(4) = -phm(4)
                    CASE(2)
                        phm(5) = M_PI - phm(5)
                        IF ( phm(5).lt.0 ) phm(5) = phm(5) + M_PI_2
                    CASE(3)
                        phm(5) = M_PI_2 - phm(5)
                    END SELECT
                    
                ELSE
                
                    CALL diffuseB( phm, hit, judge, -1 )
                    dEdiff(phcell(1), phcell(2), phcell(3)) = &
                            dEdiff(phcell(1), phcell(2), phcell(3)) + &
                            phm(6) - &
                            dEunit(phcell(1), phcell(2), phcell(3))
                    phm(6) = dEunit(phcell(1), phcell(2), phcell(3))
                    phm(7) = dVunit(phcell(1), phcell(2), phcell(3))
                    phm(8) = iCmat(phcell(1), phcell(2), phcell(3))
                    
                ENDIF

                true = 1

            CASE(2)
            
                IF ( judge.eq.1 ) THEN
                
                    phm(hit) = 0d0
                    phcell(hit) = 1
                    cellbdy(1, hit) = 0d0
                    cellbdy(2, hit) = dLclen(hit)
                    
                ELSE IF ( judge.eq.-1 ) THEN
                
                    phm(hit) = dLdomain(hit)
                    phcell(hit) = iNcell(hit)
                    cellbdy(2, hit) = dLdomain(hit)
                    cellbdy(1, hit) = dLdomain(hit) - dLclen(hit)
                    
                ENDIF

                IF (true.eq.0) true=-1

            CASE(3) ! possible only if hit=1

                j = phcell(2)
                k = phcell(3)
                
                !------------------------------------------------------
                ! WAY_DIR = 1 represents periodic injection method
                ! ......... 2 represents random injection method
                ! 
                ! The 3rd index of mlost: 1 represents the higher 
                ! temperature surface (entry of heat flux), and
                ! 2 represents lower temperature surface (heat flux 
                ! outlet.) The phonon leaves the computational 
                ! domain from outlet (inlet) surface of heat flux, will
                ! be saved into the pool of inlet (outlet) surface.
                !------------------------------------------------------
                IF ( WAY_DIR.eq.1 ) THEN
                
                    IF ( judge.eq.1 ) THEN
                    
                        mlost(j, k, 1) = mlost(j, k, 1) + 1
                        IF ( mlost(j, k, 1).gt.iNmakeup ) &
                                                    mlost(j, k, 1) = 1
                        dPpool(1, mlost(j, k, 1), j, k, 1) = dtremain
                        dPpool(2:5, mlost(j,k,1), j, k, 1) = phm(2:5)
                        dPpool(6, mlost(j, k, 1), j, k, 1) = phm(8)
                    
                    ELSE
                    
                        mlost(j, k, 2) = mlost(j, k, 2) + 1
                        IF ( mlost(j, k, 2).gt.iNmakeup ) &
                                                    mlost(j, k, 2) = 1
                        dPpool(1, mlost(j, k, 2), j, k, 2) = dtremain
                        dPpool(2:5, mlost(j, k, 2), j, k, 2) = phm(2:5)
                        dPpool(6, mlost(j, k, 2), j, k, 2) = phm(8)
                        
                    ENDIF
                    
                ENDIF
                
                 !-----------------------------------------------------
                 ! WAY_HEAT = 1 represents constant heat flux
                 ! ...........2 represents constant temperature
                 !
                 ! Constant heat flux (WAY_HEAT = 1) is not supported 
                 ! in current version. (Not completed)
                 !-----------------------------------------------------
                IF ( WAY_HEAT.eq.1 ) THEN
                    WRITE(*, *) "Constant heat flux (WAY_HEAT = 1)"//&
                                " is not supported in current "//&
                                "version. (The function is not "//&
                                "completed)"
                    WRITE(*, *) "The Program is Going to Shut Down "//&
                                "in 5 Seconds."
                    CALL SLEEP(5)
                    STOP
                    IF ( judge.eq.1 ) &
                            dElost(j, k, 2) = dElost(j, k, 2) + phm(6)
                    IF ( judge.eq.-1 ) &
                            dElost(j, k, 1) = dElost(j, k, 1) + phm(6)
                ENDIF

                phm(6) = 0
                dtremain = 0
                true = -1

            END SELECT
            
        ENDIF

    END SUBROUTINE proc_outdomain
!======================================================================
!======================================================================
    SUBROUTINE proc_createdelete
    IMPLICIT NONE
    INTEGER*4:: i, j, k, m, s, bg, ed, true
    REAL*8:: random1, tmp
    REAL*8, ALLOCATABLE:: rannum(:, :), newphn(:, :)
    !------------------------------------------------------------------
    ! This subroutine is used to create or delete phonons in elements
    ! to obey energy conservation.
    !------------------------------------------------------------------
    
    
        nadd = 0
        true = 0
        
        DO k = 1, iNcell(3)
            DO j = 1, iNcell(2)
                DO i = 1, iNcell(1)
                    tmp = 0.5d0 * dEunit(i, j, k)
                    IF ( dEdiff(i, j, k).ge.tmp ) THEN
                        
                        nadd(i, j, k) = &
                        INT( dEdiff(i, j, k) / dEunit(i, j, k) + 0.5d0)
	                
                        dEdiff(i, j, k) = &
                                          dEdiff(i, j, k) - &
                                          nadd(i, j, k) * &
                                          dEunit(i, j, k)
                        true = 1
                    ELSE IF ( dEdiff(i, j, k).lt.-tmp ) THEN
                        DO WHILE ( dEdiff(i, j, k).lt.-tmp )
                            
                            CALL RANDOM_NUMBER( random1 )
                            
                            s = iNnumcell(i, j, k)
                            bg = iNbgcell(i, j, k)
                            
                            m = &
                              MIN( INT( bg + random1 * s + 1 ), bg + s)
                            
                            IF ( phn(6, m).gt.0 ) THEN
                                dEdiff(i, j, k) = dEdiff(i, j, k) + &
                                                  phn(6, m)
                                phn(6, m) = 0
                                nadd(i, j, k) = nadd(i, j, k) - 1
                            ENDIF
                        
                        ENDDO
                        true = 1
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
    !------------------------------------------------------------------
    ! 02:58 04/22/2014
    !------------------------------------------------------------------
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

END SUBROUTINE proc_createdelete
!======================================================================
!======================================================================
    SUBROUTINE Compute_qflux( phcell, phm, vel)
    IMPLICIT NONE
    INTEGER*4:: phcell(3)
    REAL*8:: phm(7), vel(3)
    !------------------------------------------------------------------
    ! This subroutine will calculate the heat pass through the middle
    ! plane of the computational domain in x-direction.
    !------------------------------------------------------------------

        IF ( (phcell(1).eq.(iNcell(1)/2)) .and. (vel(1).gt.0d0) ) THEN 
            qflow(phcell(2), phcell(3)) = qflow(phcell(2), phcell(3)) &
                                          + phm(6)
        ELSE IF ( (phcell(1).eq.iNcell(1)/2+1) .and. &
                                                (vel(1).lt.0d0) ) THEN
            qflow(phcell(2), phcell(3)) = qflow(phcell(2), phcell(3)) &
                                          - phm(6)
        ENDIF

    END SUBROUTINE Compute_qflux
!======================================================================
!======================================================================
END MODULE mod_ADVANCE