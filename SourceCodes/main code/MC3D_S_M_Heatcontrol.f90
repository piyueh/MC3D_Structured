!**********************************************************************
!   HEADING: MC3D - Structured Grids: Heat Control Module
!   AUTHOR: MJ Huang, PY Chuang
!   PURPOSE: This module maintain a constant heat flow rate along one 
!            direction.
!   DATE : 07/10/2009
!**********************************************************************
MODULE mod_heatcontrol
USE mod_VARIABLES
USE mod_ROUTINES
USE mod_ADVANCE
IMPLICIT NONE
!**********************************************************************
CONTAINS
!======================================================================
!======================================================================
    SUBROUTINE proc_heatcontrol
    IMPLICIT NONE
    INTEGER*4:: Nmax, s, m, k, ip, j, i, npool
    INTEGER*4:: bg, ed
    REAL*8, ALLOCATABLE:: phnheat(:, :), cellheat(:,:), rannum(:,:)
    INTEGER*4, ALLOCATABLE:: im(:, :, :, :)
    REAL*8:: tmp
        
         !-------------------------------------------------------------
         ! WAY_HEAT = 1 represents constant heat flux
         ! ...........2 represents constant temperature
         !
         ! Constant heat flux (WAY_HEAT = 1) is not supported 
         ! in current version. (Not completed)
         !-------------------------------------------------------------
        SELECTCASE(WAY_HEAT)
        CASE(1)
            WRITE(*, *) "Constant heat flux (WAY_HEAT = 1)"//&
                        " is not supported in current "//&
                        "version. (The function is not "//&
                        "completed)"
            WRITE(*, *) "The Program is Going to Shut Down "//&
                        "in 5 Seconds."
            CALL SLEEP(5)
            STOP
            !CALL proc_BC( dTemp(1, 1:iNcell(2), 1:iNcell(3) ), &
            !             dTemp(iNcell(1), 1:iNcell(2), 1:iNcell(3)) )
            !
            !IF ( (sumQ.eq.0) .and. (abc.eq.1) ) THEN
            !    dEheatflux(:, :, 1) = (dElost(:, :, 1) / &
            !                          SUM( dElost(:, :, 1) ) + &
            !                          dElost(:, :, 2) / &
            !                          SUM( dElost(:, :, 2) )) / 2d0
            !                          
            !    dEheatflux(:, :, 1) = dEheatflux0 * dEheatflux(:,:,1)
            !    dEheatflux(:, :, 2) = dEheatflux(:, :, 1)
            !ELSE
            !    dEheatflux(:, :, 1) = dEheatflux0 * qctrl(:, :)
            !    dEheatflux(:,:,2) = dEheatflux(:, :, 1)
            !ENDIF
            !
            !dElost(:, :, 1) = dElost(:, :, 1) + dEheatflux(:, :, 1)
            !dElost(:, :, 2) = dElost(:, :, 2) - dEheatflux(:, :, 2)
        CASE(2)
            dElost = dElost + dEheatflux
        END SELECT
        
        tmp = MAX( 0d0, MAXVAL( dElost ) / MINVAL( dEinject ) )
        i = INT( tmp ) + 1
        ALLOCATE( cellheat(2, 3) )
        ALLOCATE( im(i, iNcell(2), iNcell(3), 2) )
        im = 0
        iNemit = 0
        
        DO s = 1, 2  ! left or right boundary

            SELECTCASE(s)
            CASE(1)
                ip = 1
            CASE(2)
                ip = iNcell(1)
            END SELECT

            DO k = 1, iNcell(3)
                DO j = 1, iNcell(2)
                    !--------------------------------------------------
                    ! iNemit:
                    !       The number of phonons that is needed to be
                    !       injected.
                    !--------------------------------------------------
                    iNemit(j, k, s) = 0

                    IF ( WAY_DIR.eq.1 ) THEN
                        npool = mlost(j, k, s)
                        IF ( dPpool(6, iNmakeup, j, k, s).gt.0) &
                                                       npool = iNmakeup
                    ENDIF

                    IF ( (WAY_DIR.eq.1) .and. (npool.gt.0) ) THEN
                        DO WHILE( dElost(j, k, s).ge. &
                          (0.5d0 * dEinject(j, k, iCmat(ip, j, k), s) )
                        
                            CALL RANDOM_NUMBER( tmp )
                            m = MIN( INT( tmp * npool ) + 1, npool)
                            iNemit(j, k, s) = iNemit(j, k, s) + 1
                            im( iNemit(j, k, s), j, k, s) = m
                            dElost(j, k, s) = dElost(j, k, s) - &
                               dEinject(j, k, dPpool(6, m, j, k, s), s)
                        
                        ENDDO
                    ELSE
                        iNemit(j, k, s) = MAX( 0, &
                            INT( (dElost(j, k, s) / &
                                dEinject(j, k, iCmat(ip, j, k), s)) + &
                                                              0.5d0 ) )
                        dElost(j, k, s) = dElost(j, k, s) - 
                                 dEinject(j, k, iCmat(ip, j, k), s) * &
                                                DBLE( iNemit(j, k, s) )
                    ENDIF
                ENDDO
            ENDDO

        ENDDO

        newiNph = SUM( iNemit )
        ALLOCATE( phnheat(iNprop, iNph+newiNph) )
        phnheat(1:iNprop, 1:iNph) = phn(1:iNprop, 1:iNph)
        DEALLOCATE( phn )

        m = iNph

        DO s = 1, 2
            DO k = 1, iNcell(3)
                DO j = 1, iNcell(2)
                    IF ( iNemit(j, k, s).gt.0 ) THEN
                        ALLOCATE( dtheat(iNemit(j, k, s)) )
                        cellheat(1, 2) = dLclen(2) * DBLE(j-1)
                        cellheat(2, 2) = dLclen(2) * DBLE(j)
                        cellheat(1, 3) = dLclen(3) * DBLE(k-1)
                        cellheat(2, 3) = dLclen(3) * DBLE(k)
                        
                        SELECTCASE(s)
                        CASE(1)
                            cellheat(1, 1) = 0d0
                            cellheat(2, 1) = dLclen(1)
                            ip=1
                            phnheat(1, m+1:m+iNemit(j, k, s)) = 0d0
                        CASE(2)
                            cellheat(1, 1) = dLdomain(1) - dLclen(1)
                            cellheat(2, 1) = dLdomain(1)
                            ip = iNcell(1)
                            phnheat(1, m+1:m+iNemit(j, k, s)) = &
                                                            dLdomain(1)
                        ENDSELECT
                       
                        IF ( WAY_DIR.eq.1 ) THEN
                            npool=mlost(j,k,s)
                            IF ( dPpool(6, iNmakeup, j, k, s).gt.0 ) &
                                                         npool=iNmakeup
                        ENDIF

                        IF ( (WAY_DIR.eq.1) .and. &
                                           (mlost(j, k, s).gt.0) ) THEN
                            
                            i = iNemit(j, k, s)
                            bg = m + 1
                            ed = m + i
                            
                            phnheat(2:5, bg:ed) = &
                                 dPpool(2:5, im(1:i, j, k, s), j, k, s)
                                 
                            phnheat(8, bg:ed) = &
                                   dPpool(6, im(1:i, j, k, s), j, k, s)
	                
                            phnheat(6, bg:ed) = &
                                   dEinject(j, k, phnheat(8, bg:ed), s)
                                   
                            phnheat(7, bg:ed) = &
                                     dVinject(j, k, iCmat(ip, j, k), s)
                                     
                            dtheat = &
                                   dPpool(1, im(1:i, j, k, s), j, k, s)
                    
                            m = ed
                        
                        ELSE
                    
                            ALLOCATE( rannum(iNemit(j, k, s), 5) )
                            CALL RANDOM_NUMBER( rannum )
                            
                            i = iNemit(j, k, s)
                            bg = m + 1
                            ed = m + i
                            
                            phnheat(2, bg:ed) = &
                              cellheat(1, 2) + rannum(:, 1) * dLclen(2)
                              
                            phnheat(3, bg:ed) = &
                              cellheat(1, 3) + rannum(:, 2) * dLclen(3)
                              
                            phnheat(4, bg:ed) = &
                                        DSIGN( DSQRT( rannum(:, 3) ), &
                                                          DBLE( 2-ip ))
                            
                            phnheat(5, bg:ed) = M_PI_2 * rannum(:, 4)
                            
                            phnheat(6, bg:ed) = &
                                     dEinject(j, k, iCmat(ip, j, k), s)
                                     
                            phnheat(7, bg:ed) = &
                                     dVinject(j, k, iCmat(ip, j, k), s)
                                     
                            phnheat(8, bg:ed) = iCmat(ip, j, k)
                            dtheat = dt * rannum(:, 5)
                            m = ed
                            DEALLOCATE( rannum )
                            
                        ENDIF

                        CALL proc_advection( ip, j, k, &
                                             cellheat, i, k, s), &
                                             phnheat(:, bg:ed), -1)

                        DEALLOCATE( dtheat )
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        iNph = iNph + newiNph
        ALLOCATE( phn(iNprop, iNph) )
        phn = phnheat

        DEALLOCATE( phnheat, cellheat, im )

    END SUBROUTINE proc_heatcontrol
!======================================================================
!======================================================================
    SUBROUTINE proc_BC(TL,TR)
    IMPLICIT NONE
    INTEGER*4:: j,k,m
    REAL*8:: TL(iNcell(2), iNcell(3))
    REAL*8:: TR(iNcell(2), iNcell(3))
    REAL*8:: dNum,dUL,dUR
    !------------------------------------------------------------------
    ! Calculate the properties of phonons injected from heat contral 
    ! boundary.
    !------------------------------------------------------------------
    
        DO m=1,2
            DO k=1,iNcell(3)
                DO j=1,iNcell(2)
                    CALL proc_energy( m, TL(j, k), dUL )
                    CALL Etable( m, 4, dUL, dVinject(j, k, m, 1) )
                    CALL Etable( m, 2, dUL, dNum )
                    dEinject(j, k, m, 1) = dUL / dNum * bundle(m)

                    CALL proc_energy( m, TR(j, k), dUR)
                    CALL Etable( m, 4, dUR, dVinject(j, k, m, 2) )
                    CALL Etable( m, 2, dUR, dNum )
                    dEinject(j, k, m, 2) = dUR / dNum * bundle(m)
                ENDDO
            ENDDO
        ENDDO

    END SUBROUTINE proc_BC
!======================================================================
!======================================================================
END MODULE mod_heatcontrol
