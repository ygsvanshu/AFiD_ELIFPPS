!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
!    FILE:     TridiagonalSolvers.F90                          !
!                                                              !
!    CONTAINS: subroutine TridiagonalSolverXM                  !
!              subroutine TridiagonalSolverXC                  !
!              subroutine TridiagonalSolverYM                  !
!              subroutine TridiagonalSolverYC                  !
!              subroutine TridiagonalSolverZM                  !
!              subroutine TridiagonalSolverZC                  !
!                                                              !
!              subroutine TridiagonalSolverXP                  !
!              subroutine TridiagonalSolverXR                  !
!              subroutine TridiagonalSolverYP                  !
!              subroutine TridiagonalSolverYR                  !
!              subroutine TridiagonalSolverZP                  !
!              subroutine TridiagonalSolverZR                  !
!                                                              !
!    PURPOSE: subroutines that are used for                    !
!             solving tridiagonal systems for diffusive        !
!             terms of the Navier-Stokes or Pressure-Poisson   !
!             equations implicitly                             !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------- FOR MOMENTUM EQUATIONS -------------------!

subroutine TridiagonalSolverXM(bcid)

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays ,only : rhx
    use boundary_arrays, only: cfbcxs,cfbcxe

    implicit none

    integer, intent(in) :: bcid

    real                :: amil(nxm),apil(nxm),acil(nxm),appi(nxm)
    real                :: rx1d(nxm)
    real                :: uflx(nxm),vflx(nxm),qflx
    integer             :: ipiv(nxm)
    integer             :: ic,jc,kc,info
    real                :: acil_b
    real                :: betadx,gvalue
    
    betadx = 0.5*al*dt/rey
        
    if (periodic(1)) then

        do jc = xstart(2),xend(2)
            do kc = xstart(3),xend(3)

                do ic = 1,nxm
                    acil_b   =  1.0/(1.0 - ac1ci(ic)*betadx)
                    amil(ic) = -am1ci(ic)*betadx*acil_b
                    acil(ic) =  1.0
                    apil(ic) = -ap1ci(ic)*betadx*acil_b
                    rx1d(ic) =  rhx(ic,jc,kc)*acil_b
                    uflx(ic) =  0.0
                    vflx(ic) =  0.0
                end do

                gvalue    = (amil(1)*apil(nxm))**0.5

                acil(1)   =  acil(1)   - gvalue
                acil(nxm) =  acil(nxm) - ((amil(1)*apil(nxm))/gvalue)

                uflx(1)   = gvalue
                uflx(nxm) = apil(nxm)

                vflx(1)   = 1.0
                vflx(nxm) = amil(1)/gvalue

                call dgttrf(nxm,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,rx1d,nxm,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,uflx,nxm,info)

                qflx = (vflx(1)*rx1d(1) + vflx(nxm)*rx1d(nxm))/(1.0 + (vflx(1)*uflx(1) + vflx(nxm)*uflx(nxm)))

                do ic = 1,nxm               
                    rhx(ic,jc,kc) = rx1d(ic) - (qflx*uflx(ic))
                enddo
            
            end do
        end do

    else

        do jc = xstart(2),xend(2)
            do kc = xstart(3),xend(3)

                do ic = 1,nxm
                    acil_b   =  1.0/(1.0 - ac1ci(ic)*betadx)
                    amil(ic) = -am1ci(ic)*betadx*acil_b
                    acil(ic) =  1.0
                    apil(ic) = -ap1ci(ic)*betadx*acil_b
                    rx1d(ic) =  rhx(ic,jc,kc)*acil_b
                end do

                ! Lower boundary
                acil(1) = acil(1) + cfbcxs(jc,kc,bcid)*amil(1)

                ! Upper boundary
                acil(nxm) = acil(nxm) + cfbcxe(jc,kc,bcid)*apil(nxm)

                ! Solve
                call dgttrf(nxm,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,rx1d,nxm,info)

                do ic = 1,nxm
                    rhx(ic,jc,kc) = rx1d(ic)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverXM

subroutine TridiagonalSolverXC(bcid)

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays ,only : rhx
    use boundary_arrays, only: cfbcxs,cfbcxe

    implicit none

    integer, intent(in) :: bcid

    real                :: amil(nxm),apil(nxm),acil(nxm),appi(nxm)
    real                :: rx1d(nxm)
    real                :: uflx(nxm),vflx(nxm),qflx
    integer             :: ipiv(nxm)
    integer             :: ic,jc,kc,info
    real                :: acil_b
    real                :: betadx,gvalue
    
    betadx = 0.5*al*dt/rey

    if (periodic(1)) then

        do jc = xstart(2),xend(2)
            do kc = xstart(3),xend(3)

                do ic = 1,nxm
                    acil_b   =  1.0/(1.0 - ac1si(ic)*betadx)
                    amil(ic) = -am1si(ic)*betadx*acil_b
                    acil(ic) =  1.0
                    apil(ic) = -ap1si(ic)*betadx*acil_b
                    rx1d(ic) =  rhx(ic,jc,kc)*acil_b
                    uflx(ic) =  0.0
                    vflx(ic) =  0.0
                end do

                gvalue    = (amil(1)*apil(nxm))**0.5

                acil(1)   =  acil(1)   - gvalue
                acil(nxm) =  acil(nxm) - ((amil(1)*apil(nxm))/gvalue)

                uflx(1)   = gvalue
                uflx(nxm) = apil(nxm)

                vflx(1)   = 1.0
                vflx(nxm) = amil(1)/gvalue

                call dgttrf(nxm,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,rx1d,nxm,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,uflx,nxm,info)

                qflx = (vflx(1)*rx1d(1) + vflx(nxm)*rx1d(nxm))/(1.0 + (vflx(1)*uflx(1) + vflx(nxm)*uflx(nxm)))

                do ic = 1,nxm               
                    rhx(ic,jc,kc) = rx1d(ic) - (qflx*uflx(ic))
                enddo

            end do
        end do

    else

        do jc = xstart(2),xend(2)
            do kc = xstart(3),xend(3)

                do ic = 2,nxm
                    acil_b   =  1.0/(1.0 - ac1si(ic)*betadx)
                    amil(ic) = -am1si(ic)*betadx*acil_b
                    acil(ic) =  1.0
                    apil(ic) = -ap1si(ic)*betadx*acil_b
                    rx1d(ic) =  rhx(ic,jc,kc)*acil_b
                end do

                ! Lower boundary
                acil(2) = acil(2) + cfbcxs(jc,kc,bcid)*amil(2)

                ! Upper boundary
                acil(nxm) = acil(nxm) + cfbcxe(jc,kc,bcid)*apil(nxm)

                ! Solve
                call dgttrf(nxm-1,amil(3:nxm),acil(2:nxm),apil(2:nxm-1),appi(2:nxm-2),ipiv,info)
                call dgttrs('N',nxm-1,1,amil(3:nxm),acil(2:nxm),apil(2:nxm-1),appi(2:nxm-2),ipiv,rx1d(2:nxm),nxm-1,info)

                do ic = 2,nxm
                    rhx(ic,jc,kc) = rx1d(ic)
                end do

                rhx(1,jc,kc) = 0.0

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverXC

subroutine TridiagonalSolverYM(bcid)

    use param
    use decomp_2d, only: ystart,yend
    use local_arrays ,only : rhy
    use boundary_arrays, only: cfbcys,cfbcye

    implicit none

    integer, intent(in) :: bcid

    real                :: amjl(nym),apjl(nym),acjl(nym),appj(nym)
    real                :: ry1d(nym)
    real                :: ufly(nym),vfly(nym),qfly
    integer             :: ipjv(nym)
    integer             :: ic,jc,kc,info
    real                :: acjl_b
    real                :: betadx,gvalue
    
    betadx = 0.5*al*dt/rey

    if (periodic(2)) then
    
    do ic = ystart(1),yend(1)
        do kc = ystart(3),yend(3)

                do jc = 1,nym
                    acjl_b   = 1.0/(1.0 - ac2cj(jc)*betadx)
                    amjl(jc) = -am2cj(jc)*betadx*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = -ap2cj(jc)*betadx*acjl_b
                    ry1d(jc) = rhy(ic,jc,kc)*acjl_b
                    ufly(jc) =  0.0
                    vfly(jc) =  0.0
                end do

                gvalue    = (amjl(1)*apjl(nym))**0.5

                acjl(1)   =  acjl(1)   - gvalue
                acjl(nym) =  acjl(nym) - ((amjl(1)*apjl(nym))/gvalue)

                ufly(1)   = gvalue
                ufly(nym) = apjl(nym)

                vfly(1)   = 1.0
                vfly(nym) = amjl(1)/gvalue

                call dgttrf(nym,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ry1d,nym,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ufly,nym,info)

                qfly = (vfly(1)*ry1d(1) + vfly(nym)*ry1d(nym))/(1.0 + (vfly(1)*ufly(1) + vfly(nym)*ufly(nym)))

                do jc = 1,nym                
                    rhy(ic,jc,kc) = ry1d(jc) - (qfly*ufly(jc))
                enddo

            end do
        end do

    else

        do ic = ystart(1),yend(1)
            do kc = ystart(3),yend(3)    

                do jc = 1,nym
                    acjl_b   = 1.0/(1.0 - ac2cj(jc)*betadx)
                    amjl(jc) = -am2cj(jc)*betadx*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = -ap2cj(jc)*betadx*acjl_b
                    ry1d(jc) = rhy(ic,jc,kc)*acjl_b
                end do

                ! Lower boundary
                acjl(1) = acjl(1) + cfbcys(ic,kc,bcid)*amjl(1)

                ! Upper boundary
                acjl(nym) = acjl(nym) + cfbcye(ic,kc,bcid)*apjl(nym)

                ! Solve
                call dgttrf(nym,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ry1d,nym,info)

                do jc = 1,nym
                    rhy(ic,jc,kc) = ry1d(jc)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverYM

subroutine TridiagonalSolverYC(bcid)

    use param
    use decomp_2d, only: ystart,yend
    use local_arrays ,only : rhy
    use boundary_arrays, only: cfbcys,cfbcye

    implicit none

    integer, intent(in) :: bcid

    real                :: amjl(nym),apjl(nym),acjl(nym),appj(nym)
    real                :: ry1d(nym)
    real                :: ufly(nym),vfly(nym),qfly
    integer             :: ipjv(nym)
    integer             :: ic,jc,kc,info
    real                :: acjl_b
    real                :: betadx,gvalue
    
    betadx = 0.5*al*dt/rey
    
    if (periodic(2)) then

        do ic = ystart(1),yend(1)
            do kc = ystart(3),yend(3)

                do jc = 1,nym
                    acjl_b   = 1.0/(1.0 - ac2sj(jc)*betadx)
                    amjl(jc) = -am2sj(jc)*betadx*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = -ap2sj(jc)*betadx*acjl_b
                    ry1d(jc) = rhy(ic,jc,kc)*acjl_b
                    ufly(jc) =  0.0
                    vfly(jc) =  0.0
                end do

                gvalue    = (amjl(1)*apjl(nym))**0.5

                acjl(1)   =  acjl(1)   - gvalue
                acjl(nym) =  acjl(nym) - ((amjl(1)*apjl(nym))/gvalue)

                ufly(1)   = gvalue
                ufly(nym) = apjl(nym)

                vfly(1)   = 1.0
                vfly(nym) = amjl(1)/gvalue

                call dgttrf(nym,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ry1d,nym,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ufly,nym,info)

                qfly = (vfly(1)*ry1d(1) + vfly(nym)*ry1d(nym))/(1.0 + (vfly(1)*ufly(1) + vfly(nym)*ufly(nym)))

                do jc = 1,nym                
                    rhy(ic,jc,kc) = ry1d(jc) - (qfly*ufly(jc))
                enddo

            end do
        end do

    else

        do ic = ystart(1),yend(1)
            do kc = ystart(3),yend(3)

                do jc = 2,nym
                    acjl_b   = 1.0/(1.0 - ac2sj(jc)*betadx)
                    amjl(jc) = -am2sj(jc)*betadx*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = -ap2sj(jc)*betadx*acjl_b
                    ry1d(jc) = rhy(ic,jc,kc)*acjl_b
                end do

                ! Lower boundary
                acjl(2) = acjl(2) + cfbcys(ic,kc,bcid)*amjl(2)

                ! Upper boundary
                acjl(nym) = acjl(nym) + cfbcye(ic,kc,bcid)*apjl(nym)

                ! Solve
                call dgttrf(nym-1,amjl(3:nym),acjl(2:nym),apjl(2:nym-1),appj(2:nym-2),ipjv,info)
                call dgttrs('N',nym-1,1,amjl(3:nym),acjl(2:nym),apjl(2:nym-1),appj(2:nym-2),ipjv,ry1d(2:nym),nym-1,info)

                do jc = 2,nym
                    rhy(ic,jc,kc) = ry1d(jc)
                end do

                rhy(ic,1,kc) = 0.0

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverYC

subroutine TridiagonalSolverZM(bcid)

    use param
    use decomp_2d, only: zstart,zend
    use local_arrays, only: rhz
    use boundary_arrays, only: cfbczs,cfbcze

    implicit none

    integer, intent(in) :: bcid
    
    real                :: amkl(nzm),apkl(nzm),ackl(nzm),appk(nzm)
    real                :: rz1d(nzm)
    real                :: uflz(nzm),vflz(nzm),qflz
    integer             :: ipkv(nzm)
    integer             :: ic,jc,kc,info
    real                :: ackl_b
    real                :: betadx,gvalue
    
    betadx = 0.5*al*dt/rey

    if (periodic(3)) then

        do ic = zstart(1),zend(1)
            do jc = zstart(2),zend(2)

                do kc = 1,nzm
                    ackl_b   =  1.0/(1.0 - ac3ck(kc)*betadx)
                    amkl(kc) = -am3ck(kc)*betadx*ackl_b
                    ackl(kc) =  1.0
                    apkl(kc) = -ap3ck(kc)*betadx*ackl_b
                    rz1d(kc) =  rhz(ic,jc,kc)*ackl_b
                    uflz(kc) =  0.0
                    vflz(kc) =  0.0
                end do

                gvalue    = (amkl(1)*apkl(nzm))**0.5

                ackl(1)   = ackl(1)   - gvalue
                ackl(nzm) = ackl(nzm) - ((amkl(1)*apkl(nzm))/gvalue)

                uflz(1)   = gvalue
                uflz(nzm) = apkl(nzm)

                vflz(1)   = 1.0
                vflz(nzm) = amkl(1)/gvalue

                call dgttrf(nzm,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,rz1d,nzm,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,uflz,nzm,info)

                qflz = (vflz(1)*rz1d(1) + vflz(nzm)*rz1d(nzm))/(1.0 + (vflz(1)*uflz(1) + vflz(nzm)*uflz(nzm)))

                do kc = 1,nzm                
                    rhz(ic,jc,kc) = rz1d(kc) - (qflz*uflz(kc))
                enddo

            end do
        end do

    else

        do ic = zstart(1),zend(1)
            do jc = zstart(2),zend(2)

                do kc = 1,nzm
                    ackl_b   =  1.0/(1.0 - ac3ck(kc)*betadx)
                    amkl(kc) = -am3ck(kc)*betadx*ackl_b
                    ackl(kc) =  1.0
                    apkl(kc) = -ap3ck(kc)*betadx*ackl_b
                    rz1d(kc) =  rhz(ic,jc,kc)*ackl_b
                end do

                ! Lower boundary
                ackl(1) = ackl(1) + cfbczs(ic,jc,bcid)*amkl(1)

                ! Upper boundary
                ackl(nzm) = ackl(nzm) + cfbcze(ic,jc,bcid)*apkl(nzm)

                ! Solve
                call dgttrf(nzm,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,rz1d,nzm,info)

                do kc = 1,nzm
                    rhz(ic,jc,kc) = rz1d(kc)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverZM

subroutine TridiagonalSolverZC(bcid)

    use param
    use decomp_2d, only: zstart,zend
    use local_arrays, only: rhz
    use boundary_arrays, only: cfbczs,cfbcze

    implicit none

    integer, intent(in) :: bcid
    
    real                :: amkl(nzm),apkl(nzm),ackl(nzm),appk(nzm)
    real                :: rz1d(nzm)
    real                :: uflz(nzm),vflz(nzm),qflz
    integer             :: ipkv(nzm)
    integer             :: ic,jc,kc,info
    real                :: ackl_b
    real                :: betadx,gvalue
    
    betadx = 0.5*al*dt/rey
    
    if (periodic(3)) then

        do ic = zstart(1),zend(1)
            do jc = zstart(2),zend(2)

                do kc = 1,nzm
                    ackl_b   =  1.0/(1.0 - ac3sk(kc)*betadx)
                    amkl(kc) = -am3sk(kc)*betadx*ackl_b
                    ackl(kc) =  1.0
                    apkl(kc) = -ap3sk(kc)*betadx*ackl_b
                    rz1d(kc) =  rhz(ic,jc,kc)*ackl_b
                    uflz(kc) =  0.0
                    vflz(kc) =  0.0
                end do

                gvalue    = (amkl(1)*apkl(nzm))**0.5

                ackl(1)   = ackl(1)   - gvalue
                ackl(nzm) = ackl(nzm) - ((amkl(1)*apkl(nzm))/gvalue)

                uflz(1)   = gvalue
                uflz(nzm) = apkl(nzm)

                vflz(1)   = 1.0
                vflz(nzm) = amkl(1)/gvalue

                call dgttrf(nzm,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,rz1d,nzm,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,uflz,nzm,info)

                qflz = (vflz(1)*rz1d(1) + vflz(nzm)*rz1d(nzm))/(1.0 + (vflz(1)*uflz(1) + vflz(nzm)*uflz(nzm)))

                do kc = 1,nzm                
                    rhz(ic,jc,kc) = rz1d(kc) - (qflz*uflz(kc))
                enddo

            end do
        end do

    else

        do ic = zstart(1),zend(1)
            do jc = zstart(2),zend(2)

                do kc = 2,nzm
                    ackl_b   =  1.0/(1.0 - ac3sk(kc)*betadx)
                    amkl(kc) = -am3sk(kc)*betadx*ackl_b
                    ackl(kc) =  1.0
                    apkl(kc) = -ap3sk(kc)*betadx*ackl_b
                    rz1d(kc) =  rhz(ic,jc,kc)*ackl_b
                end do

                ! Lower boundary
                ackl(2) = ackl(2) + cfbczs(ic,jc,bcid)*amkl(2)

                ! Upper boundary
                ackl(nzm) = ackl(nzm) + cfbcze(ic,jc,bcid)*apkl(nzm)

                ! Solve

                call dgttrf(nzm-1,amkl(3:nzm),ackl(2:nzm),apkl(2:nzm-1),appk(2:nzm-2),ipkv,info)
                call dgttrs('N',nzm-1,1,amkl(3:nzm),ackl(2:nzm),apkl(2:nzm-1),appk(2:nzm-2),ipkv,rz1d(2:nzm),nzm-1,info)

                do kc = 2,nzm
                    rhz(ic,jc,kc) = rz1d(kc)
                end do

                rhz(ic,jc,1) = 0.0

            end do
        end do

    end if
    
    return

end subroutine TridiagonalSolverZC

!--------------- FOR PRESSURE POISSON EQUATION  ---------------!

subroutine TridiagonalSolverXP(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : cx1

    implicit none

    integer, intent(in)             :: norm
    type(decomp_info), intent(in)   :: dcmp

    complex                         :: amil(nxm),apil(nxm),acil(nxm),appi(nxm)
    complex                         :: rx1d(nxm)
    complex                         :: uflx(nxm),vflx(nxm),qflx
    integer                         :: ipiv(nxm)
    integer                         :: ic,jc,kc,info
    complex                         :: acil_b,gvalue

    if (periodic(1)) then

        do kc = dcmp%xst(3),dcmp%xen(3)
            do jc = dcmp%xst(2),dcmp%xen(2)

                do ic = 1,nxm
                    acil_b   = 1.0/(acphi(ic) - ak2(jc) - ak3(kc))
                    amil(ic) = amphi(ic)*acil_b
                    acil(ic) = 1.0
                    apil(ic) = apphi(ic)*acil_b
                    rx1d(ic) = cx1(ic,jc,kc)*acil_b/real(norm)
                    uflx(ic) = 0.0
                    vflx(ic) = 0.0
                end do

                gvalue    = (amil(1)*apil(nxm))**0.5

                acil(1)   = acil(1)   - gvalue
                acil(nxm) = acil(nxm) - amil(1)*apil(nxm)/gvalue

                uflx(1)   = gvalue
                uflx(nxm) = apil(nxm)

                vflx(1)   = 1.0
                vflx(nxm) = amil(1)/gvalue

                call zgttrf(nxm,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,info)
                call zgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,rx1d,nxm,info)
                call zgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,uflx,nxm,info)

                qflx = (vflx(1)*rx1d(1) + vflx(nxm)*rx1d(nxm))/(1.0 + (vflx(1)*uflx(1) + vflx(nxm)*uflx(nxm)))

                do ic = 1,nxm               
                    cx1(ic,jc,kc) = rx1d(ic) - (qflx*uflx(ic))
                enddo

            end do 
        end do

    else

        do kc = dcmp%xst(3),dcmp%xen(3)
            do jc = dcmp%xst(2),dcmp%xen(2)

                do ic = 1,nxm
                    acil_b   = 1.0/(acphi(ic) - ak2(jc) - ak3(kc))
                    amil(ic) = amphi(ic)*acil_b
                    acil(ic) = 1.0
                    apil(ic) = apphi(ic)*acil_b
                    rx1d(ic) = cx1(ic,jc,kc)*acil_b/real(norm)
                end do

                acil(1)   = acil(1)   + amil(1)
                acil(nxm) = acil(nxm) + apil(nxm)

                call zgttrf(nxm,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,info)
                call zgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,rx1d,nxm,info)

                do ic = 1,nxm
                    cx1(ic,jc,kc) = rx1d(ic)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverXP

subroutine TridiagonalSolverXR(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : rx1

    implicit none

    type(decomp_info), intent(in)   :: dcmp
    integer, intent(in)             :: norm

    real                            :: amil(nxm),apil(nxm),acil(nxm),appi(nxm)
    real                            :: rx1d(nxm)
    real                            :: uflx(nxm),vflx(nxm),qflx
    integer                         :: ipiv(nxm)
    integer                         :: ic,jc,kc,info
    real                            :: gvalue,acil_b

    if (periodic(1)) then

        do kc = dcmp%xst(3),dcmp%xen(3)
            do jc = dcmp%xst(2),dcmp%xen(2)

                do ic = 1,nxm
                    acil_b   = 1.0/(acphi(ic) - ak2(jc) - ak3(kc))
                    amil(ic) = amphi(ic)*acil_b
                    acil(ic) = 1.0
                    apil(ic) = apphi(ic)*acil_b
                    rx1d(ic) = rx1(ic,jc,kc)*acil_b/real(norm)
                    uflx(ic) = 0.0
                    vflx(ic) = 0.0
                end do

                gvalue    = (amil(1)*apil(nxm))**0.5

                acil(1)   = acil(1)   - gvalue
                acil(nxm) = acil(nxm) - amil(1)*apil(nxm)/gvalue

                uflx(1)   = gvalue
                uflx(nxm) = apil(nxm)

                vflx(1)   = 1.0
                vflx(nxm) = amil(1)/gvalue

                call dgttrf(nxm,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,rx1d,nxm,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,uflx,nxm,info)

                qflx = (vflx(1)*rx1d(1) + vflx(nxm)*rx1d(nxm))/(1.0 + (vflx(1)*uflx(1) + vflx(nxm)*uflx(nxm)))

                do ic = 1,nxm               
                    rx1(ic,jc,kc) = rx1d(ic) - (qflx*uflx(ic))
                enddo
            
            end do
        end do

    else

        do kc = dcmp%xst(3),dcmp%xen(3)
            do jc = dcmp%xst(2),dcmp%xen(2)

                do ic = 1,nxm
                    acil_b   = 1.0/(acphi(ic) - ak2(jc) - ak3(kc))
                    amil(ic) = amphi(ic)*acil_b
                    acil(ic) = 1.0
                    apil(ic) = apphi(ic)*acil_b
                    rx1d(ic) = rx1(ic,jc,kc)*acil_b/real(norm)
                end do

                acil(1)   = acil(1)   + amil(1)
                acil(nxm) = acil(nxm) + apil(nxm)

                call dgttrf(nxm,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,info)
                call dgttrs('N',nxm,1,amil(2:nxm),acil(1:nxm),apil(1:nxm-1),appi(1:nxm-2),ipiv,rx1d,nxm,info)

                do ic = 1,nxm
                    rx1(ic,jc,kc) = rx1d(ic)
                end do
                
            end do
        end do

    end if

    return

end subroutine TridiagonalSolverXR

subroutine TridiagonalSolverYP(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : cy1

    implicit none

    type(decomp_info), intent(in)   :: dcmp
    integer, intent(in)             :: norm

    complex                         :: amjl(nym),apjl(nym),acjl(nym),appj(nym)
    complex                         :: ry1d(nym)
    complex                         :: ufly(nym),vfly(nym),qfly
    integer                         :: ipjv(nym)
    integer                         :: ic,jc,kc,info
    complex                         :: acjl_b,gvalue

    if (periodic(2)) then

        do kc = dcmp%yst(3),dcmp%yen(3)
            do ic = dcmp%yst(1),dcmp%yen(1)

                do jc = 1,nym
                    acjl_b   = 1.0/(acphj(jc) - ak1(ic) - ak3(kc))
                    amjl(jc) = amphj(jc)*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = apphj(jc)*acjl_b
                    ry1d(jc) = cy1(ic,jc,kc)*acjl_b/real(norm)
                    ufly(jc) = 0.0
                    vfly(jc) = 0.0
                end do

                gvalue    = (amjl(1)*apjl(nym))**0.5

                acjl(1)   =  acjl(1)   - gvalue
                acjl(nym) =  acjl(nym) - amjl(1)*apjl(nym)/gvalue

                ufly(1)   = gvalue
                ufly(nym) = apjl(nym)

                vfly(1)   = 1.0
                vfly(nym) = amjl(1)/gvalue

                call zgttrf(nym,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,info)
                call zgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ry1d,nym,info)
                call zgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ufly,nym,info)

                qfly = (vfly(1)*ry1d(1) + vfly(nym)*ry1d(nym))/(1.0 + (vfly(1)*ufly(1) + vfly(nym)*ufly(nym)))

                do jc = 1,nym 
                    cy1(ic,jc,kc) = ry1d(jc) - (qfly*ufly(jc))
                enddo

            end do
        end do

    else

        do kc = dcmp%yst(3),dcmp%yen(3)
            do ic = dcmp%yst(1),dcmp%yen(1)

                do jc = 1,nym
                    acjl_b   = 1.0/(acphj(jc) - ak1(ic) - ak3(kc))
                    amjl(jc) = amphj(jc)*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = apphj(jc)*acjl_b
                    ry1d(jc) = cy1(ic,jc,kc)*acjl_b/real(norm)
                end do

                acjl(1)   = acjl(1)   + amjl(1)
                acjl(nym) = acjl(nym) + apjl(nym)

                call zgttrf(nym,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,info)
                call zgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ry1d,nym,info)

                do jc = 1,nym
                    cy1(ic,jc,kc) = ry1d(jc)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverYP

subroutine TridiagonalSolverYR(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : ry1

    implicit none

    type(decomp_info), intent(in)   :: dcmp
    integer, intent(in)             :: norm

    real                            :: amjl(nym),apjl(nym),acjl(nym),appj(nym)
    real                            :: ry1d(nym)
    real                            :: ufly(nym),vfly(nym),qfly
    integer                         :: ipjv(nym)
    integer                         :: ic,jc,kc,info
    real                            :: acjl_b,gvalue

    if (periodic(2)) then

        do kc = dcmp%yst(3),dcmp%yen(3)
            do ic = dcmp%yst(1),dcmp%yen(1)

                do jc = 1,nym
                    acjl_b   = 1.0/(acphj(jc) - ak1(ic) - ak3(kc))
                    amjl(jc) = amphj(jc)*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = apphj(jc)*acjl_b
                    ry1d(jc) = ry1(ic,jc,kc)*acjl_b/real(norm)
                    ufly(jc) = 0.0
                    vfly(jc) = 0.0
                end do

                gvalue    = (amjl(1)*apjl(nym))**0.5

                acjl(1)   =  acjl(1)   - gvalue
                acjl(nym) =  acjl(nym) - amjl(1)*apjl(nym)/gvalue

                ufly(1)   = gvalue
                ufly(nym) = apjl(nym)

                vfly(1)   = 1.0
                vfly(nym) = amjl(1)/gvalue

                call dgttrf(nym,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ry1d,nym,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ufly,nym,info)

                qfly = (vfly(1)*ry1d(1) + vfly(nym)*ry1d(nym))/(1.0 + (vfly(1)*ufly(1) + vfly(nym)*ufly(nym)))
                if (qfly.eq.qfly+1) print *, "INF QFLY at ", nrank, ic, kc, (1.0 + (vfly(1)*ufly(1) + vfly(nym)*ufly(nym))), vfly(1)*ufly(1), vfly(nym)*ufly(nym)

                do jc = 1,nym
                    ry1(ic,jc,kc) = ry1d(jc) - (qfly*ufly(jc))
                enddo

            end do
        end do

    else

        do kc = dcmp%yst(3),dcmp%yen(3)
            do ic = dcmp%yst(1),dcmp%yen(1)

                do jc = 1,nym
                    acjl_b   = 1.0/(acphj(jc) - ak1(ic) - ak3(kc))
                    amjl(jc) = amphj(jc)*acjl_b
                    acjl(jc) = 1.0
                    apjl(jc) = apphj(jc)*acjl_b
                    ry1d(jc) = ry1(ic,jc,kc)*acjl_b/real(norm)
                end do

                acjl(1)   = acjl(1)   + amjl(1)
                acjl(nym) = acjl(nym) + apjl(nym)

                call dgttrf(nym,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,info)
                call dgttrs('N',nym,1,amjl(2:nym),acjl(1:nym),apjl(1:nym-1),appj(1:nym-2),ipjv,ry1d,nym,info)

                do jc = 1,nym
                    ry1(ic,jc,kc) = ry1d(jc)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverYR

subroutine TridiagonalSolverZP(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : cz1

    implicit none

    type(decomp_info), intent(in)   :: dcmp
    integer, intent(in)             :: norm
    
    complex                         :: amkl(nzm),apkl(nzm),ackl(nzm),appk(nzm)
    complex                         :: rz1d(nzm)
    complex                         :: uflz(nzm),vflz(nzm),qflz
    integer                         :: ipkv(nzm)
    integer                         :: ic,jc,kc,info
    complex                         :: ackl_b,gvalue

    if (periodic(3)) then

        do ic = dcmp%zst(1),dcmp%zen(1)
            do jc = dcmp%zst(2),dcmp%zen(2)

                do kc = 1,nzm
                    ackl_b   = 1.0/(acphk(kc) - ak2(jc) - ak1(ic))
                    amkl(kc) = amphk(kc)*ackl_b
                    ackl(kc) = 1.0
                    apkl(kc) = apphk(kc)*ackl_b
                    rz1d(kc) = cz1(ic,jc,kc)*ackl_b/real(norm)
                    uflz(kc) = 0.0
                    vflz(kc) = 0.0
                end do

                gvalue    = (amkl(1)*apkl(nzm))**0.5

                ackl(1)   =  ackl(1)   - gvalue
                ackl(nzm) =  ackl(nzm) - amkl(1)*apkl(nzm)/gvalue

                uflz(1)   = gvalue
                uflz(nzm) = apkl(nzm)

                vflz(1)   = 1.0
                vflz(nzm) = amkl(1)/gvalue

                call zgttrf(nzm,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,info)
                call zgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,rz1d,nzm,info)
                call zgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,uflz,nzm,info)

                qflz = (vflz(1)*rz1d(1) + vflz(nzm)*rz1d(nzm))/(1.0 + (vflz(1)*uflz(1) + vflz(nzm)*uflz(nzm)))

                do kc = 1,nzm                
                    cz1(ic,jc,kc) = rz1d(kc) - (qflz*uflz(kc))
                enddo

            end do
        end do

    else

        do ic = dcmp%zst(1),dcmp%zen(1)
            do jc = dcmp%zst(2),dcmp%zen(2)

                do kc = 1,nzm
                    ackl_b   = 1.0/(acphk(kc) - ak2(jc) - ak1(ic))
                    amkl(kc) = amphk(kc)*ackl_b
                    ackl(kc) = 1.0
                    apkl(kc) = apphk(kc)*ackl_b
                    rz1d(kc) = cz1(ic,jc,kc)*ackl_b/real(norm)
                end do

                ackl(1)   = ackl(1)   + amkl(1)
                ackl(nzm) = ackl(nzm) + apkl(nzm)

                call zgttrf(nzm,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,info)
                call zgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,rz1d,nzm,info)

                do kc = 1,nzm
                    cz1(ic,jc,kc) = rz1d(kc)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverZP

subroutine TridiagonalSolverZR(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : rz1

    implicit none

    type(decomp_info), intent(in)   :: dcmp
    integer, intent(in)             :: norm
    
    real                            :: amkl(nzm),apkl(nzm),ackl(nzm),appk(nzm)
    real                            :: rz1d(nzm)
    real                            :: uflz(nzm),vflz(nzm),qflz
    integer                         :: ipkv(nzm)
    integer                         :: ic,jc,kc,info
    real                            :: ackl_b,gvalue
    if (periodic(3)) then

        do ic = dcmp%zst(1),dcmp%zen(1)
            do jc = dcmp%zst(2),dcmp%zen(2)

                do kc = 1,nzm
                    ackl_b   = 1.0/(acphk(kc) - ak2(jc) - ak1(ic))
                    amkl(kc) = amphk(kc)*ackl_b
                    ackl(kc) = 1.0
                    apkl(kc) = apphk(kc)*ackl_b
                    rz1d(kc) = rz1(ic,jc,kc)*ackl_b/real(norm)
                    uflz(kc) = 0.0
                    vflz(kc) = 0.0
                end do

                gvalue    = (amkl(1)*apkl(nzm))**0.5

                ackl(1)   =  ackl(1)   - gvalue
                ackl(nzm) =  ackl(nzm) - amkl(1)*apkl(nzm)/gvalue

                uflz(1)   = gvalue
                uflz(nzm) = apkl(nzm)

                vflz(1)   = 1.0
                vflz(nzm) = amkl(1)/gvalue

                call dgttrf(nzm,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,rz1d,nzm,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,uflz,nzm,info)

                qflz = (vflz(1)*rz1d(1) + vflz(nzm)*rz1d(nzm))/(1.0 + (vflz(1)*uflz(1) + vflz(nzm)*uflz(nzm)))

                do kc = 1,nzm                
                    rz1(ic,jc,kc) = rz1d(kc) - (qflz*uflz(kc))
                enddo

            end do
        end do

    else

        do ic = dcmp%zst(1),dcmp%zen(1)
            do jc = dcmp%zst(2),dcmp%zen(2)

                do kc = 1,nzm
                    ackl_b   = 1.0/(acphk(kc) - ak2(jc) - ak1(ic))
                    amkl(kc) = amphk(kc)*ackl_b
                    ackl(kc) = 1.0
                    apkl(kc) = apphk(kc)*ackl_b
                    rz1d(kc) = rz1(ic,jc,kc)*ackl_b/real(norm)
                end do

                ackl(1)   = ackl(1)   + amkl(1)
                ackl(nzm) = ackl(nzm) + apkl(nzm)

                call dgttrf(nzm,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,info)
                call dgttrs('N',nzm,1,amkl(2:nzm),ackl(1:nzm),apkl(1:nzm-1),appk(1:nzm-2),ipkv,rz1d,nzm,info)

                do kc = 1,nzm
                    rz1(ic,jc,kc) = rz1d(kc)
                end do

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverZR