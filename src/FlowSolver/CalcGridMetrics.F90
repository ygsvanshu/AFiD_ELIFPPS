!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CalcGridMetrics.F90                            !
!    CONTAINS: subroutine CalcGridMetrics                 !
!                                                         !
!    PURPOSE: Compute the grid spacing, grid metrics      !
!    and coefficients for discrete differential operators !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcGridMetrics

    use param
    use AuxiliaryRoutines

    implicit none

    integer :: i,j,k

    !==============================================================!
    ! GRID FOR X-AXIS                                              !
    !==============================================================!

    ! CELL CENTER X-GRID
    do i = 1,nxm
        xm(i) = (xc(i) + xc(i+1))*0.5
    end do

    ! SPACING DISTANCE METRICS FOR X-GRID

    do i = 1,nxm
        dxc(i) = (xc(i+1) - xc(i))
    end do

    if (periodic(1)) then
        dxc(0)  = dxc(nxm)
        dxc(nx) = dxc(1)
    else
        dxc(0)  = dxc(1)
        dxc(nx) = dxc(nxm)
    end if

    xm(0)  = xc(1)  - 0.5*dxc(0)
    xm(nx) = xc(nx) + 0.5*dxc(nx)

    do i = 1,nx
        dxm(i) = xm(i) - xm(i-1)
    end do

    ! INTERPOLATION COEFFICIENTS FOR X-GRID
    do i = 1,nx
        ixcm(i) = 0.5*dxc(i-1)/dxm(i)
        ixcp(i) = 0.5*dxc(i  )/dxm(i)
    end do 

    ! 2ND ORDER DERIVATIVE COEFFICIENTS FOR XC-GRID
    do i = 1,nx
        am1si(i) =  2.0/(dxc(i-1)*(dxc(i) + dxc(i-1)))
        ac1si(i) = -2.0/(dxc(i)*dxc(i-1))
        ap1si(i) =  2.0/(dxc(i)*(dxc(i) + dxc(i-1)))
    end do

    ! 2ND ORDER DERIVATIVE COEFFICIENTS FOR XM-GRID
    do i = 1,nxm
        am1ci(i) =  2.0/(dxm(i)*(dxm(i+1) + dxm(i)))
        ac1ci(i) = -2.0/(dxm(i+1)*dxm(i))
        ap1ci(i) =  2.0/(dxm(i+1)*(dxm(i+1) + dxm(i)))
    end do

    ! PRESSURE SOLVER COEFFICIENTS FOR X-GRID
    do i = 1,nxm
        amphi(i) =  (1.0/dxc(i))*(1.0/dxm(i))
        acphi(i) = -(1.0/dxc(i))*((1.0/dxm(i)) + (1.0/dxm(i+1)))
        apphi(i) =  (1.0/dxc(i))*(1.0/dxm(i+1))
    end do
    
    !==============================================================!
    ! GRID FOR Y-AXIS                                              !
    !==============================================================!

    ! CELL CENTER Y-GRID
    do j = 1,nym
        ym(j) = (yc(j) + yc(j+1))*0.5
    end do

    ! SPACING DISTANCE METRICS FOR Y-GRID

    do j = 1,nym
        dyc(j) = (yc(j+1) - yc(j))
    end do

    if (periodic(2)) then
        dyc(0)  = dyc(nym)
        dyc(ny) = dyc(1)
    else
        dyc(0)  = dyc(1)
        dyc(ny) = dyc(nym)
    end if

    ym(0)  = yc(1)  - 0.5*dyc(0)
    ym(ny) = yc(ny) + 0.5*dyc(ny)

    do j = 1,ny
        dym(j) = ym(j) - ym(j-1)
    end do

    ! INTERPOLATION COEFFICIENTS FOR Y-GRID
    do j = 1,ny
        iycm(j) = 0.5*dyc(j-1)/dym(j)
        iycp(j) = 0.5*dyc(j  )/dym(j)
    end do 

    ! 2ND ORDER DERIVATIVE COEFFICIENTS FOR YC-GRID
    do j = 1,ny
        am2sj(j) =  2.0/(dyc(j-1)*(dyc(j) + dyc(j-1)))
        ac2sj(j) = -2.0/(dyc(j)*dyc(j-1))
        ap2sj(j) =  2.0/(dyc(j)*(dyc(j) + dyc(j-1)))
    end do

    ! 2ND ORDER DERIVATIVE COEFFICIENTS FOR YM-GRID
    do j = 1,nym
        am2cj(j) =  2.0/(dym(j)*(dym(j+1) + dym(j)))
        ac2cj(j) = -2.0/(dym(j+1)*dym(j))
        ap2cj(j) =  2.0/(dym(j+1)*(dym(j+1) + dym(j)))
    end do

    ! PRESSURE SOLVER COEFFICIENTS FOR Y-GRID
    do j = 1,nym
        amphj(j) =  (1.0/dyc(j))*(1.0/dym(j))
        acphj(j) = -(1.0/dyc(j))*((1.0/dym(j)) + (1.0/dym(j+1)))
        apphj(j) =  (1.0/dyc(j))*(1.0/dym(j+1))
    end do

    !==============================================================!
    ! GRID FOR Z-AXIS                                              !
    !==============================================================!

    ! CELL CENTER Z-GRID
    do k = 1,nzm
        zm(k) = (zc(k) + zc(k+1))*0.5
    end do

    ! SPACING DISTANCE METRICS FOR Z-GRID

    do k = 1,nzm
        dzc(k) = (zc(k+1) - zc(k))
    end do

    if (periodic(3)) then
        dzc(0)  = dzc(nzm)
        dzc(nz) = dzc(1)
    else
        dzc(0)  = dzc(1)
        dzc(nz) = dzc(nzm)
    end if

    zm(0)  = zc(1)  - 0.5*dzc(0)
    zm(nz) = zc(nz) + 0.5*dzc(nz)

    do k = 1,nz
        dzm(k) = zm(k) - zm(k-1)
    end do

    ! INTERPOLATION COEFFICIENTS FOR Z-GRID
    do k = 1,nz
        izcm(k) = 0.5*dzc(k-1)/dzm(k)
        izcp(k) = 0.5*dzc(k  )/dzm(k)
    end do 

    ! 2ND ORDER DERIVATIVE COEFFICIENTS FOR ZC-GRID
    do k = 1,nz
        am3sk(k) =  2.0/(dzc(k-1)*(dzc(k) + dzc(k-1)))
        ac3sk(k) = -2.0/(dzc(k)*dzc(k-1))
        ap3sk(k) =  2.0/(dzc(k)*(dzc(k) + dzc(k-1)))
    end do

    ! 2ND ORDER DERIVATIVE COEFFICIENTS FOR ZM-GRID
    do k = 1,nzm
        am3ck(k) =  2.0/(dzm(k)*(dzm(k+1) + dzm(k)))
        ac3ck(k) = -2.0/(dzm(k+1)*dzm(k))
        ap3ck(k) =  2.0/(dzm(k+1)*(dzm(k+1) + dzm(k)))
    end do

    ! PRESSURE SOLVER COEFFICIENTS FOR Z-GRID
    do k = 1,nzm
        amphk(k) =  (1.0/dzc(k))*(1.0/dzm(k))
        acphk(k) = -(1.0/dzc(k))*((1.0/dzm(k)) + (1.0/dzm(k+1)))
        apphk(k) =  (1.0/dzc(k))*(1.0/dzm(k+1))
    end do

    ! WRITE OUT GRID INFORMATION

    if (ismaster) then

        open (unit=78, file='Results/xgrid.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,15(A24,X))') 'i', 'xc(i)', 'xm(i)', 'dxc(i)', 'dxm(i)', 'ixcm(i)', 'ixcp(i)', 'am1si(i)', 'ac1si(i)', 'ap1si(i)', 'am1ci(i)', 'ac1ci(i)', 'ap1ci(i)', 'amphi(i)', 'acphi(i)', 'apphi(i)'
        do i = 1, nxm
            write (78, '(I8,X,15(ES24.16,X))') i, xc(i), xm(i), dxc(i), dxm(i), ixcm(i), ixcp(i), am1si(i), ac1si(i), ap1si(i), am1ci(i), ac1ci(i), ap1ci(i), amphi(i), acphi(i), apphi(i)
        end do
        close (78)

        open (unit=78, file='Results/ygrid.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,15(A24,X))') 'j', 'yc(j)', 'ym(j)', 'dyc(j)', 'dym(j)', 'iycm(j)', 'iycp(j)', 'am2sj(j)', 'ac2sj(j)', 'ap2sj(j)', 'am2cj(j)', 'ac2cj(j)', 'ap2cj(j)', 'amphj(j)', 'acphj(j)', 'apphj(j)'
        do j = 1, nym
            write (78, '(I8,X,15(ES24.16,X))') j, yc(j), ym(j), dyc(j), dym(j), iycm(j), iycp(j), am2sj(j), ac2sj(j), ap2sj(j), am2cj(j), ac2cj(j), ap2cj(j), amphj(j), acphj(j), apphj(j)
        end do
        close (78)

        open (unit=78, file='Results/zgrid.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,15(A24,X))') 'k', 'zc(k)', 'zm(k)', 'dzc(k)', 'dzm(k)', 'izcm(k)', 'izcp(k)', 'am3sk(k)', 'ac3sk(k)', 'ap3sk(k)', 'am3ck(k)', 'ac3ck(k)', 'ap3ck(k)', 'amphk(k)', 'acphk(k)', 'apphk(k)'
        do k = 1, nzm
            write (78, '(I8,X,15(ES24.16,X))') k, zc(k), zm(k), dzc(k), dzm(k), izcm(k), izcp(k), am3sk(k), ac3sk(k), ap3sk(k), am3ck(k), ac3ck(k), ap3ck(k), amphk(k), acphk(k), apphk(k)
        end do
        close (78)

    end if

    return

end subroutine CalcGridMetrics