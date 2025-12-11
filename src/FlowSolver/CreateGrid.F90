!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         !
!    PURPOSE: Compute the indices, grid, grid metrics     !
!     and coefficients for differentiation                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
! FOR HYPERBOLIC TANGENT CLUSTERING
subroutine CreateStretchedGrid1(nc,gc,ln)

    use param, only: ismaster,strval

    implicit none

    integer, intent(in) :: nc
    real, intent(inout) :: gc(nc)
    real, intent(in)    :: ln
    integer             :: lc

    do lc = 1,nc
        gc(lc) = (1.0 + tanh(strval*(real(2*lc-nc-1)/real(nc-1)))/tanh(strval))*0.5*ln
        if (gc(lc) .lt. 0 .or. gc(lc) .gt. ln) then
            if (ismaster) write (6,*) 'Grid is too streched'
            if (ismaster) write (6,*) 'Check stretching parameter'
            call MpiAbort
        end if
    end do

end subroutine CreateStretchedGrid1

! FOR CLIPPED CHEBYSHEV CLUSTERING
subroutine CreateStretchedGrid2(nc,gc,ln)

    use param, only: ismaster,strval,pi

    implicit none

    integer, intent(in)             :: nc
    real, intent(inout)             :: gc(nc)
    real, intent(in)                :: ln
    integer                         :: lc
    real                            :: sfac

    sfac = real(nc)/(real(nc) + 2.0*strval)
    do lc = 1,nc
        gc(lc) = 0.5*ln*(1.0 + (sin(pi*sfac*((real(lc-1)/real(nc-1)) - 0.5))/sin(0.5*pi*sfac)))
        if (gc(lc) .lt. 0 .or. gc(lc) .gt. ln) then
            if (ismaster) write (6,*) 'Grid is too streched'
            if (ismaster) write (6,*) 'Check stretching parameter'
            call MpiAbort
        end if 
    end do

end subroutine CreateStretchedGrid2

! FOR ERROR FUNCTION TYPE CLUSTERING
subroutine CreateStretchedGrid3(nc,gc,ln)

    use param, only: strval

    implicit none

    integer, intent(in)             :: nc
    real, intent(inout)             :: gc(nc)
    real, intent(in)                :: ln
    integer                         :: lc
    real                            :: gt

    do lc = 1,nc
        gt     = (real(lc - 1)/real(nc - 1)) - 0.5
        gt     = erf(gt*strval)/erf(0.5*strval)
        gc(lc) = 0.5*(1.0 + gt)*ln
    end do

end subroutine CreateStretchedGrid3

! TO READ IN GRID FROM INPUT FILE
subroutine ReadStretchedGrid(nc,gc,ln)

    use param

    implicit none

    integer, intent(in) :: nc
    real, intent(inout) :: gc(nc)
    real, intent(in)    :: ln
    integer             :: lc
    logical             :: exists
    real                :: strorg,strscl

    inquire (file='Inputs/sgrid.in', exist=exists)
    if (exists) then
        if (ismaster) write (6,*) 'Reading grid from sgrid.in'
        open (unit=78, file='Inputs/sgrid.in', status='old')
        do lc = 1,nc
            read (78,*) gc(lc)
        end do
        strscl = ln/(gc(nc) - gc(1))
        strorg = gc(1)
        do lc = 1,nc
            gc(lc) = (gc(lc) - strorg)*strscl
        end do
        close (78)
    else
        if (ismaster) write (6,*) 'Warning: Inputs/sgrid.in not found!'
        call MpiAbort
    end if

end subroutine ReadStretchedGrid

! Create grid operator coefficients, wavenumbers, metrics
subroutine CreateGrid

    use param
    use AuxiliaryRoutines

    implicit none

    integer :: i,j,k
    integer :: nxmh,nxmp
    integer :: nymh,nymp
    integer :: nzmh,nzmp

    !==============================================================!
    ! GRID FOR X-AXIS                                              !
    !==============================================================!

    if (straxs.eq.1) then
        if (strtyp.eq.1) then
            call CreateStretchedGrid1(nx,xc,xlen)
        else if (strtyp.eq.2) then
            call CreateStretchedGrid2(nx,xc,xlen)
        else if (strtyp.eq.3) then
            call CreateStretchedGrid3(nx,xc,xlen)
        else
            call ReadStretchedGrid(nx,xc,xlen)
        end if
    else
        ! UNIFORM GRID ALONG X-AXIS
        do i = 1,nx
            xc(i) = xlen*(real(i - 1)/real(nxm))
        end do
    end if

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

    ! MODIFIED WAVENUMBERS FOR X-GRID
    if (periodic(1)) then
        ! USE MODIFIED WAVENUMBERS FOR DFT

        nxmh = nxm/2 + 1
        nxmp = nxmh  + 1

        do i = 1, nxmh
            ao(i) = (i-1)*2.0*pi
        end do
        do i = nxmp, nxm
            ao(i) = -(nxm-i+1)*2.0*pi
        end do

    else
        ! USE MODIFIED WAVENUMBERS FOR DCT

        do i = 1,nxm
            ao(i) = (i-1)*pi
        end do

    end if

    do i = 1, nxm
        ak1(i) = 2.0*(1.0 - cos(ao(i)/nxm))*(float(nxm)/xlen)**2
    end do

    !==============================================================!
    ! GRID FOR Y-AXIS                                              !
    !==============================================================!                                                        

    if (straxs.eq.2) then
        if (strtyp.eq.1) then
            call CreateStretchedGrid1(ny,yc,ylen)
        else if (strtyp.eq.2) then
            call CreateStretchedGrid2(ny,yc,ylen)
        else if (strtyp.eq.3) then
            call CreateStretchedGrid3(ny,yc,ylen)
        else
            call ReadStretchedGrid(ny,yc,ylen)
        end if
    else
        ! UNIFORM GRID ALONG Y-AXIS
        do j = 1,ny
            yc(j) = ylen*(real(j-1)/real(nym))
        end do
    end if

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

    if (periodic(2)) then
        ! USE MODIFIED WAVENUMBERS FOR DFT
        
        nymh = nym/2 + 1
        nymp = nymh  + 1

        do j = 1, nymh
            ap(j) = (j-1)*2.0*pi
        end do
        do j = nymp, nym
            ap(j) = -(nym-j+1)*2.0*pi
        end do

    else
        ! USE MODIFIED WAVENUMBERS FOR DCT
        do j = 1,nym
            ap(j) = (j-1)*pi
        end do

    end if

    do j = 1, nym
        ak2(j) = 2.0*(1.0 - cos(ap(j)/nym))*(float(nym)/ylen)**2
    end do

    !==============================================================!
    ! GRID FOR Z-AXIS                                              !
    !==============================================================!

    if (straxs.eq.3) then
        if (strtyp.eq.1) then
            call CreateStretchedGrid1(nz,zc,zlen)
        else if (strtyp.eq.2) then
            call CreateStretchedGrid2(nz,zc,zlen)
        else if (strtyp.eq.3) then
            call CreateStretchedGrid3(nz,zc,zlen)
        else
            call ReadStretchedGrid(nz,zc,zlen)
        end if
    else
        ! UNIFORM GRID ALONG Z-AXIS
        do k = 1,nz
            zc(k) = zlen*(real(k-1)/real(nzm))
        end do
    end if

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

    if (periodic(3)) then
        ! USE MODIFIED WAVENUMBERS FOR DFT

        nzmh = nzm/2 + 1
        nzmp = nzmh  + 1

        do k = 1, nzmh
            aq(k) = (k-1)*2.0*pi
        end do
        do k = nzmp, nzm
            aq(k) = -(nzm-k+1)*2.0*pi
        end do

    else
        ! USE MODIFIED WAVENUMBERS FOR DCT

        do k = 1,nzm
            aq(k) = (k-1)*pi
        end do

    end if

    do k = 1, nzm
        ak3(k) = 2.0*(1.0 - cos(aq(k)/nzm))*(float(nzm)/zlen)**2
    end do

    ! WRITE OUT GRID INFORMATION

    if (ismaster) then

        open (unit=78, file='Results/xgrid.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,16(A24,X))') 'i', 'xc(i)', 'xm(i)', 'dxc(i)', 'dxm(i)', 'ixcm(i)', 'ixcp(i)', 'am1si(i)', 'ac1si(i)', 'ap1si(i)', 'am1ci(i)', 'ac1ci(i)', 'ap1ci(i)', 'amphi(i)', 'acphi(i)', 'apphi(i)', 'ak1(i)'
        do i = 1, nxm
            write (78, '(I8,X,16(ES24.16,X))') i, xc(i), xm(i), dxc(i), dxm(i), ixcm(i), ixcp(i), am1si(i), ac1si(i), ap1si(i), am1ci(i), ac1ci(i), ap1ci(i), amphi(i), acphi(i), apphi(i), ak1(i)
        end do
        close (78)

        open (unit=78, file='Results/ygrid.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,16(A24,X))') 'j', 'yc(j)', 'ym(j)', 'dyc(j)', 'dym(j)', 'iycm(j)', 'iycp(j)', 'am2sj(j)', 'ac2sj(j)', 'ap2sj(j)', 'am2cj(j)', 'ac2cj(j)', 'ap2cj(j)', 'amphj(j)', 'acphj(j)', 'apphj(j)', 'ak2(j)'
        do j = 1, nym
            write (78, '(I8,X,16(ES24.16,X))') j, yc(j), ym(j), dyc(j), dym(j), iycm(j), iycp(j), am2sj(j), ac2sj(j), ap2sj(j), am2cj(j), ac2cj(j), ap2cj(j), amphj(j), acphj(j), apphj(j), ak2(j)
        end do
        close (78)

        open (unit=78, file='Results/zgrid.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,16(A24,X))') 'k', 'zc(k)', 'zm(k)', 'dzc(k)', 'dzm(k)', 'izcm(k)', 'izcp(k)', 'am3sk(k)', 'ac3sk(k)', 'ap3sk(k)', 'am3ck(k)', 'ac3ck(k)', 'ap3ck(k)', 'amphk(k)', 'acphk(k)', 'apphk(k)', 'ak3(k)'
        do k = 1, nzm
            write (78, '(I8,X,16(ES24.16,X))') k, zc(k), zm(k), dzc(k), dzm(k), izcm(k), izcp(k), am3sk(k), ac3sk(k), ap3sk(k), am3ck(k), ac3ck(k), ap3ck(k), amphk(k), acphk(k), apphk(k), ak3(k)
        end do
        close (78)

    end if

    return

end subroutine CreateGrid