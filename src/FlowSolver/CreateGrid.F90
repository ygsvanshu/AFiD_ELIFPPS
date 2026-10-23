!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateGrid.F90                                 !
!    CONTAINS: subroutine CreateGrid                      !
!                                                         !
!    PURPOSE: Create/Read in stretched grids based on     !
!    user input selection                                 !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateGrid

    use param

    !==============================================================!
    ! GRID FOR X-AXIS                                              !
    !==============================================================!

    if (gtpx.eq.0) then
        ! UNIFORM GRID ALONG X-AXIS
        do i = 1,nx
            xc(i) = xlen*(real(i - 1)/real(nxm))
        end do
    else if (gtpx.eq.1) then
        call CreateStretchedGrid1(nx,xc,xlen,gvlx)
    else if (gtpx.eq.2) then
        call CreateStretchedGrid2(nx,xc,xlen,gvlx)
    else if (gtpx.eq.3) then
        call CreateStretchedGrid3(nx,xc,xlen,gvlx)
    else
        call ReadStretchedGrid(nx,xc,xlen,'x')
    end if

    !==============================================================!
    ! GRID FOR Y-AXIS                                              !
    !==============================================================!                                                        

    if (gtpy.eq.0) then
        ! UNIFORM GRID ALONG Y-AXIS
        do j = 1,ny
            yc(j) = ylen*(real(j-1)/real(nym))
        end do
    else if (gtpy.eq.1) then
        call CreateStretchedGrid1(ny,yc,ylen,gvly)
    else if (gtpy.eq.2) then
        call CreateStretchedGrid2(ny,yc,ylen,gvly)
    else if (gtpy.eq.3) then
        call CreateStretchedGrid3(ny,yc,ylen,gvly)
    else
        call ReadStretchedGrid(ny,yc,ylen,'y')
    end if

    !==============================================================!
    ! GRID FOR Z-AXIS                                              !
    !==============================================================!

    if (gtpz.eq.0) then
        ! UNIFORM GRID ALONG Z-AXIS
        do k = 1,nz
            zc(k) = zlen*(real(k-1)/real(nzm))
        end do
    else if (gtpz.eq.1) then
        call CreateStretchedGrid1(nz,zc,zlen,gvlz)
    else if (gtpz.eq.2) then
        call CreateStretchedGrid2(nz,zc,zlen,gvlz)
    else if (gtpz.eq.3) then
        call CreateStretchedGrid3(nz,zc,zlen,gvlz)
    else
        call ReadStretchedGrid(nz,zc,zlen,'z')
    end if

end subroutine CreateGrid