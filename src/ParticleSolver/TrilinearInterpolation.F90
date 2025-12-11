!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: TrilinearInterpolation.F90                     !
!    CONTAINS: subroutine TrilinearInterpolation          !
!                                                         !
!    PURPOSE: To compute trilinear interpolation between  !
!    Eulerian and Lagrangian frames of reference.         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetLocationCellIndex(gpos,grid,g_st,g_en,indx)

    ! IMPORTANT: Ensure that gpos>=grid(g_st) and gpos<grid(g_en)

    implicit none

    real, intent(in)        :: gpos
    real, intent(in)        :: grid(g_st:g_en)
    integer, intent(in)     :: g_st,g_en
    integer, intent(inout)  :: indx
    logical                 :: gloc

    ! For uninitialized particles or particles transported from periodic boundaries
    indx = max(min(indx,g_en-1),g_st)

    ! Compute only if particle is not still in the same cell
    if (.not.((gpos.ge.grid(indx)).and.(gpos.lt.grid(indx+1)))) then
        ! If particle has moved forward
        if (gpos.ge.grid(indx+1)) then
            ! Start searching from next coordinate index
            indx = indx + 1
            gloc = .false.
            do while ((indx.lt.g_en).and.(.not.gloc))
                if ((gpos.ge.grid(indx)).and.(gpos.lt.grid(indx+1))) then
                    ! Particle position located in the current cell
                    gloc = .true.
                else
                    ! Particle position not located in current cell, increment cell index
                    indx = indx + 1
                end if
            end do
        ! If particle has moved backward
        else if (gpos.lt.grid(indx)) then
            ! Start searching from previous coordinate index
            indx = indx - 1
            gloc = .false.
            do while ((indx.lt.g_en).and.(.not.gloc))
                if ((gpos.ge.grid(indx)).and.(gpos.lt.grid(indx+1))) then
                    ! Particle position located in the current cell
                    gloc = .true.
                else
                    ! Particle position not located in current cell, decrement cell index
                    indx = indx - 1
                end if
            end do
        end if
    end if

    return

end subroutine GetLocationCellIndex

subroutine CalcTrilinearInterpolationCoefficients(indm,indc,ipos,icfm,icfc)

    use param, only: xc,xm,yc,ym,zc,zm

    implicit none

    integer, intent(in)     :: indm(3)
    integer, intent(in)     :: indc(3)
    real,    intent(in)     :: ipos(3)
    real,    intent(out)    :: icfm(3)
    real,    intent(out)    :: icfc(3)

    icfc(1) = (xc(indc(1)+1) - ipos(1))/(xc(indc(1)+1) - xc(indc(1)))
    icfm(1) = (xm(indm(1)+1) - ipos(1))/(xm(indm(1)+1) - xm(indm(1)))

    icfc(2) = (yc(indc(2)+1) - ipos(2))/(yc(indc(2)+1) - yc(indc(2)))
    icfm(2) = (ym(indm(2)+1) - ipos(2))/(ym(indm(2)+1) - ym(indm(2)))

    icfc(3) = (zc(indc(3)+1) - ipos(3))/(zc(indc(3)+1) - zc(indc(3)))
    icfm(3) = (zm(indm(3)+1) - ipos(3))/(zm(indm(3)+1) - zm(indm(3)))

end subroutine CalcTrilinearInterpolationCoefficients

subroutine ApplyTrilinearInterpolation(indm,indc,icfm,icfc,grid,idir,qarr,qval)

    use decomp_2d, only: xstart,xend
    use param, only: lvlhalo

    implicit none

    integer, intent(in)     :: indm(3)
    integer, intent(in)     :: indc(3)
    real, intent(in)        :: icfm(3)
    real, intent(in)        :: icfc(3)
    character, intent(in)   :: grid
    integer, intent(in)     :: idir
    real, intent(inout)     :: qarr(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo)
    real, intent(inout)     :: qval

    character               :: glst(3)
    integer                 :: i,j,k
    integer                 :: ival
    integer                 :: indq(3)
    real                    :: coef(3,0:1)

    glst = (/'x','y','z'/)

    do ival = 1,3

        if (grid.eq.glst(ival)) then
            indq(ival)   = indc(ival)
            coef(ival,0) = icfc(ival)
            coef(ival,1) = 1.0 - icfc(ival)
        else
            indq(ival)   = indm(ival)
            coef(ival,0) = icfm(ival)
            coef(ival,1) = 1.0 - icfm(ival)
        end if

    end do

    if (idir.gt.0) then

        do i = 0,1
            do j = 0,1
                do k = 0,1
                    qval = qval + coef(1,i)*coef(2,j)*coef(3,k)*qarr(indq(1)+i,indq(2)+j,indq(3)+k)
                end do
            end do
        end do

    end if

    if (idir.lt.0) then

        do i = 0,1
            do j = 0,1
                do k = 0,1
                    qarr(indq(1)+i,indq(2)+j,indq(3)+k) = qarr(indq(1)+i,indq(2)+j,indq(3)+k) + coef(1,i)*coef(2,j)*coef(3,k)*qval
                end do
            end do
        end do

    end if

end subroutine ApplyTrilinearInterpolation