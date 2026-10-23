!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: StretchedGrid.F90                              !
!    CONTAINS: subroutine CreateStretchedGrid1            !
!              subroutine CreateStretchedGrid2            !
!              subroutine CreateStretchedGrid3            !
!              subroutine ReadStretchedGrid               !
!                                                         !
!    PURPOSE: Create/Read in stretched grids              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! FOR HYPERBOLIC TANGENT CLUSTERING
subroutine CreateStretchedGrid1(nc,gc,ln,sv)

    use param, only: ismaster

    implicit none

    integer, intent(in) :: nc
    real, intent(inout) :: gc(nc)
    real, intent(in)    :: ln
    real, intent(in)    :: sv
    
    integer             :: lc

    do lc = 1,nc
        gc(lc) = (1.0 + tanh(sv*(real(2*lc-nc-1)/real(nc-1)))/tanh(sv))*0.5*ln
        if (gc(lc) .lt. 0 .or. gc(lc) .gt. ln) then
            if (ismaster) write (6,*) 'Grid is too streched'
            if (ismaster) write (6,*) 'Check stretching parameter'
            call MpiAbort
        end if
    end do

end subroutine CreateStretchedGrid1

! FOR CLIPPED CHEBYSHEV CLUSTERING
subroutine CreateStretchedGrid2(nc,gc,ln,sv)

    use param, only: ismaster,pi

    implicit none

    integer, intent(in) :: nc
    real, intent(inout) :: gc(nc)
    real, intent(in)    :: ln
    real, intent(in)    :: sv

    integer             :: lc
    real                :: sfac

    sfac = real(nc)/(real(nc) + 2.0*sv)
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
subroutine CreateStretchedGrid3(nc,gc,ln,sv)

    implicit none

    integer, intent(in) :: nc
    real, intent(inout) :: gc(nc)
    real, intent(in)    :: ln
    real, intent(in)    :: sv

    integer             :: lc
    real                :: gt

    do lc = 1,nc
        gt     = (real(lc - 1)/real(nc - 1)) - 0.5
        gt     = erf(gt*sv)/erf(0.5*sv)
        gc(lc) = 0.5*(1.0 + gt)*ln
    end do

end subroutine CreateStretchedGrid3

! TO READ IN GRID FROM INPUT FILE
subroutine ReadStretchedGrid(nc,gc,ln,ax)

    use param

    implicit none

    integer, intent(in) :: nc
    real, intent(inout) :: gc(nc)
    real, intent(in)    :: ln
    integer             :: lc
    logical             :: exists
    real                :: strorg,strscl
    character           :: ax

    inquire (file='Inputs/s'//ax//'grid.in', exist=exists)
    if (exists) then
        if (ismaster) write (6,*) 'Reading grid from s'//ax//'grid.in'
        open (unit=78, file='Inputs/s'//ax//'grid.in', status='old')
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
        if (ismaster) write (6,*) 'Warning: Inputs/s'//ax//'grid.in not found!'
        call MpiAbort
    end if

end subroutine ReadStretchedGrid