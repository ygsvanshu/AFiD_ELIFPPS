!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SourceRoutines.F90                             !
!    CONTAINS: subroutine CheckIsSourceLocal              !
!                                                         !
!    PURPOSE: Subroutines that perform operations on      !
!    single source                                        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckIsSourceLocal(s,n)

    use decomp_2d, only: xstart,xend
    use param, only: nxm,nym,nzm,xc,yc,zc,xlen,ylen,zlen,periodic
    use lagrangian_point_particle

    implicit none

    type(particle_source), intent(inout)    :: s
    logical, intent(out)                    :: n
    integer                                 :: i
    integer                                 :: ngrm(3)
    real                                    :: nlen(3),ncst(3),ncen(3)

    ngrm(1) = nxm
    ngrm(2) = nym
    ngrm(3) = nzm

    nlen(1) = xlen
    nlen(2) = ylen
    nlen(3) = zlen

    ncst(1) = xc(xstart(1))
    ncst(2) = yc(xstart(2))
    ncst(3) = zc(xstart(3))

    ncen(1) = xc(xend(1)+1)
    ncen(2) = yc(xend(2)+1)
    ncen(3) = zc(xend(3)+1)

    n = .true.

    do i = 1,3
        if (periodic(i)) then
            s%src_pos(i) = modulo(s%src_pos(i),nlen(i))
            n = n.and.((s%src_pos(i).ge.ncst(i)).and.(s%src_pos(i).lt.ncen(i)))
        else
            if (xend(i).eq.ngrm(i)) then
                n = n.and.((s%src_pos(i).ge.ncst(i)).and.(s%src_pos(i).le.ncen(i)))
            else
                n = n.and.((s%src_pos(i).ge.ncst(i)).and.(s%src_pos(i).lt.ncen(i)))
            end if
        end if
    end do

    return

end subroutine CheckIsSourceLocal