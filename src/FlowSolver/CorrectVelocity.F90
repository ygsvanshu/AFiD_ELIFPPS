!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectVelocity.F90                            !
!    CONTAINS: subroutine CorrectVelocity                 !
!                                                         ! 
!    PURPOSE: Update velocities with the pressure         !
!     correction to enforce incompresibility              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectVelocity

    use param
    use local_arrays, only: vy,vx,vz,dphhalo
    use decomp_2d

    implicit none

    integer :: km,kc
    integer :: jm,jc
    integer :: im,ic

    integer :: xst1,xen1
    integer :: xst2,xen2
    integer :: xst3,xen3

    xst1 = 1
    xen1 = nx
    xst2 = xstart(2)
    xen2 = xend(2)
    xst3 = xstart(3)
    xen3 = xend(3)

    do ic = xst1,xen1
        im = ic-1
        vx(ic,xst2:xen2,xst3:xen3) = vx(ic,xst2:xen2,xst3:xen3) - (dphhalo(ic,xst2:xen2,xst3:xen3)-dphhalo(im,xst2:xen2,xst3:xen3))*al*dt/dxm(ic)
    end do

    call ImposeExternalBoundaryVX

    xst1 = 1
    xen1 = nxm
    xst2 = xstart(2)
    xen2 = xend(2)
    xst3 = xstart(3)
    xen3 = xend(3)

    if (xend(2).eq.nym) xen2 = ny

    do jc = xst2,xen2
        jm = jc-1
        vy(xst1:xen1,jc,xst3:xen3) = vy(xst1:xen1,jc,xst3:xen3) - (dphhalo(xst1:xen1,jc,xst3:xen3)-dphhalo(xst1:xen1,jm,xst3:xen3))*al*dt/dym(jc)
    end do

    call ImposeExternalBoundaryVY

    xst1 = 1
    xen1 = nxm
    xst2 = xstart(2)
    xen2 = xend(2)
    xst3 = xstart(3)
    xen3 = xend(3)

    if (xend(3).eq.nzm) xen3 = nz

    do kc = xst3,xen3
        km = kc-1
        vz(xst1:xen1,xst2:xen2,kc) = vz(xst1:xen1,xst2:xen2,kc) - (dphhalo(xst1:xen1,xst2:xen2,kc)-dphhalo(xst1:xen1,xst2:xen2,km))*al*dt/dzm(kc)
    end do

    call ImposeExternalBoundaryVZ

    return

end subroutine CorrectVelocity
