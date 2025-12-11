!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ImplicitAndUpdateVY.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVY             !
!                                                         !
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the y (horizontal) dimension        !
!     and call the implicit solver                        !
!     After this routine, the velocity field in y has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVY

    use param
    use local_arrays, only : vy,pr,ady,dfy,ruy,rhx
    use decomp_2d, only : xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp

    real    :: alre,pry
    real    :: dxxvy,dyyvy,dzzvy

    alre = al/rey
    rhx(:,:,:) = 0.0

    if (pre_diff) then

        do kc=xstart(3),xend(3)
            km=kc-1
            kp=kc+1
            do jc=xstart(2),xend(2)
                jm=jc-1
                jp=jc+1
                do ic=xstart(1),xend(1)
                    im=ic-1
                    ip=ic+1

                    pry = (pr(ic,jc,kc) - pr(ic,jm,kc))/dym(jc)

                    rhx(ic,jc,kc) = (ga*ady(ic,jc,kc) + ro*ruy(ic,jc,kc) + alre*dfy(ic,jc,kc) - al*pry)*dt
                    ruy(ic,jc,kc) = ady(ic,jc,kc)

                end do
            end do
        end do

    else

        do kc=xstart(3),xend(3)
            km=kc-1
            kp=kc+1
            do jc=xstart(2),xend(2)
                jm=jc-1
                jp=jc+1
                do ic=xstart(1),xend(1)
                    im=ic-1
                    ip=ic+1

                    dxxvy = vy(ip,jc,kc)*ap1ci(ic) + vy(ic,jc,kc)*ac1ci(ic) + vy(im,jc,kc)*am1ci(ic)
                    dyyvy = vy(ic,jp,kc)*ap2sj(jc) + vy(ic,jc,kc)*ac2sj(jc) + vy(ic,jm,kc)*am2sj(jc)
                    dzzvy = vy(ic,jc,kp)*ap3ck(kc) + vy(ic,jc,kc)*ac3ck(kc) + vy(ic,jc,km)*am3ck(kc)

                    pry = (pr(ic,jc,kc) - pr(ic,jm,kc))/dym(jc)

                    rhx(ic,jc,kc) = (ga*ady(ic,jc,kc) + ro*ruy(ic,jc,kc) + alre*(dxxvy + dyyvy + dzzvy) - al*pry)*dt

                    ruy(ic,jc,kc) = ady(ic,jc,kc)

                end do
            end do
        end do

    end if

    call AddBoundaryRHSTermsVY
    call SolveImpEqnUpdateVY
    call ImposeInternalBoundaryVY

    return

end subroutine ImplicitAndUpdateVY