!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ImplicitAndUpdateVX.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVX             !
!                                                         !
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the X (vertical) direction and call !
!     the implicit solver. After this routine, the        !
!     vertical velocity has been updated to the new       !
!     timestep                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVX

    use param
    use local_arrays, only : vx,pr,adx,dfx,rux,rhx
    use decomp_2d, only : xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp

    real    :: alre,prx
    real    :: dxxvx,dyyvx,dzzvx

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

                    prx = (pr(ic,jc,kc) - pr(im,jc,kc))/dxm(ic)

                    rhx(ic,jc,kc) = (ga*adx(ic,jc,kc) + ro*rux(ic,jc,kc) + alre*dfx(ic,jc,kc) - al*prx)*dt
                    rux(ic,jc,kc) = adx(ic,jc,kc)
                    
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

                    dxxvx = vx(ip,jc,kc)*ap1si(ic) + vx(ic,jc,kc)*ac1si(ic) + vx(im,jc,kc)*am1si(ic)
                    dyyvx = vx(ic,jp,kc)*ap2cj(jc) + vx(ic,jc,kc)*ac2cj(jc) + vx(ic,jm,kc)*am2cj(jc)
                    dzzvx = vx(ic,jc,kp)*ap3ck(kc) + vx(ic,jc,kc)*ac3ck(kc) + vx(ic,jc,km)*am3ck(kc)

                    prx = (pr(ic,jc,kc) - pr(im,jc,kc))/dxm(ic)

                    rhx(ic,jc,kc) = (ga*adx(ic,jc,kc) + ro*rux(ic,jc,kc) + alre*(dxxvx + dyyvx + dzzvx) - al*prx)*dt
                    
                    rux(ic,jc,kc) = adx(ic,jc,kc)
                    
                end do
            end do
        end do

    end if

    call AddBoundaryRHSTermsVX
    call SolveImpEqnUpdateVX
    call ImposeInternalBoundaryVX

    return

end subroutine ImplicitAndUpdateVX
