!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ImplicitAndUpdateVZ.F90                        !
!    CONTAINS: subroutine ImplicitAndUpdateVZ             !
!                                                         !
!    PURPOSE: Compute the linear terms associated to      !
!     the velocity in the z (horizontal) dimension        !
!     and call the implicit solver.                       !
!     After this routine, the velocity field in z has     !
!     been updated to the new timestep                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ImplicitAndUpdateVZ

    use param
    use local_arrays, only : vz,pr,adz,dfz,ruz,rhx
    use decomp_2d, only : xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp

    real    :: alre,prz
    real    :: dxxvz,dyyvz,dzzvz

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

                    prz = (pr(ic,jc,kc) - pr(ic,jc,km))/dzm(kc)

                    rhx(ic,jc,kc) = (ga*adz(ic,jc,kc) + ro*ruz(ic,jc,kc) + alre*dfz(ic,jc,kc) - al*prz)*dt
                    ruz(ic,jc,kc) = adz(ic,jc,kc)

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

                    dxxvz = vz(ip,jc,kc)*ap1ci(ic) + vz(ic,jc,kc)*ac1ci(ic) + vz(im,jc,kc)*am1ci(ic)
                    dyyvz = vz(ic,jp,kc)*ap2cj(jc) + vz(ic,jc,kc)*ac2cj(jc) + vz(ic,jm,kc)*am2cj(jc)
                    dzzvz = vz(ic,jc,kp)*ap3sk(kc) + vz(ic,jc,kc)*ac3sk(kc) + vz(ic,jc,km)*am3sk(kc)

                    prz = (pr(ic,jc,kc) - pr(ic,jc,km))/dzm(kc)

                    rhx(ic,jc,kc) = (ga*adz(ic,jc,kc) + ro*ruz(ic,jc,kc) + alre*(dxxvz + dyyvz + dzzvz) - al*prz)*dt
                    
                    ruz(ic,jc,kc) = adz(ic,jc,kc)

                end do
            end do
        end do

    end if

    call AddBoundaryRHSTermsVZ
    call SolveImpEqnUpdateVZ
    call ImposeInternalBoundaryVZ
    
    return

end subroutine ImplicitAndUpdateVZ