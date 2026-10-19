!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectPressure.F90                            !
!    CONTAINS: subroutine CorrectPressure                 !
!                                                         ! 
!    PURPOSE: Apply the pressure correction to the        !
!     pressure                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectPressure

    use param
    use local_arrays, only: pr,dphhalo
    use decomp_2d

    implicit none

    integer :: km,kc,kp
    integer :: jm,jc,jp
    integer :: im,ic,ip
    real    :: be,dp

    be = 0.5*al*dt/rey

    do kc=xstart(3),xend(3)
        km=kc-1
        kp=kc+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do ic=1,nxm
                ip=ic+1
                im=ic-1

                dp = dphhalo(ic,jc,kc)
                dp = dp - be*(dphhalo(im,jc,kc)*amphi(ic) + dphhalo(ic,jc,kc)*acphi(ic) + dphhalo(ip,jc,kc)*apphi(ic))
                dp = dp - be*(dphhalo(ic,jm,kc)*amphj(jc) + dphhalo(ic,jc,kc)*acphj(jc) + dphhalo(ic,jp,kc)*apphj(jc))
                dp = dp - be*(dphhalo(ic,jc,km)*amphk(kc) + dphhalo(ic,jc,kc)*acphk(kc) + dphhalo(ic,jc,kp)*apphk(kc))

                pr(ic,jc,kc) = pr(ic,jc,kc) + dp

                if (.not.periodic(1)) then
                    if (ic.eq.1)   pr(0 ,jc,kc) = pr(0 ,jc,kc) + dp + (dphhalo(0 ,jc,kc) - dphhalo(1  ,jc,kc))
                    if (ic.eq.nxm) pr(nx,jc,kc) = pr(nx,jc,kc) + dp + (dphhalo(nx,jc,kc) - dphhalo(nxm,jc,kc))
                end if

                if (.not.periodic(2)) then
                    if (jc.eq.1)   pr(ic,0 ,kc) = pr(ic,0 ,kc) + dp + (dphhalo(ic,0 ,kc) - dphhalo(ic,1  ,kc))
                    if (jc.eq.nym) pr(ic,ny,kc) = pr(ic,ny,kc) + dp + (dphhalo(ic,ny,kc) - dphhalo(ic,nym,kc))
                end if

                if (.not.periodic(3)) then
                    if (kc.eq.1)   pr(ic,jc,0 ) = pr(ic,jc,0 ) + dp + (dphhalo(ic,jc,0 ) - dphhalo(ic,jc,1  ))
                    if (kc.eq.nzm) pr(ic,jc,nz) = pr(ic,jc,nz) + dp + (dphhalo(ic,jc,nz) - dphhalo(ic,jc,nzm))
                end if

            enddo
        enddo
    enddo

    ! Halo exchange doesn't handle the pencils in x-direction, manually corrected.
    if (periodic(1)) then
        pr(0 ,:,:) = pr(nxm,:,:)
        pr(nx,:,:) = pr(1  ,:,:)
    end if

    return

end subroutine CorrectPressure
