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
    real    :: be

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

                pr(ic,jc,kc) = pr(ic,jc,kc) + dphhalo(ic,jc,kc) 
                pr(ic,jc,kc) = pr(ic,jc,kc) - be*(dphhalo(im,jc,kc)*amphi(ic) + dphhalo(ic,jc,kc)*acphi(ic) + dphhalo(ip,jc,kc)*apphi(ic))
                pr(ic,jc,kc) = pr(ic,jc,kc) - be*(dphhalo(ic,jm,kc)*amphj(jc) + dphhalo(ic,jc,kc)*acphj(jc) + dphhalo(ic,jp,kc)*apphj(jc))
                pr(ic,jc,kc) = pr(ic,jc,kc) - be*(dphhalo(ic,jc,km)*amphk(kc) + dphhalo(ic,jc,kc)*acphk(kc) + dphhalo(ic,jc,kp)*apphk(kc))

            enddo
        enddo
    enddo

    if (periodic(1)) then
        pr(0 ,:,:) = pr(nxm,:,:)
        pr(nx,:,:) = pr(1  ,:,:)
    else
        pr(0 ,:,:) = pr(1  ,:,:)
        pr(nx,:,:) = pr(nxm,:,:)
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            pr(:,0 ,:) = pr(:,1  ,:)
        end if
        if (xend(2).eq.nym) then
            pr(:,ny,:) = pr(:,nym,:)
        end if
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            pr(:,:,0 ) = pr(:,:,1  )
        end if
        if (xend(3).eq.nzm) then
            pr(:,:,nz) = pr(:,:,nzm)
        end if
    end if
    
    return

end subroutine CorrectPressure
