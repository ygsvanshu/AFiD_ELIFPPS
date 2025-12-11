
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ExplicitTermsVY.F90                            !
!    CONTAINS: subroutine ExplicitTermsVY                 !
!                                                         !
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the y (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVY

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,ady

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp

    real    :: vyvx_m,vyvx_p
    real    :: vyvy_m,vyvy_p
    real    :: vyvz_m,vyvz_p

    do kc=xstart(3),xend(3)
        km=kc-1
        kp=kc+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do ic=xstart(1),xend(1)
                im=ic-1
                ip=ic+1

                vyvx_m = ((iycm(jc)*vx(ic,jm,kc)) + (iycp(jc)*vx(ic,jc,kc))) * ((ixcm(ic)*vy(im,jc,kc)) + (ixcp(ic)*vy(ic,jc,kc)))
                vyvx_p = ((iycm(jc)*vx(ip,jm,kc)) + (iycp(jc)*vx(ip,jc,kc))) * ((ixcm(ip)*vy(ic,jc,kc)) + (ixcp(ip)*vy(ip,jc,kc)))
                
                vyvy_m = ((0.5*vy(ic,jm,kc)) + (0.5*vy(ic,jc,kc))) * ((0.5*vy(ic,jm,kc)) + (0.5*vy(ic,jc,kc)))
                vyvy_p = ((0.5*vy(ic,jc,kc)) + (0.5*vy(ic,jp,kc))) * ((0.5*vy(ic,jc,kc)) + (0.5*vy(ic,jp,kc)))
                
                vyvz_m = ((iycm(jc)*vz(ic,jm,kc)) + (iycp(jc)*vz(ic,jc,kc))) * ((izcm(kc)*vy(ic,jc,km)) + (izcp(kc)*vy(ic,jc,kc)))
                vyvz_p = ((iycm(jc)*vz(ic,jm,kp)) + (iycp(jc)*vz(ic,jc,kp))) * ((izcm(kp)*vy(ic,jc,kc)) + (izcp(kp)*vy(ic,jc,kp)))

                ady(ic,jc,kc) = 0.0
                ady(ic,jc,kc) = ady(ic,jc,kc) - (vyvx_p - vyvx_m)/dxc(ic)
                ady(ic,jc,kc) = ady(ic,jc,kc) - (vyvy_p - vyvy_m)/dym(jc)
                ady(ic,jc,kc) = ady(ic,jc,kc) - (vyvz_p - vyvz_m)/dzc(kc)

            enddo
        enddo
    enddo

    return
    
end subroutine ExplicitTermsVY
  