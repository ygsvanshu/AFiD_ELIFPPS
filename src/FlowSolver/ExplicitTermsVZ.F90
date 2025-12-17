
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ExplicitTermsVZ.F90                            !
!    CONTAINS: subroutine ExplicitTermsVZ                 !
!                                                         !
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the z (horizontal) dimension        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVZ

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,adz

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp

    real    :: vzvx_m,vzvx_p
    real    :: vzvy_m,vzvy_p
    real    :: vzvz_m,vzvz_p

    do kc=xstart(3),xend(3)
        km=kc-1
        kp=kc+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do ic=xstart(1),xend(1)
                im=ic-1
                ip=ic+1

                vzvx_m = ((izcm(kc)*vx(ic,jc,km)) + (izcp(kc)*vx(ic,jc,kc))) * ((ixcm(ic)*vz(im,jc,kc)) + (ixcp(ic)*vz(ic,jc,kc)))
                vzvx_p = ((izcm(kc)*vx(ip,jc,km)) + (izcp(kc)*vx(ip,jc,kc))) * ((ixcm(ip)*vz(ic,jc,kc)) + (ixcp(ip)*vz(ip,jc,kc)))
                
                vzvy_m = ((izcm(kc)*vy(ic,jc,km)) + (izcp(kc)*vy(ic,jc,kc))) * ((iycm(jc)*vz(ic,jm,kc)) + (iycp(jc)*vz(ic,jc,kc)))
                vzvy_p = ((izcm(kc)*vy(ic,jp,km)) + (izcp(kc)*vy(ic,jp,kc))) * ((iycm(jp)*vz(ic,jc,kc)) + (iycp(jp)*vz(ic,jp,kc)))

                vzvz_m = ((0.5*vz(ic,jc,km)) + (0.5*vz(ic,jc,kc))) * ((0.5*vz(ic,jc,km)) + (0.5*vz(ic,jc,kc)))
                vzvz_p = ((0.5*vz(ic,jc,kc)) + (0.5*vz(ic,jc,kp))) * ((0.5*vz(ic,jc,kc)) + (0.5*vz(ic,jc,kp)))

                adz(ic,jc,kc) = 0.0
                adz(ic,jc,kc) = adz(ic,jc,kc) - (vzvx_p - vzvx_m)/dxc(ic)
                adz(ic,jc,kc) = adz(ic,jc,kc) - (vzvy_p - vzvy_m)/dyc(jc)
                adz(ic,jc,kc) = adz(ic,jc,kc) - (vzvz_p - vzvz_m)/dzm(kc)

            enddo
        enddo
    enddo

    return

end subroutine ExplicitTermsVZ
  