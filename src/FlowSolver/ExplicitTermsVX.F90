!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: ExplicitTermsVX.F90                            !
!    CONTAINS: subroutine ExplicitTermsVX                 !
!                                                         ! 
!    PURPOSE: Compute the non-linear terms associated to  !
!     the velocity in the x (vertical) dimension          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExplicitTermsVX

    use param
    use local_arrays, only: vx,vy,vz,adx
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: im,ic,ip
    integer :: jm,jc,jp
    integer :: km,kc,kp

    real    :: vxvx_m,vxvx_p
    real    :: vxvy_m,vxvy_p
    real    :: vxvz_m,vxvz_p

    do kc=xstart(3),xend(3)
        km=kc-1
        kp=kc+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do ic=xstart(1),xend(1)
                im=ic-1
                ip=ic+1

                vxvx_m = ((0.5*vx(im,jc,kc)) + (0.5*vx(ic,jc,kc))) * ((0.5*vx(im,jc,kc)) + (0.5*vx(ic,jc,kc)))
                vxvx_p = ((0.5*vx(ic,jc,kc)) + (0.5*vx(ip,jc,kc))) * ((0.5*vx(ic,jc,kc)) + (0.5*vx(ip,jc,kc)))
                
                vxvy_m = ((ixcm(ic)*vy(im,jc,kc)) + (ixcp(ic)*vy(ic,jc,kc))) * ((iycm(jc)*vx(ic,jm,kc)) + (iycp(jc)*vx(ic,jc,kc)))
                vxvy_p = ((ixcm(ic)*vy(im,jp,kc)) + (ixcp(ic)*vy(ic,jp,kc))) * ((iycm(jp)*vx(ic,jc,kc)) + (iycp(jp)*vx(ic,jp,kc)))
                
                vxvz_m = ((ixcm(ic)*vz(im,jc,kc)) + (ixcp(ic)*vz(ic,jc,kc))) * ((izcm(kc)*vx(ic,jc,km)) + (izcp(kc)*vx(ic,jc,kc))) 
                vxvz_p = ((ixcm(ic)*vz(im,jc,kp)) + (ixcp(ic)*vz(ic,jc,kp))) * ((izcm(kp)*vx(ic,jc,kc)) + (izcp(kp)*vx(ic,jc,kp))) 

                adx(ic,jc,kc) = 0.0
                adx(ic,jc,kc) = adx(ic,jc,kc) - (vxvx_p - vxvx_m)/dxm(ic)
                adx(ic,jc,kc) = adx(ic,jc,kc) - (vxvy_p - vxvy_m)/dyc(jc)
                adx(ic,jc,kc) = adx(ic,jc,kc) - (vxvz_p - vxvz_m)/dzc(kc)

            enddo
        enddo
    enddo

    return
    
end subroutine ExplicitTermsVX