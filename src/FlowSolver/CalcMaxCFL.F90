!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcMaxCFL.F90                                 !
!    CONTAINS: subroutine CalcMaxCFL                      !
!                                                         ! 
!    PURPOSE: Compute the maximum value of the local CFL  !
!     stability condition for the explicit terms          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcMaxCFL(cflm)

    use param
    use local_arrays, only: vx,vy,vz
    use decomp_2d
    use mpih

    implicit none

    real,intent(out)    :: cflm
    integer             :: kc,kp
    integer             :: jc,jp
    integer             :: ic,ip
    real                :: qcf
    real                :: res_dummy

    cflm = tiny(0.0)

    do kc=xstart(3),xend(3)
        kp=kc+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            do ic=xstart(1),xend(1)
                ip=ic+1
                
                qcf  = (abs(0.5*(vz(ic,jc,kc)+vz(ic,jc,kp))/dzc(kc)) + abs(0.5*(vy(ic,jc,kc)+vy(ic,jp,kc))/dyc(jc)) + abs(0.5*(vx(ic,jc,kc)+vx(ip,jc,kc))/dxc(ic)))
                cflm = max(cflm,qcf)

            enddo
        enddo
    enddo

    call MpiAllMaxRealScalar(cflm,res_dummy)
    cflm = res_dummy

    return

end subroutine CalcMaxCFL