!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CheckDivergence.F90                            !
!    CONTAINS: subroutine CheckDivergence                 !
!                                                         ! 
!    PURPOSE: Check the maximum divergence of velocity    !
!     in the domain                                       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckDivergence(qmax,qavg)

    use param
    use local_arrays, only: vx,vy,vz
    use mpih
    use decomp_2d, only: xstart,xend

    implicit none

    real,intent(out)    :: qmax,qavg
    integer             :: kc,kp
    integer             :: jc,jp
    integer             :: ic,ip
    real                :: dqcap
    real                :: res_dummy
        
    qmax = -huge(0.0)
    qavg = 0.0

    do kc=xstart(3),xend(3)
        kp=kc+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            do ic=xstart(1),xend(1)
                ip=ic+1

                dqcap = ((vx(ip,jc,kc)-vx(ic,jc,kc))/dxc(ic) + (vy(ic,jp,kc)-vy(ic,jc,kc))/dyc(jc) + (vz(ic,jc,kp)-vz(ic,jc,kc))/dzc(kc))
                qmax  = max(abs(dqcap),qmax)
                qavg  = qavg + (dqcap*dxc(ic)*dyc(jc)*dzc(kc))
            
            enddo
        enddo
    enddo

    call MpiMaxRealScalar(qmax,res_dummy)
    qmax = res_dummy

    call MpiSumRealScalar(qavg,res_dummy)
    qavg = res_dummy

    return     

end subroutine CheckDivergence