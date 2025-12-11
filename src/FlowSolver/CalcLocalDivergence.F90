!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CalcLocalDivergence.F90                        !
!    CONTAINS: subroutine CalcLocalDivergence             !
!                                                         ! 
!    PURPOSE: Compute the divergence of the intermediate  !
!     velocity at every point for the pressure            !
!     correction step                                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcLocalDivergence

    use param
    use local_arrays, only: vx,vy,vz,dph
    use decomp_2d, only: xstart,xend

    implicit none

    integer :: kc,kp
    integer :: jc,jp
    integer :: ic,ip 

    do kc=xstart(3),xend(3)
        kp=kc+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            do ic=xstart(1),xend(1)
                ip=ic+1

                dph(ic,jc,kc) = ((vx(ip,jc,kc)-vx(ic,jc,kc))/dxc(ic) + (vy(ic,jp,kc)-vy(ic,jc,kc))/dyc(jc) + (vz(ic,jc,kp)-vz(ic,jc,kc))/dzc(kc))/(al*dt)

            enddo
        enddo
    enddo

    return
      
end subroutine CalcLocalDivergence
