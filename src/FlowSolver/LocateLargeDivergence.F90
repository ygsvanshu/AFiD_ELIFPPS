!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: LocateLargeDivergence.F90                      !
!    CONTAINS: subroutine LocateLargeDivergence           !
!                                                         !
!    PURPOSE: Debugging routine. Output the location(s)   !
!     of excessive divergence.                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine LocateLargeDivergence

    use param
    use local_arrays,only:vy,vx,vz
    use mpih
    use decomp_2d,only:xstart,xend,nrank

    implicit none

    integer :: ic,ip
    integer :: jc,jp
    integer :: kc,kp
    real    :: dqcap

    if (nrank .eq. 0) then
        write (6,"(A68)") repeat("-",68)
        write (6,"(A4,X,A,X,A6,X,A6,X,A6,X,A,X,A12,X,A12,X,A12)") "Rank","|","I","J","K","|","d(vx)/dx","d(vy)/dy","d(vz)/dz"
        write (6,"(A68)") repeat("-",68)
    end if

    do kc = xstart(3),xend(3)
        kp = kc+1
        do jc = xstart(2),xend(2)
            jp = jc+1
            do ic = 1,nxm
                ip = ic+1

                dqcap = ((vx(ip,jc,kc)-vx(ic,jc,kc))/dxc(ic) + (vy(ic,jp,kc)-vy(ic,jc,kc))/dyc(jc) + (vz(ic,jc,kp)-vz(ic,jc,kc))/dzc(kc))

                if (abs(dqcap) .gt. resid) then
                    write (6,"(I4,X,A,X,I6,X,I6,X,I6,X,A,X,E12.3,X,E12.3,X,E12.3)") nrank,'|',ic,jc,kc,'|',(vx(ip,jc,kc)-vx(ic,jc,kc))/dxc(ic),(vy(ic,jp,kc)-vy(ic,jc,kc))/dyc(jc),(vz(ic,jc,kp)-vz(ic,jc,kc))/dzc(kc)
                end if
                
            end do
        end do
    end do

    return

end subroutine LocateLargeDivergence
