!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: AddParticleForces.F90                          !
!    CONTAINS: subroutine AddParticleForces               !
!                                                         !
!    PURPOSE: Adds the computed particle forces as source !
!    terms to the explicit terms of momentum equations    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine AddParticleForces

    use param, only: dxm,dxc,dym,dyc,dzm,dzc
    use decomp_2d, only: xstart,xend
    use local_arrays, only: adx,ady,adz
    use lagrangian_point_particle, only: lpp_bdfx,lpp_bdfy,lpp_bdfz

    implicit none

    integer :: i,j,k

    do i = xstart(1),xend(1)
        do j = xstart(2),xend(2)
            do k = xstart(3),xend(3)
                adx(i,j,k) = adx(i,j,k) + lpp_bdfx(i,j,k)/(dxm(i)*dyc(j)*dzc(k))
                ady(i,j,k) = ady(i,j,k) + lpp_bdfy(i,j,k)/(dxc(i)*dym(j)*dzc(k))
                adz(i,j,k) = adz(i,j,k) + lpp_bdfz(i,j,k)/(dxc(i)*dyc(j)*dzm(k))
            end do
        end do
    end do

end subroutine AddParticleForces