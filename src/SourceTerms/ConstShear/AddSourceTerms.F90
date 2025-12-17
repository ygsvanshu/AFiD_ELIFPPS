!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: AddSourceTerms.F90                             !
!    CONTAINS: subroutine AddSourceTerms                  !
!                                                         !
!    PURPOSE: Adds source terms to advection-diffusion    !
!             equations.                                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine AddSourceTerms

    use decomp_2d, only: xstart,xend
    use param
    use local_arrays, only: vy,ady

    implicit none

    integer :: i,j,k
    real    :: vbar,res_dummy

    ! Compute the global friction velocity
    !! Calculate contribution of current pencil/process
    vbar = 0.0
    do j = xstart(2),xend(2)
        do k = xstart(3),xend(3)
            vbar = vbar + 0.25*((vy(1  ,j  ,k) - vy(0 ,j  ,k))/dxm(1))*dyc(j)*dzc(k)
            vbar = vbar + 0.25*((vy(1  ,j+1,k) - vy(0 ,j+1,k))/dxm(1))*dyc(j)*dzc(k)
            vbar = vbar + 0.25*((vy(nxm,j  ,k) - vy(nx,j  ,k))/dxm(1))*dyc(j)*dzc(k)
            vbar = vbar + 0.25*((vy(nxm,j+1,k) - vy(nx,j+1,k))/dxm(1))*dyc(j)*dzc(k)
        end do
    end do
    !! Get sum on full domain
    call MpiAllSumRealScalar(vbar,res_dummy)
    vbar = res_dummy/(ylen*zlen)
    vbar = sqrt(vbar/rey)

    ! Apply constant source term (acceleration) to remove the deficit/excess
    ! This will enforce the friction Reynolds number (Re = 0.5*U_mean*H/nu) to be statistically constant 
    do i = xstart(1),xend(1)
        do j = xstart(2),xend(2)
            do k = xstart(3),xend(3)
                ady(i,j,k) = ady(i,j,k) + (0.5 - vbar)/(al*dt)
            end do
        end do
    end do

    return

end subroutine AddSourceTerms  