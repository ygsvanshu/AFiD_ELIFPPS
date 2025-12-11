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

    ! Compute the global flux
    !! Calculate contribution of current pencil/process
    vbar = 0.0
    do i = xstart(1),xend(1)
        do j = xstart(2),xend(2)
            do k = xstart(3),xend(3)
                vbar = vbar + 0.5*(vy(i,j,k) + vy(i,j+1,k))*dxc(i)*dyc(j)*dzc(k)
            end do
        end do
    end do
    !! Get sum on full domain
    call MpiAllSumRealScalar(vbar,res_dummy)
    vbar = res_dummy/(xlen*ylen*zlen)

    ! Apply constant source term (acceleration) to remove the deficit/excess
    ! This will enforce the bulk Reynolds number (Re = U_mean*H/nu) to be statistically constant 
    do i = xstart(1),xend(1)
        do j = xstart(2),xend(2)
            do k = xstart(3),xend(3)
                ady(i,j,k) = ady(i,j,k) + (1.0 - vbar)/(al*dt)
            end do
        end do
    end do

    return

end subroutine AddSourceTerms  