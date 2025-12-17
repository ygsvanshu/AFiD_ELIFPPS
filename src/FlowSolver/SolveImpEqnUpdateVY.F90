!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdateVY.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdateVY             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for vy        !
!    and updates it to time t + (al*dt)                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdateVY

    use decomp_2d
    use param
    use local_arrays ,only : vy,rhx,rhy,rhz
    use boundary_arrays

    implicit none

    call TridiagonalSolverXM(2)
    call transpose_x_to_y(rhx,rhy)
    call TridiagonalSolverYC(2)
    if (nzm.gt.1) then
        call transpose_y_to_z(rhy,rhz)
        call TridiagonalSolverZM(2)
        call transpose_z_to_x(rhz,rhx)
    else
        call transpose_y_to_x(rhy,rhx)
    end if

    vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) = vy(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) + rhx(1:nxm,xstart(2):xend(2),xstart(3):xend(3))

    if (periodic(1)) then
        vy(0 ,xstart(2):xend(2),xstart(3):xend(3)) = vy(nxm,xstart(2):xend(2),xstart(3):xend(3))
        vy(nx,xstart(2):xend(2),xstart(3):xend(3)) = vy(1  ,xstart(2):xend(2),xstart(3):xend(3))
    end if

    return

end subroutine SolveImpEqnUpdateVY