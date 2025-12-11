!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdateVY.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdateVY             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for vx        !
!    and updates it to time t + (al*dt)                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdateVX

    use decomp_2d
    use param
    use local_arrays ,only : vx,rhx,rhy,rhz
    use boundary_arrays

    implicit none

    call TridiagonalSolverXC(1)
    call transpose_x_to_y(rhx,rhy)
    call TridiagonalSolverYM(1)
    if (nzm.gt.1) then
        call transpose_y_to_z(rhy,rhz)
        call TridiagonalSolverZM(1)
        call transpose_z_to_x(rhz,rhx)
    else
        call transpose_y_to_x(rhy,rhx)
    end if

    vx(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) = vx(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) + rhx(1:nxm,xstart(2):xend(2),xstart(3):xend(3))
    
    if (periodic(1)) then
        vx(nx,xstart(2):xend(2),xstart(3):xend(3)) = vx(1,xstart(2):xend(2),xstart(3):xend(3))
    end if

    return

end subroutine SolveImpEqnUpdateVX