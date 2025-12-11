!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: SolveImpEqnUpdateVY.F90                        !
!    CONTAINS: subroutine SolveImpEqnUpdateVY             !
!                                                         ! 
!    PURPOSE: Inverts the implicit equation for vz        !
!    and updates it to time t + (al*dt)                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolveImpEqnUpdateVZ

    use decomp_2d
    use param
    use local_arrays ,only : vz,rhx,rhy,rhz
    use boundary_arrays

    implicit none

    call TridiagonalSolverXM(3)
    call transpose_x_to_y(rhx,rhy)
    call TridiagonalSolverYM(3)
    if (nzm.gt.1) then
        call transpose_y_to_z(rhy,rhz)
        call TridiagonalSolverZC(3)
        call transpose_z_to_x(rhz,rhx)
    else
        call transpose_y_to_x(rhy,rhx)
    end if

    vz(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) = vz(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) + rhx(1:nxm,xstart(2):xend(2),xstart(3):xend(3))

    if (periodic(1)) then
        vz(0 ,xstart(2):xend(2),xstart(3):xend(3)) = vz(nxm,xstart(2):xend(2),xstart(3):xend(3))
        vz(nx,xstart(2):xend(2),xstart(3):xend(3)) = vz(1  ,xstart(2):xend(2),xstart(3):xend(3))
    end if

    return

end subroutine SolveImpEqnUpdateVZ