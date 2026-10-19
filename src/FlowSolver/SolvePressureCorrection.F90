!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SolvePressureCorrection.F90                    !
!    CONTAINS: subroutine SolvePressureCorrection         !
!              subroutine MultiplyEigenvectors(FW/BW)_X   !
!              subroutine MultiplyEigenvectors(FW/BW)_Y   !
!              subroutine MultiplyEigenvectors(FW/BW)_Z   !
!                                                         !
!    PURPOSE: Compute the pressure correction by solving  !
!    the pressure Poisson equation                        !
!                                                         !
!    The right-hand side is multiplied with the (back-    !
!    ward) eigenvectors of the second derivative in two   !
!    directions (BW: physical space -> eigenspace), a     !
!    tridiagonal system is solved in the third direction  !
!    (hdir), and the result is multiplied with the        !
!    (forward) eigenvectors (FW: eigenspace -> physical   !
!    space). Per direction, the multiplication is done    !
!    with dense eigenvector matrices (GEMM) or with real- !
!    to-real FFTs on uniform grids (tmet = TGEMM/TFFT):   !
!    the same operation, computed with a faster method.   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolvePressureCorrection

    use decomp_2d
    use param
    use local_arrays, only: rhx, rhy, rhz
    use pressure_solver

    implicit none

    if (hdir.eq.1) then
        ! Helmholtz solver in the X direction
        call transpose_x_to_y(rhx,rhy)
        call MultiplyEigenvectorsBW_Y
        call transpose_y_to_z(rhy,rhz)
        call MultiplyEigenvectorsBW_Z
        call transpose_z_to_x(rhz,rhx)
        call TridiagonalSolverXR
        call transpose_x_to_z(rhx,rhz)
        call MultiplyEigenvectorsFW_Z
        call transpose_z_to_y(rhz,rhy)
        call MultiplyEigenvectorsFW_Y
        call transpose_y_to_x(rhy,rhx)
    else if (hdir.eq.2) then
        ! Helmholtz solver in the Y direction
        call MultiplyEigenvectorsBW_X
        call transpose_x_to_z(rhx,rhz)
        call MultiplyEigenvectorsBW_Z
        call transpose_z_to_y(rhz,rhy)
        call TridiagonalSolverYR
        call transpose_y_to_z(rhy,rhz)
        call MultiplyEigenvectorsFW_Z
        call transpose_z_to_x(rhz,rhx)
        call MultiplyEigenvectorsFW_X
    else if (hdir.eq.3) then
        ! Helmholtz solver in the Z direction
        call MultiplyEigenvectorsBW_X
        call transpose_x_to_y(rhx,rhy)
        call MultiplyEigenvectorsBW_Y
        call transpose_y_to_z(rhy,rhz)
        call TridiagonalSolverZR
        call transpose_z_to_y(rhz,rhy)
        call MultiplyEigenvectorsFW_Y
        call transpose_y_to_x(rhy,rhx)
        call MultiplyEigenvectorsFW_X
    else
        if (ismaster) write(*,*) 'Error: Invalid helmholtz solver direction: ', hdir
        call MpiAbort
    end if

end subroutine SolvePressureCorrection

! X-pencil: rhx(1:nxm, xstart(2):xend(2), xstart(3):xend(3)) is viewed as
! a (nxm) x (xsize(2)*xsize(3)) matrix, so a single DGEMM transforms all lines.
! The FFT plan also transforms all lines at once.

subroutine MultiplyEigenvectorsFW_X

    use decomp_2d
    use param
    use local_arrays, only: rhx
    use pressure_solver

    implicit none

    if (tmet(1).eq.TFFT) then
        call fftw_execute_r2r(pfwx,rhx,rhx)
        rhx = rhx*(1.0/fnrm(1))
    else
        call dgemm('N','N',nxm,gszx,nxm,1.0,qfwx,nxm,rhx,nxm,0.0,gwrk,nxm)
        call dcopy(nxm*gszx,gwrk,1,rhx,1)
    end if

end subroutine MultiplyEigenvectorsFW_X

subroutine MultiplyEigenvectorsBW_X

    use decomp_2d
    use param
    use local_arrays, only: rhx
    use pressure_solver

    implicit none

    if (tmet(1).eq.TFFT) then
        call fftw_execute_r2r(pbwx,rhx,rhx)
    else
        call dgemm('N','N',nxm,gszx,nxm,1.0,qbwx,nxm,rhx,nxm,0.0,gwrk,nxm)
        call dcopy(nxm*gszx,gwrk,1,rhx,1)
    end if

end subroutine MultiplyEigenvectorsBW_X

! Y-pencil: Y is the middle index, so each kc plane rhy(:,1:nym,kc) is
! viewed as a (ysize(1)) x (nym) matrix and transformed with one DGEMM.
! The FFT plan transforms all lines of all kc planes at once.

subroutine MultiplyEigenvectorsFW_Y

    use decomp_2d, only: ystart, yend
    use param
    use local_arrays, only: rhy
    use pressure_solver

    implicit none

    integer :: kc

    if (tmet(2).eq.TFFT) then
        call fftw_execute_r2r(pfwy,rhy,rhy)
        rhy = rhy*(1.0/fnrm(2))
    else
        do kc = ystart(3), yend(3)
            call dgemm('N','T',gszy,nym,nym,1.0,rhy(ystart(1),1,kc),gszy,qfwy,nym,0.0,gwrk,gszy)
            call dcopy(gszy*nym,gwrk,1,rhy(ystart(1),1,kc),1)
        end do
    end if

end subroutine MultiplyEigenvectorsFW_Y

subroutine MultiplyEigenvectorsBW_Y

    use decomp_2d, only: ystart, yend
    use param
    use local_arrays, only: rhy
    use pressure_solver

    implicit none

    integer :: kc

    if (tmet(2).eq.TFFT) then
        call fftw_execute_r2r(pbwy,rhy,rhy)
    else
        do kc = ystart(3), yend(3)
            call dgemm('N','T',gszy,nym,nym,1.0,rhy(ystart(1),1,kc),gszy,qbwy,nym,0.0,gwrk,gszy)
            call dcopy(gszy*nym,gwrk,1,rhy(ystart(1),1,kc),1)
        end do
    end if

end subroutine MultiplyEigenvectorsBW_Y

! Z-pencil: rhz(zstart(1):zend(1), zstart(2):zend(2), 1:nzm) is viewed as
! a (zsize(1)*zsize(2)) x (nzm) matrix, so a single DGEMM transforms all lines.
! The FFT plan also transforms all lines at once.

subroutine MultiplyEigenvectorsFW_Z

    use decomp_2d
    use param
    use local_arrays, only: rhz
    use pressure_solver

    implicit none

    if (tmet(3).eq.TFFT) then
        call fftw_execute_r2r(pfwz,rhz,rhz)
        rhz = rhz*(1.0/fnrm(3))
    else
        call dgemm('N','T',gszz,nzm,nzm,1.0,rhz,gszz,qfwz,nzm,0.0,gwrk,gszz)
        call dcopy(gszz*nzm,gwrk,1,rhz,1)
    end if

end subroutine MultiplyEigenvectorsFW_Z

subroutine MultiplyEigenvectorsBW_Z

    use decomp_2d
    use param
    use local_arrays, only: rhz
    use pressure_solver

    implicit none

    if (tmet(3).eq.TFFT) then
        call fftw_execute_r2r(pbwz,rhz,rhz)
    else
        call dgemm('N','T',gszz,nzm,nzm,1.0,rhz,gszz,qbwz,nzm,0.0,gwrk,gszz)
        call dcopy(gszz*nzm,gwrk,1,rhz,1)
    end if

end subroutine MultiplyEigenvectorsBW_Z
