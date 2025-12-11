!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SolvePressureCorrection.F90                    !
!    CONTAINS: subroutine SolvePressureCorrection         !
!                                                         !
!    PURPOSE: Compute the pressure correction by solving  !
!    the pressure Poisson equation                        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SolvePressureCorrection

    use decomp_2d
    use param
    use fftw_params
    use local_arrays, only: dph
    use pressure_decomp

    implicit none

    if (allocated(rx1)) rx1(:,:,:) = 0.0
    if (allocated(cx1)) cx1(:,:,:) = 0.0
    if (allocated(ry1)) ry1(:,:,:) = 0.0
    if (allocated(cy1)) cy1(:,:,:) = 0.0
    if (allocated(rz1)) rz1(:,:,:) = 0.0
    if (allocated(cz1)) cz1(:,:,:) = 0.0

    if (straxs.eq.1) then
        ! Stretching in X axis
        ! X--->Y--->Y--->Z--->Z--->Y--->X--->X--->Y--->Z--->Z--->Y--->Y--->X
        !   T    F    T    F    T    T    S    T    T    B    T    B    T
        call transpose_x_to_y(dph, ry1, ph_info_n)
        if (periodic(2)) then
            call dfftw_execute_dft_r2c(fwd_guruplan_y, ry1, cy1)
            call transpose_y_to_z(cy1, cz1, ph_info_y)
            if (periodic(3)) then
                call dfftw_execute_dft(fwd_guruplan_z, cz1, cz1)
                call transpose_z_to_x(cz1, cx1, ph_info_y)
                call TridiagonalSolverXP(nym*nzm, ph_info_y)
                call transpose_x_to_z(cx1, cz1, ph_info_y)
                call dfftw_execute_dft(bwd_guruplan_z, cz1, cz1)
            else
                call dfftw_execute_r2r(fwd_guruplan_z, cz1%re, cz1%re)
                call dfftw_execute_r2r(fwd_guruplan_z, cz1%im, cz1%im)
                call transpose_z_to_x(cz1, cx1, ph_info_y)
                call TridiagonalSolverXP(2*nym*nzm, ph_info_y)
                call transpose_x_to_z(cx1, cz1, ph_info_y)
                call dfftw_execute_r2r(bwd_guruplan_z, cz1%re, cz1%re)
                call dfftw_execute_r2r(bwd_guruplan_z, cz1%im, cz1%im)
            end if
            call transpose_z_to_y(cz1, cy1, ph_info_y)
            call dfftw_execute_dft_c2r(bwd_guruplan_y, cy1, ry1)
        else
            call dfftw_execute_r2r(fwd_guruplan_y, ry1, ry1)
            call transpose_y_to_z(ry1, rz1, ph_info_n)
            if (periodic(3)) then
                call dfftw_execute_dft_r2c(fwd_guruplan_z, rz1, cz1)
                call transpose_z_to_x(cz1, cx1, ph_info_z)
                call TridiagonalSolverXP(2*nym*nzm, ph_info_z)
                call transpose_x_to_z(cx1, cz1, ph_info_z)
                call dfftw_execute_dft_c2r(bwd_guruplan_z, cz1, rz1)
            else
                call dfftw_execute_r2r(fwd_guruplan_z, rz1, rz1)
                call transpose_z_to_x(rz1, rx1, ph_info_n)
                call TridiagonalSolverXR(4*nym*nzm, ph_info_n)
                call transpose_x_to_z(rx1, rz1, ph_info_n)
                call dfftw_execute_r2r(bwd_guruplan_z, rz1, rz1)
            end if
            call transpose_z_to_y(rz1, ry1, ph_info_n)
            call dfftw_execute_r2r(bwd_guruplan_y, ry1, ry1)
        end if
        call transpose_y_to_x(ry1, dph, ph_info_n)
    else if (straxs.eq.2) then
        ! Stretching in Y axis
        ! X--->X--->Y--->Z--->Z--->Y--->Y--->Z--->Z--->Y--->X--->X
        !   F    T    T    F    T    S    T    B    T    T    B   
        if (periodic(1)) then
            call dfftw_execute_dft_r2c(fwd_guruplan_x, dph, cx1)
            call transpose_x_to_z(cx1, cz1, ph_info_x)
            if (periodic(3)) then
                call dfftw_execute_dft(fwd_guruplan_z, cz1, cz1)
                call transpose_z_to_y(cz1, cy1, ph_info_x)
                call TridiagonalSolverYP(nxm*nzm, ph_info_x)
                call transpose_y_to_z(cy1, cz1, ph_info_x)
                call dfftw_execute_dft(bwd_guruplan_z, cz1, cz1)
            else
                call dfftw_execute_r2r(fwd_guruplan_z, cz1%re, cz1%re)
                call dfftw_execute_r2r(fwd_guruplan_z, cz1%im, cz1%im)
                call transpose_z_to_y(cz1, cy1, ph_info_x)
                call TridiagonalSolverYP(2*nxm*nzm, ph_info_x)
                call transpose_y_to_z(cy1, cz1, ph_info_x)
                call dfftw_execute_r2r(bwd_guruplan_z, cz1%re, cz1%re)
                call dfftw_execute_r2r(bwd_guruplan_z, cz1%im, cz1%im)
            end if
            call transpose_z_to_x(cz1, cx1, ph_info_x)
            call dfftw_execute_dft_c2r(bwd_guruplan_x, cx1, dph)
        else
            call dfftw_execute_r2r(fwd_guruplan_x, dph, rx1)
            call transpose_x_to_z(rx1, rz1, ph_info_n)
            if (periodic(3)) then
                call dfftw_execute_dft_r2c(fwd_guruplan_z, rz1, cz1)
                call transpose_z_to_y(cz1, cy1, ph_info_z)
                call TridiagonalSolverYP(2*nxm*nzm, ph_info_z)
                call transpose_y_to_z(cy1, cz1, ph_info_z)
                call dfftw_execute_dft_c2r(bwd_guruplan_z, cz1, rz1)
            else
                call dfftw_execute_r2r(fwd_guruplan_z, rz1, rz1)
                call transpose_z_to_y(rz1, ry1, ph_info_n)
                call TridiagonalSolverYR(4*nxm*nzm, ph_info_n)
                call transpose_y_to_z(ry1, rz1, ph_info_n)
                call dfftw_execute_r2r(bwd_guruplan_z, rz1, rz1)
            end if
            call transpose_z_to_x(rz1, rx1, ph_info_n)
            call dfftw_execute_r2r(bwd_guruplan_x, rx1, dph)
        end if
    else if (straxs.eq.3) then
        ! Stretching in Z axis
        ! X--->X--->Y--->Y--->Z--->Z--->Y--->Y--->X--->X
        !   F    T    F    T    S    T    B    T    B  
        if (periodic(1)) then
            call dfftw_execute_dft_r2c(fwd_guruplan_x, dph, cx1)
            call transpose_x_to_y(cx1, cy1, ph_info_x)
            if (periodic(2)) then
                call dfftw_execute_dft(fwd_guruplan_y, cy1, cy1)
                call transpose_y_to_z(cy1, cz1, ph_info_x)
                call TridiagonalSolverZP(nxm*nym, ph_info_x)
                call transpose_z_to_y(cz1, cy1, ph_info_x)
                call dfftw_execute_dft(bwd_guruplan_y, cy1, cy1)
            else
                call dfftw_execute_r2r(fwd_guruplan_y, cy1%re, cy1%re)
                call dfftw_execute_r2r(fwd_guruplan_y, cy1%im, cy1%im)
                call transpose_y_to_z(cy1, cz1, ph_info_x)
                call TridiagonalSolverZP(2*nxm*nym, ph_info_x)
                call transpose_z_to_y(cz1, cy1, ph_info_x)
                call dfftw_execute_r2r(bwd_guruplan_y, cy1%re, cy1%re)
                call dfftw_execute_r2r(bwd_guruplan_y, cy1%im, cy1%im)
            end if
            call transpose_y_to_x(cy1, cx1, ph_info_x)
            call dfftw_execute_dft_c2r(bwd_guruplan_x, cx1, dph)
        else
            call dfftw_execute_r2r(fwd_guruplan_x, dph, rx1)
            call transpose_x_to_y(rx1, ry1, ph_info_n)
            if (periodic(2)) then
                call dfftw_execute_dft_r2c(fwd_guruplan_y, ry1, cy1)
                call transpose_y_to_z(cy1, cz1, ph_info_y)
                call TridiagonalSolverZP(2*nxm*nym, ph_info_y)
                call transpose_z_to_y(cz1, cy1, ph_info_y)
                call dfftw_execute_dft_c2r(bwd_guruplan_y, cy1, ry1)
            else
                call dfftw_execute_r2r(fwd_guruplan_y, ry1, ry1)
                call transpose_y_to_z(ry1, rz1, ph_info_n)
                call TridiagonalSolverZR(4*nxm*nym, ph_info_n)
                call transpose_z_to_y(rz1, ry1, ph_info_n)
                call dfftw_execute_r2r(bwd_guruplan_y, ry1, ry1)
            end if
            call transpose_y_to_x(ry1, rx1, ph_info_n)
            call dfftw_execute_r2r(bwd_guruplan_x, rx1, dph)
        end if
    else
        ! Uniform grid
        ! X--->X--->Y--->Y--->Z--->Z--->Y--->Y--->X--->X
        !   F    T    F    T    S    T    B    T    B
        if (periodic(1)) then
            call dfftw_execute_dft_r2c(fwd_guruplan_x, dph, cx1)
            call transpose_x_to_y(cx1, cy1, ph_info_x)
            if (periodic(2)) then
                call dfftw_execute_dft(fwd_guruplan_y, cy1, cy1)
                call transpose_y_to_z(cy1, cz1, ph_info_x)
                if (periodic(3)) then
                    call dfftw_execute_dft(fwd_guruplan_z, cz1, cz1)
                    call SpectralSolverZP(nxm*nym*nzm, ph_info_x)
                    call dfftw_execute_dft(bwd_guruplan_z, cz1, cz1)
                else
                    call dfftw_execute_r2r(fwd_guruplan_z, cz1%re, cz1%re)
                    call dfftw_execute_r2r(fwd_guruplan_z, cz1%im, cz1%im)
                    call SpectralSolverZP(2*nxm*nym*nzm, ph_info_x)
                    call dfftw_execute_r2r(bwd_guruplan_z, cz1%re, cz1%re)
                    call dfftw_execute_r2r(bwd_guruplan_z, cz1%im, cz1%im)
                end if
                call transpose_z_to_y(cz1, cy1, ph_info_x)
                call dfftw_execute_dft(bwd_guruplan_y, cy1, cy1)
            else
                call dfftw_execute_r2r(fwd_guruplan_y, cy1%re, cy1%re)
                call dfftw_execute_r2r(fwd_guruplan_y, cy1%im, cy1%im)
                call transpose_y_to_z(cy1, cz1, ph_info_x)
                if (periodic(3)) then
                    call dfftw_execute_dft(fwd_guruplan_z, cz1, cz1)
                    call SpectralSolverZP(2*nxm*nym*nzm, ph_info_x)
                    call dfftw_execute_dft(bwd_guruplan_z, cz1, cz1)
                else
                    call dfftw_execute_r2r(fwd_guruplan_z, cz1%re, cz1%re)
                    call dfftw_execute_r2r(fwd_guruplan_z, cz1%im, cz1%im)
                    call SpectralSolverZP(4*nxm*nym*nzm, ph_info_x)
                    call dfftw_execute_r2r(bwd_guruplan_z, cz1%re, cz1%re)
                    call dfftw_execute_r2r(bwd_guruplan_z, cz1%im, cz1%im)
                end if
                call transpose_z_to_y(cz1, cy1, ph_info_x)
                call dfftw_execute_r2r(bwd_guruplan_y, cy1%re, cy1%re)
                call dfftw_execute_r2r(bwd_guruplan_y, cy1%im, cy1%im)
            end if
            call transpose_y_to_x(cy1, cx1, ph_info_x)
            call dfftw_execute_dft_c2r(bwd_guruplan_x, cx1, dph)
        else
            call dfftw_execute_r2r(fwd_guruplan_x, dph, rx1)
            call transpose_x_to_y(rx1, ry1, ph_info_n)
            if (periodic(2)) then
                call dfftw_execute_dft_r2c(fwd_guruplan_y, ry1, cy1)
                call transpose_y_to_z(cy1, cz1, ph_info_y)
                if (periodic(3)) then
                    call dfftw_execute_dft(fwd_guruplan_z, cz1, cz1)
                    call SpectralSolverZP(2*nxm*nym*nzm, ph_info_y)
                    call dfftw_execute_dft(bwd_guruplan_z, cz1, cz1)
                else
                    call dfftw_execute_r2r(fwd_guruplan_z, cz1%re, cz1%re)
                    call dfftw_execute_r2r(fwd_guruplan_z, cz1%im, cz1%im)
                    call SpectralSolverZP(4*nxm*nym*nzm, ph_info_y)
                    call dfftw_execute_r2r(bwd_guruplan_z, cz1%re, cz1%re)
                    call dfftw_execute_r2r(bwd_guruplan_z, cz1%im, cz1%im)
                end if
                call transpose_z_to_y(cz1, cy1, ph_info_y)
                call dfftw_execute_dft_c2r(bwd_guruplan_y, cy1, ry1)
            else
                call dfftw_execute_r2r(fwd_guruplan_y, ry1, ry1)
                call transpose_y_to_z(ry1, rz1, ph_info_n)
                if (periodic(3)) then
                    call dfftw_execute_dft_r2c(fwd_guruplan_z, rz1, cz1)
                    call SpectralSolverZP(4*nxm*nym*nzm, ph_info_z)
                    call dfftw_execute_dft_c2r(bwd_guruplan_z, cz1, rz1)
                else
                    call dfftw_execute_r2r(fwd_guruplan_z, rz1, rz1)
                    call SpectralSolverZR(8*nxm*nym*nzm, ph_info_n)
                    call dfftw_execute_r2r(bwd_guruplan_z, rz1, rz1)
                end if
                call transpose_z_to_y(rz1, ry1, ph_info_n)
                call dfftw_execute_r2r(bwd_guruplan_y, ry1, ry1)
            end if
            call transpose_y_to_x(ry1, rx1, ph_info_n)
            call dfftw_execute_r2r(bwd_guruplan_x, rx1, dph)
        end if
    end if

end subroutine SolvePressureCorrection

subroutine SpectralSolverZP(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : cz1

    implicit none

    integer, intent(in)             :: norm
    type(decomp_info), intent(in)   :: dcmp

    integer                         :: ic,jc,kc

    do ic = dcmp%zst(1),dcmp%zen(1)
        do jc = dcmp%zst(2),dcmp%zen(2)
            do kc = dcmp%zst(3),dcmp%zen(3)
                if ((ic.eq.1).and.(jc.eq.1).and.(kc.eq.1)) then
                    cz1(ic,jc,kc) = 0.0
                else
                    cz1(ic,jc,kc) = -cz1(ic,jc,kc)/(real(norm)*(ak1(ic) + ak2(jc) + ak3(kc)))
                end if
            end do
        end do
    end do

end subroutine SpectralSolverZP

subroutine SpectralSolverZR(norm,dcmp)

    use param
    use decomp_2d
    use fftw_params, only : rz1

    implicit none

    integer, intent(in)             :: norm
    type(decomp_info), intent(in)   :: dcmp

    integer                         :: ic,jc,kc

    do ic = dcmp%zst(1),dcmp%zen(1)
        do jc = dcmp%zst(2),dcmp%zen(2)
            do kc = dcmp%zst(3),dcmp%zen(3)
                if ((ic.eq.1).and.(jc.eq.1).and.(kc.eq.1)) then
                    rz1(ic,jc,kc) = 0.0
                else
                    rz1(ic,jc,kc) = -rz1(ic,jc,kc)/(real(norm)*(ak1(ic) + ak2(jc) + ak3(kc)))
                end if
            end do
        end do
    end do

end subroutine SpectralSolverZR