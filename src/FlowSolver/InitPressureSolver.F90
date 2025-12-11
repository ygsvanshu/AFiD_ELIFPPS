!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: InitPressureSolver.F90                         !
!    CONTAINS: subroutine InitPressureSolver              !
!                                                         !
!    PURPOSE: Initialization routines. Compute the metric !
!     terms and modified wavenumbers for the pressure     !
!     correction                                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitPressureSolver

    use, intrinsic :: iso_c_binding
    use mpih
    use decomp_2d
    use param
    use local_arrays, only: dph
    use fftw_params
    use pressure_decomp
    use AuxiliaryRoutines, only: AllocateRealFFTArray, AllocateComplexFFTArray

    implicit none

    type(fftw_iodim), dimension(1)  :: iodim
    type(fftw_iodim), dimension(2)  :: iodim_howmany

    logical                         :: planned_x
    logical                         :: planned_y
    logical                         :: planned_z

    ! Create sizes
    ! Create arrays
    ! Plan transforms

    call decomp_info_init(nxm, nym, nzm, ph_info_n)
    call decomp_info_init(nxm/2+1, nym, nzm, ph_info_x)
    call decomp_info_init(nxm, nym/2+1, nzm, ph_info_y)
    call decomp_info_init(nxm, nym, nzm/2+1, ph_info_z)

    if (straxs.eq.1) then
        ! Stretching in X axis
        ! X--->Y--->Y--->Z--->Z--->Y--->X--->X--->Y--->Z--->Z--->Y--->Y--->X
        !   T    F    T    F    T    T    S    T    T    B    T    B    T
        call AllocateRealFFTArray(ry1, ph_info_n, 'y')
        if (periodic(2)) then
            call AllocateComplexFFTArray(cy1, ph_info_y, 'y')
            call SetFFTiodim(ph_info_n, ph_info_y, 'y', iodim, iodim_howmany)
            fwd_guruplan_y = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, ry1, cy1, FFTW_ESTIMATE)
            call SetFFTiodim(ph_info_y, ph_info_n, 'y', iodim, iodim_howmany)
            bwd_guruplan_y = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cy1, ry1, FFTW_ESTIMATE)
            call AllocateComplexFFTArray(cz1, ph_info_y, 'z')
            call AllocateComplexFFTArray(cx1, ph_info_y, 'x')
            call SetFFTiodim(ph_info_y, ph_info_y, 'z', iodim, iodim_howmany)
            if (periodic(3)) then
                fwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_FORWARD, FFTW_ESTIMATE)
                bwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_BACKWARD, FFTW_ESTIMATE)
            else
                fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            end if
        else
            call SetFFTiodim(ph_info_n, ph_info_n, 'y', iodim, iodim_howmany)
            fwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, ry1, ry1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
            bwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, ry1, ry1, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            call AllocateRealFFTArray(rz1, ph_info_n, 'z')
            if (periodic(3)) then
                call AllocateComplexFFTArray(cz1, ph_info_z, 'z')
                call AllocateComplexFFTArray(cx1, ph_info_z, 'x')
                call SetFFTiodim(ph_info_n, ph_info_z, 'z', iodim, iodim_howmany)
                fwd_guruplan_z = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, rz1, cz1, FFTW_ESTIMATE)
                call SetFFTiodim(ph_info_z, ph_info_n, 'z', iodim, iodim_howmany)
                bwd_guruplan_z = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cz1, rz1, FFTW_ESTIMATE)
            else
                call AllocateRealFFTArray(rx1, ph_info_n, 'x')
                call SetFFTiodim(ph_info_n, ph_info_n, 'z', iodim, iodim_howmany)
                fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rz1, rz1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rz1, rz1, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            end if
        end if
        ! Check if plans are created properly
        planned_y = (c_associated(fwd_guruplan_y)).and.(c_associated(bwd_guruplan_y))
        planned_z = (c_associated(fwd_guruplan_z)).and.(c_associated(bwd_guruplan_z))
        if (planned_y.and.planned_z) planned = .true.
    else if (straxs.eq.2) then
        ! Stretching in Y axis
        ! X--->X--->Y--->Z--->Z--->Y--->Y--->Z--->Z--->Y--->X--->X
        !   F    T    T    F    T    S    T    B    T    T    B   
        if (periodic(1)) then
            call AllocateComplexFFTArray(cx1, ph_info_x, 'x')
            call SetFFTiodim(ph_info_n, ph_info_x, 'x', iodim, iodim_howmany)
            fwd_guruplan_x = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, dph, cx1, FFTW_ESTIMATE)
            call SetFFTiodim(ph_info_x, ph_info_n, 'x', iodim, iodim_howmany)
            bwd_guruplan_x = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cx1, dph, FFTW_ESTIMATE)
            call AllocateComplexFFTArray(cy1, ph_info_x, 'y')
            call AllocateComplexFFTArray(cz1, ph_info_x, 'z')
            call SetFFTiodim(ph_info_x, ph_info_x, 'z', iodim, iodim_howmany)
            if (periodic(3)) then
                fwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_FORWARD, FFTW_ESTIMATE)
                bwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_BACKWARD, FFTW_ESTIMATE)
            else
                fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            end if
        else
            call AllocateRealFFTArray(rx1, ph_info_n, 'x')
            call SetFFTiodim(ph_info_n, ph_info_n, 'x', iodim, iodim_howmany)
            fwd_guruplan_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, dph, rx1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
            bwd_guruplan_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rx1, dph, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            call AllocateRealFFTArray(rz1, ph_info_n, 'z')
            if (periodic(3)) then
                call AllocateComplexFFTArray(cz1, ph_info_z, 'z')
                call SetFFTiodim(ph_info_n, ph_info_z, 'z', iodim, iodim_howmany)
                fwd_guruplan_z = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, rz1, cz1, FFTW_ESTIMATE)
                call SetFFTiodim(ph_info_z, ph_info_n, 'z', iodim, iodim_howmany)
                bwd_guruplan_z = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cz1, rz1, FFTW_ESTIMATE)
                call AllocateComplexFFTArray(cy1, ph_info_z, 'y')
            else
                call SetFFTiodim(ph_info_n, ph_info_n, 'z', iodim, iodim_howmany)
                fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rz1, rz1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rz1, rz1, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                call AllocateRealFFTArray(ry1, ph_info_n, 'y')
            end if
        end if
        ! Check if plans are created properly 
        planned_x = (c_associated(fwd_guruplan_x)).and.(c_associated(bwd_guruplan_x))
        planned_z = (c_associated(fwd_guruplan_z)).and.(c_associated(bwd_guruplan_z))
        if (planned_x.and.planned_z) planned = .true.
    else if (straxs.eq.3) then
        ! Stretching in Z axis
        ! X--->X--->Y--->Y--->Z--->Z--->Y--->Y--->X--->X
        !   F    T    F    T    S    T    B    T    B  
        if (periodic(1)) then
            call AllocateComplexFFTArray(cx1, ph_info_x, 'x')
            call SetFFTiodim(ph_info_n, ph_info_x, 'x', iodim, iodim_howmany)
            fwd_guruplan_x = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, dph, cx1, FFTW_ESTIMATE)
            call SetFFTiodim(ph_info_x, ph_info_n, 'x', iodim, iodim_howmany)
            bwd_guruplan_x = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cx1, dph, FFTW_ESTIMATE)
            call AllocateComplexFFTArray(cy1, ph_info_x, 'y')
            call SetFFTiodim(ph_info_x, ph_info_x, 'y', iodim, iodim_howmany)
            if (periodic(2)) then
                fwd_guruplan_y = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cy1, cy1, FFTW_FORWARD, FFTW_ESTIMATE)
                bwd_guruplan_y = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cy1, cy1, FFTW_BACKWARD, FFTW_ESTIMATE)
            else
                fwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cy1%re, cy1%re, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cy1%re, cy1%re, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            end if
            call AllocateComplexFFTArray(cz1, ph_info_x, 'z')
        else
            call AllocateRealFFTArray(rx1, ph_info_n, 'x')
            call SetFFTiodim(ph_info_n, ph_info_n, 'x', iodim, iodim_howmany)
            fwd_guruplan_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, dph, rx1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
            bwd_guruplan_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rx1, dph, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            call AllocateRealFFTArray(ry1, ph_info_n, 'y')
            if (periodic(2)) then
                call AllocateComplexFFTArray(cy1, ph_info_y, 'y')
                call SetFFTiodim(ph_info_n, ph_info_y, 'y', iodim, iodim_howmany)
                fwd_guruplan_y = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, ry1, cy1, FFTW_ESTIMATE)
                call SetFFTiodim(ph_info_y, ph_info_n, 'y', iodim, iodim_howmany)
                bwd_guruplan_y = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cy1, ry1, FFTW_ESTIMATE)
                call AllocateComplexFFTArray(cz1, ph_info_y, 'z')
            else
                call SetFFTiodim(ph_info_n, ph_info_n, 'y', iodim, iodim_howmany)
                fwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, ry1, ry1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, ry1, ry1, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                call AllocateRealFFTArray(rz1, ph_info_n, 'z')
            end if
        end if
        ! Check if plans are created properly 
        planned_x = (c_associated(fwd_guruplan_x)).and.(c_associated(bwd_guruplan_x))
        planned_y = (c_associated(fwd_guruplan_y)).and.(c_associated(bwd_guruplan_y))
        if (planned_x.and.planned_y) planned = .true.
    else
        ! Uniform grid
        ! X--->X--->Y--->Y--->Z--->Z--->Y--->Y--->X--->X
        !   F    T    F    T    S    T    B    T    B
        if (periodic(1)) then
            call AllocateComplexFFTArray(cx1, ph_info_x, 'x')
            call AllocateComplexFFTArray(cy1, ph_info_x, 'y')
            call AllocateComplexFFTArray(cz1, ph_info_x, 'z')
            call SetFFTiodim(ph_info_n, ph_info_x, 'x', iodim, iodim_howmany)
            fwd_guruplan_x = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, dph, cx1, FFTW_ESTIMATE)
            call SetFFTiodim(ph_info_x, ph_info_n, 'x', iodim, iodim_howmany)
            bwd_guruplan_x = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cx1, dph, FFTW_ESTIMATE)
            if (periodic(2)) then
                call SetFFTiodim(ph_info_x, ph_info_x, 'y', iodim, iodim_howmany)
                fwd_guruplan_y = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cy1, cy1, FFTW_FORWARD, FFTW_ESTIMATE)
                bwd_guruplan_y = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cy1, cy1, FFTW_BACKWARD, FFTW_ESTIMATE)
                call SetFFTiodim(ph_info_x, ph_info_x, 'z', iodim, iodim_howmany)
                if (periodic(3)) then
                    fwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_FORWARD, FFTW_ESTIMATE)
                    bwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_BACKWARD, FFTW_ESTIMATE)
                else
                    fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                    bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                end if
            else
                call SetFFTiodim(ph_info_x, ph_info_x, 'y', iodim, iodim_howmany)
                fwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cy1%re, cy1%re, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cy1%re, cy1%re, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                call SetFFTiodim(ph_info_x, ph_info_x, 'z', iodim, iodim_howmany)
                if (periodic(3)) then
                    fwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_FORWARD, FFTW_ESTIMATE)
                    bwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_BACKWARD, FFTW_ESTIMATE)
                else
                    fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                    bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                end if
            end if
        else
            call AllocateRealFFTArray(rx1, ph_info_n, 'x')
            call SetFFTiodim(ph_info_n, ph_info_n, 'x', iodim, iodim_howmany)
            fwd_guruplan_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, dph, rx1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
            bwd_guruplan_x = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rx1, dph, (/FFTW_REDFT01/), FFTW_ESTIMATE)
            call AllocateRealFFTArray(ry1, ph_info_n, 'y')
            if (periodic(2)) then
                call AllocateComplexFFTArray(cy1, ph_info_y, 'y')
                call SetFFTiodim(ph_info_n, ph_info_y, 'y', iodim, iodim_howmany)
                fwd_guruplan_y = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, ry1, cy1, FFTW_ESTIMATE)
                call SetFFTiodim(ph_info_y, ph_info_n, 'y', iodim, iodim_howmany)
                bwd_guruplan_y = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cy1, ry1, FFTW_ESTIMATE)
                call AllocateComplexFFTArray(cz1, ph_info_y, 'z')
                call SetFFTiodim(ph_info_y, ph_info_y, 'z', iodim, iodim_howmany)
                if (periodic(3)) then
                    fwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_FORWARD, FFTW_ESTIMATE)
                    bwd_guruplan_z = fftw_plan_guru_dft(1, iodim, 2, iodim_howmany, cz1, cz1, FFTW_BACKWARD, FFTW_ESTIMATE)
                else
                    fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                    bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, cz1%re, cz1%re, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                end if
            else
                call SetFFTiodim(ph_info_n, ph_info_n, 'y', iodim, iodim_howmany)
                fwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, ry1, ry1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                bwd_guruplan_y = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, ry1, ry1, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                call AllocateRealFFTArray(rz1, ph_info_n, 'z')
                if (periodic(3)) then
                    call AllocateComplexFFTArray(cz1, ph_info_z, 'z')
                    call SetFFTiodim(ph_info_n, ph_info_z, 'z', iodim, iodim_howmany)
                    fwd_guruplan_z = fftw_plan_guru_dft_r2c(1, iodim, 2, iodim_howmany, rz1, cz1, FFTW_ESTIMATE)
                    call SetFFTiodim(ph_info_z, ph_info_n, 'z', iodim, iodim_howmany)
                    bwd_guruplan_z = fftw_plan_guru_dft_c2r(1, iodim, 2, iodim_howmany, cz1, rz1, FFTW_ESTIMATE)
                else
                    call SetFFTiodim(ph_info_n, ph_info_n, 'z', iodim, iodim_howmany)
                    fwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rz1, rz1, (/FFTW_REDFT10/), FFTW_ESTIMATE)
                    bwd_guruplan_z = fftw_plan_guru_r2r(1, iodim, 2, iodim_howmany, rz1, rz1, (/FFTW_REDFT01/), FFTW_ESTIMATE)
                end if
            end if
        end if
        ! Check if plans are created properly 
        planned_x = (c_associated(fwd_guruplan_x)).and.(c_associated(bwd_guruplan_x))
        planned_y = (c_associated(fwd_guruplan_y)).and.(c_associated(bwd_guruplan_y))
        planned_z = (c_associated(fwd_guruplan_z)).and.(c_associated(bwd_guruplan_z))
        if (planned_x.and.planned_y.and.planned_z) planned = .true.
    end if

    if (.not.planned) then
        if (ismaster) print *, 'Failed to create FFTW guru plans'
        if (ismaster) print *, 'FFTW3 must be linked before MKL.'
        if (ismaster) print *, 'Please check linking order.'
        call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_ierr)
    end if

end subroutine InitPressureSolver

subroutine SetFFTiodim(info_i,info_o,pdir,iodim,iodim_howmany)

    ! use, intrinsic :: iso_c_binding
    use decomp_2d
    use param, only: nxm, nym, nzm
    use fftw_params

    implicit none

    type(decomp_info), intent(in)   :: info_i
    type(decomp_info), intent(in)   :: info_o
    character, intent(in)           :: pdir
    type(fftw_iodim), intent(inout) :: iodim(1)
    type(fftw_iodim), intent(inout) :: iodim_howmany(2)

    if (pdir.eq.'x') then

        iodim(1)%n          = nxm
        iodim(1)%is         = 1
        iodim(1)%os         = 1

        iodim_howmany(1)%n  = (info_i%xen(2) - info_i%xst(2) + 1)
        iodim_howmany(1)%is = (info_i%xen(1) - info_i%xst(1) + 1)
        iodim_howmany(1)%os = (info_o%xen(1) - info_o%xst(1) + 1)

        iodim_howmany(2)%n  = (info_i%xen(3) - info_i%xst(3) + 1)
        iodim_howmany(2)%is = (info_i%xen(1) - info_i%xst(1) + 1)*(info_i%xen(2) - info_i%xst(2) + 1)
        iodim_howmany(2)%os = (info_o%xen(1) - info_o%xst(1) + 1)*(info_o%xen(2) - info_o%xst(2) + 1)

    else if (pdir.eq.'y') then

        iodim(1)%n          = nym
        iodim(1)%is         = (info_i%yen(1) - info_i%yst(1) + 1)
        iodim(1)%os         = (info_o%yen(1) - info_o%yst(1) + 1)

        iodim_howmany(1)%n  = (info_i%yen(1) - info_i%yst(1) + 1)
        iodim_howmany(1)%is = 1
        iodim_howmany(1)%os = 1

        iodim_howmany(2)%n  = (info_i%yen(3) - info_i%yst(3) + 1)
        iodim_howmany(2)%is = (info_i%yen(1) - info_i%yst(1) + 1)*(info_i%yen(2) - info_i%yst(2) + 1)
        iodim_howmany(2)%os = (info_o%yen(1) - info_o%yst(1) + 1)*(info_o%yen(2) - info_o%yst(2) + 1)

    else if (pdir.eq.'z') then

        iodim(1)%n          = nzm
        iodim(1)%is         = (info_i%zen(1) - info_i%zst(1) + 1)*(info_i%zen(2) - info_i%zst(2) + 1)
        iodim(1)%os         = (info_o%zen(1) - info_o%zst(1) + 1)*(info_o%zen(2) - info_o%zst(2) + 1)

        iodim_howmany(1)%n  = (info_i%zen(1) - info_i%zst(1) + 1)
        iodim_howmany(1)%is = 1
        iodim_howmany(1)%os = 1

        iodim_howmany(2)%n  = (info_i%zen(2) - info_i%zst(2) + 1)
        iodim_howmany(2)%is = (info_i%zen(1) - info_i%zst(1) + 1)
        iodim_howmany(2)%os = (info_o%zen(1) - info_o%zst(1) + 1)

    end if

end subroutine SetFFTiodim

subroutine DeallocateFFTArrays

    use decomp_2d, only: decomp_info_finalize
    use fftw_params
    use pressure_decomp
    use AuxiliaryRoutines, only: DestroyRealFFTArray, DestroyComplexFFTArray

    implicit none

    call decomp_info_finalize(ph_info_n)
    call decomp_info_finalize(ph_info_x)
    call decomp_info_finalize(ph_info_y)
    call decomp_info_finalize(ph_info_z)

    call DestroyRealFFTArray(rx1)
    call DestroyRealFFTArray(ry1)
    call DestroyRealFFTArray(rz1)
    call DestroyComplexFFTArray(cx1)
    call DestroyComplexFFTArray(cy1)
    call DestroyComplexFFTArray(cz1)

    return

end subroutine DeallocateFFTArrays