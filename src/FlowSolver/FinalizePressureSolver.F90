subroutine FinalizePressureSolver

    use pressure_solver

    implicit none

    ! Deallocate Poisson solver arrays

    if (allocated(qfwx))  deallocate(qfwx)
    if (allocated(qbwx))  deallocate(qbwx)
    if (allocated(geix)) deallocate(geix)
    if (allocated(feix))   deallocate(feix)
    if (allocated(eigx))  deallocate(eigx)

    if (allocated(qfwy))  deallocate(qfwy)
    if (allocated(qbwy))  deallocate(qbwy)
    if (allocated(geiy)) deallocate(geiy)
    if (allocated(feiy))   deallocate(feiy)
    if (allocated(eigy))  deallocate(eigy)

    if (allocated(qfwz))  deallocate(qfwz)
    if (allocated(qbwz))  deallocate(qbwz)
    if (allocated(geiz)) deallocate(geiz)
    if (allocated(feiz))   deallocate(feiz)
    if (allocated(eigz))  deallocate(eigz)

    ! Deallocate GEMM buffer
    
    if (allocated(gwrk)) deallocate(gwrk)

    ! Destroy FFTW plans

    if (c_associated(pfwx)) call fftw_destroy_plan(pfwx)
    if (c_associated(pbwx)) call fftw_destroy_plan(pbwx)
    if (c_associated(pfwy)) call fftw_destroy_plan(pfwy)
    if (c_associated(pbwy)) call fftw_destroy_plan(pbwy)
    if (c_associated(pfwz)) call fftw_destroy_plan(pfwz)
    if (c_associated(pbwz)) call fftw_destroy_plan(pbwz)

end subroutine FinalizePressureSolver