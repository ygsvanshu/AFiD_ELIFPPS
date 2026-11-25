!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CopyPressure.F90                               !
!    CONTAINS: subroutine CopyPressure                    !
!                                                         !
!    PURPOSE: Copies the solved pressure correction       !
!    from array rhx which does not have halo cells        !
!    to array dphhalo which has halo cells                !
!    in order to be able to compute gradients, and        !
!    imposes the pressure correction boundary conditions. !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CopyPressureHalo

    use decomp_2d
    use param
    use local_arrays, only: rhx, dphhalo

    implicit none

    ! This copy can be avoided by changing transpose_x_to_y_real and
    ! transpose_y_to_x_real so these routines can handle arrays with
    ! halo. This copy is a defacto array temporary. Using inferred size
    ! arrays in the transpose calls results in 5 more of these, and more
    ! memory usage. Time spent on this copy is 0.1% for 65^3 grid.

    dphhalo(:,:,:) = 0.0
    
    dphhalo(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) = rhx(1:nxm,xstart(2):xend(2),xstart(3):xend(3))

    call ImposeBoundaryConditionPC

end subroutine CopyPressureHalo