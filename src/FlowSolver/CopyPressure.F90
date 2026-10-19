!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CopyPressure.F90                               !
!    CONTAINS: subroutine CopyPressure                    !
!                                                         !
!    PURPOSE: Copies the solved pressure correction       !
!    from array rhx which does not have halo cells        !
!    to array dphhalo which has halo cells                !
!    in order to be able to compute gradients.            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CopyPressureHalo

    use decomp_2d
    use param
    use local_arrays, only: rhx, dphhalo
    use boundary_arrays, only: DIRICHLET, NEUMANN
    use pressure_solver, only: pcbcxs, pcbcxe, pcbcys, pcbcye, pcbczs, pcbcze

    implicit none

    ! This copy can be avoided by changing transpose_x_to_y_real and
    ! transpose_y_to_x_real so these routines can handle arrays with
    ! halo. This copy is a defacto array temporary. Using inferred size
    ! arrays in the transpose calls results in 5 more of these, and more
    ! memory usage. Time spent on this copy is 0.1% for 65^3 grid.

    dphhalo(:,:,:) = 0.0
    
    dphhalo(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) = rhx(1:nxm,xstart(2):xend(2),xstart(3):xend(3))

    if (periodic(1)) then
        dphhalo(0 ,:,:) = dphhalo(nxm,:,:)
        dphhalo(nx,:,:) = dphhalo(1  ,:,:)
    else
        if (pcbcxs.eq.DIRICHLET) dphhalo(0 ,:,:) = -dphhalo(1  ,:,:)
        if (pcbcxs.eq.NEUMANN)   dphhalo(0 ,:,:) =  dphhalo(1  ,:,:)
        if (pcbcxe.eq.DIRICHLET) dphhalo(nx,:,:) = -dphhalo(nxm,:,:)
        if (pcbcxe.eq.NEUMANN)   dphhalo(nx,:,:) =  dphhalo(nxm,:,:)
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            if (pcbcys.eq.DIRICHLET) dphhalo(:,0,:) = -dphhalo(:,1,:)
            if (pcbcys.eq.NEUMANN)   dphhalo(:,0,:) =  dphhalo(:,1,:)
        end if
        if (xend(2).eq.nym) then
            if (pcbcye.eq.DIRICHLET) dphhalo(:,ny,:) = -dphhalo(:,nym,:)
            if (pcbcye.eq.NEUMANN)   dphhalo(:,ny,:) =  dphhalo(:,nym,:)
        end if
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            if (pcbczs.eq.DIRICHLET) dphhalo(:,:,0) = -dphhalo(:,:,1)
            if (pcbczs.eq.NEUMANN)   dphhalo(:,:,0) =  dphhalo(:,:,1)
        end if
        if (xend(3).eq.nzm) then
            if (pcbcze.eq.DIRICHLET) dphhalo(:,:,nz) = -dphhalo(:,:,nzm)
            if (pcbcze.eq.NEUMANN)   dphhalo(:,:,nz) =  dphhalo(:,:,nzm)
        end if
    end if

end subroutine CopyPressureHalo