!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CopyPressure.F90                               !
!    CONTAINS: subroutine CopyPressure                    !
!                                                         !
!    PURPOSE: Copies the solved pressure correction       !
!    from array dph which does not have halo cells        !
!    to array dphhalo which has halo cells                !
!    in order to be able to compute gradients.            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CopyPressureHalo

    use param
    use local_arrays, only: dph, dphhalo
    use decomp_2d

    implicit none

    ! This copy can be avoided by changing transpose_x_to_y_real and
    ! transpose_y_to_x_real so these routines can handle arrays with
    ! halo. This copy is a defacto array temporary. Using inferred size
    ! arrays in the transpose calls results in 5 more of these, and more
    ! memory usage. Time spent on this copy is 0.1% for 65^3 grid.

    dphhalo(:,:,:) = 0.0
    
    dphhalo(1:nxm,xstart(2):xend(2),xstart(3):xend(3)) = dph(1:nxm,xstart(2):xend(2),xstart(3):xend(3))

    if (periodic(1)) then
        dphhalo(0 ,:,:) = dphhalo(nxm,:,:)
        dphhalo(nx,:,:) = dphhalo(1  ,:,:)
    else
        dphhalo(0 ,:,:) = dphhalo(1  ,:,:)
        dphhalo(nx,:,:) = dphhalo(nxm,:,:)
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            dphhalo(:,0,:) = dphhalo(:,1,:)
        end if
        if (xend(2).eq.nym) then
            dphhalo(:,ny,:) = dphhalo(:,nym,:)
        end if
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            dphhalo(:,:,0) = dphhalo(:,:,1)
        end if
        if (xend(3).eq.nzm) then
            dphhalo(:,:,nz) = dphhalo(:,:,nzm)
        end if
    end if

end subroutine CopyPressureHalo