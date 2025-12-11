!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         !
!    PURPOSE: Initialization routine. Sets initial        !
!    conditions for pressure and velocity components.     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateInitialConditions

    use decomp_2d, only: xstart,xend
    use param
    use local_arrays, only: vx,vy,vz,pr
    
    implicit none

    vx(:,:,:) = 0.0
    vy(:,:,:) = 0.0
    vz(:,:,:) = 0.0
    pr(:,:,:) = 0.0

    return

end subroutine CreateInitialConditions