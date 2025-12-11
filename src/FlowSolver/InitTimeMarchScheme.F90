!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: InitTimeMarchScheme.F90                        !
!    CONTAINS: subroutine InitTimeMarchScheme             !
!                                                         !
!    PURPOSE: Initialize the time-marching constants for  !
!     the integrator                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitTimeMarchScheme

    use param

    implicit none

    if (nsst .gt. 1) then
        gam(1) =   8.0/15.0
        gam(2) =   5.0/12.0
        gam(3) =   3.0/4.0
        rom(1) =   0.0
        rom(2) = -17.0/60.0
        rom(3) =  -5.0/12.0
    else
        gam(1) =   1.5
        gam(2) =   0.0
        gam(3) =   0.0
        rom(1) =  -0.5
        rom(2) =   0.0
        rom(3) =   0.0
    end if

    do ns = 1, nsst
        alm(ns) = (gam(ns) + rom(ns))
    end do

    return

end subroutine InitTimeMarchScheme