!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE:     ParticleRoutines.F90                       !
!    CONTAINS: Subroutine PremultipliedDragCoefficient    !
!                                                         !
!    PURPOSE:  Calculates the drag coefficient multiplied !
!    by the particle Reynolds number computed  using      !
!    slip velocity and particle diameter.                 !
!    Premultiplied in this way to ensure that the drag    !
!    force doesn't encounter singularities (Inf/NaN)      !
!    during drag forcecomputation if slip/Reynolds number !
!    of the particle becomes zero.                        !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PremultipliedDragCoefficient(prey,pcfd)

    use lagrangian_point_particle

    implicit none

    real, intent(in)    :: prey     ! Particle Reynolds number
    real, intent(out)   :: pcfd     ! Premultiplied drag coefficient

    if (lpp_dmod.eq.STOKES) then
        ! Stokes 
        pcfd = 24.0
    else if (lpp_dmod.eq.SCHNAU) then
        ! Schiller-Naumann 
        if (prey.le.1000.0) then
            pcfd = 24.0*(1.0 + (0.15*(prey**0.687)))
        else
            pcfd = 0.44*prey
        end if
    end if

    return

end subroutine PremultipliedDragCoefficient