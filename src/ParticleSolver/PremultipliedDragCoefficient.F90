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

subroutine PremultipliedDragCoefficient(p,pcfd)

    use lagrangian_point_particle

    implicit none

    type(particle_data), intent(in) :: p        ! Particle data
    real, intent(out)               :: pcfd     ! Premultiplied drag coefficient

    if (lpp_dmod.eq.STOKES) then
        ! Stokes 
        pcfd = 24.0
    else if (lpp_dmod.eq.SCHNAU) then
        ! Schiller-Naumann 
        if (p%lpp_rey.le.1000.0) then
            pcfd = 24.0*(1.0 + (0.15*(p%lpp_rey**0.687)))
        else
            pcfd = 0.44*p%lpp_rey
        end if
    else if (lpp_dmod.eq.MORALE) then
        ! Morsi-Alexander
        if (p%lpp_rey.le.0.1) then
            pcfd = 24.0
        else if ((p%lpp_rey.gt.0.1).and.(p%lpp_rey.le.1.0)) then
            pcfd = 3.690*p%lpp_rey + 22.73 + 0.0903/p%lpp_rey
        else if ((p%lpp_rey.gt.1.0).and.(p%lpp_rey.le.10.0)) then
            pcfd = 1.222*p%lpp_rey + 29.1667 - 3.8889/p%lpp_rey
        else if ((p%lpp_rey.gt.10.0).and.(p%lpp_rey.le.100.0)) then
            pcfd = 0.6167*p%lpp_rey + 46.50 - 116.67/p%lpp_rey
        else if ((p%lpp_rey.gt.100.0).and.(p%lpp_rey.le.1000.0)) then
            pcfd = 0.3644*p%lpp_rey + 98.33 - 2778/p%lpp_rey
        else if ((p%lpp_rey.gt.1000.0).and.(p%lpp_rey.le.5000.0)) then
            pcfd = 0.357*p%lpp_rey + 148.62 - 47500/p%lpp_rey
        else if ((p%lpp_rey.gt.5000.0).and.(p%lpp_rey.le.10000.0)) then
            pcfd = 0.46*p%lpp_rey - 490.546 - 578700/p%lpp_rey
        else if (p%lpp_rey.gt.10000.0) then
            pcfd = 0.5191*p%lpp_rey - 1662.5 - 5416700/p%lpp_rey
        end if
    end if

end subroutine PremultipliedDragCoefficient