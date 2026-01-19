!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CalcMaxParticleTimestep.F90                    !
!    CONTAINS: subroutine CalcMaxPDT                      !
!                                                         !
!    PURPOSE: Computes the maximum timestep for stable    !
!    timestepping of particles                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcMaxPDT(instPDT)

    use mpih
    use param, only: al,dt,time,rey,dxc,dyc,dzc
    use lagrangian_point_particle

    implicit none

    real, intent(inout) :: instPDT
    integer             :: nlpp,nsrc
    real                :: pdtx,pdty,pdtz
    real                :: res_dummy

    ! Initialize the inverse of max timestep as smallest floating point value (i.e., the max timestep is initialized to be very large)
    instPDT = tiny(0.0)
    ! Loop over all the active particles
    do nlpp = 1,lpp_actv
        !! Ensure that no particle can skip grid cells in a timestep
        !! i.e. The distance covered by a particle in a single timestep cannot be greater than the grid dimensions 
        pdtx    = abs(lpp_list(nlpp)%lpp_vel(1))/(dxc(lpp_list(nlpp)%grc_idx(1)))
        pdty    = abs(lpp_list(nlpp)%lpp_vel(2))/(dyc(lpp_list(nlpp)%grc_idx(2)))
        pdtz    = abs(lpp_list(nlpp)%lpp_vel(3))/(dzc(lpp_list(nlpp)%grc_idx(3)))
        instPDT = max(instPDT,((pdtx + pdty + pdtz)/lpp_clim))
        !! Ensure that the timestep is at most a fixed fraction of the particle response timescale
        instPDT = max(instPDT,((18.0/(lpp_list(nlpp)%lpp_den*lpp_list(nlpp)%lpp_dia*lpp_list(nlpp)%lpp_dia*rey))/lpp_tlim))
    end do

    call MpiAllMaxRealScalar(instPDT,res_dummy)
    instPDT = res_dummy

end subroutine CalcMaxPDT