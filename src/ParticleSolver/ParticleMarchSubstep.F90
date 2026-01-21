!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ParticleMarchSubstep.F90                       !
!    CONTAINS: subroutine ParticleMarchSubstep            !
!                                                         !
!    PURPOSE: The main subroutine that advances particles !
!    in time at each substep                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ParticleMarchSubstep

    use mpih
    use param, only: time
    use lagrangian_point_particle

    implicit none

    logical                         :: lbound,gbound
    integer                         :: np
    real                            :: nfrc(3,3)
    real                            :: dlft
    type(particle_exit)             :: tpex

    ! Reset previous body force calculation
    lpp_bdfx(:,:,:) = 0.0
    lpp_bdfy(:,:,:) = 0.0
    lpp_bdfz(:,:,:) = 0.0

    ! Spawn new particles if there are any sources within pencil/process
    if (src_size.gt.0) call SpawnNewParticles

    ! Transfer particles that have moved in/out of other pencils/processes and in/out of the current pencil/process
    call TransferParticles

    ! Perform halo cell exchange on velocity curvature terms for interpolation during slip correction 
    if (lpp_scor) call CalcCurvatureTerms

    !! Preemptively check if there's enough allocated space to write exit data and extend if needed
    if ((pex_save).and.(time.ge.pex_ssta)) call ExtendParticleExitBuffer(pex_actv+lpp_actv)

    ! Loop over all existing active particles
    np = 1
    do while(np.le.lpp_actv)
        !! Check if particle is exclusive to the pencil/process
        call CheckIsParticleLocal(lpp_list(np),lbound)
        if (lbound) then
            !!! Update particle
            call UpdateParticleGridIndices(lpp_list(np))                ! Update cell-indices of the particles at their new positions
            call InitParticleAcceleration(lpp_list(np))                 ! Initialize particle acceleration
            call AddParticleAccelerationDrag(lpp_list(np),nfrc)         ! Add fluid drag forces to particle acceleration   
            call AddParticleAccelerationGravity(lpp_list(np))           ! Add particle acceleration due to particle apparent weight   
            call UpdateParticleVelocity(lpp_list(np))                   ! Update velocity of particles using 3rd order Runge-Kutta time-stepping scheme
            call UpdateParticlePosition(lpp_list(np))                   ! Update position of particles using (2nd order) Crank-Nicolson time-stepping scheme
            call UpdateParticleLifeTime(lpp_list(np))                   ! Update life time of particles (time since injection)
            call CalculateParticleSubStepTime(lpp_list(np),dlft)        ! Calculate time spent by particle in the domain during the substep
            call CheckIsParticleGlobal(lpp_list(np),gbound)             ! Check if particle is going to exit the simulation domain
            if (.not.gbound) then
                call CalculateParticleExit(lpp_list(np),tpex)           ! Calculate the exit in a temporary event buffer primarily to know the exit time.
                call UpdateParticleSubStepTime(lpp_list(np),tpex,dlft)  ! Update time spent by particle in the domain during the substep
                if ((pex_save).and.(time.ge.pex_ssta)) then
                    lpp_exit = lpp_exit + 1                             ! Update count of exit events
                    pex_actv = pex_actv + 1                             ! Update total number of exit events are to be written
                    pex_list(pex_actv)  = tpex                          ! Update the new exit event are to be written
                end if
            end if
            call ApplyParticleDragForce(lpp_list(np),nfrc,dlft)         ! Apply the drag force from the particle on the Eulerian field
            !!! Proceed to the next particle
            np = np + 1
        else
            !!! Particle is out of bounds and should be deactivated.
            !!! Copy the last active particle unless not already the last particle
            if (np.lt.lpp_actv) lpp_list(np) = lpp_list(lpp_actv)
            !!! Decrement the number of active particle count in lpp_list, deactivating the last (duplicate) particle
            lpp_actv = lpp_actv - 1
        end if
    end do

    ! Apply particle forces on the fluid
    !! Exchange halo forces from particles in the boundary cells of the pencil/process
    call UpdateHaloForces(lpp_bdfx)
    call UpdateHaloForces(lpp_bdfy)
    call UpdateHaloForces(lpp_bdfz)
    
    return

end subroutine ParticleMarchSubstep