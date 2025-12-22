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
    use param, only: ns,time
    use lagrangian_point_particle

    implicit none

    logical                         :: lbound,gbound
    integer                         :: np

    ! Set the number of spawned and exited particles in the substep to zero
    sub_exit = 0

    ! Reset previous body force calculation
    lpp_bdfx(:,:,:) = 0.0
    lpp_bdfy(:,:,:) = 0.0
    lpp_bdfz(:,:,:) = 0.0

    ! Spawn new particles at the beginning of the timestep if there are any sources within pencil/process
    if ((ns.eq.1).and.(src_size.gt.0)) call SpawnNewParticles

    ! Transfer particles that have moved in/out of other pencils/processes and in/out of the current pencil/process
    call TransferParticles

    ! Perform halo cell exchange on velocity curvature terms for interpolation during slip correction 
    if (lpp_scor) call CalcCurvatureTerms

    ! Loop over all existing active particles
    np = 1
    do while(np.le.lpp_actv)
        !! Check if particle is exclusive to the pencil/process
        call CheckIsParticleLocal(lpp_list(np),lbound)
        if (lbound) then
            !!! Update particle
            call UpdateParticleGridIndices(lpp_list(np))    ! Update cell-indices of the particles at their new positions
            call UpdateParticleAcceleration(lpp_list(np))   ! Add fluid-particle and particle-particle forces to particle acceleration
            call UpdateParticleVelocity(lpp_list(np))       ! Update velocity of particles using 3rd order Runge-Kutta time-stepping scheme
            call UpdateParticlePosition(lpp_list(np))       ! Update position of particles using (2nd order) Crank-Nicolson time-stepping scheme
            call UpdateParticleLifeTime(lpp_list(np))       ! Update life time of particles (time since injection)
            !!! Check if particle is going to exit the simulation domain and needs to be written to file
            call CheckIsParticleGlobal(lpp_list(np),gbound)
            if (.not.gbound) sub_exit = sub_exit + 1
            !!! Proceed to the next particle
            np = np + 1
        else
            !! Particle is out of bounds and should be deactivated.
            !!! Exchange with the last active particle unless not already the last particle
            if (np.lt.lpp_actv) lpp_list(np) = lpp_list(lpp_actv)
            !!! Decrement the number of active particle count in lpp_list
            lpp_actv = lpp_actv - 1
        end if
    end do

    ! If there are exited particles in the current substep, save the exit events
    ! Only if saving particle exit events is enabled and flow time is greater than particle exit event save start time
    if ((pex_save).and.(time.ge.pex_ssta).and.(sub_exit.gt.0)) then
        !! Increment the total number of exited particles
        lpp_exit = lpp_exit + sub_exit
        !! Check there's enough allocated space to write exit data and extend if needed (with one space extra to be safe)
        call ExtendParticleExitBuffer(pex_actv+sub_exit+1)
        !! Loop over all active particles, calculate particle exit events, and store them
        do np = 1,lpp_actv
            !!! Check if particle has exited the simulation domain
            call CheckIsParticleGlobal(lpp_list(np),gbound)
            if (.not.gbound) then
                !!! Calculate the exit information and store them
                pex_actv = pex_actv + 1
                call CalculateParticleExit(lpp_list(np),pex_list(pex_actv))
            end if
        end do
    end if

    ! Apply particle forces on the fluid
    !! Exchange halo forces from particles in the boundary cells of the pencil/process
    call UpdateHaloForces(lpp_bdfx)
    call UpdateHaloForces(lpp_bdfy)
    call UpdateHaloForces(lpp_bdfz)
    
    return


end subroutine ParticleMarchSubstep
