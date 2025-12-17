!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: TimeMarcher.F90                                !
!    CONTAINS: subroutine TimeMarcher                     !
!                                                         !
!    PURPOSE: Main time integrating routine, which calls  !
!     other subroutines for calculating the Navier-Stokes !
!     equations and advancing velocity and temperature in !
!     time                                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine TimeMarcher

    use mpih
    use decomp_2d
    use param
    use local_arrays
    use lagrangian_point_particle

    implicit none

    do ns=1,nsst

        al=alm(ns)
        ga=gam(ns)
        ro=rom(ns)

        time = time + (al*dt)

        call CalcBoundaryConditionVX
        call CalcBoundaryConditionVY
        call CalcBoundaryConditionVZ

        call ExplicitTermsVX
        call ExplicitTermsVY
        call ExplicitTermsVZ

        call AddSourceTerms

        if (particle) then
            ! COUPLE LAGRANGIAN POINT PARTICLES
            call ParticleMarchSubstep
            call AddParticleForces
            ! PrintPencilParticleCount ! FOR DEBUGGING
        end if

        call ImplicitAndUpdateVX
        call ImplicitAndUpdateVY
        call ImplicitAndUpdateVZ

        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)

        call CorrectGlobalDivergence
        call CalcLocalDivergence
        call SolvePressureCorrection

        call CopyPressureHalo
        call update_halo(dphhalo,lvlhalo)

        call CorrectPressure
        call update_halo(pr,lvlhalo)
        
        call CorrectVelocity

        call update_halo(vx,lvlhalo)
        call update_halo(vy,lvlhalo)
        call update_halo(vz,lvlhalo)

    enddo

    return

end subroutine TimeMarcher