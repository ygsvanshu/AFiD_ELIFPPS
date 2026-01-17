!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SpawnNewParticles.F90                          !
!    CONTAINS: subroutine SpawnNewParticles               !
!                                                         !
!    PURPOSE: To spawn new particles and initialize them  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SpawnNewParticles

    use param, only: al,dt,time,rey
    use lagrangian_point_particle

    implicit none

    integer :: spwn,indx
    integer :: nsrc,nlpp
    real    :: nprv,nnxt
    real    :: ntot
    real    :: tsub

    ! Find the time just before the beginning of substep
    tsub = time - al*dt

    ! Initialize source spawn counts to zero
    src_spwn(:) = 0

    ! Count the total number of spawned particles to be added to main particle list
    do nsrc = 1,src_size
        !! Check for non-zero frequency
        if (src_list(nsrc)%src_frq.gt.0.0) then
            nprv = (tsub - src_list(nsrc)%src_sta)*src_list(nsrc)%src_frq                   ! Fractional number of spawn event triggers before beginning of substep
            nnxt = (time - src_list(nsrc)%src_sta)*src_list(nsrc)%src_frq                   ! Fractional number of spawn event triggers until end of substep
            ntot = (src_list(nsrc)%src_end - src_list(nsrc)%src_sta)*src_list(nsrc)%src_frq ! Fractional number of spawn event triggers until source end time
            spwn = ceiling(nprv)                                                            ! Integer number of spawn event triggers already from source 
            do while (((real(spwn).ge.nprv).and.(real(spwn).lt.nnxt)).and.((real(spwn).ge.0).and.(real(spwn).le.ntot)))
                src_spwn(nsrc) = src_spwn(nsrc) + 1
                spwn = spwn + 1
            end do
        end if
    end do

    ! Update the total spawned particle count
    lpp_spwn = lpp_spwn + sum(src_spwn)

    ! Find index position at which spawned particles must be added
    indx = lpp_actv

    ! Update the total active particle count
    lpp_actv = lpp_actv + sum(src_spwn)

    ! If needed, extend the particle list
    call ExtendParticleListBuffer(lpp_actv)

    ! Loop over sources and initialize spawned particles
    do nsrc = 1,src_size
        !! Check for non-zero frequency and non-zero spawned particles
        if ((src_list(nsrc)%src_frq.gt.0.0).and.(src_spwn(nsrc).gt.0)) then
            !!! Get the first spawn trigger instance after the start of substep 
            spwn = ceiling((tsub - src_list(nsrc)%src_sta)*src_list(nsrc)%src_frq)
            do nlpp = 1,src_spwn(nsrc)
                !!!! Update the index at which particle is spawned
                indx = indx + 1
                !!!! Copy the source index
                lpp_list(indx)%src_idx    = src_list(nsrc)%src_idx
                !!!! Initialize the [xc,yc,zc] grid indices to the source grid indices
                lpp_list(indx)%grc_idx(1) = src_list(nsrc)%grc_idx(1)
                lpp_list(indx)%grc_idx(2) = src_list(nsrc)%grc_idx(2)
                lpp_list(indx)%grc_idx(3) = src_list(nsrc)%grc_idx(3)
                !!!! Initialize the [xm,ym,zm] grid indices to the source grid indices
                lpp_list(indx)%grm_idx(1) = src_list(nsrc)%grm_idx(1)
                lpp_list(indx)%grm_idx(2) = src_list(nsrc)%grm_idx(2)
                lpp_list(indx)%grm_idx(3) = src_list(nsrc)%grm_idx(3)
                !!!! Set the lifetime of particle at injection to a suitable negative value 
                !!!! Important to keep it negative or at least zero as this helps ParticleMarchSubstep to detect the spawn and use truncated RK3 
                lpp_list(indx)%lpp_lft    = tsub - src_list(nsrc)%src_sta - ((spwn + nlpp - 1)/src_list(nsrc)%src_frq)
                !!!! Initialize particle diameter
                lpp_list(indx)%lpp_dia    = src_list(nsrc)%src_dia
                !!!! Initialize particle density
                lpp_list(indx)%lpp_den    = src_list(nsrc)%src_den
                !!!! Initialize particle Reynolds number
                lpp_list(indx)%lpp_rey    = src_list(nsrc)%src_dia*norm2(src_list(nsrc)%src_vel)*rey
                !!!! Initialize particle position to source position
                lpp_list(indx)%lpp_pos(1) = src_list(nsrc)%src_pos(1)
                lpp_list(indx)%lpp_pos(2) = src_list(nsrc)%src_pos(2) 
                lpp_list(indx)%lpp_pos(3) = src_list(nsrc)%src_pos(3) 
                !!!! Initialize particle velocity to injection velocity
                lpp_list(indx)%lpp_vel(1) = src_list(nsrc)%src_vel(1)
                lpp_list(indx)%lpp_vel(2) = src_list(nsrc)%src_vel(2)
                lpp_list(indx)%lpp_vel(3) = src_list(nsrc)%src_vel(3)
                !!!! Initialize old acceleration to zero (will be overwritten during ParticleMarchSubstep)
                lpp_list(indx)%acc_old(1) = 0.0
                lpp_list(indx)%acc_old(2) = 0.0
                lpp_list(indx)%acc_old(3) = 0.0
                !!!! Initialize current acceleration to zero (will be overwritten during ParticleMarchSubstep)
                lpp_list(indx)%acc_now(1) = 0.0
                lpp_list(indx)%acc_now(2) = 0.0
                lpp_list(indx)%acc_now(3) = 0.0
            end do
        end if
    end do

end subroutine SpawnNewParticles