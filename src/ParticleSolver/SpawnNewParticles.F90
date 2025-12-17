!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SpawnNewParticles.F90                          !
!    CONTAINS: subroutine SpawnNewParticles               !
!                                                         !
!    PURPOSE: To spawn new particles and initialize them  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SpawnNewParticles

    use mpih
    use param, only: al,dt,time,pi,rey
    use local_arrays, only: vx,vy,vz
    use lagrangian_point_particle

    implicit none

    logical             :: spwn
    integer             :: nsrc,ncnt,nlpp,indx
    real                :: tsta,tend
    real                :: cffc(3)
    real                :: cffm(3)
    real                :: slvx,slvy,slvz,slip

    ! Set the tolerance interval for triggering spawning of particles
    tsta = time - al*dt - lpp_stol
    tend = time - al*dt + lpp_stol

    ! Count the total number of spawned particles to be added to main particle list
    ncnt = 0
    do nsrc = 1,src_size
        !! Initialize boolean spawn as true
        spwn = .true.
        !! Check if it's before the source end time
        spwn = spwn.and.(tsta.le.src_list(nsrc)%src_end)
        !! Check if the next spawn event is within tolerance of the current time
        spwn = spwn.and.(src_list(nsrc)%src_sta.ge.tsta)
        spwn = spwn.and.(src_list(nsrc)%src_sta.le.tend)
        !! Check if the frequency is a non-zero positive value
        spwn = spwn.and.(src_list(nsrc)%src_frq.gt.0.0)
        !! If spawn conditions are met, increment spawn count
        if (spwn) ncnt = ncnt + 1
    end do

    ! Update the total spawned particle count
    lpp_spwn = lpp_spwn + ncnt

    ! Find index position at which spawned particles must be added
    indx = lpp_actv

    ! Calcualte total size needed to append to the current particle list
    ncnt = ncnt + indx
    ! If needed, extend the particle list
    call ExtendParticleListBuffer(ncnt+1)

    ! Loop over sources and initialize spawned particles
    do nsrc = 1,src_size

        !! Initialize boolean spawn as true
        spwn = .true.
        !! Check if it's before the source end time
        spwn = spwn.and.(tsta.le.src_list(nsrc)%src_end)
        !! Check if the next spawn event is within tolerance of the current time
        spwn = spwn.and.(src_list(nsrc)%src_sta.ge.tsta)
        spwn = spwn.and.(src_list(nsrc)%src_sta.le.tend)
        !! Check if the frequency is a non-zero positive value
        spwn = spwn.and.(src_list(nsrc)%src_frq.gt.0.0)
        !! If spawn conditions are met, increment spawn count
        if (spwn) then

            !!! Update the index at which particle is spawned
            indx = indx + 1
            !!! Copy the source index
            lpp_list(indx)%src_idx    = src_list(nsrc)%src_idx
            !!! Initialize the [xc,yc,zc] grid indices to the source grid indices
            lpp_list(indx)%grc_idx(1) = src_list(nsrc)%grc_idx(1)
            lpp_list(indx)%grc_idx(2) = src_list(nsrc)%grc_idx(2)
            lpp_list(indx)%grc_idx(3) = src_list(nsrc)%grc_idx(3)
            !!! Initialize the [xm,ym,zm] grid indices to the source grid indices
            lpp_list(indx)%grm_idx(1) = src_list(nsrc)%grm_idx(1)
            lpp_list(indx)%grm_idx(2) = src_list(nsrc)%grm_idx(2)
            lpp_list(indx)%grm_idx(3) = src_list(nsrc)%grm_idx(3)
            !!! Set the lifetime of particle at injection to zero 
            lpp_list(indx)%lpp_lft    = 0.0
            !!! Initialize particle diameter and density
            lpp_list(indx)%lpp_dia    = src_list(nsrc)%src_dia
            lpp_list(indx)%lpp_den    = src_list(nsrc)%src_den
            !!! Initialize particle position to source position
            lpp_list(indx)%lpp_pos(1) = src_list(nsrc)%src_pos(1)
            lpp_list(indx)%lpp_pos(2) = src_list(nsrc)%src_pos(2) 
            lpp_list(indx)%lpp_pos(3) = src_list(nsrc)%src_pos(3) 
            !!! Initialize particle velocity to injection velocity
            lpp_list(indx)%lpp_vel(1) = src_list(nsrc)%src_vel(1)
            lpp_list(indx)%lpp_vel(2) = src_list(nsrc)%src_vel(2)
            lpp_list(indx)%lpp_vel(3) = src_list(nsrc)%src_vel(3)
            !!! Initialize old acceleration to zero
            lpp_list(indx)%acc_old(1) = 0.0
            lpp_list(indx)%acc_old(2) = 0.0
            lpp_list(indx)%acc_old(3) = 0.0
            !!! Initialize current acceleration to zero
            lpp_list(indx)%acc_now(1) = 0.0
            lpp_list(indx)%acc_now(2) = 0.0
            lpp_list(indx)%acc_now(3) = 0.0
            !!! Initialize the particle reynolds number using the slip (no slip correction needed at initialization)
            !!!! Compute interpolation coefficients
            call CalcTrilinearInterpolationCoefficients(lpp_list(indx),cffm,cffc)
            slvx = 0.0
            slvy = 0.0
            slvz = 0.0
            !!!! Apply trilinear interpolation with previously calculated coefficients to get fluid velocity
            call ApplyTrilinearInterpolation(lpp_list(indx),cffm,cffc,'x',1,vx,slvx)
            call ApplyTrilinearInterpolation(lpp_list(indx),cffm,cffc,'y',1,vy,slvy)
            call ApplyTrilinearInterpolation(lpp_list(indx),cffm,cffc,'z',1,vz,slvz)
            !!!! Subtract fluid velocity from source velocity to get slip velocity
            slvx = src_list(nsrc)%src_vel(1) - slvx
            slvy = src_list(nsrc)%src_vel(2) - slvy
            slvz = src_list(nsrc)%src_vel(3) - slvz
            slip = sqrt((slvx**2.0) + (slvy**2.0) + (slvz**2.0))
            lpp_list(indx)%lpp_rey    = lpp_list(indx)%lpp_dia*rey*slip
            !!! Frequency is a non-zero positive number, so update the next spawn time instant
            src_list(nsrc)%src_sta = src_list(nsrc)%src_sta + (1.0/src_list(nsrc)%src_frq)

        end if

    end do

    ! Update the number of active particles to include the newly spawned particles
    lpp_actv = ncnt

    return

end subroutine SpawnNewParticles