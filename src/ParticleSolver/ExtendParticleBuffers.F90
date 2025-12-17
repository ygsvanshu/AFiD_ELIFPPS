!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ExtendParticleBuffers.F90                      ! 
!    CONTAINS: subroutines ExtendParticleListBuffer       !
!    and ExtendParticleSendBuffer                         !
!                                                         !
!    PURPOSE: Computes the maximum timestep for stable    !
!    timestepping of particles                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ExtendParticleListBuffer(bfsz)

    use mpih
    use lagrangian_point_particle

    implicit none

    integer, intent(in)                             :: bfsz
    integer                                         :: ierr,np
    type(particle_data), allocatable, dimension(:)  :: lpp_buff

    ! Check if particle list buffer is already allocated
    if (allocated(lpp_list)) then
        ! Check if particle list buffer has insufficient size
        if (size(lpp_list).lt.bfsz) then
            ! Allocate the temporary particle list buffer
            allocate(lpp_buff(bfsz),stat=ierr)
            ! If errors are encountered, it's probably because we've run out of memory
            if (ierr.ne.0) then
                write(*,*) "Dynamic memory allocation failed for temporary particle list buffer at rank ",mpi_rank
                call MpiAbort
            end if
            ! Copy existing force list data
            do np = 1,lpp_actv
                lpp_buff(np) = lpp_list(np)
            end do
            ! Move allocation from temporary force list to the main force list
            call move_alloc(lpp_buff,lpp_list)
        end if
    else
        ! Particle list buffer is not allocated yet, simply allocate to required size
        allocate(lpp_list(bfsz),stat=ierr)
        if (ierr.ne.0) then
            write(*,*) "Dynamic memory allocation failed for particle list buffer at rank ",mpi_rank
            call MpiAbort
        end if
    end if

    return

end subroutine ExtendParticleListBuffer

subroutine ExtendParticleSendBuffer(bfsz,dirc)

    use mpih
    use lagrangian_point_particle

    implicit none

    integer, intent(in)     :: bfsz
    character, intent(in)   :: dirc
    integer                 :: ierr

    if (dirc.eq."p") then
        if (allocated(bfm_send)) then
            if (size(bfm_send).lt.bfsz) then
                deallocate(bfm_send)
                allocate(bfm_send(bfsz),stat=ierr)
                if (ierr.ne.0) then
                    write(*,*) "Dynamic memory allocation failed for particle send (prev) buffer at rank ",mpi_rank
                    call MpiAbort
                end if
            end if
        else
            allocate(bfm_send(bfsz),stat=ierr)
            if (ierr.ne.0) then
                write(*,*) "Dynamic memory allocation failed for particle send (prev) buffer at rank ",mpi_rank
                call MpiAbort
            end if
        end if
    else if (dirc.eq."n") then
        if (allocated(bfp_send)) then
            if (size(bfp_send).lt.bfsz) then
                deallocate(bfp_send)
                allocate(bfp_send(bfsz),stat=ierr)
                if (ierr.ne.0) then
                    write(*,*) "Dynamic memory allocation failed for particle send (next) buffer at rank ",mpi_rank
                    call MpiAbort
                end if
            end if
        else
            allocate(bfp_send(bfsz),stat=ierr)
            if (ierr.ne.0) then
                write(*,*) "Dynamic memory allocation failed for particle send (next) buffer at rank ",mpi_rank
                call MpiAbort
            end if
        end if
    end if

    return

end subroutine ExtendParticleSendBuffer

subroutine ExtendParticleExitBuffer(bfsz)

    use mpih
    use lagrangian_point_particle

    implicit none

    integer, intent(in)                             :: bfsz
    integer                                         :: ierr,np
    type(particle_exit), allocatable, dimension(:)  :: pex_buff

    ! Check if particle list buffer is already allocated
    if (allocated(pex_list)) then
        ! Check if particle list buffer has insufficient size
        if (size(pex_list).lt.bfsz) then
            ! Allocate the temporary particle list buffer
            allocate(pex_buff(bfsz),stat=ierr)
            ! If errors are encountered, it's probably because we've run out of memory
            if (ierr.ne.0) then
                write(*,*) "Dynamic memory allocation failed for temporary particle exit buffer at rank ",mpi_rank
                call MpiAbort
            end if
            ! Copy existing force list data
            do np = 1,pex_actv
                pex_buff(np) = pex_list(np)
            end do
            ! Move allocation from temporary force list to the main force list
            call move_alloc(pex_buff,pex_list)
        end if
    else
        ! Particle list buffer is not allocated yet, simply allocate to required size
        allocate(pex_list(bfsz),stat=ierr)
        if (ierr.ne.0) then
            write(*,*) "Dynamic memory allocation failed for particle exit buffer at rank ",mpi_rank
            call MpiAbort
        end if
    end if

    return

end subroutine ExtendParticleExitBuffer