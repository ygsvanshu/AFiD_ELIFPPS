!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ReadParticleSources.F90                        !
!    CONTAINS: subroutine ReadParticleSources             !
!                                                         !
!    PURPOSE: To read particle sources during the start   !
!    of a simulation run.                                 !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadParticleSources

    use, intrinsic :: iso_fortran_env, only : iostat_end
    use mpih
    use decomp_2d, only: xstart,xend,nrank
    use param
    use lagrangian_point_particle

    implicit none

    logical                 :: exists
    logical                 :: nloc
    character*200           :: filename,dsetname
    character*4             :: charsrun
    integer                 :: ierror,numsrc
    integer                 :: srun
    integer                 :: psrr,pstt,psst,psen
    type(particle_source)   :: nsrc

    filename = trim("Inputs/particle_sources.in")
    inquire(file=filename, exist=exists)

    ! Check if the particle source input file exists
    if (exists) then

        !! Read to find out the total number of local sources
        open(99, file=trim(filename))
        ierror = 0
        do while (ierror.ne.iostat_end)
            read(99, *, iostat=ierror) nsrc%src_sta,nsrc%src_end,nsrc%src_frq,nsrc%src_dia,nsrc%src_den,nsrc%src_pos(1),nsrc%src_pos(2),nsrc%src_pos(3),nsrc%src_vel(1),nsrc%src_vel(2),nsrc%src_vel(3)
            if (ierror.eq.0) then
                !! Check if the source is local to the pencil/process
                call CheckIsSourceLocal(nsrc,nloc)
                if (nloc) then
                    src_size = src_size + 1
                end if
            end if
        end do
        close(99)

        ! Allocate required amount of space in the source buffer
        allocate(src_list(src_size))

        ! Actually read the sources
        open(99, file=trim(filename))
        ierror = 0
        numsrc = 1
        do while (ierror.ne.iostat_end)
            read(99, *, iostat=ierror) nsrc%src_sta,nsrc%src_end,nsrc%src_frq,nsrc%src_dia,nsrc%src_den,nsrc%src_pos(1),nsrc%src_pos(2),nsrc%src_pos(3),nsrc%src_vel(1),nsrc%src_vel(2),nsrc%src_vel(3)
            if (ierror.eq.0) then
                !! Check if the source is local to the pencil/process
                call CheckIsSourceLocal(nsrc,nloc)
                if (nloc) then
                    !!! Update the nozzle index
                    nsrc%src_idx = numsrc
                    !!! Find and update cell indices for c-grid (nodes)
                    call GetLocationCellIndex(nsrc%src_pos(1),xc(xstart(1)  :xend(1)+1),xstart(1)  ,xend(1)+1,nsrc%grc_idx(1))
                    call GetLocationCellIndex(nsrc%src_pos(2),yc(xstart(2)  :xend(2)+1),xstart(2)  ,xend(2)+1,nsrc%grc_idx(2))
                    call GetLocationCellIndex(nsrc%src_pos(3),zc(xstart(3)  :xend(3)+1),xstart(3)  ,xend(3)+1,nsrc%grc_idx(3))
                    !!! Find and update cell indices for m-grid (cell-center)
                    call GetLocationCellIndex(nsrc%src_pos(1),xm(xstart(1)-1:xend(1)+1),xstart(1)-1,xend(1)+1,nsrc%grm_idx(1))
                    call GetLocationCellIndex(nsrc%src_pos(2),ym(xstart(2)-1:xend(2)+1),xstart(2)-1,xend(2)+1,nsrc%grm_idx(2))
                    call GetLocationCellIndex(nsrc%src_pos(3),zm(xstart(3)-1:xend(3)+1),xstart(3)-1,xend(3)+1,nsrc%grm_idx(3))
                    !!! Set the current source in the list
                    src_list(numsrc) = nsrc
                    !!! Incriment the index of source list
                    numsrc = numsrc + 1
                end if
            end if
        end do
        close(99)

    end if

    ! Call an MPI_ALLREDUCE to get the total particle source count from every pencil/process
    call MPI_ALLREDUCE(src_size,src_ntot,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)

    return

end subroutine ReadParticleSources 