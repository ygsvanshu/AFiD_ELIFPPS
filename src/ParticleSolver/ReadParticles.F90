!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ReadParticles.F90                              !
!    CONTAINS: subroutine ReadParticles                   !
!                                                         !
!    PURPOSE: To read saved particles from a previous run !
!    while continuing the run.                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadParticles

    use mpih
    use decomp_2d
    use param, only: ismaster,periodic,nx,ny,nz,nxm,nym,nzm,xc,yc,zc,xm,ym,zm
    use lagrangian_point_particle

    implicit none

    character*200                                   :: dsetname,filename
    logical                                         :: locl,exists
    integer                                         :: lnum,rnum
    integer                                         :: mloc,tint
    integer                                         :: lppr,lppt
    integer                                         :: lpst,lpen
    type(particle_data)                             :: tmpp
    type(particle_data), allocatable, dimension(:)  :: rlpp
    type(particle_data), allocatable, dimension(:)  :: tlpp
    logical, allocatable, dimension(:)              :: mask
    integer, allocatable, dimension(:)              :: trnk
    integer, dimension(mpi_size)                    :: scnt,sdsp,rcnt,rdsp
    real   , dimension(mpi_size)                    :: xst2,xen2,xst3,xen3

    filename = trim("continua.h5")

    ! Particle file may not exist even if fluid continua exists
    ! If no particle file is found, no particles are read in
    inquire(file=filename,exist=exists)
    call MPI_ALLREDUCE(MPI_IN_PLACE,exists,1,MPI_LOGICAL,MPI_LOR,mpi_comm,mpi_ierr)
    if (.not.exists) then 
        write(*,*) 'Continuation file not found'
        call MpiAbort
    end if

    ! Read in the particle counts
    if (ismaster) then
        dsetname = trim("/lpp_num")
        call HdfSerialReadIntScalar(filename,dsetname,lppt)
        dsetname = trim("/lpp_spwn")
        call HdfSerialReadIntScalar(filename,dsetname,tot_spwn)
        dsetname = trim("/lpp_exit")
        call HdfSerialReadIntScalar(filename,dsetname,tot_exit)
    end if
    !! Call MPI Barrier
    call MPI_BARRIER(mpi_comm,mpi_ierr)
    !! Broadcast total number of particles to all pencils/processes
    call MPI_BCAST(lppt,1,MPI_INTEGER,0,mpi_comm,mpi_ierr)
    call MPI_BCAST(tot_spwn,1,MPI_INTEGER,0,mpi_comm,mpi_ierr)
    call MPI_BCAST(tot_exit,1,MPI_INTEGER,0,mpi_comm,mpi_ierr)

    ! Read only if there's a non-zero number of particles
    if (lppt.gt.0) then

        !! Distribute the total number of particles into local pencils/processes
        lppr = int(lppt/mpi_size)
        if (mpi_rank.lt.modulo(lppt,mpi_size)) lppr = lppr + 1

        !! Get the start and end indices for the local pencil/process
        call MPI_SCAN(lppr,lpen,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
        lpst = lpen - lppr + 1

        !! Allocate temporary read buffer
        allocate(rlpp(lppr))

        dsetname = trim("")
        call HdfParallelReadParticle1D(filename,dsetname,rlpp(1:lppr),lpst,lpen,1,lppt,mpi_comm)

        lnum = 1
        do while(lnum.le.lppr)
            !!! Check if particle is within the computational domain
            call CheckIsParticleGlobal(rlpp(lnum),locl)
            if (locl) then
                !!!! Proceed to the next particle
                lnum = lnum + 1
            else
                !!!! Particle is out of bounds and should be not be read in.
                !!!! Copy the last active particle unless not already the last particle
                if (lnum.lt.lppr) rlpp(lnum) = rlpp(lppr)
                !!!! Decrement the number of active particle count in temporary lpp list, last particle is deactivated
                lppr = lppr - 1
            end if
        end do

        !! Allocate temporary arrays
        allocate(tlpp(lppr))
        allocate(trnk(lppr))
        allocate(mask(lppr))
        
        !! Copy over valid particles
        do lnum = 1,lppr
            tlpp(lnum) = rlpp(lnum)
        end do

        !! Deallocate the read buffer
        deallocate(rlpp)

        !! Set the decomp information
        xst2(mpi_rank+1) = yc(xstart(2))
        xen2(mpi_rank+1) = yc(xend(2)+1)
        xst3(mpi_rank+1) = zc(xstart(3))
        xen3(mpi_rank+1) = zc(xend(3)+1)

        !! Call MPI_Allgather in place to get decomp information from all pencils/processes
        call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,xst2,1,MPI_DOUBLE_PRECISION,mpi_comm,mpi_ierr)
        call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,xen2,1,MPI_DOUBLE_PRECISION,mpi_comm,mpi_ierr)
        call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,xst3,1,MPI_DOUBLE_PRECISION,mpi_comm,mpi_ierr)
        call MPI_ALLGATHER(MPI_IN_PLACE,1,MPI_DOUBLE_PRECISION,xen3,1,MPI_DOUBLE_PRECISION,mpi_comm,mpi_ierr)

        !! Find local pencil/process rank for each read particle
        do lnum = 1,lppr
            ! Initialize cell indices for c-grid (nodes)
            tlpp(lnum)%grc_idx(1) = int(nxm/2)
            tlpp(lnum)%grc_idx(2) = int(nym/2)
            tlpp(lnum)%grc_idx(3) = int(nzm/2)
            ! Find and update cell indices for c-grid (nodes)
            call GetLocationCellIndex(tlpp(lnum)%lpp_pos(1),xc(1:nx),1,nx,tlpp(lnum)%grc_idx(1))
            call GetLocationCellIndex(tlpp(lnum)%lpp_pos(2),yc(1:ny),1,ny,tlpp(lnum)%grc_idx(2))
            call GetLocationCellIndex(tlpp(lnum)%lpp_pos(3),zc(1:nz),1,nz,tlpp(lnum)%grc_idx(3))
            ! Initialize cell indices for m-grid (cell-center)
            tlpp(lnum)%grm_idx(1) = int(nxm/2)
            tlpp(lnum)%grm_idx(2) = int(nym/2)
            tlpp(lnum)%grm_idx(3) = int(nzm/2)
            ! Find and update cell indices for m-grid (cell-center)
            call GetLocationCellIndex(tlpp(lnum)%lpp_pos(1),xm(0:nx),0,nx,tlpp(lnum)%grm_idx(1))
            call GetLocationCellIndex(tlpp(lnum)%lpp_pos(2),ym(0:ny),0,ny,tlpp(lnum)%grm_idx(2))
            call GetLocationCellIndex(tlpp(lnum)%lpp_pos(3),zm(0:nz),0,nz,tlpp(lnum)%grm_idx(3))
            rnum = 1
            locl = .false.
            do while ((rnum.le.mpi_size).and.(.not.locl))
                locl        = .true.
                if (periodic(2)) then
                    locl    = locl.and.(tlpp(lnum)%lpp_pos(2).ge.xst2(rnum))
                else
                    locl    = locl.and.(tlpp(lnum)%lpp_pos(2).gt.xst2(rnum))
                end if
                locl        = locl.and.(tlpp(lnum)%lpp_pos(2).lt.xen2(rnum))
                if (periodic(3)) then
                    locl    = locl.and.(tlpp(lnum)%lpp_pos(3).ge.xst3(rnum))
                else
                    locl    = locl.and.(tlpp(lnum)%lpp_pos(3).gt.xst3(rnum))
                end if
                locl        = locl.and.(tlpp(lnum)%lpp_pos(3).lt.xen3(rnum))
                trnk(lnum)  = rnum
                rnum        = rnum + 1
            end do
        end do

        !! Rearrange particles and sort them according to locality of pencil/process in the order of increasing rank
        mask(:) = .true.
        do lnum = 1,lppr
            mloc = minloc(trnk,dim=1,mask=mask)
            if (mloc.gt.lnum) then
                !!!! Rearrange ranks
                tint       = trnk(lnum)
                trnk(lnum) = trnk(mloc)
                trnk(mloc) = tint
                !!!! Rearrange particles
                tmpp       = tlpp(lnum)
                tlpp(lnum) = tlpp(mloc)
                tlpp(mloc) = tmpp
            end if
            mask(lnum) = .false.
        end do

        !! Get the send counts and displacements
        scnt(:) = 0
        do lnum = 1,lppr
            scnt(trnk(lnum)) = scnt(trnk(lnum)) + 1
        end do
        sdsp(1) = 0
        do rnum = 2,mpi_size
            sdsp(rnum) = sdsp(rnum-1) + scnt(rnum-1)
        end do

        !! Call MPI_Alltoall to get the receive count of particles in the correct pencil/process
        call MPI_ALLTOALL(scnt,1,MPI_INTEGER,rcnt,1,MPI_INTEGER,mpi_comm,mpi_ierr)

        !! Get total number of particles received by current pencil/process
        lppr = sum(rcnt)

        !! Evaluate displacements for receiving
        rdsp(1) = 0
        do rnum = 2,mpi_size
            rdsp(rnum) = rdsp(rnum-1) + rcnt(rnum-1)
        end do

        !! Check if particle list needs to be extended (no need to copy existing uninitialised dummy particles)
        if (lppr.gt.size(lpp_list)) then
            if (allocated(lpp_list)) deallocate(lpp_list)
            allocate(lpp_list(lppr))
        end if

        !! Update active particles
        lpp_actv = lppr

        !! Call MPI_Alltoallv to exchange the particles and send them to the correct pencil/process where they locally belong
        call MPI_ALLTOALLV(tlpp,scnt,sdsp,mpi_pdat,lpp_list,rcnt,rdsp,mpi_pdat,mpi_comm,mpi_ierr)

        !! Dellocate temporary arrays
        deallocate(tlpp)
        deallocate(trnk)
        deallocate(mask)

        do lnum = 1,lpp_actv
            call CheckIsParticleLocal(lpp_list(lnum),locl)
            if (.not.locl) then
                write(*,*) 'Error reading in particles from continua.h5'
                call MpiAbort
            end if
        end do

    end if

    return
    
end subroutine ReadParticles