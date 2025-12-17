program AFiD

    use mpih
    use decomp_2d
    use param
    use local_arrays, only: vx,vy,vz,pr
    use lagrangian_point_particle

    implicit none

    integer         :: errorcode, ierror
    real            :: instCFL,instPDT,dmax,davg
    real            :: ti(2),tin(3)
    real            :: ts
    integer         :: prow = 0, pcol = 0
    character(100)  :: arg
    logical         :: exists

    call ReadSolverInputs
    call ReadParticleInputs

    if (command_argument_count().eq.2) then
        call get_command_argument(1,arg)
        read (arg,'(i10)') prow
        call get_command_argument(2,arg)
        read (arg,'(i10)') pcol
    end if

    call decomp_2d_init(nxm,nym,nzm,prow,pcol,periodic)
        
    ts = MPI_WTIME()
    tin(1) = MPI_WTIME()

    call MpiBarrier

    call HdfStart

    if (nrank.eq.master) ismaster = .true.

    ! Get MPI communicator and ranks
    mpi_comm = DECOMP_2D_COMM_CART_X
    call MPI_COMM_RANK(DECOMP_2D_COMM_CART_X,mpi_rank,mpi_ierr)
    call MPI_COMM_SIZE(DECOMP_2D_COMM_CART_X,mpi_size,mpi_ierr)

    ! Get MPI cartesian topology data
    call MPI_CART_GET(DECOMP_2D_COMM_CART_X,2,mpi_dims,mpi_perc,mpi_pidx,mpi_ierr)

    ! Get MPI ranks for neighbours
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X,0,1,mpi_nbrm,mpi_nbrp,mpi_ierr)
    call MPI_CART_SHIFT(DECOMP_2D_COMM_CART_X,1,1,mpi_nbcm,mpi_nbcp,mpi_ierr)

    ! Get MPI subcommunicators for writing movie slices
    mpi_xcut = DECOMP_2D_COMM_CART_X
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.false.,.true./), mpi_ycut, ierror)
    call MPI_CART_SUB(DECOMP_2D_COMM_CART_X, (/.true.,.false./), mpi_zcut, ierror)

    ! Create results directory
    if (ismaster) call system('mkdir -p Results')

    call InitTimeMarchScheme
    call InitVariables
    call CreateGrid
    call InitPressureSolver
    if (save1d) call InitProfiles
    if (save2d) call InitMovie

    if (particle) call InitParticleSolver
    if (particle) call ReadParticleSources
    if (particle) call InitParticleExitFile
    if (particle.and.lpp_save) call InitParticleHistoryFile

    cvel = limitCFL*max(maxval(dxc),maxval(dyc),maxval(dzc))/dtmax

    call PrintCaseInfo
    if (particle) call PrintParticleCaseInfo

    if (nread) then
        if (ismaster) write (6,*) 'Reading initial condition from file'
        call ReadFlowField
        if (particle) call ReadParticles
    else
        if (ismaster) write (6,*) 'Creating initial condition'
        ntime = 0
        time = 0.0
        call CreateInitialConditions
        if (ismaster) write (6,*) 'Writing initial condition to file'
        call WriteFlowField(.false.)
        if (particle) call WriteParticles(.false.)
    end if

    call update_halo(vx,lvlhalo)
    call update_halo(vy,lvlhalo)
    call update_halo(vz,lvlhalo)
    call update_halo(pr,lvlhalo)

    call CheckDivergence(dmax,davg)
    call GlobalQuantities

    tin(2) = MPI_WTIME()

    call PrintStepInfo(tin(2)-tin(1),tin(2)-tin(1),dmax,davg,0.0)

    !  ********* starts the time dependent calculation ***
    errorcode = 0 !EP set errocode to 0 (OK)

    do ntime = 1,ntst

        ti(1) = MPI_WTIME()

        call CalcMaxCFL(instCFL)
        if (particle) call CalcMaxPDT(instPDT)

        if (vardt) then
            if (instCFL.lt.(limitCFL/dtmax)) then 
                dt = dtmax
            else
                dt = limitCFL/instCFL
            end if
            if (instPDT.gt.(1.0/dt)) dt = 1.0/instPDT
            if (dt.gt.dtmax) dt=dtmax
            if (dt.lt.dtmin) errorcode = 166
        else
            instCFL = instCFL*dt
            if (instCFL.gt.limitCFL) errorcode = 165
        end if

        call TimeMarcher

        ti(2) = MPI_WTIME()

        if (mod(ntime,nout).eq.0) then

            call GlobalQuantities
            if ((vmax(1).gt.limitVel).or.(vmax(2).gt.limitVel).or.(vmax(3).gt.limitVel)) errorcode = 167

            call CalcMaxCFL(instCFL)
            if (.not.vardt) instCFL = instCFL*dt

            call CheckDivergence(dmax,davg)
            if (abs(dmax).gt.resid) errorcode = 169

            if (particle) then
                call GlobalParticleStatistics
                if ((lpp_vmax(1).gt.limitVel).or.(lpp_vmax(2).gt.limitVel).or.(lpp_vmax(3).gt.limitVel)) errorcode = 168
            end if

            call PrintStepInfo(ti(2)-ti(1),ti(2)-tin(2),dmax,davg,instCFL)

        end if

        if ((save1d).and.(time.ge.tsta1d).and.(mod(time,freq1d).lt.dt)) call WriteStatProfiles
        if ((save2d).and.(time.ge.tsta2d).and.(mod(time,freq2d).lt.dt)) call WriteMovieSlices
        if ((save3d).and.(time.ge.tsta3d).and.(mod(time,freq3d).lt.dt)) call WriteFlowField(.false.)

        if (particle) then
            if ((save3d).and.(time.ge.tsta3d).and.(mod(time,freq3d).lt.dt)) call WriteParticles(.false.)
            if ((lpp_save).and.(time.ge.lpp_ssta).and.(mod(time,lpp_sfrq).lt.dt)) call WriteParticleHistory
            if ((pex_save).and.(time.ge.pex_ssta).and.(mod(time,pex_sfrq).lt.dt)) call WriteExitedParticles
        end if

        inquire(file=abortfile,exist=exists)
        if (exists) errorcode = 222

        if (time.gt.tmax) errorcode = 333

        if ((ti(2) - tin(1)).gt.walltimemax) errorcode = 334

        if (ntime.eq.ntst) errorcode = 555

        call MpiBcastInt(errorcode)

        ! Conditional exits
        if (errorcode .ne. 0) then
            
            ! dt too small
            if (errorcode.eq.166) call QuitRoutine(tin,.false.,errorcode)
            ! cfl too high
            if (errorcode.eq.165) call QuitRoutine(tin,.false.,errorcode)
            ! velocities diverged
            if (errorcode.eq.167) call QuitRoutine(tin,.false.,errorcode)
            ! particle velocities diverged
            if (errorcode.eq.168) call QuitRoutine(tin,.false.,errorcode)
            ! mass not conserved
            if (errorcode.eq.169) call QuitRoutine(tin,.false.,errorcode)
            ! Abort file detected, no error; normal quit
            if (errorcode.eq.222) call QuitRoutine(tin,.true.,errorcode)
            ! Physical time exceeded tmax, no error; normal quit
            if (errorcode.eq.333) call QuitRoutine(tin,.true.,errorcode)
            ! walltime exceeded walltimemax, no error; normal quit
            if (errorcode.eq.334) call QuitRoutine(tin,.true.,errorcode)
            ! maximum number of timesteps reached, no error; normal quit
            if (errorcode.eq.555) call QuitRoutine(tin,.true.,errorcode)
            ! already finalized
            errorcode = 100

            exit

        end if

    end do

    call QuitRoutine(tin,.true.,errorcode)

end program AFiD