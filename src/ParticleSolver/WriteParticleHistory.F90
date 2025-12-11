!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: WriteParticleHistory.F90                       !
!    CONTAINS: subroutine InitParticleHistoryFile         !
!              subroutine WriteParticleHistory            !
!                                                         !
!    PURPOSE: To write particle snapshot/movie/history    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitParticleHistoryFile

    use mpih
    use param, only: ismaster
    use lagrangian_point_particle

    implicit none

    logical         :: exists
    integer         :: srun,ssnp
    integer         :: psrr,pstt,psst,psen
    character*200   :: filename,dsetname
    character*4     :: charsrun

    ! Check if particle history file already exists
    filename = trim("Results/particle_history.h5")
    inquire(file=filename,exist=exists)
    call MPI_ALLREDUCE(MPI_IN_PLACE,exists,1,MPI_LOGICAL,MPI_LOR,mpi_comm,mpi_ierr)
    if (.not.exists) then
        !! Particle history doesn't exist, initialize runs and snapshots to zero
        srun = 0
        ssnp = 0
    else
        !! Particle history already exists, so read in the number of previous runs and steps
        if (ismaster) then
            dsetname = trim("/runs")
            call HdfSerialReadIntScalar(filename,dsetname,srun)
            dsetname = trim("/steps")
            call HdfSerialReadIntScalar(filename,dsetname,ssnp)
        end if
    end if
    !! MPI Barrier and MPI Bcast to broadcast the number of previous runs
    call MpiBarrier
    call MpiBcastInt(srun)
    call MpiBcastInt(ssnp)
    srun = srun + 1
    ssnp = ssnp + 1
    write(charsrun,"(I4.4)") srun
    !! Write a bunch of info about the run
    if (ismaster) then
        dsetname = trim("/runs")
        call HdfSerialWriteIntScalar(filename,dsetname,srun)
        dsetname = trim("/info/"//charsrun//"/stt_indx")
        call HdfSerialWriteIntScalar(filename,dsetname,ssnp)
        dsetname = trim("/info/"//charsrun//"/lpp_dmod")
        call HdfSerialWriteIntScalar(filename,dsetname,lpp_dmod)
        dsetname = trim("/info/"//charsrun//"/lpp_scor")
        call HdfSerialWriteIntScalar(filename,dsetname,merge(1,0,lpp_scor))
        dsetname = trim("/info/"//charsrun//"/lpp_grav")
        call HdfSerialWriteReal1D(filename,dsetname,lpp_grav,1,3)
        dsetname = trim("/info/"//charsrun//"/e2l_mult")
        call HdfSerialWriteRealScalar(filename,dsetname,e2l_mult)
        dsetname = trim("/info/"//charsrun//"/l2e_mult")
        call HdfSerialWriteRealScalar(filename,dsetname,l2e_mult)
    end if

    ! Set the start number of new snapshots 
    lpp_snap = ssnp

    ! Get the number of particle sources in pencil/process and the total number of active particles 
    psrr = src_size
    pstt = src_ntot
    
    ! Get the start and end indices for the local pencil/process
    call MPI_SCAN(psrr,psen,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
    psst = psen - psrr + 1

    ! Call an MPI Barrier just in case. Since it's just once at initialization, no performance hit
    call MpiBarrier
    dsetname = trim("/info/"//charsrun//"/lpp_srcs")
    call HdfParallelWriteSource1D(filename,dsetname,src_list,psst,psen,1,pstt,mpi_comm)

end subroutine InitParticleHistoryFile

subroutine WriteParticleHistory

    use mpih
    use param
    use lagrangian_point_particle

    implicit none

    character*200                       :: dsetname,filename
    character*8                         :: charsnap
    integer                             :: lppr,lppt
    integer                             :: lpst,lpen

    ! Get the number of active particles in pencil/process and the total number of active particles 
    lppr = lpp_actv
    call MPI_ALLREDUCE(lppr,lppt,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
    
    ! Get the start and end indices for the local pencil/process
    call MPI_SCAN(lppr,lpen,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
    lpst = lpen - lppr + 1

    ! Writing to particle history
    filename = trim("Results/particle_history.h5")
    !! Figure out writing the snapshot info
    dsetname = trim("/steps")
    !! File must exist from initialization to no need to check, just write the snapshot info
    if (ismaster) then
        call HdfSerialWriteIntScalar(filename,dsetname,lpp_snap)
    end if
    !! Create the dsetname for history file
    write(charsnap,"(I8.8)") lpp_snap

    ! Increment the particle history snapshot count
    lpp_snap = lpp_snap + 1

    !! Write data global time and timestep
    if (ismaster) then
        dsetname = trim("/data/"//charsnap//"/time")
        call HdfSerialWriteRealScalar(filename,dsetname,time)
        dsetname = trim("/data/"//charsnap//"/dt")
        call HdfSerialWriteRealScalar(filename,dsetname,dt)
    end if

    !! Write the particle data
    dsetname = trim("/data/"//charsnap)
    call HdfParallelWriteParticle1D(filename,dsetname,lpp_list(1:lppr),lpst,lpen,1,lppt,mpi_comm)
    
    return

end subroutine WriteParticleHistory