!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: WriteExitedParticles.F90                       !
!    CONTAINS: subroutine InitParticleExitFile            !
!    CONTAINS: subroutine WriteExitedParticles            !
!              subroutine PrintExitedParticles            !
!                                                         !
!    PURPOSE: To write the details of exited particles to !
!    output file.                                         ! 
!                                                         !
!    It is expected that the routine PrintExitedParticles !
!    will only be used for debugging purposes for a small !
!    number of particles as it writes out exit events to  !
!    a text file. For large number of particles, use the  !
!    subroutine WriteExitedParticles which writes binary  !
!    HDF5 files.                                          !        
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitParticleExitFile

    use mpih
    use param, only: ismaster
    use lagrangian_point_particle
    
    implicit none

    logical         :: exists
    integer         :: srun,ssnp
    integer         :: psrr,pstt,psst,psen
    character*200   :: filename,dsetname
    character*4     :: charsrun

    ! Initialize the particle exit file if it doesn't exist
    filename = trim("Results/particle_exit.h5")
    inquire(file=filename,exist=exists)
    call MPI_ALLREDUCE(MPI_IN_PLACE,exists,1,MPI_LOGICAL,MPI_LOR,mpi_comm,mpi_ierr)
    if (.not.exists) then
        !! Exit file doesn't exist, create the exit file
        dsetname = trim("/data")
        call HdfParallelCreateParticleExit1D(filename,dsetname,mpi_comm)
        !! Set this to be the first run
        srun = 0
        ssnp = 0
    else
        !! Exit file already exists, so read in the number of previous runs and number of exits
        if (ismaster) then
            dsetname = trim("/runs")
            call HdfSerialReadIntScalar(filename,dsetname,srun)
            dsetname = trim("/data/pex_num")
            call HdfSerialReadIntScalar(filename,dsetname,ssnp)
        end if
    end if
    ! MPI Barrier and MPI Bcast to broadcast the number of previous runs, and number of exits
    call MpiBarrier
    call MpiBcastInt(srun)
    call MpiBcastInt(ssnp)
    srun = srun + 1
    ssnp = ssnp + 1
    write(charsrun,"(I4.4)") srun
    ! Write a bunch of info about the run
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

end subroutine InitParticleExitFile

subroutine WriteExitedParticles

    use mpih
    use lagrangian_point_particle

    implicit none

    character(len=200)                      :: filename,dsetname
    integer                                 :: ernk,esta,eend,etot
    integer                                 :: nlpp,elpp
    logical                                 :: gbound

    ! Get the number of exit events in pencil/process and the total number of active particles 
    ernk = pex_actv
    call MPI_ALLREDUCE(ernk,etot,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)

    ! Get the start and end indices of exit events for the local pencil/process
    call MPI_SCAN(ernk,eend,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
    esta = eend - ernk + 1

    ! Only participate in the write process if there are any exit events (globally) to write
    if (etot.gt.0) then

        !! Set the file name and dataset path
        filename = trim("Results/particle_exit.h5")
        dsetname = trim("/data/")

        !! Append the exited particles to the HDF5 file
        call HdfParallelAppendParticleExit1D(filename,dsetname,pex_list(1:pex_actv),esta,eend,1,etot,mpi_comm)

    end if

    ! Reset the pending exit events 
    pex_actv = 0
    
end subroutine WriteExitedParticles