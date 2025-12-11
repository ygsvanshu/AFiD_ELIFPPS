!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: WriteParticles.F90                             !
!    CONTAINS: subroutine WriteParticles                  !
!                                                         !
!    PURPOSE: To write the current snapshot of particles  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteParticles(end)

    use mpih
    use param
    use lagrangian_point_particle

    implicit none

    logical, intent(in)                 :: end
    logical                             :: exists
    character*200                       :: dsetname,filename
    character*8                         :: chartime
    integer                             :: ninttime
    integer                             :: int_dummy
    integer                             :: lppr,lppt
    integer                             :: lpst,lpen

    ! Get the number of active particles in pencil/process and the total number of active particles 
    lppr = lpp_actv
    call MPI_ALLREDUCE(lppr,lppt,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
    
    ! Get the start and end indices for the local pencil/process
    call MPI_SCAN(lppr,lpen,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
    lpst = lpen - lppr + 1

    !! Check if writing at the end or intermediate snapshots
    if (end) then
        !!! Writing at the end, no need to keep track of snapshot number
        filename = trim("continua.h5")
    else
        ninttime = nint(time)
        write (chartime,"(I8.8)") ninttime
        filename = trim(trim('Results/continua_')//trim(chartime)//'.h5')
    end if

    ! Get the total counts of spawned and exited particles
    !! Call MPI_Reduce to get the total counts
    call MPI_REDUCE(lpp_spwn,int_dummy,1,MPI_INTEGER,MPI_SUM,0,mpi_comm,mpi_ierr)
    tot_spwn = tot_spwn + int_dummy
    call MPI_REDUCE(lpp_exit,int_dummy,1,MPI_INTEGER,MPI_SUM,0,mpi_comm,mpi_ierr)
    tot_exit = tot_exit + int_dummy
    call MPI_REDUCE((lpp_actv-sub_exit),int_dummy,1,MPI_INTEGER,MPI_SUM,0,mpi_comm,mpi_ierr)
    tot_actv = int_dummy
    !! Reset the count of spawned and exited particles
    lpp_spwn = 0
    lpp_exit = 0

    !! Write data global time and timestep
    if (ismaster) then
        dsetname = "/lpp_dmod"
        call HdfSerialWriteIntScalar(filename,dsetname,lpp_dmod)
        dsetname = "/lpp_scor"
        call HdfSerialWriteIntScalar(filename,dsetname,merge(1,0,lpp_scor))
        dsetname = "/lpp_grav"
        call HdfSerialWriteReal1D(filename,dsetname,lpp_grav,1,3)
        dsetname = "/e2l_mult"
        call HdfSerialWriteRealScalar(filename,dsetname,e2l_mult)
        dsetname = "/l2e_mult"
        call HdfSerialWriteRealScalar(filename,dsetname,l2e_mult)
        dsetname = "/lpp_spwn"
        call HdfSerialWriteIntScalar(filename,dsetname,tot_spwn)
        dsetname = "/lpp_exit"
        call HdfSerialWriteIntScalar(filename,dsetname,tot_exit)
    end if

    !! Write the particle data
    dsetname = ""
    call HdfParallelWriteParticle1D(filename,dsetname,lpp_list(1:lppr),lpst,lpen,1,lppt,mpi_comm)
    
    return

end subroutine WriteParticles