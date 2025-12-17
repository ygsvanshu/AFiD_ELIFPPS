!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: WriteFlowField.F90                             !
!    CONTAINS: subroutine WriteFlowField                  !
!                                                         !
!    PURPOSE: Write down the full flow snapshot for       !
!     restarting the simulation at a later date           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine WriteFlowField(end)
    
    use mpih
    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: vx,vy,vz,pr
    
    implicit none

    logical, intent(in) :: end
    character*200       :: filename,dsetname
    character*8         :: citime
    integer             :: st1,en1,st2,en2,st3,en3
    integer             :: int_time
    integer             :: perx,pery,perz

    if (periodic(1)) then
        perx = 1
    else
        perx = 0
    end if

    if (periodic(2)) then
        pery = 1
    else
        pery = 0
    end if

    if (periodic(3)) then
        perz = 1
    else
        perz = 0
    end if

    if (end) then
        filename = trim('continua.h5')
    else
        int_time = nint(time)
        write (citime,"(I8.8)") int_time
        filename = trim(trim('Results/continua_')//trim(citime)//'.h5')
    end if

    if (ismaster) then

        dsetname = trim('/nxm')
        call HdfSerialWriteIntScalar(filename,dsetname,nxm)
        dsetname = trim('/nym')
        call HdfSerialWriteIntScalar(filename,dsetname,nym)
        dsetname = trim('/nzm')
        call HdfSerialWriteIntScalar(filename,dsetname,nzm)

        dsetname = trim('/straxs')
        call HdfSerialWriteIntScalar(filename,dsetname,straxs)

        dsetname = trim('/perx')
        call HdfSerialWriteIntScalar(filename,dsetname,perx)
        dsetname = trim('/pery')
        call HdfSerialWriteIntScalar(filename,dsetname,pery)
        dsetname = trim('/perz')
        call HdfSerialWriteIntScalar(filename,dsetname,perz)

        dsetname = trim('/xlen')
        call HdfSerialWriteRealScalar(filename,dsetname,xlen)
        dsetname = trim('/ylen')
        call HdfSerialWriteRealScalar(filename,dsetname,ylen)
        dsetname = trim('/zlen')
        call HdfSerialWriteRealScalar(filename,dsetname,zlen)

        dsetname = trim('/time')
        call HdfSerialWriteRealScalar(filename,dsetname,time)

        dsetname = trim('/dt')
        call HdfSerialWriteRealScalar(filename,dsetname,dt)

        dsetname = trim('/rey')
        call HdfSerialWriteRealScalar(filename,dsetname,rey)

        dsetname = trim('/xc')
        call HdfSerialWriteReal1D(filename,dsetname,xc,1,nx)
        dsetname = trim('/xm')
        call HdfSerialWriteReal1D(filename,dsetname,xm,0,nx)
        dsetname = trim('/yc')
        call HdfSerialWriteReal1D(filename,dsetname,yc,1,ny)
        dsetname = trim('/ym')
        call HdfSerialWriteReal1D(filename,dsetname,ym,0,ny)
        dsetname = trim('/zc')
        call HdfSerialWriteReal1D(filename,dsetname,zc,1,nz)
        dsetname = trim('/zm')
        call HdfSerialWriteReal1D(filename,dsetname,zm,0,nz)

    end if

    !--------------------------------- X-VELOCITY ---------------------------------!
    
    dsetname = trim('vx')

    st1 = 1
    en1 = nx

    st2 = xstart(2)
    en2 = xend(2)

    if (xstart(2).eq.1) st2 = 0
    if (xend(2).eq.nym) en2 = ny

    st3 = xstart(3)
    en3 = xend(3)

    if (xstart(3).eq.1) st3 = 0
    if (xend(3).eq.nzm) en3 = nz

    call HdfParallelWriteReal3D(filename,dsetname,vx(st1:en1,st2:en2,st3:en3),st1,en1,st2,en2,st3,en3,1,nx,0,ny,0,nz,mpi_comm)
    
    !--------------------------------- Y-VELOCITY ---------------------------------!

    dsetname = trim('/vy')

    st1 = 0
    en1 = nx

    st2 = xstart(2)
    en2 = xend(2)

    if (xend(2).eq.nym) en2 = ny

    st3 = xstart(3)
    en3 = xend(3)

    if (xstart(3).eq.1) st3 = 0
    if (xend(3).eq.nzm) en3 = nz

    call HdfParallelWriteReal3D(filename,dsetname,vy(st1:en1,st2:en2,st3:en3),st1,en1,st2,en2,st3,en3,0,nx,1,ny,0,nz,mpi_comm)
    
    !--------------------------------- Z-VELOCITY ---------------------------------!

    dsetname = trim('/vz')

    st1 = 0
    en1 = nx

    st2 = xstart(2)
    en2 = xend(2)

    if (xstart(2).eq.1) st2 = 0
    if (xend(2).eq.nym) en2 = ny

    st3 = xstart(3)
    en3 = xend(3)

    if (xend(3).eq.nzm) en3 = nz

    call HdfParallelWriteReal3D(filename,dsetname,vz(st1:en1,st2:en2,st3:en3),st1,en1,st2,en2,st3,en3,0,nx,0,ny,1,nz,mpi_comm)
    
    !---------------------------------- PRESSURE ----------------------------------!

    dsetname = trim('/pr')

    st1 = 0
    en1 = nx

    st2 = xstart(2)
    en2 = xend(2)

    if (xstart(2).eq.1) st2 = 0
    if (xend(2).eq.nym) en2 = ny

    st3 = xstart(3)
    en3 = xend(3)

    if (xstart(3).eq.1) st3 = 0
    if (xend(3).eq.nzm) en3 = nz

    call HdfParallelWriteReal3D(filename,dsetname,pr(st1:en1,st2:en2,st3:en3),st1,en1,st2,en2,st3,en3,0,nx,0,ny,0,nz,mpi_comm)
    
end subroutine WriteFlowField

