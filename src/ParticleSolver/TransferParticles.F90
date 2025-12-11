!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: TransferParticles.F90                          !
!    CONTAINS: subroutine TransferParticles               !
!                                                         !
!    PURPOSE: To transfer particles between MPI processes !
!    when they crossover between pencils                  !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine TransferParticles

    use mpih
    use decomp_2d, only: xstart,xend
    use param, only: xc,yc,zc,xlen,ylen,zlen,periodic
    use lagrangian_point_particle

    implicit none

    integer :: np,nb,nl
    integer :: ierror
    integer :: tag_nbrm,tag_nbrp
    integer :: tag_nbcm,tag_nbcp
    integer :: requests(4)
    integer :: rqstatus(MPI_STATUS_SIZE,4)
    integer :: idx1,idx2
    integer :: lpp_size
    integer :: snd_prev,snd_next
    integer :: rcv_prev,rcv_next
    real    :: nlen(3)

    !-------- Row transfer --------!
    
    ! Initialize send and receive counts to zero
    snd_prev = 0
    snd_next = 0
    rcv_prev = 0
    rcv_next = 0
    
    ! For sending to previous row 
    if ((mpi_nbrm.ne.MPI_PROC_NULL).and.(lpp_actv.gt.0)) then 
        !! Check if bfm_send needs to be expanded
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(2).lt.yc(xstart(2))) snd_prev = snd_prev + 1
        end do
        call ExtendParticleSendBuffer(snd_prev,'p')
        !! Populate bfm_send buffer
        nb = 1
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(2).lt.yc(xstart(2))) then
                ! call CopyParticle(bfm_send(nb),lpp_list(np))
                bfm_send(nb) = lpp_list(np)
                nb = nb + 1
            end if
        end do
    end if
    
    ! For sending to next row 
    if ((mpi_nbrp.ne.MPI_PROC_NULL).and.(lpp_actv.gt.0)) then
        !! Check if bfp_send needs to be expanded
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(2).ge.yc(xend(2)+1)) snd_next = snd_next + 1
        end do
        call ExtendParticleSendBuffer(snd_next,'n')
        !! Populate bfp_send buffer
        nb = 1
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(2).ge.yc(xend(2)+1)) then
                ! call CopyParticle(bfp_send(nb),lpp_list(np))
                bfp_send(nb) = lpp_list(np)
                nb = nb + 1
            end if
        end do
    end if

    ! Compute tags
    tag_nbrm = mpi_pidx(1)
    if ((mpi_pidx(1)==mpi_dims(1)-1).and.periodic(2)) then
        tag_nbrp = 0
    else
        tag_nbrp = mpi_pidx(1) + 1
    end if
    
    ! Compute recv counts
    !! Receive from previous row
    call MPI_IRECV(rcv_prev,1,MPI_INTEGER,mpi_nbrm,tag_nbrm,mpi_comm,requests(1),ierror)
    !! Receive from next row
    call MPI_IRECV(rcv_next,1,MPI_INTEGER,mpi_nbrp,tag_nbrp,mpi_comm,requests(2),ierror)
    !! Send to previous row
    call MPI_ISSEND(snd_prev,1,MPI_INTEGER,mpi_nbrm,tag_nbrm,mpi_comm,requests(3),ierror)
    !! Send to next row
    call MPI_ISSEND(snd_next,1,MPI_INTEGER,mpi_nbrp,tag_nbrp,mpi_comm,requests(4),ierror)
    !! Wait for all requests to be processed
    call MPI_WAITALL(4, requests, rqstatus, ierror)

    ! Check if lpp_list needs to be expanded
    lpp_size = lpp_actv + rcv_prev + rcv_next
    call ExtendParticleListBuffer(lpp_size+1)

    ! Compute indices for particle list
    idx1 = lpp_actv + 1
    idx2 = lpp_actv + rcv_prev + 1

    ! print *, "Rank = ",mpi_rank," Send size = ",size(bfp_send)," Send count = ",snd_next," Recv size = ",size(lpp_list)," Recv count = ",rcv_prev

    ! Halo exchange particles in the row direction
    !! Receive from previous row
    call MPI_IRECV(lpp_list(idx1),rcv_prev,mpi_pdat,mpi_nbrm,tag_nbrm,mpi_comm,requests(1),ierror)
    !! Receive from next row
    call MPI_IRECV(lpp_list(idx2),rcv_next,mpi_pdat,mpi_nbrp,tag_nbrp,mpi_comm,requests(2),ierror)
    !! Send to previous row
    call MPI_ISSEND(bfm_send,snd_prev,mpi_pdat,mpi_nbrm,tag_nbrm,mpi_comm,requests(3),ierror)
    !! Send to next row
    call MPI_ISSEND(bfp_send,snd_next,mpi_pdat,mpi_nbrp,tag_nbrp,mpi_comm,requests(4),ierror)
    !! Wait for all requests to be processed
    call MPI_WAITALL(4, requests, rqstatus, ierror)

    ! Update active particles
    lpp_actv = lpp_actv + rcv_prev + rcv_next

    !------ Column transfer -------!

    ! Initialize send and receive counts to zero
    snd_prev = 0
    snd_next = 0
    rcv_prev = 0
    rcv_next = 0
    
    ! For sending to previous column 
    if ((mpi_nbcm.ne.MPI_PROC_NULL).and.(lpp_actv.gt.0)) then
        !! Check if bfm_send needs to be expanded
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(3).lt.zc(xstart(3))) snd_prev = snd_prev + 1
        end do
        call ExtendParticleSendBuffer(snd_prev,'p')
        !! Populate bfm_send buffer
        nb = 1
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(3).lt.zc(xstart(3))) then
                ! call CopyParticle(bfm_send(nb),lpp_list(np))
                bfm_send(nb) = lpp_list(np)
                nb = nb + 1
            end if
        end do
    end if
    
    ! For sending to next column 
    if ((mpi_nbcp.ne.MPI_PROC_NULL).and.(lpp_actv.gt.0)) then
        !! Check if bfp_send needs to be expanded
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(3).ge.zc(xend(3)+1)) snd_next = snd_next + 1
        end do
        call ExtendParticleSendBuffer(snd_next,'n')
        !! Populate bfp_send buffer
        nb = 1
        do np = 1,lpp_actv
            if (lpp_list(np)%lpp_pos(3).ge.zc(xend(3)+1)) then
                ! call CopyParticle(bfp_send(nb),lpp_list(np))
                bfp_send(nb) = lpp_list(np)
                nb = nb + 1
            end if
        end do
    end if

    ! Compute tags
    tag_nbcm = mpi_pidx(2)
    if ((mpi_pidx(2)==mpi_dims(2)-1).and.periodic(3)) then
        tag_nbcp = 0
    else
        tag_nbcp = mpi_pidx(2) + 1
    end if
    
    ! Compute recv counts
    !! Receive from previous column
    call MPI_IRECV(rcv_prev,1,MPI_INTEGER,mpi_nbcm,tag_nbcm,mpi_comm,requests(1),ierror)
    !! Receive from next column
    call MPI_IRECV(rcv_next,1,MPI_INTEGER,mpi_nbcp,tag_nbcp,mpi_comm,requests(2),ierror)
    !! Send to previous column
    call MPI_ISSEND(snd_prev,1,MPI_INTEGER,mpi_nbcm,tag_nbcm,mpi_comm,requests(3),ierror)
    !! Send to next column
    call MPI_ISSEND(snd_next,1,MPI_INTEGER,mpi_nbcp,tag_nbcp,mpi_comm,requests(4),ierror)
    !! Wait for all requests to be processed
    call MPI_WAITALL(4, requests, rqstatus, ierror)

    ! Check if lpp_list needs to be expanded
    lpp_size = lpp_actv + rcv_prev + rcv_next
    call ExtendParticleListBuffer(lpp_size+1)

    ! Compute indices for particle list
    idx1 = lpp_actv + 1
    idx2 = lpp_actv + rcv_prev + 1

    ! Halo exchange particles in the column direction
    !! Receive from previous column
    call MPI_IRECV(lpp_list(idx1),rcv_prev,mpi_pdat,mpi_nbcm,tag_nbcm,mpi_comm,requests(1),ierror)
    !! Receive from next column
    call MPI_IRECV(lpp_list(idx2),rcv_next,mpi_pdat,mpi_nbcp,tag_nbcp,mpi_comm,requests(2),ierror)
    !! Send to previous column
    call MPI_ISSEND(bfm_send,snd_prev,mpi_pdat,mpi_nbcm,tag_nbcm,mpi_comm,requests(3),ierror)
    !! Send to next column
    call MPI_ISSEND(bfp_send,snd_next,mpi_pdat,mpi_nbcp,tag_nbcp,mpi_comm,requests(4),ierror)
    !! Wait for all requests to be processed
    call MPI_WAITALL(4, requests, rqstatus, ierror)

    ! Update active particles
    lpp_actv = lpp_actv + rcv_prev + rcv_next

    !-------- Periodic boundary treatment --------!

    ! Note down the domain bounds
    nlen(1) = xlen
    nlen(2) = ylen
    nlen(3) = zlen

    ! Ensure periodic behaviour of particles
    do np = 1,lpp_actv
        do nl = 1,3
            if (periodic(nl)) lpp_list(np)%lpp_pos(nl) = modulo(lpp_list(np)%lpp_pos(nl),nlen(nl))
        end do
    end do

    return

end subroutine TransferParticles