!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: UpdateHaloForces.F90                           !
!    CONTAINS: subroutine UpdateHaloForces                !
!                                                         !
!    PURPOSE: To transfer forces applied on halo cells    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine UpdateHaloForces(qval)

    use mpih
    use decomp_2d
    use param, only: periodic,lvlhalo
    use lagrangian_point_particle
    
    implicit none

    real, intent(inout) :: qval(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo)

    integer             :: ilvl
    integer             :: snd_hsta,snd_hend
    integer             :: rcv_hsta,rcv_hend
    integer             :: hcnt
    integer             :: requests(4)
    integer             :: rqstatus(MPI_STATUS_SIZE,4)
    integer             :: tag_nbrm,tag_nbrp
    integer             :: tag_nbcm,tag_nbcp

    !Calculate tags
    tag_nbrm = mpi_pidx(1)
    if ((mpi_pidx(1)==mpi_dims(1)-1).and.periodic(2)) then
        tag_nbrp = 0
    else
        tag_nbrp = mpi_pidx(1) + 1
    end if

    ! Start and end indices for sending
    snd_hsta = xstart(2)-lvlhalo-1
    snd_hend = xend(2)
    ! Start and end indices for receiving
    rcv_hsta = xstart(2)-1
    rcv_hend = xend(2)-lvlhalo

    if(tag_nbrp.ne.tag_nbrm) then
        ! Compute count
        hcnt = ((xend(1)+lvlhalo) - (xstart(1)-lvlhalo) + 1)*(lvlhalo)*((xend(3)+lvlhalo) - (xstart(3)-lvlhalo) + 1)
        ! Copy data to send buffer
        do ilvl=1,lvlhalo
            snd_nbrm(:,ilvl,:) = qval(:,snd_hsta+ilvl,:)
            snd_nbrp(:,ilvl,:) = qval(:,snd_hend+ilvl,:)
            rcv_nbrm(:,ilvl,:) = 0.0
            rcv_nbrp(:,ilvl,:) = 0.0
        end do
        ! Receive from previous row
        call MPI_IRECV(rcv_nbrm,hcnt,MPI_DOUBLE_PRECISION,mpi_nbrm,tag_nbrm,mpi_comm,requests(1),mpi_ierr)
        ! Receive from next row
        call MPI_IRECV(rcv_nbrp,hcnt,MPI_DOUBLE_PRECISION,mpi_nbrp,tag_nbrp,mpi_comm,requests(2),mpi_ierr)
        ! Send to previous row
        call MPI_ISSEND(snd_nbrm,hcnt,MPI_DOUBLE_PRECISION,mpi_nbrm,tag_nbrm,mpi_comm,requests(3),mpi_ierr)
        ! Send to next row
        call MPI_ISSEND(snd_nbrp,hcnt,MPI_DOUBLE_PRECISION,mpi_nbrp,tag_nbrp,mpi_comm,requests(4),mpi_ierr)
        ! Wait for all requests
        call MPI_WAITALL(4,requests,rqstatus,mpi_ierr)
        ! Add data from receive buffer
        do ilvl=1,lvlhalo
            qval(:,rcv_hsta+ilvl,:) = qval(:,rcv_hsta+ilvl,:) + rcv_nbrm(:,ilvl,:)
            qval(:,rcv_hend+ilvl,:) = qval(:,rcv_hend+ilvl,:) + rcv_nbrp(:,ilvl,:) 
        end do
    else
        ! If periodic and single pencil in 2nd dimension, copy halo cells for ease of indexing
        do ilvl=1,lvlhalo
            qval(:,rcv_hsta+ilvl,:) = qval(:,rcv_hsta+ilvl,:) + qval(:,snd_hend+ilvl,:)
            qval(:,rcv_hend+ilvl,:) = qval(:,rcv_hend+ilvl,:) + qval(:,snd_hsta+ilvl,:)
        end do
    endif

    !Calculate tags
    tag_nbcm = mpi_pidx(2)
    if ((mpi_pidx(2)==mpi_dims(2)-1).and.periodic(3)) then
        tag_nbcp = 0
    else
        tag_nbcp = mpi_pidx(2) + 1
    end if

    ! Start and end indices for sending
    snd_hsta = xstart(3)-lvlhalo-1
    snd_hend = xend(3)
    ! Start and end indices for receiving
    rcv_hsta = xstart(3)-1
    rcv_hend = xend(3)-lvlhalo

    if(tag_nbcp.ne.tag_nbcm) then
        ! Compute count
        hcnt = ((xend(1)+lvlhalo) - (xstart(1)-lvlhalo) + 1)*((xend(2)+lvlhalo) - (xstart(2)-lvlhalo) + 1)*(lvlhalo)
        ! Copy data to send buffer
        do ilvl=1,lvlhalo
            snd_nbcm(:,:,ilvl) = qval(:,:,snd_hsta+ilvl)
            snd_nbcp(:,:,ilvl) = qval(:,:,snd_hend+ilvl)
            rcv_nbcm(:,:,ilvl) = 0.0
            rcv_nbcp(:,:,ilvl) = 0.0
        end do
        ! Receive from previous row
        call MPI_IRECV(rcv_nbcm,hcnt,MPI_DOUBLE_PRECISION,mpi_nbcm,tag_nbcm,mpi_comm,requests(1),mpi_ierr)
        ! Receive from next row
        call MPI_IRECV(rcv_nbcp,hcnt,MPI_DOUBLE_PRECISION,mpi_nbcp,tag_nbcp,mpi_comm,requests(2),mpi_ierr)
        ! Send to previous row
        call MPI_ISSEND(snd_nbcm,hcnt,MPI_DOUBLE_PRECISION,mpi_nbcm,tag_nbcm,mpi_comm,requests(3),mpi_ierr)
        ! Send to next row
        call MPI_ISSEND(snd_nbcp,hcnt,MPI_DOUBLE_PRECISION,mpi_nbcp,tag_nbcp,mpi_comm,requests(4),mpi_ierr)
        ! Wait for all requests
        call MPI_WAITALL(4,requests,rqstatus,mpi_ierr)
        ! Add data from receive buffer
        do ilvl=1,lvlhalo
            qval(:,:,rcv_hsta+ilvl) = qval(:,:,rcv_hsta+ilvl) + rcv_nbcm(:,:,ilvl)
            qval(:,:,rcv_hend+ilvl) = qval(:,:,rcv_hend+ilvl) + rcv_nbcp(:,:,ilvl) 
        end do
    else
        ! If periodic and single pencil in 2nd dimension, copy halo cells for ease of indexing
        do ilvl=1,lvlhalo
            qval(:,:,rcv_hsta+ilvl) = qval(:,:,rcv_hsta+ilvl) + qval(:,:,snd_hend+ilvl)
            qval(:,:,rcv_hend+ilvl) = qval(:,:,rcv_hend+ilvl) + qval(:,:,snd_hsta+ilvl)
        end do
    endif

    return

end subroutine UpdateHaloForces