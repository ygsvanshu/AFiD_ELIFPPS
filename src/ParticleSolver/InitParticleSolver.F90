!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: InitParticleSolver.F90                         !
!    CONTAINS: subroutine InitParticleSolver              !
!                                                         !
!    PURPOSE: Subroutine to initialize MPI information,   !
!    MPI types for derived datatypes, and global arrays   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitParticleSolver

    use mpih
    use param
    use decomp_2d
    use lagrangian_point_particle

    implicit none

    integer                                                     :: hdf_error
    integer                                                     :: bcnt,bnum
    integer, allocatable, dimension(:)                          :: blen,btyp
    integer(KIND=MPI_ADDRESS_KIND)                              :: bref,bext,bblb
    integer(KIND=MPI_ADDRESS_KIND), allocatable, dimension(:)   :: bdsp
    type(particle_data)                                         :: tlpp
    character*200                                               :: filename

    ! Check lvlhalo is at least 1
    if (lvlhalo.lt.1) then
        if (ismaster) write (6,*) 'Error: number of halo layer cells must be greater than zero'
        call MpiAbort
    end if

    ! Make derived MPI datatype for particles
    !! Set block count
    bcnt = 11
    !! Allocate arrays for datatype definition
    allocate(blen(bcnt))
    allocate(bdsp(bcnt))
    allocate(btyp(bcnt))
    !! Set block lengths
    blen( 1) = 1
    blen( 2) = 3
    blen( 3) = 3
    blen( 4) = 1
    blen( 5) = 1
    blen( 6) = 1
    blen( 7) = 1
    blen( 8) = 3
    blen( 9) = 3
    blen(10) = 3
    blen(11) = 3
    !! Get block addresses from the test particle
    call MPI_GET_ADDRESS(tlpp%src_idx,bdsp( 1),mpi_ierr)
    call MPI_GET_ADDRESS(tlpp%grc_idx,bdsp( 2),mpi_ierr) 
    call MPI_GET_ADDRESS(tlpp%grm_idx,bdsp( 3),mpi_ierr) 
    call MPI_GET_ADDRESS(tlpp%lpp_lft,bdsp( 4),mpi_ierr)
    call MPI_GET_ADDRESS(tlpp%lpp_dia,bdsp( 5),mpi_ierr)
    call MPI_GET_ADDRESS(tlpp%lpp_den,bdsp( 6),mpi_ierr)
    call MPI_GET_ADDRESS(tlpp%lpp_rey,bdsp( 7),mpi_ierr) 
    call MPI_GET_ADDRESS(tlpp%lpp_pos,bdsp( 8),mpi_ierr) 
    call MPI_GET_ADDRESS(tlpp%lpp_vel,bdsp( 9),mpi_ierr) 
    call MPI_GET_ADDRESS(tlpp%acc_old,bdsp(10),mpi_ierr) 
    call MPI_GET_ADDRESS(tlpp%acc_now,bdsp(11),mpi_ierr)
    !! Get block displacements from block address reference
    call MPI_GET_ADDRESS(tlpp,bref,mpi_ierr)
    do bnum = 1,bcnt
        bdsp(bnum) = bdsp(bnum) - bref
    end do
    !! Set block datatypes
    btyp( 1) = MPI_INTEGER
    btyp( 2) = MPI_INTEGER
    btyp( 3) = MPI_INTEGER
    btyp( 4) = MPI_DOUBLE_PRECISION
    btyp( 5) = MPI_DOUBLE_PRECISION
    btyp( 6) = MPI_DOUBLE_PRECISION
    btyp( 7) = MPI_DOUBLE_PRECISION
    btyp( 8) = MPI_DOUBLE_PRECISION
    btyp( 9) = MPI_DOUBLE_PRECISION
    btyp(10) = MPI_DOUBLE_PRECISION
    btyp(11) = MPI_DOUBLE_PRECISION
    !! Create the MPI data structure and commit it
    call MPI_TYPE_CREATE_STRUCT(bcnt,blen,bdsp,btyp,mpi_pdat,mpi_ierr)
    call MPI_TYPE_COMMIT(mpi_pdat,mpi_ierr)
    !! Deallocate arrays
    deallocate(blen)
    deallocate(bdsp)
    deallocate(btyp)
    
    ! Allocate body force arrays
    if (.not.allocated(lpp_bdfx)) allocate(lpp_bdfx(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
    if (.not.allocated(lpp_bdfy)) allocate(lpp_bdfy(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
    if (.not.allocated(lpp_bdfz)) allocate(lpp_bdfz(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))

    ! Allocate arrays for curvature terms
    if (.not.allocated(lpp_d2vx)) allocate(lpp_d2vx(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
    if (.not.allocated(lpp_d2vy)) allocate(lpp_d2vy(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
    if (.not.allocated(lpp_d2vz)) allocate(lpp_d2vz(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))

    ! Allocate halo force buffers
    if (.not.allocated(snd_nbrm)) allocate(snd_nbrm(xstart(1)-lvlhalo:xend(1)+lvlhalo,1:lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
    if (.not.allocated(snd_nbrp)) allocate(snd_nbrp(xstart(1)-lvlhalo:xend(1)+lvlhalo,1:lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))

    if (.not.allocated(rcv_nbrm)) allocate(rcv_nbrm(xstart(1)-lvlhalo:xend(1)+lvlhalo,1:lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))
    if (.not.allocated(rcv_nbrp)) allocate(rcv_nbrp(xstart(1)-lvlhalo:xend(1)+lvlhalo,1:lvlhalo,xstart(3)-lvlhalo:xend(3)+lvlhalo))

    if (.not.allocated(snd_nbcm)) allocate(snd_nbcm(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,1:lvlhalo))
    if (.not.allocated(snd_nbcp)) allocate(snd_nbcp(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,1:lvlhalo))

    if (.not.allocated(rcv_nbcm)) allocate(rcv_nbcm(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,1:lvlhalo))
    if (.not.allocated(rcv_nbcp)) allocate(rcv_nbcp(xstart(1)-lvlhalo:xend(1)+lvlhalo,xstart(2)-lvlhalo:xend(2)+lvlhalo,1:lvlhalo))

    ! Allocate seed for particle list, particle buffers, and particle exit events
    if (.not.allocated(lpp_list)) allocate(lpp_list(1))
    if (.not.allocated(pex_list)) allocate(pex_list(1))
    
    if (.not.allocated(bfm_send)) allocate(bfm_send(1))
    if (.not.allocated(bfp_send)) allocate(bfp_send(1))

    ! Initialize counts 
    src_size = 0
    src_ntot = 0

    lpp_actv = 0
    pex_actv = 0

    lpp_spwn = 0
    lpp_exit = 0

    tot_spwn = 0
    tot_exit = 0

    lpp_snap = 0

    ! Initialize the slip correction model with data if slip correction is on
    if (lpp_scor) call InitSlipCorrectionData

    return

end subroutine InitParticleSolver