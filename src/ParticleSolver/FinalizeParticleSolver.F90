!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: FinalizeParticleSolver.F90                     ! 
!    CONTAINS: subroutine FinalizeParticleSolver          !
!                                                         !
!    PURPOSE: Deallocate global arrays and free MPI type  !
!    of the derived particle datatype                     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine FinalizeParticleSolver !{DONE}!

    use mpih
    use lagrangian_point_particle

    implicit none

    if (allocated(lpp_list)) deallocate(lpp_list)
    if (allocated(src_list)) deallocate(src_list)
    if (allocated(pex_list)) deallocate(pex_list)

    if (allocated(bfm_send)) deallocate(bfm_send)
    if (allocated(bfp_send)) deallocate(bfp_send)

    if (allocated(snd_nbrm)) deallocate(snd_nbrm)
    if (allocated(snd_nbrp)) deallocate(snd_nbrp)
    if (allocated(snd_nbcm)) deallocate(snd_nbcm)
    if (allocated(snd_nbcp)) deallocate(snd_nbcp)
    
    if (allocated(rcv_nbrm)) deallocate(rcv_nbrm)
    if (allocated(rcv_nbrp)) deallocate(rcv_nbrp)
    if (allocated(rcv_nbcm)) deallocate(rcv_nbcm)
    if (allocated(rcv_nbcp)) deallocate(rcv_nbcp)

    if (allocated(lpp_bdfx)) deallocate(lpp_bdfx)
    if (allocated(lpp_bdfy)) deallocate(lpp_bdfy)
    if (allocated(lpp_bdfz)) deallocate(lpp_bdfz)

    if (allocated(lpp_d2vx)) deallocate(lpp_d2vx)
    if (allocated(lpp_d2vy)) deallocate(lpp_d2vy)
    if (allocated(lpp_d2vz)) deallocate(lpp_d2vz)

    if (allocated(dat_rdia)) deallocate(dat_rdia)
    if (allocated(dat_aspc)) deallocate(dat_aspc)
    if (allocated(dat_reyn)) deallocate(dat_reyn)
    if (allocated(dat_coef)) deallocate(dat_coef)
    
    call MPI_TYPE_FREE(mpi_pdat,mpi_ierr)

end subroutine FinalizeParticleSolver