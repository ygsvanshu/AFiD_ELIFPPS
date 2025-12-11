!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: PrintParticleCaseInfo.F90                      !
!    CONTAINS: subroutine PrintParticleCaseInfo           !
!                                                         !
!    PURPOSE: PrintParticleCaseInfo prints information    !
!             about global parameters of particle solver  !
!             at the beginning of the simulation          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PrintParticleCaseInfo

    use param, only: ismaster
    use lagrangian_point_particle

    implicit none

    character*32    :: dopt,sopt

    if (lpp_dmod.eq.STOKES) dopt = trim("Stokes")
    if (lpp_dmod.eq.SCHNAU) dopt = trim("Schiller-Naumann")

    if (lpp_scor) then
        sopt = trim("Enabled")
    else
        sopt = trim("Disabled")
    end if

    if (ismaster) then

        write (6,'(A77)') '============================================================================='
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '                     Lagrangian point particle tracking                      '
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** DRAG MODEL ======================================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,A37)') '     Drag coefficient model           = ',dopt
        write (6,'(A40,A37)') '     Slip correction                  = ',sopt
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** PARTICLE COUPLING ================================================= ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,3F8.3)') '     Eulerian -> Lagrangian multiplier = ',e2l_mult
        write (6,'(A40,3F8.3)') '     Lagrangian -> Eulerian multiplier = ',l2e_mult
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** PARTICLE TIMESTEP LIMITERS ======================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,ES10.3)') '     Spawn tolerance                   = ',lpp_stol
        write (6,'(A40,3F8.3)') '     Particle CFL analogue             = ',lpp_clim
        write (6,'(A40,3F8.3)') '     Particle response time ratio      = ',lpp_tlim
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** PARTICLE SOURCES ================================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,I8)') '     Number of particle sources read   = ',src_ntot
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '============================================================================='
        write (6,'(A77)') '                                                                             '

    end if

end subroutine PrintParticleCaseInfo