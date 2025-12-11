!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: QuitRoutine.F90                                !
!    CONTAINS: subroutine QuitRoutine, NotifyError        !
!                                                         !
!    PURPOSE: Routines to exit the program and write the  !
!     data if necessary                                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine QuitRoutine(tin,normalexit,errorcode)

    use mpih
    use hdf5
    use decomp_2d, only: decomp_2d_finalize
    use param
    use lagrangian_point_particle

    implicit none

    logical,intent(in)  :: normalexit
    integer             :: errorcode
    real                :: tin(3)

    if (ismaster) write (6,*) ''

    if (errorcode .ne. 100) then ! skip if already finalized

        tin(3) = MPI_WTIME()
        if (ismaster) then
            call NotifyError(errorcode)
        end if

        if (normalexit) then
            if (ismaster) write (6,'(a,f10.2,a)') 'Total Iteration Time = ',tin(3) - tin(2),' sec.'
            call WriteFlowField(.true.)
            call WriteParticles(.true.)
            call WriteExitedParticles
        else
            call MPI_ABORT(MPI_COMM_WORLD,1)
        end if

        call FinalizeParticleSolver
        call DeallocateVariables
        call DeallocateFFTArrays
        call HdfClose
        call decomp_2d_finalize

    end if

end subroutine QuitRoutine

subroutine NotifyError(errorcode)

    use param
    
    implicit none

    integer,intent(in) :: errorcode

    if (errorcode .eq. 166) then
        write (*,*) "dt too small, DT = ", dt
    else if (errorcode .eq. 165) then
        write (*,*) "cfl too large"
    else if (errorcode .eq. 167) then
        write (*,*) "velocities diverged"
    else if (errorcode .eq. 168) then
        write (*,*) "particle velocities diverged"
    else if (errorcode .eq. 169) then
        write (*,*) "too large local residue for mass conservation"
        write (*,*) "probably the matrix in SolvePressureCorrection becomes singular"
        write (*,*) "try changing boundary conditions or grid stretching"
        call LocateLargeDivergence
    else if (errorcode .eq. 222) then
        write (*,*) "abort file detected, manual abort triggered"
        write (*,*) "continuation updated"
    else if (errorcode .eq. 333) then
        write (*,*) "time greater than tmax"
        write (*,*) "continuation updated"
    else if (errorcode .eq. 334) then
        write (*,*) "walltime greater than walltimemax"
        write (*,*) "continuation updated"
    else
        write (*,*) "maximum number of timesteps reached"
        write (*,*) "continuation updated"
    end if

    return

end subroutine NotifyError