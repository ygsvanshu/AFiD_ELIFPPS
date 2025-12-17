!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: MpiAuxRoutines.F90                             !
!    CONTAINS: subroutines MPI*                           !
!                                                         ! 
!    PURPOSE: Wrappers for MPI Routines                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine MpiBcastInt(n)

    use mpih
    implicit none
    integer, intent(in) :: n
    call MPI_BCAST(n,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiBcastInt

subroutine MpiBcastReal(n)

    use mpih
    implicit none
    real, intent(in) :: n
    call MPI_BCAST(n,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiBcastReal

subroutine MpiBarrier

    use mpih
    implicit none
    call MPI_BARRIER(MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiBarrier

subroutine MpiMinRealScalar(var,res)

    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_REDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_MIN,0,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiMinRealScalar

subroutine MpiMaxRealScalar(var,res)

    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_REDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiMaxRealScalar

subroutine MpiSumRealScalar(var,res)

    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_REDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiSumRealScalar

subroutine MpiSumIntScalar(var,res)

    use mpih
    implicit none
    integer, intent(in)  :: var
    integer, intent(out) :: res
    call MPI_REDUCE(var,res,1,MPI_INTEGER,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiSumIntScalar

subroutine MpiAllMinRealScalar(var,res)

    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_ALLREDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiAllMinRealScalar
    
subroutine MpiAllMaxRealScalar(var,res)

    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_ALLREDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiAllMaxRealScalar

subroutine MpiAllSumRealScalar(var,res)

    use mpih
    implicit none
    real, intent(in)  :: var
    real, intent(out) :: res
    call MPI_ALLREDUCE(var,res,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiAllSumRealScalar

subroutine MpiAllMaxIntScalar(var,res)

    use mpih
    implicit none
    integer, intent(in)  :: var
    integer, intent(out) :: res
    call MPI_ALLREDUCE(var,res,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiAllMaxIntScalar

subroutine MpiAllSumIntScalar(var,res)

    use mpih
    implicit none
    integer, intent(in)  :: var
    integer, intent(out) :: res
    call MPI_ALLREDUCE(var,res,1,MPI_INTEGER,MPI_SUM,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiAllSumIntScalar

subroutine MpiSumReal1D(var,res,st,en)

    use mpih
    implicit none
    integer, intent(in) :: st,en
    real, intent(in),  dimension(st:en) :: var
    real, intent(out), dimension(st:en) :: res
    call MPI_REDUCE(var,res,en-st+1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,mpi_ierr)
    return

end subroutine MpiSumReal1D

subroutine MpiAbort

    use mpih
    implicit none
    call MPI_ABORT(MPI_COMM_WORLD,1,mpi_ierr)
    return
    
end subroutine MpiAbort
