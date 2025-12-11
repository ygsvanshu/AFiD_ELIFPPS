!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ReadInputFile.F90                              !
!    CONTAINS: subroutine ReadInputFile                   !
!                                                         !
!    PURPOSE: Read parameters from sovler.in file         !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadSolverInputs

    use mpih
    use param

    implicit none

    logical             :: exists
    integer             :: narg,ierror
    integer,parameter   :: nimp = 100
    character*200       :: filename,line
    character*100       :: ss(nimp)
    character*1         :: stringdummy1

    filename = "Inputs/solver.in"

    inquire(file=filename,exist=exists)
    if (.not.exists) then 
        write(*,*) 'Solver inputs not found'
        call MpiAbort
    end if

    open (unit=15,file=filename,status='old')

    ierror = 0
    do while(ierror.eq.0)

        read (15,'(a200)',iostat=ierror) line

        if (line(1:3) == '101') then
            !###### READ DATA ######
            call scan_string(line,1,ss,narg)
            stringdummy1 = ss(1)
            if ('y' == stringdummy1) then
                nread = .true.
            elseif ('n' == stringdummy1) then
                nread = .false.
            else
                write (*,*) "ERROR: Input value of parameter NREAD not valid"
                write (*,*) "       ===> valid values 'n' or 'y'            "
                call stop_config
            end if

        else if (line(1:3) == '201') then
            !###### GRID DIMENSION ######
            call scan_string(line,3,ss,narg)
            read (ss(1),*) nxm
            read (ss(2),*) nym
            read (ss(3),*) nzm
            nx = nxm + 1
            ny = nym + 1
            nz = nzm + 1

        else if (line(1:3) == '202') then
            !###### AXIS EXTENTS ######
            call scan_string(line,3,ss,narg)
            read (ss(1),*) xlen
            read (ss(2),*) ylen
            read (ss(3),*) zlen

        else if (line(1:3) == '203') then
            !###### STRETCHING DIRECTION ######
            call scan_string(line,1,ss,narg)
            stringdummy1 = ss(1)
            if ('n' == stringdummy1) then
                straxs = 0
            elseif ('x' == stringdummy1) then
                straxs = 1
            elseif ('y' == stringdummy1) then
                straxs = 2
            elseif ('z' == stringdummy1) then
                straxs = 3
            else
                write (*,*) "ERROR: Input value of parameter STRAXS not valid"
                write (*,*) "       ===> valid values 'n', 'x', 'y', or 'z'  "
                call stop_config
            end if

        else if (line(1:3) == '204') then
            !###### STRETCHING TYPE AND PARAMETER ######
            call scan_string(line,2,ss,narg)
            read (ss(1),*) strtyp
            read (ss(2),*) strval

        else if (line(1:3) == '205') then
            !###### PERIODICITY ######
            call scan_string(line,3,ss,narg)
            stringdummy1 = ss(1)
            if ('y' == stringdummy1) then
                periodic(1) = .true.
            elseif ('n' == stringdummy1) then
                periodic(1) = .false.
            else
                write (*,*) "ERROR: Input value of parameter PERX not valid"
                write (*,*) "       ===> valid values 'n' or 'y'            "
                call stop_config
            end if
            stringdummy1 = ss(2)
            if ('y' == stringdummy1) then
                periodic(2) = .true.
            elseif ('n' == stringdummy1) then
                periodic(2) = .false.
            else
                write (*,*) "ERROR: Input value of parameter PERY not valid"
                write (*,*) "       ===> valid values 'n' or 'y'            "
                call stop_config
            end if
            stringdummy1 = ss(3)
            if ('y' == stringdummy1) then
                periodic(3) = .true.
            elseif ('n' == stringdummy1) then
                periodic(3) = .false.
            else
                write (*,*) "ERROR: Input value of parameter PERZ not valid"
                write (*,*) "       ===> valid values 'n' or 'y'            "
                call stop_config
            end if

        else if (line(1:3) == '301') then
            !###### MAX STEPS ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) ntst

        else if (line(1:3) == '302') then
            !###### MAX TIME ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) tmax

        else if (line(1:3) == '303') then
            !###### MAX WALLCLOCK TIME ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) walltimemax

        else if (line(1:3) == '304') then
            !###### ABORT FILENAME ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) abortfile

        else if (line(1:3) == '401') then
            !###### NUMERICAL SCHEME ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) nsst

        else if (line(1:3) == '402') then
            !###### ADAPTIVE TIME STEP ######
            call scan_string(line,1,ss,narg)
            stringdummy1 = ss(1)
            if ('y' == stringdummy1) then
                vardt = .true.
            elseif ('n' == stringdummy1) then
                vardt = .false.
            else
                write (*,*) "ERROR: Input value of parameter VARDT not valid "
                write (*,*) "       ===> valid values 'n' or 'y'             "
                call stop_config
            end if

        else if (line(1:3) == '403') then
            !###### FIRST TIME STEP ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) dt

        else if (line(1:3) == '404') then
            !###### MIN MAX STEP SIZE ######
            call scan_string(line,2,ss,narg)
            read (ss(1),*) dtmin
            read (ss(2),*) dtmax

        else if (line(1:3) == '405') then
            !###### TIMESTEP STABILITY ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) limitCFL

        else if (line(1:3) == '406') then
            !###### RESID ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) resid

        else if (line(1:3) == '407') then
            !##### MAX VELOCITY ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) limitVel

        else if (line(1:3) == '408') then
            !###### INITIAL PERTURBATION ######
            call scan_string(line,2,ss,narg)
            read (ss(1),*) epsnum
            read (ss(2),*) eps

        else if (line(1:3) == '501') then
            !###### REYNOLDS NUMBER ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) rey

        else if (line(1:3) == '601') then
            !###### OUTPUT INTERVAL ######
            call scan_string(line,1,ss,narg)
            read (ss(1),*) nout

        else if (line(1:3) == '602') then
            !###### 1D STATISTICAL PROFILES ######
            call scan_string(line,3,ss,narg)
            stringdummy1 = ss(1)
            if ('y' == stringdummy1) then
                save1d = .true.
            elseif ('n' == stringdummy1) then
                save1d = .false.
            else
                write (*,*) "ERROR: Input value of parameter 1DON not valid   "
                write (*,*) "       ===> valid values 'n' or 'y'              "
                call stop_config
            end if
            read (ss(2),*) tsta1d
            read (ss(3),*) freq1d

        else if (line(1:3) == '603') then
            !###### 2D MOVIE SLICES ######
            call scan_string(line,3,ss,narg)
            stringdummy1 = ss(1)
            if ('y' == stringdummy1) then
                save2d = .true.
            elseif ('n' == stringdummy1) then
                save2d = .false.
            else
                write (*,*) "ERROR: Input value of parameter 2DON not valid   "
                write (*,*) "       ===> valid values 'n' or 'y'              "
                call stop_config
            end if
            read (ss(2),*) tsta2d
            read (ss(3),*) freq2d

        else if (line(1:3) == '604') then
            !###### 3D FLOW FIELD SNAPSHOTS ######
            call scan_string(line,3,ss,narg)
            stringdummy1 = ss(1)
            if ('y' == stringdummy1) then
                save3d = .true.
            elseif ('n' == stringdummy1) then
                save3d = .false.
            else
                write (*,*) "ERROR: Input value of parameter 3DON not valid   "
                write (*,*) "       ===> valid values 'n' or 'y'              "
                call stop_config
            end if
            read (ss(2),*) tsta3d
            read (ss(3),*) freq3d
        end if

    end do

    close (15)

end subroutine ReadSolverInputs

subroutine scan_string(S,N,SS,M)

    implicit none

    integer I,J,N,M,APO(2*N)
    character*200 S; character*100 SS(N)
    character*1,parameter :: PRIME = "'"

    ! Given the input string s, this routine determines the sub-strings ss of s
    ! which are delimited by consecutive pairs of primes ('). n (input) is the
    ! expected number of substrings, and m (output) is the effective number
    ! found. Apo is a pointer to the primes. The execution is stopped whenever
    ! i) n=0, ii) m=0, iii) m/=n, or iv) if m is odd. The maximum lenght of the
    ! substrings is 20, that of the input string is 200.  == GS Nov 17 2007 ==

    ! ----- Exits for n=0
    if (n == 0) then
        write (*,*) "SCAN_STRING has nothing to do"
        call stop_config
    end if

    ! ----- Looking for primes (', not ") within the input string
    i = 0
    do j = 1,200
        if (s(j:j) == prime) then
            i = i + 1; apo(i) = j
        end if
    end do
    m = i/2

    if (i == 0) then
        ! ----- Exits if no primes are found
        write (*,*) "SCAN_STRING has found NO primes"
        write (*,'(a100)') s
        call stop_config

    elseif (mod(i,2) /= 0) then
        ! ----- Exits if an odd number of primes is found
        write (*,*) "SCAN_STRING has found an odd number of primes:",i
        call stop_config

    elseif (m /= n) then
        ! ----- Exits if m/=n, otherwise determines the substrings
        write (*,*) "SCAN_STRING has found ",m," substrings"
        write (*,*) "SCAN_STRING  expected ",n," substrings"
        call stop_config
    else
        do i = 1,n
            ss(i) = s(apo(2*i - 1) + 1:apo(2*i) - 1)
        end do
    end if

end subroutine scan_string

subroutine stop_config

    use mpih
    implicit none
    write (*,*) 'STOP_CONFIG: The program will STOP! -----------------------------'
    call MPI_ABORT(MPI_COMM_WORLD,1,mpi_ierr)
    STOP 1
    
end subroutine stop_config