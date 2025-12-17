subroutine PrintCaseInfo

    use mpih
    use param

    implicit none

    character :: strchr
    
    if (straxs.eq.1) then
        strchr = 'x'
    else if (straxs.eq.2) then
        strchr = 'y'
    else if (straxs.eq.3) then
        strchr = 'z'
    else
        strchr = 'n'
    end if

    if (ismaster) then

        ! "Nothing wrong with a little bit of artistic liberty :3 " - Vanshu

        write (6,'(A77)') '============================================================================='
        write (6,'(A77)') '                           __________        _____                           '
        write (6,'(A77)') '                          /    _____/  _    |  __ \                          '
        write (6,'(A77)') '                         / /| |       (_)   | |  \ \                         '
        write (6,'(A77)') '========================/ /=| |=============| |===\ \========================'
        write (6,'(A77)') '=======================/ ___   ___/===| |===| |===/ /========================'
        write (6,'(A77)') '                      / /   | |       | |   | |__/ /                         '
        write (6,'(A77)') '                     /_/    |_|       |_|   |_____/                          '
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '============================================================================='
        write (6,'(A77)') '          _          _     _  _       _   _      _           _   _           '
        write (6,'(A77)') '         |_) |_| \/ (_` | /  (_`     / \ |_     |_ |  | | | | \ (_`          '
        write (6,'(A77)') '         |   | | /  ._) | \_ ._)     \_/ |      |  |_ |_| | |_/ ._)          '
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '============================================================================='
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '                            Navier-Stokes Solver                             '
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** PERIODICITY  ====================================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,L8)') '     Periodic in X                    = ',periodic(1)
        write (6,'(A40,L8)') '     Periodic in Y                    = ',periodic(2)
        write (6,'(A40,L8)') '     Periodic in Z                    = ',periodic(3)
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** 3D CELL DIMENSIONS ================================================ ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,F8.4)') '     Domain size in X                 = ',xlen
        write (6,'(A40,F8.4)') '     Domain size in Y                 = ',ylen
        write (6,'(A40,F8.4)') '     Domain size in Z                 = ',zlen
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** GRID RESOLUTION =================================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,I8)') '     Number of grid points in X       = ',nxm
        write (6,'(A40,I8)') '     Number of grid points in Y       = ',nym
        write (6,'(A40,I8)') '     Number of grid points in Z       = ',nzm
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** GRID STRETCHING =================================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,A8)') '     Grid stretching axis             = ',strchr
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** PARALLELIZATION =================================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,I8)') '     Number of MPI processes          = ',mpi_size
        write (6,'(A40,I8)') '     Number of MPI pencils in Y       = ',mpi_dims(1)
        write (6,'(A40,I8)') '     Number of MPI pencils in Z       = ',mpi_dims(2)
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** PHYSICAL PARAMETERS =============================================== ****'
        write (6,'(A77)') '                                                                             '
        write (6,'(A40,F10.3)') '     Reynolds Number                  = ',rey
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '**** TIME-STEPPING ===================================================== ****'
        write (6,'(A77)') '                                                                             '

        if (vardt) then
            write (6,'(A40,F8.4)') '     Variable dt, Fixed CFL           = ',limitCFL
            write (6,'(A40,ES10.3)') '     Variable dt, maximum dt          = ',dtmax
            write (6,'(A40,ES10.3)') '     Variable dt, minimum dt          = ',dtmin
        else
            write (6,'(A40,F8.4)') '     Fixed dt, dt                     = ',dtmax
            write (6,'(A40,F8.4)') '     Fixed dt, maximum CFL            = ',limitCFL
        end if
        write (6,'(A40,I15)') '     Maximum number of timesteps      = ',ntst

        if (nsst .gt. 1) then
            write (6,'(A62)') '     Time-stepping scheme             =  III order Runge-Kutta'
            write (6,'(A40,3F8.3)') '     gam                              = ',(gam(ns),ns=1,nsst)
            write (6,'(A40,3F8.3)') '     ro                               = ',(rom(ns),ns=1,nsst)
        else
            write (6,'(A55)') '     Time-stepping scheme             =  Adams-Bashfort'
            write (6,'(A40,F8.3)') '     gam                              = ',gam(1)
            write (6,'(A40,F8.3)') '     ro                               = ',rom(1)
        end if
        write (6,'(A77)') '                                                                             '
        write (6,'(A77)') '============================================================================='
        write (6,'(A77)') '                                                                             '

    end if

end subroutine PrintCaseInfo

subroutine PrintStepInfo(wct,ela,dmax,davg,icfl)

    use decomp_2d, only: nproc
    use param

    implicit none

    real,intent(in) :: wct,ela
    real,intent(in) :: dmax,davg,icfl
    real            :: eta,res,rut,prf
    integer         :: hr1,mn1,sc1
    integer         :: hr2,mn2,sc2

    call MpiMaxRealScalar(wct,res)
    rut = res

    res = ela

    if (ismaster) then

        prf = (rut*nproc*1.0e6)/(nxm*nym*nzm)

        hr1 = int(res/3600)
        res = res - real(3600*hr1)
        mn1 = int(res/60)
        res = res - real(60*mn1)
        sc1 = int(res)

        if ((time.gt.0.0).and.(ntime.gt.0)) then
            eta = min( ela*((tmax - time)/time), ela*(real(ntst - ntime)/real(ntime)) )
        else
            eta = walltimemax
        end if

        hr2 = int(eta/3600)
        eta = eta - real(3600*hr2)
        mn2 = int(eta/60)
        eta = eta - real(60*mn2)
        sc2 = int(eta)

        write(6,'(A)') '-----------+------------+------------+------------+------------+------------+------------+------------+------------+-----------+------------+------------+------------'
        write(6,'(A)') '  Timestep |  FlowTime  |     DT     | Global Div |  Max Div   |  vmax(1)   |  vmax(2)   |  vmax(3)   | CFL number | CPU Time  |   CPU DT   |  CPU Perf  |     ETA    '
        write(6,'(I10,A,8(ES10.3,A),I3.2,2(A,I2.2),A,ES10.3,A,ES10.3,A,I3.2,2(A,I2.2))') ntime,' | ',time,' | ',dt,' | ',davg,' | ',dmax,' | ',vmax(1),' | ',vmax(2),' | ',vmax(3),' | ',icfl*dt,' | ',hr1,':',mn1,':',sc1,' | ',rut,' | ',prf,' | ',hr2,':',mn2,':',sc2

    end if


end subroutine PrintStepInfo