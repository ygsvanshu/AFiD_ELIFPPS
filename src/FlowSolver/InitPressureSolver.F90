!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: InitPressureSolver.F90                         !
!    CONTAINS: subroutine InitPressureSolver              !
!              subroutine InitPoissonOperator             !
!              subroutine CreateWavenumbers               !
!              subroutine InitFFTOperator                 !
!              subroutine SetPressureSolver               !
!              subroutine CheckPressureSolver             !
!              subroutine ChooseBestPressureSolver        !
!                                                         !
!    PURPOSE: Initialization routines for the pressure    !
!     solver. The pressure Poisson equation is solved by  !
!     multiplying with the eigenvectors of the discrete   !
!     second derivative in two directions and solving a   !
!     tridiagonal (Helmholtz) system in the third one.    !
!     The multiplication with the eigenvectors is done    !
!     with dense eigenvector matrices (GEMM, any grid) or !
!     with real-to-real FFTs (uniform grids only). The    !
!     fastest valid combination is chosen at start-up.    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitPressureSolver

    use mpih
    use decomp_2d, only: xstart,xend,ystart,yend,ysize,zstart,zend
    use param
    use local_arrays, only: rhx,rhy,rhz
    use pressure_solver

    implicit none

    integer :: rbfx,rbfy,rbfz
    integer :: i,j,k

    ! Allocate Poisson solver arrays
    if (.not.allocated(qfwx)) allocate(qfwx(1:nxm,1:nxm))
    if (.not.allocated(qbwx)) allocate(qbwx(1:nxm,1:nxm))
    if (.not.allocated(geix)) allocate(geix(1:nxm))
    if (.not.allocated(feix)) allocate(feix(1:nxm))
    if (.not.allocated(eigx)) allocate(eigx(1:nxm))

    if (.not.allocated(qfwy)) allocate(qfwy(1:nym,1:nym))
    if (.not.allocated(qbwy)) allocate(qbwy(1:nym,1:nym))
    if (.not.allocated(geiy)) allocate(geiy(1:nym))
    if (.not.allocated(feiy)) allocate(feiy(1:nym))
    if (.not.allocated(eigy)) allocate(eigy(1:nym))

    if (.not.allocated(qfwz)) allocate(qfwz(1:nzm,1:nzm))
    if (.not.allocated(qbwz)) allocate(qbwz(1:nzm,1:nzm))
    if (.not.allocated(geiz)) allocate(geiz(1:nzm))
    if (.not.allocated(feiz)) allocate(feiz(1:nzm))
    if (.not.allocated(eigz)) allocate(eigz(1:nzm))

    ! Compute the size of LWORK required for DSTEVD/DSYEVD
    ! In X-direction
    if (periodic(1)) then
        rbfx = 1 + 6*nxm + 2*nxm**2
    else
        rbfx = 1 + 4*nxm + nxm**2
    end if
    ! In Y-direction
    if (periodic(2)) then
        rbfy = 1 + 6*nym + 2*nym**2
    else
        rbfy = 1 + 4*nym + nym**2
    end if
    ! In Z-direction
    if (periodic(3)) then
        rbfz = 1 + 6*nzm + 2*nzm**2
    else
        rbfz = 1 + 4*nzm + nzm**2
    end if
    ! Determine the size of the temporary buffers based on the problem size
    niwk = 3 + 5*max(nxm,nym,nzm)
    nrwk = max(rbfx,rbfy,rbfz)
    ! Allocate temporary buffers 
    if (.not.allocated(iwrk)) allocate(iwrk(niwk))
    if (.not.allocated(rwrk)) allocate(rwrk(nrwk))

    ! Initialize the eigenvector matrices (GEMM) and their eigenvalues
    call InitPoissonOperator(nxm,periodic(1),pcbcxs,pcbcxe,dxc(1:nxm),amphi,acphi,apphi,qfwx,qbwx,geix)
    call InitPoissonOperator(nym,periodic(2),pcbcys,pcbcye,dyc(1:nym),amphj,acphj,apphj,qfwy,qbwy,geiy)
    call InitPoissonOperator(nzm,periodic(3),pcbczs,pcbcze,dzc(1:nzm),amphk,acphk,apphk,qfwz,qbwz,geiz)

    ! Deallocate temporary buffers
    if (allocated(iwrk)) deallocate(iwrk)
    if (allocated(rwrk)) deallocate(rwrk)

    ! Modified wavenumbers (eigenvalues -ak for the FFT, only used on uniform grids)
    call CreateWavenumbers(nxm,xlen,periodic(1),pcbcxs,pcbcxe,feix)
    call CreateWavenumbers(nym,ylen,periodic(2),pcbcys,pcbcye,feiy)
    call CreateWavenumbers(nzm,zlen,periodic(3),pcbczs,pcbcze,feiz)

    ! Write out the eigenvalues of both methods
    if (ismaster) then

        open (unit=78, file='Results/eigx.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,2(A24,X))') 'i', 'geix(i)', 'feix(i)'
        do i = 1, nxm
            write (78, '(I8,X,2(ES24.16,X))') i, geix(i), feix(i)
        end do
        close (78)

        open (unit=78, file='Results/eigy.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,2(A24,X))') 'j', 'geiy(j)', 'feiy(j)'
        do j = 1, nym
            write (78, '(I8,X,2(ES24.16,X))') j, geiy(j), feiy(j)
        end do
        close (78)

        open (unit=78, file='Results/eigz.out', status='unknown', position='rewind', access='sequential')
        write (78, '(A8,X,2(A24,X))') 'k', 'geiz(k)', 'feiz(k)'
        do k = 1, nzm
            write (78, '(I8,X,2(ES24.16,X))') k, geiz(k), feiz(k)
        end do
        close (78)

    end if

    ! Allocate GEMM buffer
    gszx = (xend(2)-xstart(2)+1)*(xend(3)-xstart(3)+1)  ! Solved batched without transpose
    gszy = (yend(1)-ystart(1)+1)                        ! Solved per k-index with transpose
    gszz = (zend(1)-zstart(1)+1)*(zend(2)-zstart(2)+1)  ! Solved batched with transpose
    if (.not.allocated(gwrk)) allocate(gwrk(max(nxm*gszx,nym*gszy,nzm*gszz)))

    ! Initialize the FFT transforms in the directions with a uniform grid. Each
    ! plan transforms all lines of a pencil in place. Arguments after the grid
    ! spacing: memory distance between two points of a line, and two loops over
    ! the lines (number of lines, memory distance between lines). Planning with
    ! FFTW_MEASURE overwrites rhx/rhy/rhz, which do not hold any data yet.
    call InitFourierOperator(nxm,periodic(1),pcbcxs,pcbcxe,dxc(1:nxm),1,   (/gszx,1/),       (/nxm,0/),     rhx,pfwx,pbwx,fnrm(1),fftok(1))
    call InitFourierOperator(nym,periodic(2),pcbcys,pcbcye,dyc(1:nym),gszy,(/gszy,ysize(3)/),(/1,gszy*nym/),rhy,pfwy,pbwy,fnrm(2),fftok(2))
    call InitFourierOperator(nzm,periodic(3),pcbczs,pcbcze,dzc(1:nzm),gszz,(/gszz,1/),       (/1,0/),       rhz,pfwz,pbwz,fnrm(3),fftok(3))

    ! The FFT is only used in a direction if it is possible on all ranks
    call MPI_Allreduce(MPI_IN_PLACE,fftok,3,MPI_LOGICAL,MPI_LAND,MPI_COMM_WORLD,mpi_ierr)

end subroutine InitPressureSolver

subroutine InitPoissonOperator(ngm,perc,bcst,bcen,dgc,amphg,acphg,apphg,qfwd,qbwd,eign)

    use param, only: ismaster
    use boundary_arrays, only: DIRICHLET, NEUMANN
    use pressure_solver, only: niwk, nrwk, iwrk, rwrk

    implicit none

    integer, intent(in) :: ngm
    logical, intent(in) :: perc
    integer, intent(in) :: bcst,bcen
    real, intent(in)    :: dgc(1:ngm)
    real, intent(in)    :: amphg(1:ngm)
    real, intent(in)    :: acphg(1:ngm)
    real, intent(in)    :: apphg(1:ngm)
    real, intent(inout) :: qfwd(1:ngm,1:ngm)
    real, intent(inout) :: qbwd(1:ngm,1:ngm)
    real, intent(inout) :: eign(1:ngm)

    integer             :: gm,gc,gp,nc
    real                :: fm,fc,fp
    integer             :: info
    real                :: uppr(1:ngm)
    real                :: diag(1:ngm)

    ! Initialize to zeros
    qfwd(1:ngm,1:ngm) = 0.0
    qbwd(1:ngm,1:ngm) = 0.0
    eign(1:ngm) = 0.0

    ! Construct the tridiagonal discrete second derivative operators
    do gc = 1,ngm
        gm = gc - 1
        gp = gc + 1
    
        fm = amphg(gc)
        fc = acphg(gc)
        fp = apphg(gc)
        
        if (gc.eq.1) then
            if (perc) then
                qfwd(gc,ngm) = qfwd(gc,ngm) + fm
            else
                if (bcst == DIRICHLET) fc = fc - fm
                if (bcst == NEUMANN)   fc = fc + fm
            end if
        else
            qfwd(gc,gm) = qfwd(gc,gm) + fm
        end if
        if (gc.eq.ngm) then
            if (perc) then
                qfwd(gc,1) = qfwd(gc,1) + fp
            else
                if (bcen == DIRICHLET) fc = fc - fp
                if (bcen == NEUMANN)   fc = fc + fp
            end if
        else
            qfwd(gc,gp) = qfwd(gc,gp) + fp
        end if
        qfwd(gc,gc) = qfwd(gc,gc) + fc
        diag(gc) = sqrt(1.0/dgc(gc))
    end do

    ! Symmetrize the discrete second derivative operators
    do gc = 1,ngm
        do nc = 1,ngm
            qfwd(gc,nc) = qfwd(gc,nc)*(diag(nc)/diag(gc))
        end do
    end do
    
    if (perc) then
        ! Initialize the eigenvector and eigenvalue operators for the Poisson solver
        call dsyevd('V','U',ngm,qfwd,ngm,eign,rwrk,nrwk,iwrk,niwk,info)
    else
        ! Keep the diagonal elements for the eigenvalue solver
        do gc = 1,ngm
            eign(gc) = qfwd(gc,gc)
        end do
        ! Keep the upper diagonal elements for the eigenvalue solver
        do gc = 1,ngm-1
            gp = gc + 1
            uppr(gc) = qfwd(gc,gp)
        end do
        ! Initialize the eigenvector and eigenvalue operators for the Poisson solver without periodicity
        call dstevd('V',ngm,eign,uppr,qfwd,ngm,rwrk,nrwk,iwrk,niwk,info)
    end if

    ! Compute the eigenvector operators for the Poisson solver along (i) axis
    qbwd = transpose(qfwd)
    do gc = 1,ngm
        do nc = 1,ngm
            qfwd(gc,nc) = qfwd(gc,nc)*diag(gc)
            qbwd(gc,nc) = qbwd(gc,nc)/diag(nc)
        end do
    end do

end subroutine InitPoissonOperator

! MODIFIED WAVENUMBERS OF A UNIFORM GRID
! The FFT based multiplication with the eigenvectors of the second
! derivative uses real-to-real FFTW transforms, chosen by the boundary
! conditions (see InitFFTOperator). The eigenvalues are -akc, stored in
! the same order as the coefficients produced by FFTW:
!   periodic            : R2HC/HC2R (halfcomplex, modes m and nc-m share a value)
!   Neumann-Neumann     : REDFT10/REDFT01
!   Dirichlet-Dirichlet : RODFT10/RODFT01
!   mixed               : REDFT11 or RODFT11
! On a stretched grid the values are not used (the FFT is not possible).
subroutine CreateWavenumbers(nc,ln,perc,bcst,bcen,akc)

    use param, only: pi
    use boundary_arrays, only: DIRICHLET, NEUMANN

    implicit none

    integer, intent(in) :: nc
    real, intent(in)    :: ln
    logical, intent(in) :: perc
    integer, intent(in) :: bcst,bcen
    real, intent(inout) :: akc(nc)

    integer             :: lc
    real                :: ao

    do lc = 1,nc
        if (perc) then
            ao = 2.0*pi*real(lc-1)
        else if ((bcst.eq.NEUMANN).and.(bcen.eq.NEUMANN)) then
            ao = pi*real(lc-1)
        else if ((bcst.eq.DIRICHLET).and.(bcen.eq.DIRICHLET)) then
            ao = pi*real(lc)
        else
            ao = pi*(real(lc)-0.5)
        end if
        akc(lc) = 2.0*(1.0 - cos(ao/real(nc)))*(real(nc)/ln)**2
    end do

end subroutine CreateWavenumbers

! Create the FFTW plans to multiply with the eigenvectors of the second
! derivative along one direction, if the grid in that direction is uniform.
! The real-to-real transform kind follows from the boundary conditions, see
! CreateWavenumbers for the matching eigenvalues (-ak).
subroutine InitFourierOperator(ngm,perc,bcst,bcen,dgc,ist,hmn,hms,rr,pfwd,pbwd,fnrmd,fftokd)

    use boundary_arrays, only: DIRICHLET, NEUMANN
    use pressure_solver

    implicit none

    integer, intent(in)     :: ngm              ! Number of cells along the line
    logical, intent(in)     :: perc             ! Periodicity
    integer, intent(in)     :: bcst,bcen        ! Boundary condition at start/end
    real, intent(in)        :: dgc(1:ngm)       ! Grid spacing along the line
    integer, intent(in)     :: ist              ! Memory distance between two points of a line
    integer, intent(in)     :: hmn(2),hms(2)    ! Two loops over the lines: number of lines, memory distance
    real, intent(inout)     :: rr(*)            ! Array that is transformed in place
    type(C_PTR), intent(out):: pfwd,pbwd        ! Plans back to physical space / to eigenspace
    real, intent(out)       :: fnrmd            ! Normalisation of the two transforms together
    logical, intent(out)    :: fftokd           ! True if the FFT can be used

    integer(C_FFTW_R2R_KIND):: kfw(1),kbw(1)    ! FFTW real-to-real transform kinds
    type(fftw_iodim)        :: iodim(1)         ! Transform along the line
    type(fftw_iodim)        :: iodim_howmany(2) ! Two loops over the lines

    pfwd   = C_NULL_PTR
    pbwd   = C_NULL_PTR
    fnrmd  = 1.0
    fftokd = .false.

    ! The FFT needs a uniform grid with at least two cells
    if (ngm.lt.2) return
    if (maxval(abs(dgc(1:ngm) - dgc(1))).gt.(1.0e-10*dgc(1))) return

    ! Transform kinds and normalisation belonging to the boundary conditions
    if (perc) then
        kbw   = FFTW_R2HC
        kfw   = FFTW_HC2R
        fnrmd = real(ngm)
    else if ((bcst.eq.NEUMANN).and.(bcen.eq.NEUMANN)) then
        kbw   = FFTW_REDFT10
        kfw   = FFTW_REDFT01
        fnrmd = 2.0*real(ngm)
    else if ((bcst.eq.DIRICHLET).and.(bcen.eq.DIRICHLET)) then
        kbw   = FFTW_RODFT10
        kfw   = FFTW_RODFT01
        fnrmd = 2.0*real(ngm)
    else if (bcst.eq.NEUMANN) then
        kbw   = FFTW_REDFT11
        kfw   = FFTW_REDFT11
        fnrmd = 2.0*real(ngm)
    else
        kbw   = FFTW_RODFT11
        kfw   = FFTW_RODFT11
        fnrmd = 2.0*real(ngm)
    end if

    ! Sizes and memory distances (in place, so input and output are the same)
    iodim(1)%n          = ngm
    iodim(1)%is         = ist
    iodim(1)%os         = ist

    iodim_howmany(1)%n  = hmn(1)
    iodim_howmany(1)%is = hms(1)
    iodim_howmany(1)%os = hms(1)

    iodim_howmany(2)%n  = hmn(2)
    iodim_howmany(2)%is = hms(2)
    iodim_howmany(2)%os = hms(2)

    ! Plan both transforms (all lines of the pencil at once). The transforms are
    ! in place, so rr is both input and output; compilers may warn about this,
    ! which is harmless for the FFTW C library.
    pbwd = fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,rr,rr,kbw,fftflg)
    pfwd = fftw_plan_guru_r2r(1,iodim,2,iodim_howmany,rr,rr,kfw,fftflg)

    ! A plan is null if FFTW could not create it (e.g. FFTW wrappers without real-to-real support)
    fftokd = c_associated(pbwd).and.c_associated(pfwd)

end subroutine InitFourierOperator

! Select the Helmholtz direction and the method to multiply with the
! eigenvectors in each direction, together with the matching eigenvalues
subroutine SetPressureSolver(hd,tx,ty,tz)

    use param
    use boundary_arrays, only: NEUMANN
    use pressure_solver

    implicit none

    integer, intent(in) :: hd           ! Direction of the tridiagonal (Helmholtz) solve
    integer, intent(in) :: tx,ty,tz     ! Method (TGEMM or TFFT) in X, Y and Z

    hdir = hd
    tmet = (/tx,ty,tz/)

    ! Eigenvalues belonging to the selected methods
    if (tx.eq.TFFT) then
        eigx(1:nxm) = -feix(1:nxm)
    else
        eigx(1:nxm) = geix(1:nxm)
    end if
    if (ty.eq.TFFT) then
        eigy(1:nym) = -feiy(1:nym)
    else
        eigy(1:nym) = geiy(1:nym)
    end if
    if (tz.eq.TFFT) then
        eigz(1:nzm) = -feiz(1:nzm)
    else
        eigz(1:nzm) = geiz(1:nzm)
    end if

    ! Zero-mode index for each direction (only exists if periodic or Neumann-Neumann)
    izmx = 0
    izmy = 0
    izmz = 0
    if (periodic(1).or.((pcbcxs.eq.NEUMANN).and.(pcbcxe.eq.NEUMANN))) izmx = maxloc(eigx,1)
    if (periodic(2).or.((pcbcys.eq.NEUMANN).and.(pcbcye.eq.NEUMANN))) izmy = maxloc(eigy,1)
    if (periodic(3).or.((pcbczs.eq.NEUMANN).and.(pcbcze.eq.NEUMANN))) izmz = maxloc(eigz,1)
    if (izmx.gt.0) eigx(izmx) = 0.0
    if (izmy.gt.0) eigy(izmy) = 0.0
    if (izmz.gt.0) eigz(izmz) = 0.0

end subroutine SetPressureSolver

! Check a combination of methods against the one using eigenvector matrices
! (GEMM) in all directions, for the same Helmholtz direction and the same
! right-hand side. Both solve the same discrete equations, so the solutions
! must agree up to round-off. This catches any problem with the FFTs (e.g.
! an FFTW library that does not do real-to-real transforms correctly).
subroutine CheckPressureSolver(hd,tx,ty,tz,rref,ok)

    use mpih
    use decomp_2d, only: xstart,xend
    use param
    use local_arrays, only: rhx
    use pressure_solver

    implicit none

    integer, intent(in)  :: hd,tx,ty,tz
    real, intent(in)     :: rref(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))
    logical, intent(out) :: ok

    real, allocatable    :: rsol(:,:,:)
    real                 :: errl(2),errg(2)

    ok = .true.
    if ((tx.eq.TGEMM).and.(ty.eq.TGEMM).and.(tz.eq.TGEMM)) return

    allocate(rsol(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))

    ! Reference solution with eigenvector matrices in all directions
    call SetPressureSolver(hd,TGEMM,TGEMM,TGEMM)
    rhx = rref
    call SolvePressureCorrection
    rsol = rhx

    ! Solution with the combination to check
    call SetPressureSolver(hd,tx,ty,tz)
    rhx = rref
    call SolvePressureCorrection

    ! Largest difference relative to the largest value, over all ranks
    errl(1) = maxval(abs(rhx - rsol))
    errl(2) = maxval(abs(rsol))
    call MPI_Allreduce(errl,errg,2,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_ierr)
    ok = (errg(1).le.(1.0e-6*errg(2)))

    deallocate(rsol)

end subroutine CheckPressureSolver

subroutine ChooseBestPressureSolver(verbose)

    use mpih
    use decomp_2d, only: xstart,xend
    use param
    use local_arrays, only: rhx
    use pressure_solver

    implicit none

    logical, intent(in) :: verbose
    integer             :: hd,tx,ty,tz,inum
    integer             :: bhd,btx,bty,btz          ! Selected combination
    integer             :: ngl(3),tcmb(3)
    real                :: t_start,t_run,t_best
    logical             :: ok
    character(4)        :: tnam(3)
    character(1)        :: dnam(3) = (/'X','Y','Z'/)
    real, allocatable   :: rref(:,:,:)

    ngl = (/nxm,nym,nzm/)

    ! Default combination: FFT wherever possible, and the tridiagonal (Helmholtz)
    ! solve along the largest direction without FFT, or else along the largest
    ! direction. The tridiagonal solve needs at least 3 cells.
    bhd = maxloc(ngl + merge(sum(ngl),0,.not.fftok), 1, mask=(ngl.ge.3))
    if (bhd.eq.0) then
        if (ismaster) write (6,*) 'Error: the pressure solver needs at least 3 cells in one direction'
        call MpiAbort
    end if
    tcmb      = merge(TFFT,TGEMM,fftok)
    tcmb(bhd) = TGEMM
    btx       = tcmb(1)
    bty       = tcmb(2)
    btz       = tcmb(3)

    ! Random right-hand side for checking and timing the combinations
    allocate(rref(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    call random_number(rref)

    if (itpr.le.0) then

        ! No speed test: only check the default combination, fall back to GEMM if it fails
        call CheckPressureSolver(bhd,btx,bty,btz,rref,ok)
        if (.not.ok) then
            btx = TGEMM
            bty = TGEMM
            btz = TGEMM
        end if

    else

        if (ismaster.and.verbose) then
            write (6,'(A77)') '**** Starting pressure solver speed test =============================== ****'
            write (6,'(A77)') '                                                                             '
        end if

        t_best = huge(1.0)

        ! Loop over all combinations of Helmholtz direction and methods
        do hd = 1,3
            do tx = TGEMM,TFFT
                do ty = TGEMM,TFFT
                    do tz = TGEMM,TFFT

                        tcmb = (/tx,ty,tz/)

                        ! The Helmholtz direction needs at least 3 cells and is not transformed
                        if (ngl(hd).lt.3) cycle
                        if (tcmb(hd).ne.TGEMM) cycle
                        ! The FFT can only be used where it is possible
                        if (any((tcmb.eq.TFFT).and.(.not.fftok))) cycle

                        tnam(1:3) = merge(' FFT','GEMM',tcmb.eq.TFFT)
                        tnam(hd)  = '   -'

                        ! Check the combination against the eigenvector matrices
                        call CheckPressureSolver(hd,tx,ty,tz,rref,ok)
                        if (.not.ok) then
                            if (ismaster.and.verbose) write (6,'(A,A1,A,A4,A,A4,A,A4,A)') '     Helmholtz in ',dnam(hd),'   X:',tnam(1),'  Y:',tnam(2),'  Z:',tnam(3),'   failed check, skipped'
                            cycle
                        end if

                        ! Mean time of a pressure solve (the same on all ranks)
                        call SetPressureSolver(hd,tx,ty,tz)
                        t_run = 0.0
                        do inum = 1,itpr
                            rhx = rref
                            call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)
                            t_start = MPI_WTIME()
                            call SolvePressureCorrection
                            call MPI_Barrier(MPI_COMM_WORLD, mpi_ierr)
                            t_run = t_run + (MPI_WTIME() - t_start)
                        end do
                        t_run = t_run/real(itpr)
                        call MPI_Allreduce(MPI_IN_PLACE,t_run,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,mpi_ierr)

                        if (ismaster.and.verbose) write (6,'(A,A1,A,A4,A,A4,A,A4,A,ES10.3,A)') '     Helmholtz in ',dnam(hd),'   X:',tnam(1),'  Y:',tnam(2),'  Z:',tnam(3),'   time:',t_run,' s'

                        if (t_run.lt.t_best) then
                            t_best = t_run
                            bhd    = hd
                            btx    = tx
                            bty    = ty
                            btz    = tz
                        end if

                    end do
                end do
            end do
        end do

        if (ismaster.and.verbose) then
            write (6,'(A77)') '                                                                             '
            write (6,'(A77)') '**** Finished pressure solver speed test =============================== ****'
            write (6,'(A77)') '                                                                             '
        end if

    end if

    call SetPressureSolver(bhd,btx,bty,btz)

    ! Eigenvector matrices that are not used anymore (FFT directions and the
    ! Helmholtz direction) are freed: they are stored on every rank and can
    ! be as large as the flow field arrays. After this, the combination can
    ! not be changed anymore.
    if ((tmet(1).eq.TFFT).or.(hdir.eq.1)) deallocate(qfwx,qbwx)
    if ((tmet(2).eq.TFFT).or.(hdir.eq.2)) deallocate(qfwy,qbwy)
    if ((tmet(3).eq.TFFT).or.(hdir.eq.3)) deallocate(qfwz,qbwz)

    deallocate(rref)

    return

end subroutine ChooseBestPressureSolver