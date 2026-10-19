!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: Param.F90                                      !
!    CONTAINS: module param                               !
!              module local_arrays                        !
!              module boundary_arrays                     !
!              module stat_arrays                         !
!              module movie_arrays                        !
!              module body_force                          !
!              module mpih                                !
!              module tridiagonal_buffer                  !
!              module fftw_params                         !
!              module pressure_solver                     !
!                                                         ! 
!    PURPOSE: All modules required for the flow solver    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declaration of global variables
module param

    implicit none
        
    ! Read from input file solver.in

    ! 1XX -----------------------
    logical         :: nread                                        ! Flag indicating if the input file has been read
    ! 2XX -----------------------
    integer         :: nxm, nym, nzm                                ! Number of grid cells in x, y, z directions
    real            :: xlen,ylen,zlen                               ! Domain lengths in x, y, z directions
    integer         :: gtpx,gtpy,gtpz                               ! Grid stretching type in x, y, z directions 
    real            :: gvlx,gvly,gvlz                               ! Grid stretching parameter value in x, y, z directions
    logical         :: periodic(3) = (/.false.,.false.,.false./)    ! Periodic boundary conditions in x, y, z directions
    ! 3XX -----------------------
    integer         :: ntst                                         ! Maximum number of time steps  
    real            :: tmax                                         ! Maximum simulation time
    real            :: walltimemax                                  ! Maximum allowed wallclock time
    character*200   :: abortfile                                    ! Filename for abort condition
    ! 4XX -----------------------
    integer         :: nsst 
    logical         :: vardt = .true.                               ! Flag for adaptive time stepping
    real            :: dt                                           ! First time step (or constant time step if vardt = .false.)
    real            :: dtmin,dtmax                                  ! Minimum and maximum allowed time step sizes
    real            :: limitCFL                                     ! CFL limit for adaptive time stepping
    real            :: resid                                        ! Maximum allowed divergence
    real            :: limitVel                                     ! Maximum allowed velocity
    integer         :: epsnum                                       ! Number of initial perturbations
    real            :: eps                                          ! Initial perturbation magnitude
    integer         :: itpr = 10                                    ! Pressure solver speed test iterations
    ! 5XX -----------------------
    real            :: rey                                          ! Reynolds number
    real            :: pra                                          ! Prandtl number
    ! 6XX -----------------------
    integer         :: nout                                         ! Output interval
    real            :: csav                                         ! Time interval for saving checkpoint data (continua.h5)
    logical         :: save1d,save2d,save3d                         ! Flags, start times, and save intervals for saving 1D, 2D, 3D data        
    real            :: tsta1d,tsta2d,tsta3d                         ! Flags, start times, and save intervals for saving 1D, 2D, 3D data        
    real            :: freq1d,freq2d,freq3d                         ! Flags, start times, and save intervals for saving 1D, 2D, 3D data

    ! Other global variables

    integer                 :: nx, ny, nz                           ! Number of grid points in x, y, z directions
    integer                 :: ntime,ns                             ! Number of time steps and sub-steps

    real                    :: time,tsta                            ! Current simulation time and start time
    real                    :: pi = 2.0*asin(1.0)                   ! Constant pi
    real                    :: al,ga,ro                             ! Timestepping parameters: alpha, gamma, rho  

    real                    :: cvel = 1.0                           ! Non-reflecting (SOMMERFELD) boundary characteristic velocity   

    real, dimension(1:3)    :: vmax                                 ! Maximum velocity in x, y, z directions
    real, dimension(1:3)    :: vavg                                 ! Average velocity in x, y, z directions
    real, dimension(1:3)    :: vrms                                 ! Root mean square velocity in x, y, z directions
    real                    :: vmag                                 ! Magnitude of velocity
    real                    :: hmax                                 ! Maximum temperature
    real                    :: havg                                 ! Average temperature
    real                    :: hrms                                 ! Root mean square temperature

    real, dimension(1:3)    :: gam,rom,alm                          ! Timestepping parameter arrays
              
    logical                 :: ismaster = .false.                   ! Flag indicating if the current process is the master process
    integer                 :: lvlhalo = 1                          ! Level of halo cells
    integer                 :: mpi_xcut,mpi_ycut,mpi_zcut           ! MPI domain decomposition cuts in x, y, z directions  

    logical                 :: particle = .false.                   ! If particle solver should run
    logical                 :: pre_diff = .false.                   ! If diffusive terms is precalculated for particle solver

    ! Grid parameters 

    real, allocatable, dimension(:)     :: xc,xm
    real, allocatable, dimension(:)     :: yc,ym
    real, allocatable, dimension(:)     :: zc,zm

    real, allocatable, dimension(:)     :: dxc,dxm
    real, allocatable, dimension(:)     :: dyc,dym
    real, allocatable, dimension(:)     :: dzc,dzm

    real, allocatable, dimension(:)     :: ixcm,ixcp
    real, allocatable, dimension(:)     :: iycm,iycp
    real, allocatable, dimension(:)     :: izcm,izcp

    real, allocatable, dimension(:)     :: ap1si,ac1si,am1si
    real, allocatable, dimension(:)     :: ap1ci,ac1ci,am1ci
    real, allocatable, dimension(:)     :: amphi,acphi,apphi

    real, allocatable, dimension(:)     :: ap2sj,ac2sj,am2sj
    real, allocatable, dimension(:)     :: ap2cj,ac2cj,am2cj
    real, allocatable, dimension(:)     :: amphj,acphj,apphj

    real, allocatable, dimension(:)     :: ap3sk,ac3sk,am3sk
    real, allocatable, dimension(:)     :: ap3ck,ac3ck,am3ck
    real, allocatable, dimension(:)     :: amphk,acphk,apphk
        
end module param

module local_arrays

    implicit none

    real,allocatable,dimension(:,:,:)   :: vx,vy,vz,pr
    real,allocatable,dimension(:,:,:)   :: adx,ady,adz
    real,allocatable,dimension(:,:,:)   :: dfx,dfy,dfz
    real,allocatable,dimension(:,:,:)   :: rhx,rhy,rhz
    real,allocatable,dimension(:,:,:)   :: rux,ruy,ruz
    real,allocatable,dimension(:,:,:)   :: dph,dphhalo

end module local_arrays

module boundary_arrays

    implicit none

    ! Generic boundary condition types
    integer, parameter                      :: DIRICHLET  = 1
    integer, parameter                      :: NEUMANN    = 2
    integer, parameter                      :: SOMMERFELD = 3
    ! Arrays to store boundary types (for default pencil)
    integer,allocatable,dimension(:,:,:)    :: btypxs,btypxe
    integer,allocatable,dimension(:,:,:)    :: btypys,btypye
    integer,allocatable,dimension(:,:,:)    :: btypzs,btypze
    ! Arrays to store boundary values (for default pencil)
    real,allocatable,dimension(:,:,:)       :: bvalxs,bvalxe
    real,allocatable,dimension(:,:,:)       :: bvalys,bvalye
    real,allocatable,dimension(:,:,:)       :: bvalzs,bvalze
    ! Arrays to store coefficients for implicit solver
    real,allocatable,dimension(:,:,:)       :: cfbcxs,cfbcxe
    real,allocatable,dimension(:,:,:)       :: cfbcys,cfbcye
    real,allocatable,dimension(:,:,:)       :: cfbczs,cfbcze
    ! Arrays to store magnitude of global flux correction allowed
    real,allocatable,dimension(:,:)         :: glpcxs,glpcxe
    real,allocatable,dimension(:,:)         :: glpcys,glpcye
    real,allocatable,dimension(:,:)         :: glpczs,glpcze

end module boundary_arrays

module stat_arrays

    implicit none

    integer                             :: pnum
    integer                             :: snpp

    real,allocatable,dimension(:)       :: vx_m1_xp,vx_m2_xp
    real,allocatable,dimension(:)       :: vx_m1_yp,vx_m2_yp
    real,allocatable,dimension(:)       :: vx_m1_zp,vx_m2_zp

    real,allocatable,dimension(:)       :: vy_m1_xp,vy_m2_xp
    real,allocatable,dimension(:)       :: vy_m1_yp,vy_m2_yp
    real,allocatable,dimension(:)       :: vy_m1_zp,vy_m2_zp

    real,allocatable,dimension(:)       :: vz_m1_xp,vz_m2_xp
    real,allocatable,dimension(:)       :: vz_m1_yp,vz_m2_yp
    real,allocatable,dimension(:)       :: vz_m1_zp,vz_m2_zp

end module stat_arrays

module movie_arrays

    implicit none

    integer                             :: mnum
    integer                             :: snpm
    logical,allocatable,dimension(:,:)  :: sqnx,sqny,sqnz
    real,allocatable,dimension(:)       :: slcx,slcy,slcz
    integer,allocatable,dimension(:,:)  :: sidx,sidy,sidz
    real,allocatable,dimension(:,:)     :: scfx,scfy,scfz
    real,allocatable,dimension(:,:)     :: mslx,msly,mslz

end module movie_arrays

module mpih

    implicit none

    include 'mpif.h'

    integer, parameter  :: master = 0   ! For param logical variable ismaster

    integer             :: mpi_ierr     ! To redirect all MPI error codes

    integer             :: mpi_rank     ! MPI rank
    integer             :: mpi_comm     ! MPI comm 
    integer             :: mpi_size     ! MPI comm size
    integer             :: mpi_dims(2)  ! MPI cartesian topology dimensions
    logical             :: mpi_perc(2)  ! MPI cartesian topology periodicity
    integer             :: mpi_pidx(2)  ! MPI cartesian topology process coordinate indices
    integer             :: mpi_nbrm     ! MPI rank at cartesian topology neighbour previous row
    integer             :: mpi_nbrp     ! MPI rank at cartesian topology neighbour next row
    integer             :: mpi_nbcm     ! MPI rank at cartesian topology neighbour previous column
    integer             :: mpi_nbcp     ! MPI rank at cartesian topology neighbour next column

    integer             :: comm_xct     ! MPI subcommunicator for writing X movie slices
    integer             :: comm_yct     ! MPI subcommunicator for writing Y movie slices
    integer             :: comm_zct     ! MPI subcommunicator for writing Z movie slices

end module mpih

module tridiagonal_buffer

    implicit none

    integer                         :: tdsz             ! Size of Thomas work buffers
    real, allocatable, dimension(:) :: twk1, twk2     ! Inverse pivots / Sherman-Morrison vector

end module tridiagonal_buffer

! FFTW interface (Fortran 2003, calls the FFTW C library directly).
! fftw3.f03 is installed with FFTW, next to fftw3.h.
module fftw_params

    use, intrinsic :: iso_c_binding

    implicit none

    include 'fftw3.f03'

end module fftw_params

module pressure_solver

    use fftw_params

    implicit none

    integer                             :: hdir     ! The direction of helmholtz solver (1=X, 2=Y, 3=Z)

    real, allocatable, dimension(:,:)   :: qfwx, qfwy, qfwz ! Forward eigenvector arrays
    real, allocatable, dimension(:,:)   :: qbwx, qbwy, qbwz ! Backward eigenvector arrays
    real, allocatable, dimension(:)     :: geix,geiy,geiz! Eigenvalues belonging to the eigenvector arrays
    real, allocatable, dimension(:)     :: feix, feiy, feiz    ! Modified wavenumbers: -ak are the eigenvalues for the FFT
    real, allocatable, dimension(:)     :: eigx, eigy, eigz ! Eigenvalues of the selected method (used by the tridiagonal solvers)

    integer                             :: niwk, nrwk       ! Buffer sizes (integer/real)
    integer, allocatable, dimension(:)  :: iwrk             ! Integer buffer for eigenvalue factorization
    real, allocatable, dimension(:)     :: rwrk             ! Real buffer for eigenvalue factorization

    integer                             :: gszx, gszy, gszz ! Sizes of GEMM buffers in X, Y, Z directions
    real, allocatable, dimension(:)     :: gwrk             ! Real buffer for batched GEMM

    integer                             :: pcbcxs, pcbcxe   ! Pressure correction boundary condition start/end in X
    integer                             :: pcbcys, pcbcye   ! Pressure correction boundary condition start/end in Y
    integer                             :: pcbczs, pcbcze   ! Pressure correction boundary condition start/end in Z

    integer                             :: izmx = 0         ! Zero-mode index for x direction (0 = none)
    integer                             :: izmy = 0         ! Zero-mode index for y direction (0 = none)
    integer                             :: izmz = 0         ! Zero-mode index for z direction (0 = none)

    ! Methods to multiply with the eigenvectors of the pressure operator: dense
    ! eigenvector matrices (GEMM, any grid) or real-to-real FFTs (uniform grids)
    integer, parameter                  :: TGEMM = 1, TFFT = 2
    integer                             :: tmet(3)          ! Method per direction (X, Y, Z)
    logical                             :: fftok(3)         ! FFT possible per direction (uniform grid, plans created)
    real                                :: fnrm(3)          ! FFT normalisation per direction
    type(C_PTR)                         :: pfwx = C_NULL_PTR, pfwy = C_NULL_PTR, pfwz = C_NULL_PTR ! FFTW plans: from eigenspace back to physical space
    type(C_PTR)                         :: pbwx = C_NULL_PTR, pbwy = C_NULL_PTR, pbwz = C_NULL_PTR ! FFTW plans: from physical space to eigenspace

    ! FFTW planner flag (change here and recompile, e.g. to FFTW_ESTIMATE)
    integer(C_INT), parameter           :: fftflg = FFTW_MEASURE

end module pressure_solver
