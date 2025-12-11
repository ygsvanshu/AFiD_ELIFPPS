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
!              module pressure_decomp                     !
!              module fftw_params                         !
!                                                         ! 
!    PURPOSE: All modules required for the flow solver    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Declaration of global variables
module param

    implicit none
        
    ! Read from input file solver.in

    ! 1XX -----------------------
    logical         :: nread 
    ! 2XX -----------------------
    integer         :: nxm, nym, nzm
    real            :: xlen,ylen,zlen
    integer         :: straxs
    integer         :: strtyp
    real            :: strval
    logical         :: periodic(3) = (/.false.,.false.,.false./)
    ! 3XX -----------------------
    integer         :: ntst
    real            :: tmax
    real            :: walltimemax
    character*200   :: abortfile
    ! 4XX -----------------------
    integer         :: nsst 
    logical         :: vardt = .true.
    real            :: dt
    real            :: dtmin,dtmax
    real            :: limitCFL
    real            :: resid
    real            :: limitVel
    integer         :: epsnum
    real            :: eps
    ! 5XX -----------------------
    real            :: rey
    ! 6XX -----------------------
    integer         :: nout
    logical         :: save1d,save2d,save3d
    real            :: tsta1d,tsta2d,tsta3d
    real            :: freq1d,freq2d,freq3d

    ! Other global variables

    integer                 :: nx, ny, nz
    integer                 :: ntime,ns

    real                    :: time
    real                    :: pi = 2.0*asin(1.0)
    real                    :: al,ga,ro

    real                    :: cvel

    real, dimension(1:3)    :: vmax
    real, dimension(1:3)    :: vavg
    real, dimension(1:3)    :: vrms
    real                    :: vmag

    real, dimension(1:3)    :: gam,rom,alm
              
    logical                 :: ismaster = .false.
    integer                 :: lvlhalo = 1
    integer                 :: mpi_xcut,mpi_ycut,mpi_zcut

    logical                 :: particle = .false.
    logical                 :: pre_diff = .false.

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

    ! Variables for FFTW and Poisson solver

    real, allocatable, dimension(:)     :: ak1,ao
    real, allocatable, dimension(:)     :: ak2,ap
    real, allocatable, dimension(:)     :: ak3,aq
        
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
    ! Arrays for boundary condition coefficients
    real,dimension(3,3)                     :: coefxs,coefxe
    real,dimension(3,3)                     :: coefys,coefye
    real,dimension(3,3)                     :: coefzs,coefze
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
    ! Arrays to store whether global flux balance is allowed
    logical,allocatable,dimension(:,:)      :: glpcxs,glpcxe
    logical,allocatable,dimension(:,:)      :: glpcys,glpcye
    logical,allocatable,dimension(:,:)      :: glpczs,glpcze

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

module pressure_decomp

    use decomp_2d

    implicit none

    type(decomp_info)                   :: ph_info_n
    type(decomp_info)                   :: ph_info_x
    type(decomp_info)                   :: ph_info_y
    type(decomp_info)                   :: ph_info_z

end module pressure_decomp

module fftw_params

    use iso_c_binding

    type, bind(C) :: fftw_iodim
        integer(C_INT) n, is, os
    end type fftw_iodim

    integer, parameter :: C_FFTW_R2R_KIND =  C_INT32_T
    integer, parameter :: FFTW_PATIENT    =  32
    integer, parameter :: FFTW_MEASURE    =  0
    integer, parameter :: FFTW_ESTIMATE   =  64   
    integer, parameter :: FFTW_FORWARD    = -1   
    integer, parameter :: FFTW_BACKWARD   =  1
    integer, parameter :: FFTW_REDFT01    =  4
    integer, parameter :: FFTW_REDFT10    =  5
    integer, parameter :: FFTW_R2HC       =  0
    integer, parameter :: FFTW_HC2R       =  1

    interface

        type(C_PTR) function fftw_plan_guru_dft(rank,dims,howmany_rank,howmany_dims,in,out,sign,flags) bind(C, name='fftw_plan_guru_dft')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: sign
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft

        type(C_PTR) function fftw_plan_guru_dft_r2c(rank,dims,howmany_rank,howmany_dims,in,out,flags) bind(C, name='fftw_plan_guru_dft_r2c')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_r2c
        
        type(C_PTR) function fftw_plan_guru_dft_c2r(rank,dims,howmany_rank,howmany_dims,in,out,flags) bind(C, name='fftw_plan_guru_dft_c2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            complex(C_DOUBLE_COMPLEX), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_INT), value :: flags
        end function fftw_plan_guru_dft_c2r

        type(C_PTR) function fftw_plan_guru_r2r(rank,dims,howmany_rank,howmany_dims,in,out,kind,flags) bind(C, name='fftw_plan_guru_r2r')
            import
            integer(C_INT), value :: rank
            type(fftw_iodim), dimension(*), intent(in) :: dims
            integer(C_INT), value :: howmany_rank
            type(fftw_iodim), dimension(*), intent(in) :: howmany_dims
            real(C_DOUBLE), dimension(*), intent(out) :: in
            real(C_DOUBLE), dimension(*), intent(out) :: out
            integer(C_FFTW_R2R_KIND), dimension(*), intent(in) :: kind
            integer(C_INT), value :: flags
        end function fftw_plan_guru_r2r

    end interface

    type(C_PTR) :: fwd_guruplan_x,fwd_guruplan_y,fwd_guruplan_z
    type(C_PTR) :: bwd_guruplan_x,bwd_guruplan_y,bwd_guruplan_z
    
    logical :: planned=.false.

    real   ,allocatable,dimension(:,:,:) :: rx1,ry1,rz1
    complex,allocatable,dimension(:,:,:) :: cx1,cy1,cz1

end module fftw_params