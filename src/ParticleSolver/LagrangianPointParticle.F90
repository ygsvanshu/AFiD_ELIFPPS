!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: LagrangianPointParticle.F90                    !
!    CONTAINS: module lagrangian_point_particle           !
!                                                         !
!    PURPOSE: The main module for the lagrangian point    !
!    particle model. Contains derived datatypes for       !
!    particles and sources, as well as global variables   !
!    and arrays                                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module lagrangian_point_particle

    implicit none

    ! Main derived datatype for the particle data
    type particle_data
        integer             :: src_idx      ! Index of source from which particle has been injected
        integer             :: grc_idx(3)   ! Grid cell indices on the [xm,ym,zm] grid
        integer             :: grm_idx(3)   ! Grid cell indices on the [xc,yc,zc] grid
        real                :: lpp_lft      ! Life time of particle since injection
        real                :: lpp_dia      ! Diameter of particle
        real                :: lpp_den      ! Density of particle
        real                :: lpp_rey      ! Reynolds number computed from slip velocity and particle diameter
        real                :: lpp_pos(3)   ! Position at next step
        real                :: lpp_vel(3)   ! Velocity at next step
        real                :: acc_old(3)   ! Acceleration at previous step
        real                :: acc_now(3)   ! Acceleration at current step
    end type particle_data

    ! Derived datatype for the particle source
    type particle_source
        integer             :: src_idx      ! Source index (in the order as entered in the file "particle_soures.in")
        integer             :: grc_idx(3)   ! Grid cell indices on the [xm,ym,zm] grid
        integer             :: grm_idx(3)   ! Grid cell indices on the [xc,yc,zc] grid
        real                :: src_sta      ! Injection start time instance of the source
        real                :: src_end      ! Injection stop time instance of the source
        real                :: src_frq      ! Injection frequency of the source
        real                :: src_dia      ! Diameter of injected particles
        real                :: src_den      ! Density of injected particles
        real                :: src_pos(3)   ! Position of the source
        real                :: src_vel(3)   ! Velocity of injection at the source
    end type particle_source

    ! Derived datatype for the particle exit event
    type particle_exit
        integer             :: src_idx      ! Source index (in the order as entered in the file "particle_soures.in")
        character           :: pex_pln      ! Plane of exit (X/Y/Z)
        real                :: pex_lft      ! Life time of the particle at exit
        real                :: pex_eft      ! Flow time at particle exit
        real                :: pex_pos(3)   ! Position of the particle at exit
        real                :: pex_vel(3)   ! Velocity of the particle at exit
    end type particle_exit

    ! Drag coefficient models

    integer, parameter                                  :: STOKES  =  0 ! Stokes drag model (for Reynolds number << 1)
    integer, parameter                                  :: SCHNAU  =  1 ! Schiller-Naumann drag model (for All Reynolds numbers)

    ! Global input parameters

    integer                                             :: lpp_dmod     ! Drag model
    logical                                             :: lpp_scor     ! Slip correction
    real                                                :: lpp_grav(3)  ! Acceleration due to gravity
    real                                                :: e2l_mult     ! Multiplier for Eulerian to Langrangian coupling
    real                                                :: l2e_mult     ! Multiplier for Langrangian to Eulerian coupling
    real                                                :: lpp_stol     ! Spawn tolerance interval 
    real                                                :: lpp_clim     ! Limit on the CFL analogue (to limit timestep)
    real                                                :: lpp_tlim     ! Limit on the characteristic response time (to limit timestep)
    logical                                             :: lpp_save     ! Save particle snapshots
    real                                                :: lpp_ssta     ! Start time for saving particle snapshots 
    real                                                :: lpp_sfrq     ! Save time interval for particle snapshots
    logical                                             :: pex_save     ! Save particle exit events
    real                                                :: pex_ssta     ! Start time for saving particle exit events 
    real                                                :: pex_sfrq     ! Save time interval for particle exit events


    ! MPI related module parameters

    integer                                             :: mpi_pdat     ! MPI type for particle data

    ! Statistics

    integer                                             :: tot_spwn     ! Total number of spawned particles since the beginning of the simulation (in all pencils/processes)
    integer                                             :: tot_actv     ! Total number of active particles at the current time (in all pencils/processes)
    integer                                             :: tot_exit     ! Total number of exited particles since the beginning of the simulation (in all pencils/processes)
    real                                                :: lpp_vmax(3)  ! Maximum particle velocity along each cardinal direction
    real                                                :: lpp_vavg(3)  ! Average particle velocity along each cardinal direction
    real                                                :: lpp_vrms(3)  ! RMS particle velocity along each cardinal direction
    real                                                :: lpp_bmax(3)  ! Maximum body force along each cardinal direction
    real                                                :: lpp_bavg(3)  ! Average body force along each cardinal direction
    real                                                :: lpp_brms(3)  ! RMS body force along each cardinal direction

    ! Counts for sources, particles, spawn events, and exit events for current pencil/process

    integer                                             :: src_ntot     ! Total number of sources 
    integer                                             :: src_size     ! Size of buffers required for storing the particle sources for current pencil/process
    integer                                             :: lpp_actv     ! Number of active particles in the buffer lpp_list for current pencil/process
    integer                                             :: pex_actv     ! Number of active exit events in the buffer pex_list for current pencil/process
    integer                                             :: lpp_spwn     ! Number of spawned particles for current pencil/process
    integer                                             :: lpp_exit     ! Number of exited particles for current pencil/process
    integer                                             :: sub_exit     ! Number of exited particles at current substep for current pencil/process
    integer                                             :: lpp_snap     ! The number of particle history snapshots 

    ! Source and particle arrays

    type(particle_data), allocatable, dimension(:)      :: lpp_list     ! Main buffer for all particles stored per pencil/process
    type(particle_data), allocatable, dimension(:)      :: bfm_send     ! Send buffer for particles on the start side 
    type(particle_data), allocatable, dimension(:)      :: bfp_send     ! Send buffer for particles on the end side 

    type(particle_source), allocatable, dimension(:)    :: src_list     ! Buffer to store all particle sources that can add to the pencil/process

    type(particle_exit), allocatable, dimension(:)      :: pex_list     ! Main buffer for all exit events stored per pencil/process

    ! Send and recieve buffers for halo force exchange

    real,allocatable,dimension(:,:,:)                   :: snd_nbrm     ! Send buffer for halo force update at cartesian topology neighbour previous row
    real,allocatable,dimension(:,:,:)                   :: snd_nbrp     ! Send buffer for halo force update at cartesian topology neighbour next row
    real,allocatable,dimension(:,:,:)                   :: snd_nbcm     ! Send buffer for halo force update at cartesian topology neighbour previous column
    real,allocatable,dimension(:,:,:)                   :: snd_nbcp     ! Send buffer for halo force update at cartesian topology neighbour next column

    real,allocatable,dimension(:,:,:)                   :: rcv_nbrm     ! Receive buffer for halo force update at cartesian topology neighbour previous row
    real,allocatable,dimension(:,:,:)                   :: rcv_nbrp     ! Receive buffer for halo force update at cartesian topology neighbour next row
    real,allocatable,dimension(:,:,:)                   :: rcv_nbcm     ! Receive buffer for halo force update at cartesian topology neighbour previous column
    real,allocatable,dimension(:,:,:)                   :: rcv_nbcp     ! Receive buffer for halo force update at cartesian topology neighbour next column

    ! Body force arrays

    real,allocatable,dimension(:,:,:)                   :: lpp_bdfx     ! Body force array in X direction
    real,allocatable,dimension(:,:,:)                   :: lpp_bdfy     ! Body force array in Y direction
    real,allocatable,dimension(:,:,:)                   :: lpp_bdfz     ! Body force array in Z direction

    ! Curvature terms

    real,allocatable,dimension(:,:,:)                   :: lpp_d2vx     ! Velocity curvature array in X direction (d2vx/dx2 + d2vx/dy2 + d2vx/dz2)
    real,allocatable,dimension(:,:,:)                   :: lpp_d2vy     ! Velocity curvature array in Y direction (d2vy/dx2 + d2vy/dy2 + d2vy/dz2)
    real,allocatable,dimension(:,:,:)                   :: lpp_d2vz     ! Velocity curvature array in Z direction (d2vz/dx2 + d2vz/dy2 + d2vz/dz2)

    ! Slip correction data arrays

    integer                                             :: num_aspc     ! Number of data points for aspect ratio (ratio of particle diameter to grid-cell dimension)
    integer                                             :: num_reyn     ! Number of data points for particle Reynolds number
    integer                                             :: num_data     ! Number of data points for the slip correction coefficients and timescales
    real, allocatable, dimension(:)                     :: dat_aspc     ! Slip correction data coordinates for aspect ratio
    real, allocatable, dimension(:)                     :: dat_reyn     ! Slip correction data coordinates for Reynolds number
    real, allocatable, dimension(:,:,:,:,:)             :: dat_coef     ! Slip correction coefficients

end module lagrangian_point_particle