!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE:     ParticleRoutines.F90                       !
!                                                         !
!    CONTAINS: Subroutine CheckIsParticleGlobal           !
!              Subroutine CheckIsParticleLocal            !
!              Subroutine UpdateParticleGridIndices       !
!              Subroutine UpdateParticleAcceleration      !
!              Subroutine UpdateParticleVelocity          !
!              Subroutine UpdateParticlePosition          !
!              Subroutine UpdateParticleLifetime          !
!              Subroutine CalculateParticleExit           !
!              Subroutine CalcApplyParticleDrag           !
!                                                         !
!    PURPOSE:  Subroutines that perform operations on     !
!              single particles                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CheckIsParticleGlobal(p,n) 

    use param, only: xlen,ylen,zlen,periodic
    use lagrangian_point_particle, only: particle_data

    implicit none

    type(particle_data), intent(inout)  :: p
    logical, intent(out)                :: n
    integer                             :: i
    real                                :: nlen(3)

    nlen(1) = xlen
    nlen(2) = ylen
    nlen(3) = zlen

    n = .true.

    do i = 1,3
        n = n.and.(((p%lpp_pos(i).gt.0.0).and.(p%lpp_pos(i).lt.nlen(i))).or.(periodic(i)))
    end do

    return

end subroutine CheckIsParticleGlobal

subroutine CheckIsParticleLocal(p,n) 

    use decomp_2d, only: xstart,xend
    use param, only: xc,yc,zc,xlen,ylen,zlen,periodic
    use lagrangian_point_particle, only: particle_data

    implicit none

    type(particle_data), intent(in) :: p
    logical, intent(out)            :: n
    integer                         :: i
    real                            :: nlen(3),ncst(3),ncen(3)

    nlen(1) = xlen
    nlen(2) = ylen
    nlen(3) = zlen

    ncst(1) = xc(xstart(1))
    ncst(2) = yc(xstart(2))
    ncst(3) = zc(xstart(3))

    ncen(1) = xc(xend(1)+1)
    ncen(2) = yc(xend(2)+1)
    ncen(3) = zc(xend(3)+1)

    n = .true.

    do i = 1,3
        if ((xstart(i).eq.1).and.(.not.periodic(i))) then
            n = n.and.((p%lpp_pos(i).gt.ncst(i)).and.(p%lpp_pos(i).lt.ncen(i)))
        else
            n = n.and.((p%lpp_pos(i).ge.ncst(i)).and.(p%lpp_pos(i).lt.ncen(i)))
        end if
    end do

    return

end subroutine CheckIsParticleLocal

subroutine UpdateParticleGridIndices(p) 

    use decomp_2d, only: xstart,xend
    use param, only: xc,xm,yc,ym,zc,zm
    use lagrangian_point_particle, only: particle_data

    implicit none

    type(particle_data), intent(inout)  :: p

    ! Find and update cell indices for c-grid (nodes)
    call GetLocationCellIndex(p%lpp_pos(1),xc(xstart(1)  :xend(1)+1),xstart(1)  ,xend(1)+1,p%grc_idx(1))
    call GetLocationCellIndex(p%lpp_pos(2),yc(xstart(2)  :xend(2)+1),xstart(2)  ,xend(2)+1,p%grc_idx(2))
    call GetLocationCellIndex(p%lpp_pos(3),zc(xstart(3)  :xend(3)+1),xstart(3)  ,xend(3)+1,p%grc_idx(3))
    ! Find and update cell indices for m-grid (cell-center)
    call GetLocationCellIndex(p%lpp_pos(1),xm(xstart(1)-1:xend(1)+1),xstart(1)-1,xend(1)+1,p%grm_idx(1))
    call GetLocationCellIndex(p%lpp_pos(2),ym(xstart(2)-1:xend(2)+1),xstart(2)-1,xend(2)+1,p%grm_idx(2))
    call GetLocationCellIndex(p%lpp_pos(3),zm(xstart(3)-1:xend(3)+1),xstart(3)-1,xend(3)+1,p%grm_idx(3))

    return

end subroutine UpdateParticleGridIndices 

subroutine UpdateParticleAcceleration(p) 

    use param
    use local_arrays, only: vx,vy,vz
    use lagrangian_point_particle, only: particle_data,lpp_grav

    implicit none

    type(particle_data), intent(inout)  :: p

    integer                             :: im,jm,km
    integer                             :: ip,jp,kp

    real                                :: cfxm,cfym,cfzm
    real                                :: cfxp,cfyp,cfzp
    real                                :: slvx,slvy,slvz
    real                                :: svel,srey,pcfd,cvol

    ! Store old acceleration
    p%acc_old(1) = p%acc_now(1)
    p%acc_old(2) = p%acc_now(2)
    p%acc_old(3) = p%acc_now(3)

    ! Initialize current acceleration to zero
    p%acc_now(1) = 0.0
    p%acc_now(2) = 0.0
    p%acc_now(3) = 0.0

    ! Add fluid drag
    call CalcApplyParticleDrag(p)

    ! Add gravitational forces
    p%acc_now(1) = p%acc_now(1) + (1.0 - 1.0/p%lpp_den)*lpp_grav(1)
    p%acc_now(2) = p%acc_now(2) + (1.0 - 1.0/p%lpp_den)*lpp_grav(2)
    p%acc_now(3) = p%acc_now(3) + (1.0 - 1.0/p%lpp_den)*lpp_grav(3)

    ! Add electrostatic forces if necesssary and update ... an eulerian electric field? idk ¯\_(ツ)_/¯ 

    return

end subroutine UpdateParticleAcceleration

subroutine UpdateParticleVelocity(p) 

    use param, only: dt,ga,ro
    use lagrangian_point_particle, only: particle_data

    implicit none

    type(particle_data), intent(inout)  :: p
    integer                             :: i

    ! Update velocity at next sub step using RK-3 time-stepping scheme
    do i = 1,3
        p%lpp_vel(i) = p%lpp_vel(i) + dt*(ga*p%acc_now(i) + ro*p%acc_old(i))
    end do

    return

end subroutine UpdateParticleVelocity

subroutine UpdateParticlePosition(p) 

    use param, only: dt,al,ga,ro
    use lagrangian_point_particle, only: particle_data

    implicit none

    type(particle_data), intent(inout)  :: p
    integer                             :: i

    ! Update position at next sub-step assuming constant acceleration during sub-step (equivalent to Crank-Nicolson)
    do i = 1,3
        p%lpp_pos(i) = p%lpp_pos(i) + al*dt*(p%lpp_vel(i) - 0.5*dt*(ga*p%acc_now(i) + ro*p%acc_old(i)))
    end do

    return

end subroutine UpdateParticlePosition

subroutine UpdateParticleLifeTime(p) 

    use param, only: dt,al
    use lagrangian_point_particle, only: particle_data

    implicit none

    type(particle_data), intent(inout)  :: p

    ! Update particle life time
    p%lpp_lft = p%lpp_lft + al*dt

end subroutine UpdateParticleLifeTime

subroutine CalculateParticleExit(p,q)

    use param
    use lagrangian_point_particle, only: particle_data,particle_exit

    implicit none

    type(particle_data), intent(in)     :: p
    type(particle_exit), intent(out)    :: q
    integer                             :: ndir,edir
    real                                :: pold(3),pnew(3),vold(3),vnew(3),accl(3)
    real                                :: sdet,den1,den2,stme,tcrs
    real                                :: nlen(3)
    character(len=4)                    :: caxs

    ! For compact processing over 3D space
    nlen(1) = xlen
    nlen(2) = ylen
    nlen(3) = zlen

    ! Useful for converting the integer direction from CalcParticleExit to a character value
    caxs = "XYZ!"
    ! Set default such that in case there's a problem the exit plane shows up as "!"
    edir = 4

    ! Initialize the shortest cross time values to full sub-step
    stme = al*dt
    tcrs = al*dt
    ! Initialize the particle lifetime at exit to the value just before the sub-step advance
    q%pex_lft = p%lpp_lft - stme
    ! Initialize the flow time at exit to the value just before the sub-step advance
    q%pex_eft = time - stme

    do ndir = 1,3
        ! Unpack acceleration, old velocity, and old position
        accl(ndir) = (ga*p%acc_now(ndir) + ro*p%acc_old(ndir))/al
        vnew(ndir) = p%lpp_vel(ndir)
        vold(ndir) = p%lpp_vel(ndir) - al*dt*accl(ndir)
        pnew(ndir) = p%lpp_pos(ndir)
        pold(ndir) = p%lpp_pos(ndir) - al*dt*(p%lpp_vel(ndir) - 0.5*al*dt*accl(ndir))
        ! Only check crossing if not periodic in that direction
        if (.not.periodic(ndir)) then
            !! Check crossing at lower end of axis
            if (((pnew(ndir) - 0.0)*(pold(ndir) - 0.0)).le.0.0) then
                sdet = (vold(ndir)**2.0) + (2.0*accl(ndir)*(0.0 - pold(ndir)))
                if (sdet.ge.0.0) then
                    den1 = (vold(ndir) + sqrt(sdet))
                    den2 = (vold(ndir) - sqrt(sdet))
                    if (abs(den1).gt.abs(den2)) then
                        tcrs = 2.0*(0.0 - pold(ndir))/den1
                    else
                        tcrs = 2.0*(0.0 - pold(ndir))/den2
                    end if
                end if
            end if
            !! Check crossing at higher end of axis
            if (((pnew(ndir) - nlen(ndir))*(pold(ndir) - nlen(ndir))).le.0.0) then
                sdet = (vold(ndir)**2.0) + (2.0*accl(ndir)*(nlen(ndir) - pold(ndir)))
                if (sdet.ge.0.0) then
                    den1 = (vold(ndir) + sqrt(sdet))
                    den2 = (vold(ndir) - sqrt(sdet))
                    if (abs(den1).gt.abs(den2)) then
                        tcrs = 2.0*(nlen(ndir) - pold(ndir))/den1
                    else
                        tcrs = 2.0*(nlen(ndir) - pold(ndir))/den2
                    end if
                end if
            end if
            !! Crossing with shortest time must be the first exit plane
            if (tcrs.lt.stme) then
                stme = tcrs
                edir = ndir
            end if
        end if
    end do

    ! Set the source index
    q%src_idx = p%src_idx
    ! Get the axis of exit
    q%pex_pln = caxs(edir:edir)
    ! Add the shortest cross time to the pre sub-step particle lifetime to get full particle lifetime
    q%pex_lft = q%pex_lft + stme
    ! Add the shortest cross time to the pre sub-step flow time to get exit flow time
    q%pex_eft = q%pex_eft + stme

    ! Update the exit position and velocity from the shortest cross time
    do ndir = 1,3
        q%pex_vel(ndir) = vold(ndir) + accl(ndir)*stme
        q%pex_pos(ndir) = pold(ndir) + stme*(vold(ndir) + 0.5*stme*accl(ndir))
    end do
    
    return

end subroutine CalculateParticleExit

subroutine CalcApplyParticleDrag(p)

    use decomp_2d
    use param
    use local_arrays, only: vx,vy,vz
    use lagrangian_point_particle

    implicit none

    type(particle_data), intent(inout)  :: p
    
    integer                             :: i,j,k

    real                                :: cffc(3)
    real                                :: cffm(3)
    
    real                                :: slip(3)
    real                                :: curv(3)
    real                                :: coef(3)
    real                                :: forc(3)
    real                                :: pvol,pcfd

    ! Compute interpolation coefficients
    call CalcTrilinearInterpolationCoefficients(p,cffm,cffc)

    ! Compute slip velocity
    !! Initialize slip velocity to zero
    slip(:) = 0.0
    !! Apply trilinear interpolation with previously calculated coefficients to get fluid velocity
    call ApplyTrilinearInterpolation(p,cffm,cffc,'x',1,vx,slip(1))
    call ApplyTrilinearInterpolation(p,cffm,cffc,'y',1,vy,slip(2))
    call ApplyTrilinearInterpolation(p,cffm,cffc,'z',1,vz,slip(3))

    !! If slip correction is enabled, correct the slip velocity to account for self-induced velocity
    if (lpp_scor) then
        
        !!! Initialize curvatures to zero
        curv(:) = 0.0
        !!! Compute the curvature at the particle using trilinear interpolation
        call ApplyTrilinearInterpolation(p,cffm,cffc,'x',1,lpp_d2vx,curv(1))
        call ApplyTrilinearInterpolation(p,cffm,cffc,'y',1,lpp_d2vy,curv(2))
        call ApplyTrilinearInterpolation(p,cffm,cffc,'z',1,lpp_d2vz,curv(3))
        !!! Compute the correction coefficient using trilinear interpolation on lookup data array 
        !!! Note: Reynolds number calculated at previous substep is being used to avoid complicated implicit calculation
        call CalcSlipCorrectionCoefficient(p,'x',coef(1))
        call CalcSlipCorrectionCoefficient(p,'y',coef(2))
        call CalcSlipCorrectionCoefficient(p,'z',coef(3))
        !!! Compute the slip correction using coefficients and curvatures
        slip(1) = slip(1) + coef(1)*curv(1)
        slip(2) = slip(2) + coef(2)*curv(2)
        slip(3) = slip(3) + coef(3)*curv(3)

    end if

    !! Subtract fluid velocity from particle velocity to get slip velocity
    slip(1) = p%lpp_vel(1) - slip(1)
    slip(2) = p%lpp_vel(2) - slip(2)
    slip(3) = p%lpp_vel(3) - slip(3)

    ! Compute particle volume
    pvol = pi*(p%lpp_dia**3.0)/6.0
    
    ! Compute paricle Reynolds number using slip velocity
    p%lpp_rey = rey*p%lpp_dia*norm2(slip)

    ! Calculate drag force and acceleration
    !! Compute premultiplied drag coefficient using empirical drag models (Stokes/Schiller-Naumann/Morsi-Alexander)
    call PremultipliedDragCoefficient(p,pcfd)
    !! Compute drag force using premultiplied drag coefficient
    forc(1) = (pi*slip(1)*pcfd*p%lpp_dia)/(8.0*rey)
    forc(2) = (pi*slip(2)*pcfd*p%lpp_dia)/(8.0*rey)
    forc(3) = (pi*slip(3)*pcfd*p%lpp_dia)/(8.0*rey)

    ! Compute particle accelerations
    p%acc_now(1) = p%acc_now(1) - e2l_mult*forc(1)/(p%lpp_den*pvol)
    p%acc_now(2) = p%acc_now(2) - e2l_mult*forc(2)/(p%lpp_den*pvol)
    p%acc_now(3) = p%acc_now(3) - e2l_mult*forc(3)/(p%lpp_den*pvol)

    ! Apply opposite drag force on fluid (from Newton's third law)
    !! Apply trilinear interpolation with previously calculated coefficients to get body forces
    call ApplyTrilinearInterpolation(p,cffm,cffc,'x',-1,lpp_bdfx,(l2e_mult*forc(1)))
    call ApplyTrilinearInterpolation(p,cffm,cffc,'y',-1,lpp_bdfy,(l2e_mult*forc(2)))
    call ApplyTrilinearInterpolation(p,cffm,cffc,'z',-1,lpp_bdfz,(l2e_mult*forc(3)))
    
    return

end subroutine CalcApplyParticleDrag