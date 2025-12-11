!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: SlipCorrectionRoutines.F90                     !
!                                                         !
!    CONTAINS: subroutine InitSlipCorrectionData          !
!              subroutine CalcSlipCorrectionCoefficient   !
!              subroutine CalcCurvatureTerms              !
!                                                         !
!    PURPOSE: Subroutine InitSlipCorrection reads in and  !
!    initializes the coefficients and rise times used to  !
!    correct the self-induced velocity in the             ! 
!    slip velocity computed by interpolation.             !
!                                                         !
!    Subroutine CalcSlipCorrectionCoefficient computes    !
!    coefficients that relate the velocity curvature to   !
!    the slip correction using multilinear interpolation  !
!    of the slip correction data.                         !
!                                                         !
!    Subroutine CalcCurvatureTerms computes the velocity  !
!    curvature flow field for all three components of     !
!    fluid velocity.                                      !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitSlipCorrectionData

    use, intrinsic :: iso_fortran_env, only : iostat_end
    use param, only: ismaster
    use lagrangian_point_particle

    implicit none

    logical             :: exists
    integer             :: istt,nlen
    integer             :: ierror,nlines
    integer             :: indl(5)
    character(len=500)  :: afidpath,datapath
    real                :: valc(8)
    
    ! Get the path variable for the solver source
    call get_environment_variable("AFID_PATH",value=afidpath,status=istt,trim_name=.true.)
    if (istt.ne.0) then
        if (ismaster) write (6,*) 'Error: Environment variable AFID_PATH not set or not found'
        call MpiAbort
    end if
    ! Check if the solver path ends with slash. Correspondingly update the path to slip correction data
    nlen = len_trim(afidpath)
    if (afidpath(nlen:nlen).eq.'/') then
        datapath = trim(afidpath)//"data/slip_correction.dat"
    else
        datapath = trim(afidpath)//"/data/slip_correction.dat"
    end if
    ! Check that the data file exists
    inquire(file=datapath, exist=exists)
    if (exists) then
        !! Open the file
        open(99, file=trim(datapath))
        !! Read the size of the slip correction data aspect ratios
        nlines = 1
        ierror = 0
        do while ((ierror.ne.iostat_end).and.(nlines.le.1))
            read(99, *, iostat=ierror) num_aspc
            if (ierror.eq.0) nlines = nlines + 1
        end do
        !! Read the size of the slip correction data Reynolds numbers
        nlines = 1
        ierror = 0
        do while ((ierror.ne.iostat_end).and.(nlines.le.1))
            read(99, *, iostat=ierror) num_reyn
            if (ierror.eq.0) nlines = nlines + 1
        end do
        !! Set the fixed number of data values and drag models
        num_data = 8 !!!! (8 coefficients at corners of the cell and 8 rise-time values)
        !! Initialize the slip correction arrays
        if (.not.allocated(dat_aspc)) allocate(dat_aspc(num_aspc))
        if (.not.allocated(dat_reyn)) allocate(dat_reyn(num_reyn))
        if (.not.allocated(dat_coef)) allocate(dat_coef(num_aspc,num_aspc,num_aspc,num_reyn,num_data))
        !! Read the aspect ratios of the slip correction data
        nlines = 1
        ierror = 0
        do while ((ierror.ne.iostat_end).and.(nlines.le.num_aspc))
            read(99, *, iostat=ierror) dat_aspc(nlines)
            if (ierror.eq.0) nlines = nlines + 1
        end do
        !! Read the Reynolds numbers of the slip correction data
        nlines = 1
        ierror = 0
        do while ((ierror.ne.iostat_end).and.(nlines.le.num_reyn))
            read(99, *, iostat=ierror) dat_reyn(nlines)
            if (ierror.eq.0) nlines = nlines + 1
        end do
        !! Read the slip correction coefficients and rise time
        ierror = 0
        nlines = 0
        do while (ierror.ne.iostat_end)
            read(99, *, iostat=ierror) indl,valc
            !!! Check for valid data for the given drag model
            if ((ierror.eq.0).and.(indl(1).eq.lpp_dmod)) then
                dat_coef(indl(2),indl(3),indl(4),indl(5),1:8) = valc
                nlines = nlines + 1
            end if
        end do
        !! Close file
        close(99)
        !! If all required data has not been read, abort as the data may be incomplete/corrupted/incorrect/incompatible 
        if (nlines.ne.(num_aspc*num_aspc*num_aspc*num_reyn)) then
            if (ismaster) write (6,*) 'Error: Slip correction data not incomplete/corrupted at ',datapath
            call MpiAbort
        end if
    else
        !! Slip correection data doesn't exist, abort solver.
        if (ismaster) write (6,*) 'Error: Slip correction data not found at ',datapath
        call MpiAbort
    end if

end subroutine InitSlipCorrectionData

subroutine CalcSlipCorrectionCoefficient(p,grid,coef)

    use param, only: xc,xm,yc,ym,zc,zm,rey
    use lagrangian_point_particle

    implicit none

    type(particle_data), intent(inout)  :: p
    character, intent(in)               :: grid
    real,      intent(out)              :: coef

    integer                             :: i,j,k,l,m
    integer                             :: icrd(4)
    real                                :: aspc(3)
    real                                :: aspf(3)
    real                                :: fcrd(4,0:1)
    real                                :: cval,mtau
    real                                :: aspm(8)
    real                                :: reyn

    ! Compute the grid interpolation multipliers/coefficients depending on which face the velocity component lies on
    if (grid.eq.'x') then

        aspc(1) = p%lpp_dia/(xc(p%grc_idx(1)+1) - xc(p%grc_idx(1)))
        aspc(2) = p%lpp_dia/(ym(p%grm_idx(2)+1) - ym(p%grm_idx(2)))
        aspc(3) = p%lpp_dia/(zm(p%grm_idx(3)+1) - zm(p%grm_idx(3)))

        aspf(1) = abs(p%lpp_pos(1) - (0.5*(xc(p%grc_idx(1)+1) + xc(p%grc_idx(1)))))/(0.5*(xc(p%grc_idx(1)+1) - xc(p%grc_idx(1))))
        aspf(2) = abs(p%lpp_pos(2) - (0.5*(ym(p%grm_idx(2)+1) + ym(p%grm_idx(2)))))/(0.5*(ym(p%grm_idx(2)+1) - ym(p%grm_idx(2))))
        aspf(3) = abs(p%lpp_pos(3) - (0.5*(zm(p%grm_idx(3)+1) + zm(p%grm_idx(3)))))/(0.5*(zm(p%grm_idx(3)+1) - zm(p%grm_idx(3))))

    else if (grid.eq.'y') then

        aspc(1) = p%lpp_dia/(yc(p%grc_idx(2)+1) - yc(p%grc_idx(2)))
        aspc(2) = p%lpp_dia/(zm(p%grm_idx(3)+1) - zm(p%grm_idx(3)))
        aspc(3) = p%lpp_dia/(xm(p%grm_idx(1)+1) - xm(p%grm_idx(1)))

        aspf(1) = abs(p%lpp_pos(2) - (0.5*(yc(p%grc_idx(2)+1) + yc(p%grc_idx(2)))))/(0.5*(yc(p%grc_idx(2)+1) - yc(p%grc_idx(2))))
        aspf(2) = abs(p%lpp_pos(3) - (0.5*(zm(p%grm_idx(3)+1) + zm(p%grm_idx(3)))))/(0.5*(zm(p%grm_idx(3)+1) - zm(p%grm_idx(3))))
        aspf(3) = abs(p%lpp_pos(1) - (0.5*(xm(p%grm_idx(1)+1) + xm(p%grm_idx(1)))))/(0.5*(xm(p%grm_idx(1)+1) - xm(p%grm_idx(1))))

    else if (grid.eq.'z') then

        aspc(1) = p%lpp_dia/(zc(p%grc_idx(3)+1) - zc(p%grc_idx(3)))
        aspc(2) = p%lpp_dia/(xm(p%grm_idx(1)+1) - xm(p%grm_idx(1)))
        aspc(3) = p%lpp_dia/(ym(p%grm_idx(2)+1) - ym(p%grm_idx(2)))

        aspf(1) = abs(p%lpp_pos(3) - (0.5*(zc(p%grc_idx(3)+1) + zc(p%grc_idx(3)))))/(0.5*(zc(p%grc_idx(3)+1) - zc(p%grc_idx(3))))
        aspf(2) = abs(p%lpp_pos(1) - (0.5*(xm(p%grm_idx(1)+1) + xm(p%grm_idx(1)))))/(0.5*(xm(p%grm_idx(1)+1) - xm(p%grm_idx(1))))
        aspf(3) = abs(p%lpp_pos(2) - (0.5*(ym(p%grm_idx(2)+1) + ym(p%grm_idx(2)))))/(0.5*(ym(p%grm_idx(2)+1) - ym(p%grm_idx(2))))

    end if

    ! Compute the grid interpolation multipliers/coefficients for the relevant sub-cell 
    aspm(1) = aspf(1)*aspf(2)*aspf(3)
    aspm(2) = (1.0 - aspf(1))*aspf(2)*aspf(3)
    aspm(3) = aspf(1)*(1.0 - aspf(2))*aspf(3)
    aspm(4) = (1.0 - aspf(1))*(1.0 - aspf(2))*aspf(3)
    aspm(5) = aspf(1)*aspf(2)*(1.0 - aspf(3))
    aspm(6) = (1.0 - aspf(1))*aspf(2)*(1.0 - aspf(3))
    aspm(7) = aspf(1)*(1.0 - aspf(2))*(1.0 - aspf(3))
    aspm(8) = (1.0 - aspf(1))*(1.0 - aspf(2))*(1.0 - aspf(3))

    ! Find out the aspect ratio indices encompassing the current aspect ratio of the particle in its cell
    do i = 1,3
        !! Check if aspect ratio lies within data coordinates
        if (aspc(i).lt.dat_aspc(1)) then
            !!! Extrapolate on lower end
            icrd(i) = 1
        else if (aspc(i).ge.dat_aspc(num_aspc)) then
            !!! Extrapolate on higher end
            icrd(i) = num_aspc-1
        else
            !!! Find the interval for interpolation
            do j = 1,num_aspc-1
                icrd(i) = j
                if ((dat_aspc(j).le.aspc(i)).and.(dat_aspc(j+1).gt.aspc(i))) exit
            end do
        endif
        !! Calculate the multipliers/coefficients for aspect ratio interpolation
        fcrd(i,0) = (dat_aspc(icrd(i)+1) - aspc(i))/(dat_aspc(icrd(i)+1) - dat_aspc(icrd(i)))
        fcrd(i,1) = (aspc(i)  -  dat_aspc(icrd(i)))/(dat_aspc(icrd(i)+1) - dat_aspc(icrd(i)))
    end do

    ! Check if aspect ratio lies within data coordinates
    ! Reynolds number can't be negative and minimum bound is zero (Stokes flow)
    if (p%lpp_rey.ge.dat_reyn(num_reyn)) then
        !! Extrapolate on higher end
        icrd(4) = num_reyn-1
    else
        !! Find the interval for interpolation
        do i = 1,num_reyn-1
            icrd(4) = i
            if ((dat_reyn(i).le.p%lpp_rey).and.(dat_reyn(i+1).gt.p%lpp_rey)) exit
        end do
    endif
    ! Calculate the multipliers/coefficients for Reynolds number interpolation
    fcrd(4,0) = (dat_reyn(icrd(4)+1) - p%lpp_rey)/(dat_reyn(icrd(4)+1) - dat_reyn(icrd(4)))
    fcrd(4,1) = (p%lpp_rey  -  dat_reyn(icrd(4)))/(dat_reyn(icrd(4)+1) - dat_reyn(icrd(4)))

    ! Calculate the final slip correction coefficient through multilinear interpolation
    coef = 0.0
    do i = 0,1
        do j = 0,1
            do k = 0,1
                do l = 0,1
                    cval = fcrd(1,i)*fcrd(2,j)*fcrd(3,k)*fcrd(4,l)
                    do m = 1,8
                        coef = coef + cval*aspm(m)*dat_coef(icrd(1)+i,icrd(2)+j,icrd(3)+k,icrd(4)+l,m)
                    end do
                end do
            end do
        end do
    end do

    coef = max(coef,0.0)

end subroutine CalcSlipCorrectionCoefficient

subroutine CalcCurvatureTerms

    use decomp_2d, only: xstart,xend,update_halo
    use param
    use local_arrays, only: vx,vy,vz,dfx,dfy,dfz
    use lagrangian_point_particle, only: lpp_d2vx,lpp_d2vy,lpp_d2vz

    implicit none

    integer :: im,jm,km
    integer :: ic,jc,kc
    integer :: ip,jp,kp

    real    :: dxxv,dyyv,dzzv

    ! Loop over all grid points in pencil/process
    do kc=xstart(3),xend(3)
        km=kc-1
        kp=kc+1
        do jc=xstart(2),xend(2)
            jm=jc-1
            jp=jc+1
            do ic=xstart(1),xend(1)
                im=ic-1
                ip=ic+1

                ! Compute (d2(vx)/dx2), (d2(vx)/dy2), (d2(vx)/dz2)
                dxxv = vx(ip,jc,kc)*ap1si(ic) + vx(ic,jc,kc)*ac1si(ic) + vx(im,jc,kc)*am1si(ic)
                dyyv = vx(ic,jp,kc)*ap2cj(jc) + vx(ic,jc,kc)*ac2cj(jc) + vx(ic,jm,kc)*am2cj(jc)
                dzzv = vx(ic,jc,kp)*ap3ck(kc) + vx(ic,jc,kc)*ac3ck(kc) + vx(ic,jc,km)*am3ck(kc)
                ! Save (d2(vx)/dx2) + (d2(vx)/dy2) + (d2(vx)/dz2) for advection-diffusion
                dfx(ic,jc,kc) = dxxv + dyyv + dzzv
                ! Save (d2(vx)/dx2)(dx*dx) + (d2(vx)/dy2)(dy*dy) + (d2(vx)/dz2)(dz*dz) for slip correction
                lpp_d2vx(ic,jc,kc) = dxxv*dxm(ic)*dxm(ic) + dyyv*dyc(jc)*dyc(jc) + dzzv*dzc(kc)*dzc(kc)

                ! Compute (d2(vy)/dx2), (d2(vy)/dy2), (d2(vy)/dz2)
                dxxv = vy(ip,jc,kc)*ap1ci(ic) + vy(ic,jc,kc)*ac1ci(ic) + vy(im,jc,kc)*am1ci(ic)
                dyyv = vy(ic,jp,kc)*ap2sj(jc) + vy(ic,jc,kc)*ac2sj(jc) + vy(ic,jm,kc)*am2sj(jc)
                dzzv = vy(ic,jc,kp)*ap3ck(kc) + vy(ic,jc,kc)*ac3ck(kc) + vy(ic,jc,km)*am3ck(kc)
                ! Save (d2(vy)/dx2) + (d2(vy)/dy2) + (d2(vy)/dz2) for advection-diffusion
                dfy(ic,jc,kc) = dxxv + dyyv + dzzv
                ! Save (d2(vy)/dx2)(dx*dx) + (d2(vy)/dy2)(dy*dy) + (d2(vy)/dz2)(dz*dz) for slip correction
                lpp_d2vy(ic,jc,kc) = dxxv*dxc(ic)*dxc(ic) + dyyv*dym(jc)*dym(jc) + dzzv*dzc(kc)*dzc(kc)

                ! Compute (d2(vz)/dx2), (d2(vz)/dy2), (d2(vz)/dz2)
                dxxv = vz(ip,jc,kc)*ap1ci(ic) + vz(ic,jc,kc)*ac1ci(ic) + vz(im,jc,kc)*am1ci(ic)
                dyyv = vz(ic,jp,kc)*ap2cj(jc) + vz(ic,jc,kc)*ac2cj(jc) + vz(ic,jm,kc)*am2cj(jc)
                dzzv = vz(ic,jc,kp)*ap3sk(kc) + vz(ic,jc,kc)*ac3sk(kc) + vz(ic,jc,km)*am3sk(kc)
                ! Save (d2(vz)/dx2) + (d2(vz)/dy2) + (d2(vz)/dz2) for advection-diffusion
                dfz(ic,jc,kc) = dxxv + dyyv + dzzv
                ! Save (d2(vz)/dx2)(dx*dx) + (d2(vz)/dy2)(dy*dy) + (d2(vz)/dz2)(dz*dz) for slip correction
                lpp_d2vz(ic,jc,kc) = dxxv*dxc(ic)*dxc(ic) + dyyv*dyc(jc)*dyc(jc) + dzzv*dzm(kc)*dzm(kc)

            end do
        end do
    end do

    ! Boundary values for the curvature terms under the assumption gradient of curvature at boundaries is zero
    !! For x-direction
    if (periodic(1)) then
        !!! For vx curvature
        lpp_d2vx(0 ,:,:) = lpp_d2vx(nxm  ,:,:)
        lpp_d2vx(nx,:,:) = lpp_d2vx(1    ,:,:)
        !!! For vy curvature
        lpp_d2vy(0 ,:,:) = lpp_d2vy(nxm  ,:,:)
        lpp_d2vy(nx,:,:) = lpp_d2vy(1    ,:,:)
        !!! For vz curvature
        lpp_d2vz(0 ,:,:) = lpp_d2vz(nxm  ,:,:)
        lpp_d2vz(nx,:,:) = lpp_d2vz(1    ,:,:)
    else
        !!! For vx curvature
        lpp_d2vx(0 ,:,:) = lpp_d2vx(2    ,:,:)
        lpp_d2vx(nx,:,:) = lpp_d2vx(nxm-1,:,:)
        !!! For vy curvature
        lpp_d2vy(0 ,:,:) = lpp_d2vy(1    ,:,:)
        lpp_d2vy(nx,:,:) = lpp_d2vy(nxm  ,:,:)
        !!! For vz curvature
        lpp_d2vz(0 ,:,:) = lpp_d2vz(1    ,:,:)
        lpp_d2vz(nx,:,:) = lpp_d2vz(nxm  ,:,:)
    end if
    !! For y-direction
    if (.not.periodic(2)) then
        !!! At y = 0
        if (xstart(2).eq.1) then
            !!!! For vx curvature
            lpp_d2vx(:,0,:) = lpp_d2vx(:,1,:)
            !!!! For vy curvature
            lpp_d2vy(:,0,:) = lpp_d2vy(:,2,:)
            !!!! For vz curvature
            lpp_d2vz(:,0,:) = lpp_d2vz(:,1,:)
        end if
        !!! At y = ylen
        if (xend(2).eq.nym) then
            !!!! For vx curvature
            lpp_d2vx(:,ny,:) = lpp_d2vx(:,nym  ,:)
            !!!! For vy curvature
            lpp_d2vy(:,ny,:) = lpp_d2vy(:,nym-1,:)
            !!!! For vz curvature
            lpp_d2vz(:,ny,:) = lpp_d2vz(:,nym  ,:)
        end if
    end if
    !! For z-direction
    if (.not.periodic(3)) then
        !!! At z = 0
        if (xstart(3).eq.1) then
            !!!! For vx curvature
            lpp_d2vx(:,:,0) = lpp_d2vx(:,:,1)
            !!!! For vy curvature
            lpp_d2vy(:,:,0) = lpp_d2vy(:,:,1)
            !!!! For vz curvature
            lpp_d2vz(:,:,0) = lpp_d2vz(:,:,2)
        end if
        !!! At z = zlen
        if (xend(3).eq.nzm) then
            !!!! For vx curvature
            lpp_d2vx(:,:,nz) = lpp_d2vx(:,:,nzm  )
            !!!! For vy curvature
            lpp_d2vy(:,:,nz) = lpp_d2vy(:,:,nzm  )
            !!!! For vz curvature
            lpp_d2vz(:,:,nz) = lpp_d2vz(:,:,nzm-1)
        end if
    end if

    ! Ensure the curvatures in halo cells are also updated 
    call update_halo(lpp_d2vx,lvlhalo)
    call update_halo(lpp_d2vy,lvlhalo)
    call update_halo(lpp_d2vz,lvlhalo)

end subroutine CalcCurvatureTerms