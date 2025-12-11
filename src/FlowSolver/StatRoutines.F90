!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutines InitProfiles and WriteProfiles !
!                                                         !
!    PURPOSE: Calculates and writes out statistical       !
!    profiles.                                            !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitProfiles

    use param
    use stat_arrays

    implicit none

    integer         :: rnum
    logical         :: exists
    character*4     :: charpnum
    character*200   :: filename,dsetname

    pnum = 1
    snpp = 1
    
    if (nread) then

        if (ismaster) then

            filename = trim("Results/stat_profiles_x.h5")
            inquire(file=filename, exist=exists)
            if (exists) then
                filename = filename
                dsetname = trim("/runs")
                call HdfSerialReadIntScalar(filename,dsetname,rnum)
                pnum = max(pnum,rnum)
                dsetname = trim("/steps")
                call HdfSerialReadIntScalar(filename,dsetname,rnum)
                snpp = max(snpp,rnum)
            end if
            
            filename = trim("Results/stat_profiles_y.h5")
            inquire(file=filename, exist=exists)
            if (exists) then
                filename = filename
                dsetname = trim("/runs")
                call HdfSerialReadIntScalar(filename,dsetname,rnum)
                pnum = max(pnum,rnum)
                dsetname = trim("/steps")
                call HdfSerialReadIntScalar(filename,dsetname,rnum)
                snpp = max(snpp,rnum)
            end if

            filename = trim("Results/stat_profiles_z.h5")
            inquire(file=filename, exist=exists)
            if (exists) then
                filename = filename
                dsetname = trim("/runs")
                call HdfSerialReadIntScalar(filename,dsetname,rnum)
                pnum = max(pnum,rnum)
                dsetname = trim("/steps")
                call HdfSerialReadIntScalar(filename,dsetname,rnum)
                snpp = max(snpp,rnum)
            end if

            pnum = pnum + 1
            snpp = snpp + 1

        end if

        call MpiBarrier
        call MpiBcastInt(pnum)
        call MpiBcastInt(snpp)
        
    else

        if (ismaster) then

            filename = trim("Results/stat_profiles_x.h5")
            inquire(file=filename, exist=exists)
            if (exists) call HdfClean(filename)
            filename = trim("Results/stat_profiles_y.h5")
            inquire(file=filename, exist=exists)
            if (exists) call HdfClean(filename)
            filename = trim("Results/stat_profiles_z.h5")
            inquire(file=filename, exist=exists)
            if (exists) call HdfClean(filename)

        end if

    end if

    if (ismaster) then

        write(charpnum,"(I4.4)") pnum

        filename = trim("Results/stat_profiles_x.h5")
        dsetname = trim("/runs")
        call HdfSerialWriteIntScalar(filename,dsetname,pnum)
        dsetname = trim("/info/"//charpnum//"/st")
        call HdfSerialWriteIntScalar(filename,dsetname,snpp)
        dsetname = trim("/info/"//charpnum//"/xc")
        call HdfSerialWriteReal1D(filename,dsetname,xc,1,nx)
        dsetname = trim("/info/"//charpnum//"/xm")
        call HdfSerialWriteReal1D(filename,dsetname,xm,0,nx)

        filename = trim("Results/stat_profiles_y.h5")
        dsetname = trim("/runs")
        call HdfSerialWriteIntScalar(filename,dsetname,pnum)
        dsetname = trim("/info/"//charpnum//"/st")
        call HdfSerialWriteIntScalar(filename,dsetname,snpp)
        dsetname = trim("/info/"//charpnum//"/yc")
        call HdfSerialWriteReal1D(filename,dsetname,yc,1,ny)
        dsetname = trim("/info/"//charpnum//"/ym")
        call HdfSerialWriteReal1D(filename,dsetname,ym,0,ny)

        filename = trim("Results/stat_profiles_z.h5")
        dsetname = trim("/runs")
        call HdfSerialWriteIntScalar(filename,dsetname,pnum)
        dsetname = trim("/info/"//charpnum//"/st")
        call HdfSerialWriteIntScalar(filename,dsetname,snpp)
        dsetname = trim("/info/"//charpnum//"/zc")
        call HdfSerialWriteReal1D(filename,dsetname,zc,1,nz)
        dsetname = trim("/info/"//charpnum//"/zm")
        call HdfSerialWriteReal1D(filename,dsetname,zm,0,nz)

    end if

    return

end subroutine InitProfiles

subroutine WriteStatProfiles

    use decomp_2d
    use param
    use local_arrays, only: vx,vy,vz
    use stat_arrays

    implicit none

    integer         :: ic,jc,kc
    integer         :: ip,jp,kp
    integer         :: st1,en1,st2,en2,st3,en3
    character*8     :: charsnap
    character*200   :: filename,dsetname
    real            :: rprof_xc(1:nx),rprof_xm(0:nx)
    real            :: rprof_yc(1:ny),rprof_ym(0:ny)
    real            :: rprof_zc(1:nz),rprof_zm(0:nz)

    !! VX STATISTICS !!

    vx_m1_xp(:) = 0.0
    vx_m2_xp(:) = 0.0
    vx_m1_yp(:) = 0.0
    vx_m2_yp(:) = 0.0
    vx_m1_zp(:) = 0.0
    vx_m2_zp(:) = 0.0

    st1 = 1
    en1 = nx
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)

    do kc = st3,en3
        do jc = st2,en2
            do ic = st1,en1
                vx_m1_xp(ic) = vx_m1_xp(ic) + (vx(ic,jc,kc)     )*dyc(jc)*dzc(kc)/(ylen*zlen)
                vx_m2_xp(ic) = vx_m2_xp(ic) + (vx(ic,jc,kc)**2.0)*dyc(jc)*dzc(kc)/(ylen*zlen)
            end do
        end do
    end do

    st1 = 1
    en1 = nxm
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)
    if (st2.eq.1  ) st2 = 0
    if (en2.eq.nym) en2 = ny

    do kc = st3,en3
        do jc = st2,en2
            do ic = st1,en1
                ip = ic+1
                vx_m1_yp(jc) = vx_m1_yp(jc) + 0.5*(vx(ic,jc,kc)      + vx(ip,jc,kc)     )*dxc(ic)*dzc(kc)/(xlen*zlen)
                vx_m2_yp(jc) = vx_m2_yp(jc) + 0.5*(vx(ic,jc,kc)**2.0 + vx(ip,jc,kc)**2.0)*dxc(ic)*dzc(kc)/(xlen*zlen)
            end do
        end do
    end do

    st1 = 1
    en1 = nxm
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)
    if (st3.eq.1  ) st3 = 0
    if (en3.eq.nzm) en3 = nz

    do kc = st3,en3
        do jc = st2,en2
            do ic = st1,en1
                ip = ic+1
                vx_m1_zp(kc) = vx_m1_zp(kc) + 0.5*(vx(ic,jc,kc)      + vx(ip,jc,kc)     )*dxc(ic)*dyc(jc)/(xlen*ylen)
                vx_m2_zp(kc) = vx_m2_zp(kc) + 0.5*(vx(ic,jc,kc)**2.0 + vx(ip,jc,kc)**2.0)*dxc(ic)*dyc(jc)/(xlen*ylen)
            end do
        end do
    end do

    call MpiSumReal1D(vx_m1_xp,rprof_xc,1,nx)
    vx_m1_xp(1:nx) = rprof_xc(1:nx)

    call MpiSumReal1D(vx_m2_xp,rprof_xc,1,nx)
    vx_m2_xp(1:nx) = rprof_xc(1:nx)

    call MpiSumReal1D(vx_m1_yp,rprof_ym,0,ny)
    vx_m1_yp(0:ny) = rprof_ym(0:ny)

    call MpiSumReal1D(vx_m2_yp,rprof_ym,0,ny)
    vx_m2_yp(0:ny) = rprof_ym(0:ny)

    call MpiSumReal1D(vx_m1_zp,rprof_zm,0,nz)
    vx_m1_zp(0:nz) = rprof_zm(0:nz)

    call MpiSumReal1D(vx_m2_zp,rprof_zm,0,nz)
    vx_m2_zp(0:nz) = rprof_zm(0:nz)

    ! VY STATISTICS !!

    vy_m1_xp(:) = 0.0
    vy_m2_xp(:) = 0.0
    vy_m1_yp(:) = 0.0
    vy_m2_yp(:) = 0.0
    vy_m1_zp(:) = 0.0
    vy_m2_zp(:) = 0.0

    st1 = 0
    en1 = nx
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)

    do kc = st3,en3
        do jc = st2,en2
            jp = jc+1
            do ic = st1,en1
                vy_m1_xp(ic) = vy_m1_xp(ic) + 0.5*(vy(ic,jc,kc)      + vy(ic,jp,kc)     )*dyc(jc)*dzc(kc)/(ylen*zlen)
                vy_m2_xp(ic) = vy_m2_xp(ic) + 0.5*(vy(ic,jc,kc)**2.0 + vy(ic,jp,kc)**2.0)*dyc(jc)*dzc(kc)/(ylen*zlen)
            end do
        end do
    end do

    st1 = 1
    en1 = nxm
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)
    if (en2.eq.nym) en2 = ny

    do kc = st3,en3
        do jc = st2,en2
            do ic = st1,en1
                vy_m1_yp(jc) = vy_m1_yp(jc) + (vy(ic,jc,kc)     )*dxc(ic)*dzc(kc)/(xlen*zlen)
                vy_m2_yp(jc) = vy_m2_yp(jc) + (vy(ic,jc,kc)**2.0)*dxc(ic)*dzc(kc)/(xlen*zlen)
            end do
        end do
    end do

    st1 = 1
    en1 = nxm
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)
    if (st3.eq.1  ) st3 = 0
    if (en3.eq.nzm) en3 = nz

    do kc = st3,en3
        do jc = st2,en2
            jp = jc+1
            do ic = st1,en1
                vy_m1_zp(kc) = vy_m1_zp(kc) + 0.5*(vy(ic,jc,kc)      + vy(ic,jp,kc)     )*dxc(ic)*dyc(jc)/(xlen*ylen)
                vy_m2_zp(kc) = vy_m2_zp(kc) + 0.5*(vy(ic,jc,kc)**2.0 + vy(ic,jp,kc)**2.0)*dxc(ic)*dyc(jc)/(xlen*ylen)
            end do
        end do
    end do

    call MpiSumReal1D(vy_m1_xp,rprof_xm,0,nx)
    vy_m1_xp(0:nx) = rprof_xm(0:nx)

    call MpiSumReal1D(vy_m2_xp,rprof_xm,0,nx)
    vy_m2_xp(0:nx) = rprof_xm(0:nx)

    call MpiSumReal1D(vy_m1_yp,rprof_yc,1,ny)
    vy_m1_yp(1:ny) = rprof_yc(1:ny)

    call MpiSumReal1D(vy_m2_yp,rprof_yc,1,ny)
    vy_m2_yp(1:ny) = rprof_yc(1:ny)

    call MpiSumReal1D(vy_m1_zp,rprof_zm,0,nz)
    vy_m1_zp(0:nz) = rprof_zm(0:nz)

    call MpiSumReal1D(vy_m2_zp,rprof_zm,0,nz)
    vy_m2_zp(0:nz) = rprof_zm(0:nz)

    ! VZ STATISTICS !!

    vz_m1_xp(:) = 0.0
    vz_m2_xp(:) = 0.0
    vz_m1_yp(:) = 0.0
    vz_m2_yp(:) = 0.0
    vz_m1_zp(:) = 0.0
    vz_m2_zp(:) = 0.0

    st1 = 0
    en1 = nx
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)

    do kc = st3,en3
        kp = kc+1
        do jc = st2,en2
            do ic = st1,en1
                vz_m1_xp(ic) = vz_m1_xp(ic) + 0.5*(vz(ic,jc,kc)      + vz(ic,jc,kp)     )*dyc(jc)*dzc(kc)/(ylen*zlen)
                vz_m2_xp(ic) = vz_m2_xp(ic) + 0.5*(vz(ic,jc,kc)**2.0 + vz(ic,jc,kp)**2.0)*dyc(jc)*dzc(kc)/(ylen*zlen)
            end do
        end do
    end do

    st1 = 1
    en1 = nxm
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)
    if (st2.eq.1  ) st2 = 0
    if (en2.eq.nym) en2 = ny

    do kc = st3,en3
        kp = kc+1
        do jc = st2,en2
            do ic = st1,en1
                vz_m1_yp(jc) = vz_m1_yp(jc) + 0.5*(vz(ic,jc,kc)      + vz(ic,jc,kp)     )*dxc(ic)*dzc(kc)/(xlen*zlen)
                vz_m2_yp(jc) = vz_m2_yp(jc) + 0.5*(vz(ic,jc,kc)**2.0 + vz(ic,jc,kp)**2.0)*dxc(ic)*dzc(kc)/(xlen*zlen)
            end do
        end do
    end do

    st1 = 1
    en1 = nxm
    st2 = xstart(2)
    en2 = xend(2)
    st3 = xstart(3)
    en3 = xend(3)
    if (en3.eq.nzm) en3 = nz

    do kc = st3,en3
        do jc = st2,en2
            do ic = st1,en1
                vz_m1_zp(kc) = vz_m1_zp(kc) + (vz(ic,jc,kc))*dxc(ic)*dyc(jc)/(xlen*ylen)
                vz_m2_zp(kc) = vz_m2_zp(kc) + (vz(ic,jc,kc))*dxc(ic)*dyc(jc)/(xlen*ylen)
            end do
        end do
    end do

    call MpiSumReal1D(vz_m1_xp,rprof_xm,0,nx)
    vz_m1_xp(0:nx) = rprof_xm(0:nx)

    call MpiSumReal1D(vz_m2_xp,rprof_xm,0,nx)
    vz_m2_xp(0:nx) = rprof_xm(0:nx)

    call MpiSumReal1D(vz_m1_yp,rprof_ym,0,ny)
    vz_m1_yp(0:ny) = rprof_ym(0:ny)

    call MpiSumReal1D(vz_m2_yp,rprof_ym,0,ny)
    vz_m2_yp(0:ny) = rprof_ym(0:ny)

    call MpiSumReal1D(vz_m1_zp,rprof_zc,1,nz)
    vz_m1_zp(1:nz) = rprof_zc(1:nz)

    call MpiSumReal1D(vz_m2_zp,rprof_zc,1,nz)
    vz_m2_zp(1:nz) = rprof_zc(1:nz)

    if (ismaster) then

        write(charsnap,"(I8.8)") snpp

        !! WRITE X-PROFILES !!

        filename = trim("Results/stat_profiles_x.h5")
        
        dsetname = trim("/steps")
        call HdfSerialWriteIntScalar(filename,dsetname,snpp)
        dsetname = trim("/data/"//charsnap//"/time")
        call HdfSerialWriteRealScalar(filename,dsetname,time)
        dsetname = trim("/data/"//charsnap//"/dt")
        call HdfSerialWriteRealScalar(filename,dsetname,dt)
        dsetname = trim("/data/"//charsnap//"/vx_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vx_m1_xp,1,nx)
        dsetname = trim("/data/"//charsnap//"/vx_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vx_m2_xp,1,nx)
        dsetname = trim("/data/"//charsnap//"/vy_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vy_m1_xp,0,nx)
        dsetname = trim("/data/"//charsnap//"/vy_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vy_m2_xp,0,nx)
        dsetname = trim("/data/"//charsnap//"/vz_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vz_m1_xp,0,nx)
        dsetname = trim("/data/"//charsnap//"/vz_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vz_m2_xp,0,nx)

        !! WRITE Y-PROFILES !!

        filename = trim("Results/stat_profiles_y.h5")

        dsetname = trim("/steps")
        call HdfSerialWriteIntScalar(filename,dsetname,snpp)
        dsetname = trim("/data/"//charsnap//"/time")
        call HdfSerialWriteRealScalar(filename,dsetname,time)
        dsetname = trim("/data/"//charsnap//"/dt")
        call HdfSerialWriteRealScalar(filename,dsetname,dt)
        dsetname = trim("/data/"//charsnap//"/vx_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vx_m1_yp,0,ny)
        dsetname = trim("/data/"//charsnap//"/vx_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vx_m2_yp,0,ny)
        dsetname = trim("/data/"//charsnap//"/vy_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vy_m1_yp,1,ny)
        dsetname = trim("/data/"//charsnap//"/vy_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vy_m2_yp,1,ny)
        dsetname = trim("/data/"//charsnap//"/vz_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vz_m1_yp,0,ny)
        dsetname = trim("/data/"//charsnap//"/vz_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vz_m2_yp,0,ny)

        !! WRITE Z-PROFILES !!

        filename = trim("Results/stat_profiles_z.h5")
        
        dsetname = trim("/steps")
        call HdfSerialWriteIntScalar(filename,dsetname,snpp)
        dsetname = trim("/data/"//charsnap//"/time")
        call HdfSerialWriteRealScalar(filename,dsetname,time)
        dsetname = trim("/data/"//charsnap//"/dt")
        call HdfSerialWriteRealScalar(filename,dsetname,dt)
        dsetname = trim("/data/"//charsnap//"/vx_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vx_m1_zp,0,nz)
        dsetname = trim("/data/"//charsnap//"/vx_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vx_m2_zp,0,nz)
        dsetname = trim("/data/"//charsnap//"/vy_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vy_m1_zp,0,nz)
        dsetname = trim("/data/"//charsnap//"/vy_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vy_m2_zp,0,nz)
        dsetname = trim("/data/"//charsnap//"/vz_m1")
        call HdfSerialWriteReal1D(filename,dsetname,vz_m1_zp,1,nz)
        dsetname = trim("/data/"//charsnap//"/vz_m2")
        call HdfSerialWriteReal1D(filename,dsetname,vz_m2_zp,1,nz)

    end if

    snpp = snpp + 1

end subroutine WriteStatProfiles