!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: StatRoutines.F90                               !
!    CONTAINS: subroutine GetLines                        !
!              subroutine GetInterpolation                !
!              subroutine InitMovie                       !
!              subroutine WriteMovieSlices                !
!                                                         !
!    PURPOSE: Calculates and writes out interpolated      !
!    2D slices for movie.                                 !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GetLines(filename,nlines)

    use, intrinsic :: iso_fortran_env, only : iostat_end

    implicit none

    character*50, intent(in)    :: filename
    integer, intent(out)        :: nlines
    integer                     :: ierror
    real                        :: locval
    logical                     :: savevx,savevy,savevz
    logical                     :: savepr
    logical                     :: saveox,saveoy,saveoz

    open(99, file=trim(filename))

    nlines = 0
    ierror = 0

    do while (ierror.ne.iostat_end)
        read(99, *, iostat=ierror) locval, savevx, savevy, savevz, savepr, saveox, saveoy, saveoz
        if (ierror.eq.0) nlines = nlines + 1
    end do

    close(99)

end subroutine GetLines

subroutine GetInterpolation(grid,stid,enid,sloc,l_id,h_id,l_cf,h_cf)

    implicit none

    real, intent(in)        :: grid(stid:enid)
    integer, intent(in)     :: stid
    integer, intent(in)     :: enid
    real, intent(in)        :: sloc
    integer, intent(out)    :: l_id
    integer, intent(out)    :: h_id
    real, intent(out)       :: l_cf
    real, intent(out)       :: h_cf
    integer                 :: lc,lp

    l_id = enid
    h_id = enid
    l_cf = 0.5
    h_cf = 0.5

    do lc = stid,enid-1
        lp = lc+1
        if ((grid(lc).le.sloc).and.(grid(lp).gt.sloc)) then
            l_id = lc
            h_id = lp
            l_cf = (grid(lp) - sloc)/(grid(lp) - grid(lc))
            h_cf = (sloc - grid(lc))/(grid(lp) - grid(lc))
        end if
    end do

end subroutine GetInterpolation

subroutine InitMovie

    use param
    use movie_arrays

    implicit none

    logical         :: exists
    logical         :: savevx,savevy,savevz
    logical         :: savepr
    logical         :: saveox,saveoy,saveoz
    integer         :: rnum,nlines,l
    integer         :: ierror
    real            :: locval
    character*4     :: charmnum
    character*200   :: filename,dsetname

    mnum = 1
    snpm = 1

    if (ismaster) then

        filename = trim("Results/movie_slices_x.h5")
        inquire(file=filename, exist=exists)
        if (exists) then
            dsetname = trim("/runs")
            call HdfSerialReadIntScalar(filename,dsetname,rnum)
            mnum = max(mnum,rnum)
            dsetname = trim("/steps")
            call HdfSerialReadIntScalar(filename,dsetname,rnum)
            snpm = max(snpm,rnum)
        end if
        
        filename = trim("Results/movie_slices_y.h5")
        inquire(file=filename, exist=exists)
        if (exists) then
            filename = filename
            dsetname = trim("/runs")
            call HdfSerialReadIntScalar(filename,dsetname,rnum)
            mnum = max(mnum,rnum)
            dsetname = trim("/steps")
            call HdfSerialReadIntScalar(filename,dsetname,rnum)
            snpm = max(snpm,rnum)
        end if

        filename = trim("Results/movie_slices_z.h5")
        inquire(file=filename, exist=exists)
        if (exists) then
            filename = filename
            dsetname = trim("/runs")
            call HdfSerialReadIntScalar(filename,dsetname,rnum)
            mnum = max(mnum,rnum)
            dsetname = trim("/steps")
            call HdfSerialReadIntScalar(filename,dsetname,rnum)
            snpm = max(snpm,rnum)
        end if

        mnum = mnum + 1
        snpm = snpm + 1

    end if

    call MpiBarrier
    call MpiBcastInt(mnum)
    call MpiBcastInt(snpm)

    filename = trim("Inputs/movie_slices_x.in")
    inquire(file=filename, exist=exists)

    if (exists) then

        call GetLines(filename,nlines)

        if (nlines.gt.0) then

            if (.not.allocated(slcx)) allocate(slcx(1:nlines))
            if (.not.allocated(sidx)) allocate(sidx(1:nlines,1:4))
            if (.not.allocated(scfx)) allocate(scfx(1:nlines,1:4))
            if (.not.allocated(sqnx)) allocate(sqnx(1:nlines,1:7))

            open(99, file=trim(filename))
            l = 1
            do while (l.le.nlines)
                read(99, *, iostat=ierror) locval, savevx, savevy, savevz, savepr, saveox, saveoy, saveoz
                if (ierror.eq.0) then
                    slcx(l)   = max(xc(1),min(xc(nx),locval))
                    sqnx(l,1) = savevx
                    sqnx(l,2) = savevy
                    sqnx(l,3) = savevz
                    sqnx(l,4) = savepr
                    sqnx(l,5) = saveox
                    sqnx(l,6) = saveoy
                    sqnx(l,7) = saveoz
                    call GetInterpolation(xc,1,nx,slcx(l),sidx(l,1),sidx(l,2),scfx(l,1),scfx(l,2))
                    call GetInterpolation(xm,0,nx,slcx(l),sidx(l,3),sidx(l,4),scfx(l,3),scfx(l,4))
                    l = l + 1
                end if
            end do
            close(99)

            write(charmnum,"(I4.4)") mnum

            if (ismaster) then

                filename = trim("Results/movie_slices_x.h5")

                dsetname = trim("/runs")
                call HdfSerialWriteIntScalar(filename,dsetname,mnum)
                dsetname = trim("/info/"//charmnum//"/st")
                call HdfSerialWriteIntScalar(filename,dsetname,snpm)
                dsetname = trim("/info/"//charmnum//"/xs")
                call HdfSerialWriteReal1D(filename,dsetname,slcx,1,nlines)
                dsetname = trim("/info/"//charmnum//"/yc")
                call HdfSerialWriteReal1D(filename,dsetname,yc,1,ny)
                dsetname = trim("/info/"//charmnum//"/ym")
                call HdfSerialWriteReal1D(filename,dsetname,ym,0,ny)
                dsetname = trim("/info/"//charmnum//"/zc")
                call HdfSerialWriteReal1D(filename,dsetname,zc,1,nz)
                dsetname = trim("/info/"//charmnum//"/zm")
                call HdfSerialWriteReal1D(filename,dsetname,zm,0,nz)

            end if

        end if

    end if
 
    filename = trim("Inputs/movie_slices_y.in")
    inquire(file=filename, exist=exists)

    if (exists) then

        call GetLines(filename,nlines)

        if (nlines.gt.0) then

            if (.not.allocated(slcy)) allocate(slcy(1:nlines))
            if (.not.allocated(sidy)) allocate(sidy(1:nlines,1:4))
            if (.not.allocated(scfy)) allocate(scfy(1:nlines,1:4))
            if (.not.allocated(sqny)) allocate(sqny(1:nlines,1:7))

            open(99, file=trim(filename))
            l = 1
            do while(l.le.nlines)
                read(99, *, iostat=ierror) locval, savevx, savevy, savevz, savepr, saveox, saveoy, saveoz
                if (ierror.eq.0) then
                    slcy(l)   = max(yc(1),min(yc(ny),locval))
                    sqny(l,1) = savevx
                    sqny(l,2) = savevy
                    sqny(l,3) = savevz
                    sqny(l,4) = savepr
                    sqny(l,5) = saveox
                    sqny(l,6) = saveoy
                    sqny(l,7) = saveoz
                    call GetInterpolation(yc,1,ny,slcy(l),sidy(l,1),sidy(l,2),scfy(l,1),scfy(l,2))
                    call GetInterpolation(ym,0,ny,slcy(l),sidy(l,3),sidy(l,4),scfy(l,3),scfy(l,4))
                    l = l + 1
                end if
            end do
            close(99)
            
            write(charmnum,"(I4.4)") mnum

            if (ismaster) then

                filename = trim("Results/movie_slices_y.h5")

                dsetname = trim("/runs")
                call HdfSerialWriteIntScalar(filename,dsetname,mnum)
                dsetname = trim("/info/"//charmnum//"/st")
                call HdfSerialWriteIntScalar(filename,dsetname,snpm)
                dsetname = trim("/info/"//charmnum//"/xc")
                call HdfSerialWriteReal1D(filename,dsetname,xc,1,nx)
                dsetname = trim("/info/"//charmnum//"/xm")
                call HdfSerialWriteReal1D(filename,dsetname,xm,0,nx)
                dsetname = trim("/info/"//charmnum//"/ys")
                call HdfSerialWriteReal1D(filename,dsetname,slcy,1,nlines)
                dsetname = trim("/info/"//charmnum//"/zc")
                call HdfSerialWriteReal1D(filename,dsetname,zc,1,nz)
                dsetname = trim("/info/"//charmnum//"/zm")
                call HdfSerialWriteReal1D(filename,dsetname,zm,0,nz)

            end if

        end if

    end if

    filename = trim("Inputs/movie_slices_z.in")
    inquire(file=filename, exist=exists)

    if (exists) then

        call GetLines(filename,nlines)

        if (nlines.gt.0) then

            if (.not.allocated(slcz)) allocate(slcz(1:nlines))
            if (.not.allocated(sidz)) allocate(sidz(1:nlines,1:4))
            if (.not.allocated(scfz)) allocate(scfz(1:nlines,1:4))
            if (.not.allocated(sqnz)) allocate(sqnz(1:nlines,1:7))

            open(99, file=trim(filename))
            l = 1
            do while (l.le.nlines)
                read(99, *, iostat=ierror) locval, savevx, savevy, savevz, savepr, saveox, saveoy, saveoz
                if (ierror.eq.0) then
                    slcz(l)   = max(zc(1),min(zc(nz),locval))
                    sqnz(l,1) = savevx
                    sqnz(l,2) = savevy
                    sqnz(l,3) = savevz
                    sqnz(l,4) = savepr
                    sqnz(l,5) = saveox
                    sqnz(l,6) = saveoy
                    sqnz(l,7) = saveoz
                    call GetInterpolation(zc,1,nz,slcz(l),sidz(l,1),sidz(l,2),scfz(l,1),scfz(l,2))
                    call GetInterpolation(zm,0,nz,slcz(l),sidz(l,3),sidz(l,4),scfz(l,3),scfz(l,4))
                    l = l + 1
                end if
            end do
            close(99)
            
            write(charmnum,"(I4.4)") mnum

            if (ismaster) then

                filename = trim("Results/movie_slices_z.h5")

                dsetname = trim("/runs")
                call HdfSerialWriteIntScalar(filename,dsetname,mnum)
                dsetname = trim("/info/"//charmnum//"/st")
                call HdfSerialWriteIntScalar(filename,dsetname,snpm)
                dsetname = trim("/info/"//charmnum//"/xc")
                call HdfSerialWriteReal1D(filename,dsetname,xc,1,nx)
                dsetname = trim("/info/"//charmnum//"/xm")
                call HdfSerialWriteReal1D(filename,dsetname,xm,0,nx)
                dsetname = trim("/info/"//charmnum//"/yc")
                call HdfSerialWriteReal1D(filename,dsetname,yc,1,ny)
                dsetname = trim("/info/"//charmnum//"/ym")
                call HdfSerialWriteReal1D(filename,dsetname,ym,0,ny)
                dsetname = trim("/info/"//charmnum//"/zs")
                call HdfSerialWriteReal1D(filename,dsetname,slcz,1,nlines)

            end if

        end if

    end if

end subroutine InitMovie

subroutine WriteMovieSlices

    use decomp_2d
    use param
    use local_arrays, only: vx,vy,vz,pr
    use movie_arrays

    implicit none
    
    integer         :: lc
    integer         :: im,jm,km
    integer         :: ic,jc,kc
    integer         :: ip,jp,kp
    integer         :: st1,en1,st2,en2,st3,en3
    real            :: fc,fp
    character*4     :: charindx
    character*8     :: charsnap
    character*200   :: dsetname,filename

    write(charsnap,"(I8.8)") snpm

    ! For X-Slices
    if (allocated(slcx)) then

        filename = trim("Results/movie_slices_x.h5")
        if (ismaster) then
            dsetname = trim("/steps")
            call HdfSerialWriteIntScalar(filename,dsetname,snpm)
            dsetname = trim("/data/"//charsnap//"/time")
            call HdfSerialWriteRealScalar(filename,dsetname,time)
            dsetname = trim("/data/"//charsnap//"/dt")
            call HdfSerialWriteRealScalar(filename,dsetname,dt)
        end if

        do lc = 1,size(slcx)

            write(charindx,"(I4.4)") lc

            ! For X-Velocity

            if (sqnx(lc,1)) then

                st2 = xstart(2)
                en2 = xend(2)

                if (xstart(2).eq.1) st2 = 0
                if (xend(2).eq.nym) en2 = ny

                st3 = xstart(3)
                en3 = xend(3)

                if (xstart(3).eq.1) st3 = 0
                if (xend(3).eq.nzm) en3 = nz
                        
                fc = scfx(lc,1)
                fp = scfx(lc,2)
                ic = sidx(lc,1)
                ip = sidx(lc,2)

                do kc = st3,en3
                    do jc = st2,en2
                        mslx(jc,kc) = (fc*vx(ic,jc,kc) + fp*vx(ip,jc,kc))
                    end do
                end do

                dsetname = trim("/data/"//charsnap//"/"//charindx//"/vx")
                call HdfParallelWriteReal2D(filename,dsetname,mslx(st2:en2,st3:en3),st2,en2,st3,en3,0,ny,0,nz,mpi_xcut)

            end if

            ! For Y-Velocity

            if (sqnx(lc,2)) then

                st2 = xstart(2)
                en2 = xend(2)

                if (xend(2).eq.nym) en2 = ny

                st3 = xstart(3)
                en3 = xend(3)

                if (xstart(3).eq.1) st3 = 0
                if (xend(3).eq.nzm) en3 = nz
                
                fc = scfx(lc,3)
                fp = scfx(lc,4)
                ic = sidx(lc,3)
                ip = sidx(lc,4)

                do kc = st3,en3
                    do jc = st2,en2
                        mslx(jc,kc) = (fc*vy(ic,jc,kc) + fp*vy(ip,jc,kc))
                    end do
                end do

                dsetname = trim("/data/"//charsnap//"/"//charindx//"/vy")
                call HdfParallelWriteReal2D(filename,dsetname,mslx(st2:en2,st3:en3),st2,en2,st3,en3,1,ny,0,nz,mpi_xcut)

            end if

            ! For Z-Velocity

            if (sqnx(lc,3)) then

                st2 = xstart(2)
                en2 = xend(2)

                if (xstart(2).eq.1) st2 = 0
                if (xend(2).eq.nym) en2 = ny

                st3 = xstart(3)
                en3 = xend(3)

                if (xend(3).eq.nzm) en3 = nz

                fc = scfx(lc,3)
                fp = scfx(lc,4)
                ic = sidx(lc,3)
                ip = sidx(lc,4)

                do kc = st3,en3
                    do jc = st2,en2
                        mslx(jc,kc) = (fc*vz(ic,jc,kc) + fp*vz(ip,jc,kc))
                    end do
                end do

                dsetname = trim("/data/"//charsnap//"/"//charindx//"/vz")
                call HdfParallelWriteReal2D(filename,dsetname,mslx(st2:en2,st3:en3),st2,en2,st3,en3,0,ny,1,nz,mpi_xcut)

            end if

            ! For Pressure

            if (sqnx(lc,4)) then

                st2 = xstart(2)
                en2 = xend(2)

                if (xstart(2).eq.1) st2 = 0
                if (xend(2).eq.nym) en2 = ny

                st3 = xstart(3)
                en3 = xend(3)

                if (xstart(3).eq.1) st3 = 0
                if (xend(3).eq.nzm) en3 = nz

                fc = scfx(lc,3)
                fp = scfx(lc,4)
                ic = sidx(lc,3)
                ip = sidx(lc,4)

                do kc = st3,en3
                    do jc = st2,en2
                        mslx(jc,kc) = (fc*pr(ic,jc,kc) + fp*pr(ip,jc,kc))
                    end do
                end do

                dsetname = trim("/data/"//charsnap//"/"//charindx//"/pr")
                call HdfParallelWriteReal2D(filename,dsetname,mslx(st2:en2,st3:en3),st2,en2,st3,en3,0,ny,0,nz,mpi_xcut)

            end if

            ! For X-Vorticity

            if (sqnx(lc,5)) then

                st2 = xstart(2)
                en2 = xend(2)

                if (xend(2).eq.nym) en2 = ny

                st3 = xstart(3)
                en3 = xend(3)

                if (xend(3).eq.nzm) en3 = nz

                fc = scfx(lc,3)
                fp = scfx(lc,4)
                ic = sidx(lc,3)
                ip = sidx(lc,4)

                do kc = st3,en3
                    km = kc-1
                    kp = kc+1
                    do jc = st2,en2
                        jm = jc-1
                        jp = jc+1
                        
                        mslx(jc,kc) = 0.0

                        mslx(jc,kc) = mslx(jc,kc) + (fc*(vz(ic,jc,kc) - vz(ic,jm,kc))/dym(jc))
                        mslx(jc,kc) = mslx(jc,kc) - (fc*(vy(ic,jc,kc) - vy(ic,jc,km))/dzm(kc))
                        mslx(jc,kc) = mslx(jc,kc) + (fp*(vz(ip,jc,kc) - vz(ip,jm,kc))/dym(jc))
                        mslx(jc,kc) = mslx(jc,kc) - (fp*(vy(ip,jc,kc) - vy(ip,jc,km))/dzm(kc))

                    end do
                end do

                dsetname = trim("/data/"//charsnap//"/"//charindx//"/ox")
                call HdfParallelWriteReal2D(filename,dsetname,mslx(st2:en2,st3:en3),st2,en2,st3,en3,1,ny,1,nz,mpi_xcut)

            end if

            ! For Y-Vorticity

            if (sqnx(lc,6)) then

                st2 = xstart(2)
                en2 = xend(2)

                if (xstart(2).eq.1) st2 = 0
                if (xend(2).eq.nym) en2 = ny

                st3 = xstart(3)
                en3 = xend(3)

                if (xend(3).eq.nzm) en3 = nz

                fc = scfx(lc,1)
                fp = scfx(lc,2)
                ic = sidx(lc,1)
                ip = sidx(lc,2)

                im = ic-1
                do kc = st3,en3
                    km = kc-1
                    kp = kc+1
                    do jc = st2,en2
                        jm = jc-1
                        jp = jc+1
                            
                        mslx(jc,kc) = 0.0

                        mslx(jc,kc) = mslx(jc,kc) + (fc*(vx(ic,jc,kc) - vx(ic,jc,km))/dzm(kc))
                        mslx(jc,kc) = mslx(jc,kc) - (fc*(vz(ic,jc,kc) - vz(im,jc,kc))/dxm(ic))
                        mslx(jc,kc) = mslx(jc,kc) + (fp*(vx(ip,jc,kc) - vx(ip,jc,km))/dzm(kc))
                        mslx(jc,kc) = mslx(jc,kc) - (fp*(vz(ip,jc,kc) - vz(ic,jc,kc))/dxm(ip))

                    end do
                end do

                dsetname = trim("/data/"//charsnap//"/"//charindx//"/oy")
                call HdfParallelWriteReal2D(filename,dsetname,mslx(st2:en2,st3:en3),st2,en2,st3,en3,0,ny,1,nz,mpi_xcut)

            end if

            ! For Z-Vorticity

            if (sqnx(lc,7)) then

                st2 = xstart(2)
                en2 = xend(2)

                if (xend(2).eq.nym) en2 = ny

                st3 = xstart(3)
                en3 = xend(3)

                if (xstart(3).eq.1) st3 = 0
                if (xend(3).eq.nzm) en3 = nz

                fc = scfx(lc,1)
                fp = scfx(lc,2)
                ic = sidx(lc,1)
                ip = sidx(lc,2)

                im = ic-1
                do kc = st3,en3
                    km = kc-1
                    kp = kc+1
                    do jc = st2,en2
                        jm = jc-1
                        jp = jc+1

                        mslx(jc,kc) = 0.0
                        
                        mslx(jc,kc) = mslx(jc,kc) + (fc*(vy(ic,jc,kc) - vy(im,jc,kc))/dxm(ic))
                        mslx(jc,kc) = mslx(jc,kc) - (fc*(vx(ic,jc,kc) - vx(ic,jm,kc))/dym(jc))
                        mslx(jc,kc) = mslx(jc,kc) + (fp*(vy(ip,jc,kc) - vy(ic,jc,kc))/dxm(ip))
                        mslx(jc,kc) = mslx(jc,kc) - (fp*(vx(ip,jc,kc) - vx(ip,jm,kc))/dym(jc))

                    end do
                end do

                dsetname = trim("/data/"//charsnap//"/"//charindx//"/oz")
                call HdfParallelWriteReal2D(filename,dsetname,mslx(st2:en2,st3:en3),st2,en2,st3,en3,1,ny,0,nz,mpi_xcut)

            end if

        end do

    end if

    ! For Y-slices
    if (allocated(slcy)) then

        filename = trim("Results/movie_slices_y.h5")
        if (ismaster) then
            dsetname = trim("/steps")
            call HdfSerialWriteIntScalar(filename,dsetname,snpm)
            dsetname = trim("/data/"//charsnap//"/time")
            call HdfSerialWriteRealScalar(filename,dsetname,time)
            dsetname = trim("/data/"//charsnap//"/dt")
            call HdfSerialWriteRealScalar(filename,dsetname,dt)
        end if

        call MpiBarrier

        do lc = 1,size(slcy)

            write(charindx,"(I4.4)") lc
            
            ! Assume that lvlhalo is at least 1
            if (((yc(xstart(2)).le.slcy(lc)).and.(yc(xend(2)+1).gt.slcy(lc))).or.((slcy(lc).eq.ylen).and.(xend(2).eq.nym))) then

                ! For X-Velocity

                if (sqny(lc,1)) then

                    st1 = 1
                    en1 = nx

                    st3 = xstart(3)
                    en3 = xend(3)

                    if (xstart(3).eq.1) st3 = 0
                    if (xend(3).eq.nzm) en3 = nz

                    fc = scfy(lc,3)
                    fp = scfy(lc,4)
                    jc = sidy(lc,3)
                    jp = sidy(lc,4)

                    do kc = st3,en3
                        do ic = st1,en1
                            msly(ic,kc) = (fc*vx(ic,jc,kc) + fp*vx(ic,jp,kc))
                        end do
                    end do

                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/vx")
                    call HdfParallelWriteReal2D(filename,dsetname,msly(st1:en1,st3:en3),st1,en1,st3,en3,1,nx,0,nz,mpi_ycut)

                end if

                ! For Y-Velocity

                if (sqny(lc,2)) then

                    st1 = 0
                    en1 = nx

                    st3 = xstart(3)
                    en3 = xend(3)

                    if (xstart(3).eq.1) st3 = 0
                    if (xend(3).eq.nzm) en3 = nz
                
                    fc = scfy(lc,1)
                    fp = scfy(lc,2)
                    jc = sidy(lc,1)
                    jp = sidy(lc,2)

                    do kc = st3,en3
                        do ic = st1,en1
                            msly(ic,kc) = (fc*vy(ic,jc,kc) + fp*vy(ic,jp,kc))
                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/vy")
                    call HdfParallelWriteReal2D(filename,dsetname,msly(st1:en1,st3:en3),st1,en1,st3,en3,0,nx,0,nz,mpi_ycut)

                end if

                ! For Z-Velocity

                if (sqny(lc,3)) then

                    st1 = 0
                    en1 = nx

                    st3 = xstart(3)
                    en3 = xend(3)

                    if (xend(3).eq.nzm) en3 = nz
                
                    fc = scfy(lc,3)
                    fp = scfy(lc,4)
                    jc = sidy(lc,3)
                    jp = sidy(lc,4)

                    do kc = st3,en3
                        do ic = st1,en1
                            msly(ic,kc) = (fc*vz(ic,jc,kc) + fp*vz(ic,jp,kc))
                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/vz")
                    call HdfParallelWriteReal2D(filename,dsetname,msly(st1:en1,st3:en3),st1,en1,st3,en3,0,nx,1,nz,mpi_ycut)

                end if

                ! For Pressure

                if (sqny(lc,4)) then

                    st1 = 0
                    en1 = nx

                    st3 = xstart(3)
                    en3 = xend(3)

                    if (xstart(3).eq.1) st3 = 0
                    if (xend(3).eq.nzm) en3 = nz
                    
                    fc = scfy(lc,3)
                    fp = scfy(lc,4)
                    jc = sidy(lc,3)
                    jp = sidy(lc,4)

                    do kc = st3,en3
                        do ic = st1,en1
                            msly(ic,kc) = (fc*pr(ic,jc,kc) + fp*pr(ic,jp,kc))
                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/pr")
                    call HdfParallelWriteReal2D(filename,dsetname,msly(st1:en1,st3:en3),st1,en1,st3,en3,0,nx,0,nz,mpi_ycut)

                end if

                ! For X-Vorticity

                if (sqny(lc,5)) then

                    st1 = 0
                    en1 = nx

                    st3 = xstart(3)
                    en3 = xend(3)

                    if (xend(3).eq.nzm) en3 = nz
                    
                    fc = scfy(lc,1)
                    fp = scfy(lc,2)
                    jc = sidy(lc,1)
                    jp = sidy(lc,2)

                    jm = jc-1
                    do kc = st3,en3
                        km = kc-1
                        kp = kc+1
                        do ic = st1,en1
                            im = ic-1
                            ip = ic+1

                            msly(ic,kc) = 0.0

                            msly(ic,kc) = msly(ic,kc) + (fc*(vz(ic,jc,kc) - vz(ic,jm,kc))/dym(jc))
                            msly(ic,kc) = msly(ic,kc) - (fc*(vy(ic,jc,kc) - vy(ic,jc,km))/dzm(kc))
                            msly(ic,kc) = msly(ic,kc) + (fp*(vz(ic,jp,kc) - vz(ic,jc,kc))/dym(jp))
                            msly(ic,kc) = msly(ic,kc) - (fp*(vy(ic,jp,kc) - vy(ic,jp,km))/dzm(kc))

                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/ox")
                    call HdfParallelWriteReal2D(filename,dsetname,msly(st1:en1,st3:en3),st1,en1,st3,en3,0,nx,1,nz,mpi_ycut)

                end if

                ! For Y-Vorticity

                if (sqny(lc,6)) then

                    st1 = 1
                    en1 = nx

                    st3 = xstart(3)
                    en3 = xend(3)

                    if (xend(3).eq.nzm) en3 = nz
                    
                    fc = scfy(lc,3)
                    fp = scfy(lc,4)
                    jc = sidy(lc,3)
                    jp = sidy(lc,4)

                    do kc = st3,en3
                        km = kc-1
                        kp = kc+1
                        do ic = st1,en1
                            im = ic-1
                            ip = ic+1

                            msly(ic,kc) = 0.0

                            msly(ic,kc) = msly(ic,kc) + (fc*(vx(ic,jc,kc) - vx(ic,jc,km))/dzm(kc))
                            msly(ic,kc) = msly(ic,kc) - (fc*(vz(ic,jc,kc) - vz(im,jc,kc))/dxm(ic))
                            msly(ic,kc) = msly(ic,kc) + (fp*(vx(ic,jp,kc) - vx(ic,jp,km))/dzm(kc))
                            msly(ic,kc) = msly(ic,kc) - (fp*(vz(ic,jp,kc) - vz(im,jp,kc))/dxm(ic))

                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/oy")
                    call HdfParallelWriteReal2D(filename,dsetname,msly(st1:en1,st3:en3),st1,en1,st3,en3,1,nx,1,nz,mpi_ycut)

                end if

                ! For Z-Vorticity

                if (sqny(lc,7)) then

                    st1 = 1
                    en1 = nx

                    st3 = xstart(3)
                    en3 = xend(3)

                    if (xstart(3).eq.1) st3 = 0
                    if (xend(3).eq.nzm) en3 = nz
                    
                    fc = scfy(lc,1)
                    fp = scfy(lc,2)
                    jc = sidy(lc,1)
                    jp = sidy(lc,2)

                    jm = jc-1
                    do kc = st3,en3
                        km = kc-1
                        kp = kc+1
                        do ic = st1,en1
                            im = ic-1
                            ip = ic+1

                            msly(ic,kc) = 0.0

                            msly(ic,kc) = msly(ic,kc) + (fc*(vy(ic,jc,kc) - vy(im,jc,kc))/dxm(ic))
                            msly(ic,kc) = msly(ic,kc) - (fc*(vx(ic,jc,kc) - vx(ic,jm,kc))/dym(jc))
                            msly(ic,kc) = msly(ic,kc) + (fp*(vy(ic,jp,kc) - vy(im,jp,kc))/dxm(ic))
                            msly(ic,kc) = msly(ic,kc) - (fp*(vx(ic,jp,kc) - vx(ic,jc,kc))/dym(jp))

                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/oz")
                    call HdfParallelWriteReal2D(filename,dsetname,msly(st1:en1,st3:en3),st1,en1,st3,en3,1,nx,0,nz,mpi_ycut)

                end if

            end if
            
            call MpiBarrier

        end do

    end if

    ! For Z-slices
    if (allocated(slcz)) then

        filename = trim("Results/movie_slices_z.h5")
        if (ismaster) then
            dsetname = trim("/steps")
            call HdfSerialWriteIntScalar(filename,dsetname,snpm)
            dsetname = trim("/data/"//charsnap//"/time")
            call HdfSerialWriteRealScalar(filename,dsetname,time)
            dsetname = trim("/data/"//charsnap//"/dt")
            call HdfSerialWriteRealScalar(filename,dsetname,dt)
        end if

        call MpiBarrier

        do lc = 1,size(slcz)

            write(charindx,"(I4.4)") lc

            ! Assume that lvlhalo is at least 1
            if (((zc(xstart(3)).le.slcz(lc)).and.(zc(xend(3)+1).gt.slcz(lc))).or.((slcz(lc).eq.zlen).and.(xend(3).eq.nzm))) then

                ! For X-Velocity

                if (sqnz(lc,1)) then

                    st1 = 1
                    en1 = nx

                    st2 = xstart(2)
                    en2 = xend(2)

                    if (xstart(2).eq.1) st2 = 0
                    if (xend(2).eq.nym) en2 = ny
                
                    fc = scfz(lc,3)
                    fp = scfz(lc,4)
                    kc = sidz(lc,3)
                    kp = sidz(lc,4)
                    
                    do jc = st2,en2
                        do ic = st1,en1
                            mslz(ic,jc) = (fc*vx(ic,jc,kc) + fp*vx(ic,jc,kc))
                        end do
                    end do

                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/vx")
                    call HdfParallelWriteReal2D(filename,dsetname,mslz(st1:en1,st2:en2),st1,en1,st2,en2,1,nx,0,ny,mpi_zcut)

                end if

                ! For Y-Velocity

                if (sqnz(lc,2)) then

                    st1 = 0
                    en1 = nx

                    st2 = xstart(2)
                    en2 = xend(2)

                    if (xend(2).eq.nym) en2 = ny
                    
                    fc = scfz(lc,3)
                    fp = scfz(lc,4)
                    kc = sidz(lc,3)
                    kp = sidz(lc,4)

                    do jc = st2,en2
                        do ic = st1,en1
                            mslz(ic,jc) = (fc*vy(ic,jc,kc) + fp*vy(ic,jc,kp))
                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/vy")
                    call HdfParallelWriteReal2D(filename,dsetname,mslz(st1:en1,st2:en2),st1,en1,st2,en2,0,nx,1,ny,mpi_zcut)

                end if

                ! For Z-Velocity

                if (sqnz(lc,3)) then

                    st1 = 0
                    en1 = nx

                    st2 = xstart(2)
                    en2 = xend(2)

                    if (xstart(2).eq.1) st2 = 0
                    if (xend(2).eq.nym) en2 = ny
                    
                    fc = scfz(lc,1)
                    fp = scfz(lc,2)
                    kc = sidz(lc,1)
                    kp = sidz(lc,2)

                    do jc = st2,en2
                        do ic = st1,en1
                            mslz(ic,jc) = (fc*vz(ic,jc,kc) + fp*vz(ic,jc,kp))
                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/vz")
                    call HdfParallelWriteReal2D(filename,dsetname,mslz(st1:en1,st2:en2),st1,en1,st2,en2,0,nx,0,ny,mpi_zcut)

                end if

                ! For Pressure

                if (sqnz(lc,4)) then

                    st1 = 0
                    en1 = nx

                    st2 = xstart(2)
                    en2 = xend(2)

                    if (xstart(2).eq.1) st2 = 0
                    if (xend(2).eq.nym) en2 = ny
                    
                    fc = scfz(lc,3)
                    fp = scfz(lc,4)
                    kc = sidz(lc,3)
                    kp = sidz(lc,4)

                    do jc = st2,en2
                        do ic = st1,en1
                            mslz(ic,jc) = (fc*pr(ic,jc,kc) + fp*pr(ic,jc,kp))
                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/pr")
                    call HdfParallelWriteReal2D(filename,dsetname,mslz(st1:en1,st2:en2),st1,en1,st2,en2,0,nx,0,ny,mpi_zcut)

                end if

                ! For X-Vorticity

                if (sqnz(lc,5)) then

                    st1 = 0
                    en1 = nx

                    st2 = xstart(2)
                    en2 = xend(2)

                    if (xend(2).eq.nym) en2 = ny
                    
                    fc = scfz(lc,1)
                    fp = scfz(lc,2)
                    kc = sidz(lc,1)
                    kp = sidz(lc,2)

                    km = kc-1
                    do jc = st2,en2
                        jm = jc-1
                        jp = jc+1
                        do ic = st1,en1
                            im = ic-1
                            ip = ic+1

                            mslz(ic,jc) = 0.0
                            mslz(ic,jc) = mslz(ic,jc) + (fc*(vz(ic,jc,kc) - vz(ic,jm,kc))/dym(jc))
                            mslz(ic,jc) = mslz(ic,jc) - (fc*(vy(ic,jc,kc) - vy(ic,jc,km))/dzm(kc))
                            mslz(ic,jc) = mslz(ic,jc) + (fp*(vz(ic,jc,kp) - vz(ic,jm,kp))/dym(jc))
                            mslz(ic,jc) = mslz(ic,jc) - (fp*(vy(ic,jc,kp) - vy(ic,jc,kc))/dzm(kp))

                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/ox")
                    call HdfParallelWriteReal2D(filename,dsetname,mslz(st1:en1,st2:en2),st1,en1,st2,en2,0,nx,1,ny,mpi_zcut)

                end if

                ! For Y-Vorticity

                if (sqnz(lc,6)) then

                    st1 = 1
                    en1 = nx

                    st2 = xstart(2)
                    en2 = xend(2)

                    if (xstart(2).eq.1) st2 = 0
                    if (xend(2).eq.nym) en2 = ny
                    
                    fc = scfz(lc,1)
                    fp = scfz(lc,2)
                    kc = sidz(lc,1)
                    kp = sidz(lc,2)

                    km = kc-1
                    do jc = st2,en2
                        jm = jc-1
                        jp = jc+1
                        do ic = st1,en1
                            im = ic-1
                            ip = ic+1

                            mslz(ic,jc) = 0.0
                            mslz(ic,jc) = mslz(ic,jc) + (fc*(vx(ic,jc,kc) - vx(ic,jc,km))/dzm(kc))
                            mslz(ic,jc) = mslz(ic,jc) - (fc*(vz(ic,jc,kc) - vz(im,jc,kc))/dxm(ic))
                            mslz(ic,jc) = mslz(ic,jc) + (fp*(vx(ic,jc,kp) - vx(ic,jc,kc))/dzm(kp))
                            mslz(ic,jc) = mslz(ic,jc) - (fp*(vz(ic,jc,kp) - vz(im,jc,kp))/dxm(ic))

                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/oy")
                    call HdfParallelWriteReal2D(filename,dsetname,mslz(st1:en1,st2:en2),st1,en1,st2,en2,1,nx,0,ny,mpi_zcut)

                end if

                ! For Z-Vorticity

                if (sqnz(lc,7)) then

                    st1 = 1
                    en1 = nx

                    st2 = xstart(2)
                    en2 = xend(2)

                    if (xend(2).eq.nym) en2 = ny
                    
                    fc = scfz(lc,3)
                    fp = scfz(lc,4)
                    kc = sidz(lc,3)
                    kp = sidz(lc,4)

                    do jc = st2,en2
                        jm = jc-1
                        jp = jc+1
                        do ic = st1,en1
                            im = ic-1
                            ip = ic+1

                            mslz(ic,jc) = 0.0
                            mslz(ic,jc) = mslz(ic,jc) + (fc*(vy(ic,jc,kc) - vy(im,jc,kc))/dxm(ic))
                            mslz(ic,jc) = mslz(ic,jc) - (fc*(vx(ic,jc,kc) - vx(ic,jm,kc))/dym(jc))
                            mslz(ic,jc) = mslz(ic,jc) + (fp*(vy(ic,jc,kp) - vy(im,jc,kp))/dxm(ic))
                            mslz(ic,jc) = mslz(ic,jc) - (fp*(vx(ic,jc,kp) - vx(ic,jm,kp))/dym(jc))

                        end do
                    end do
                    
                    dsetname = trim("/data/"//charsnap//"/"//charindx//"/oz")
                    call HdfParallelWriteReal2D(filename,dsetname,mslz(st1:en1,st2:en2),st1,en1,st2,en2,1,nx,1,ny,mpi_zcut)

                end if

            end if

            call MpiBarrier

        end do

    end if

    snpm = snpm + 1

    return

end subroutine WriteMovieSlices