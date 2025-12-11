!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: CorrectGlobalDivergence.F90                    !
!    CONTAINS: subroutine CorrectGlobalDivergence         !
!                                                         ! 
!    PURPOSE: Calculate the global divergence of velocity !
!    in the domain and correct it by applying boundary    !
!    fluxes.                                              !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CorrectGlobalDivergence

    use mpih
    use decomp_2d, only: xstart,xend
    use param
    use local_arrays, only: vx,vy,vz
    use boundary_arrays, only: glpcxs,glpcxe,glpcys,glpcye,glpczs,glpcze

    implicit none

    real                :: gdiv,gcar
    integer             :: gcnt
    integer             :: kc,kp
    integer             :: jc,jp
    integer             :: ic,ip
    real                :: res_dummy
    integer             :: ies_dummy

    gdiv = 0.0
    gcar = 0.0
        
    do kc=xstart(3),xend(3)
        kp=kc+1
        do jc=xstart(2),xend(2)
            jp=jc+1
            do ic=xstart(1),xend(1)
                ip=ic+1
            
                gdiv = gdiv + ((vx(ip,jc,kc)-vx(ic,jc,kc))*dyc(jc)*dzc(kc) + (vy(ic,jp,kc)-vy(ic,jc,kc))*dxc(ic)*dzc(kc) + (vz(ic,jc,kp)-vz(ic,jc,kc))*dxc(ic)*dyc(jc))
            
            enddo
        enddo
    enddo

    call MpiAllSumRealScalar(gdiv,res_dummy)
    gdiv = res_dummy

    ! if (ismaster) print *, "PRE ---- ", "GDIV = ", gdiv, "GCAR = ",gcar 

    gcnt = 0

    if (.not.periodic(1)) then
        do kc=xstart(3),xend(3)
            do jc=xstart(2),xend(2)
                if (glpcxs(jc,kc)) then
                    gcar = gcar + dyc(jc)*dzc(kc)
                    gcnt = gcnt + 1
                end if
                if (glpcxe(jc,kc)) then 
                    gcar = gcar + dyc(jc)*dzc(kc)
                    gcnt = gcnt + 1
                end if
            enddo
        enddo
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc=xstart(3),xend(3)
                do ic=xstart(1),xend(1)
                    if (glpcys(ic,kc)) then
                        gcar = gcar + dxc(ic)*dzc(kc)
                        gcnt = gcnt + 1
                    end if
                enddo
            enddo
        end if
        if (xend(2).eq.nym) then
            do kc=xstart(3),xend(3)
                do ic=xstart(1),xend(1)
                    if (glpcye(ic,kc)) then
                        gcar = gcar + dxc(ic)*dzc(kc)
                        gcnt = gcnt + 1
                    end if
                enddo
            enddo
        endif
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc=xstart(2),xend(2)
                do ic=xstart(1),xend(1)
                    if (glpczs(ic,jc)) then
                        gcar = gcar + dxc(ic)*dyc(jc)
                        gcnt = gcnt + 1
                    end if
                enddo
            enddo
        end if
        if (xend(3).eq.nzm) then
            do jc=xstart(2),xend(2)
                do ic=xstart(1),xend(1)
                    if (glpcze(ic,jc)) then
                        gcar = gcar + dxc(ic)*dyc(jc)
                        gcnt = gcnt + 1
                    end if
                enddo
            enddo
        endif
    end if

    call MpiAllSumRealScalar(gcar,res_dummy)
    gcar = res_dummy

    call MpiAllSumIntScalar(gcnt,ies_dummy)
    gcnt = ies_dummy

    if (gcnt.gt.0) then

        if (.not.periodic(1)) then
            do kc=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    if (glpcxs(jc,kc)) vx(1 ,jc,kc) = vx(1 ,jc,kc) + gdiv/gcar
                    if (glpcxe(jc,kc)) vx(nx,jc,kc) = vx(nx,jc,kc) - gdiv/gcar
                enddo
            enddo
        end if

        if (.not.periodic(2)) then
            if (xstart(2).eq.1) then
                do kc=xstart(3),xend(3)
                    do ic=xstart(1),xend(1)
                        if (glpcys(ic,kc)) vy(ic,1 ,kc) = vy(ic,1 ,kc) + gdiv/gcar
                    enddo
                enddo
            end if
            if (xend(2).eq.nym) then
                do kc=xstart(3),xend(3)
                    do ic=xstart(1),xend(1)
                        if (glpcye(ic,kc)) vy(ic,ny,kc) = vy(ic,ny,kc) - gdiv/gcar
                    enddo
                enddo
            endif
        end if

        if (.not.periodic(3)) then
            if (xstart(3).eq.1) then
                do jc=xstart(2),xend(2)
                    do ic=xstart(1),xend(1)
                        if (glpczs(ic,jc)) vz(ic,jc,1 ) = vz(ic,jc,1 ) + gdiv/gcar
                    enddo
                enddo
            end if
            if (xend(3).eq.nzm) then
                do jc=xstart(2),xend(2)
                    do ic=xstart(1),xend(1)
                        if (glpcze(ic,jc)) vz(ic,jc,nz) = vz(ic,jc,nz) - gdiv/gcar
                    enddo
                enddo
            endif
        end if

    else

        if (gdiv.gt.resid) then
            if (ismaster) print *, 'Global divergence is non-zero'
            if (ismaster) print *, 'Please check boundary conditions'
            call MPI_ABORT(MPI_COMM_WORLD, 1, mpi_ierr)
        end if

    end if

    return

end subroutine CorrectGlobalDivergence