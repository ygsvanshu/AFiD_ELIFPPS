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
    use boundary_arrays, only: DIRICHLET, NEUMANN, glpcxs, glpcxe, glpcys, glpcye, glpczs, glpcze
    use pressure_solver, only: pcbcxs, pcbcxe, pcbcys, pcbcye, pcbczs, pcbcze

    implicit none

    real                :: gdiv,gcar
    integer             :: gcnt
    integer             :: kc,kp
    integer             :: jc,jp
    integer             :: ic,ip
    real                :: res_dummy
    integer             :: ies_dummy
    logical             :: correction

    ! Check if the global divergence correction should be applied
    ! (Only if every non-periodic direction has Neumann-Neumann pressure correction
    ! boundary conditions, and at least one non-periodic direction exists)
    correction = .not.all(periodic)
    if (.not.periodic(1)) correction = correction.and.(pcbcxs.eq.NEUMANN).and.(pcbcxe.eq.NEUMANN)
    if (.not.periodic(2)) correction = correction.and.(pcbcys.eq.NEUMANN).and.(pcbcye.eq.NEUMANN)
    if (.not.periodic(3)) correction = correction.and.(pcbczs.eq.NEUMANN).and.(pcbcze.eq.NEUMANN)

    if (correction) then

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

        if (.not.periodic(1)) then
            do kc=xstart(3),xend(3)
                do jc=xstart(2),xend(2)
                    gcar = gcar + dyc(jc)*dzc(kc)*glpcxs(jc,kc)
                    gcar = gcar + dyc(jc)*dzc(kc)*glpcxe(jc,kc)
                enddo
            enddo
        end if

        if (.not.periodic(2)) then
            if (xstart(2).eq.1) then
                do kc=xstart(3),xend(3)
                    do ic=xstart(1),xend(1)
                        gcar = gcar + dxc(ic)*dzc(kc)*glpcys(ic,kc)
                    enddo
                enddo
            end if
            if (xend(2).eq.nym) then
                do kc=xstart(3),xend(3)
                    do ic=xstart(1),xend(1)
                        gcar = gcar + dxc(ic)*dzc(kc)*glpcye(ic,kc)
                    enddo
                enddo
            endif
        end if

        if (.not.periodic(3)) then
            if (xstart(3).eq.1) then
                do jc=xstart(2),xend(2)
                    do ic=xstart(1),xend(1)
                        gcar = gcar + dxc(ic)*dyc(jc)*glpczs(ic,jc)
                    enddo
                enddo
            end if
            if (xend(3).eq.nzm) then
                do jc=xstart(2),xend(2)
                    do ic=xstart(1),xend(1)
                        gcar = gcar + dxc(ic)*dyc(jc)*glpcze(ic,jc)
                    enddo
                enddo
            endif
        end if

        call MpiAllSumRealScalar(gcar,res_dummy)
        gcar = res_dummy

        ! Check if the division blows up
        if (abs(gcar).gt.(gdiv/limitVel)) then

            if (.not.periodic(1)) then
                do kc=xstart(3),xend(3)
                    do jc=xstart(2),xend(2)
                        vx(1 ,jc,kc) = vx(1 ,jc,kc) + (gdiv/gcar)*glpcxs(jc,kc)
                        vx(nx,jc,kc) = vx(nx,jc,kc) - (gdiv/gcar)*glpcxe(jc,kc)
                    enddo
                enddo
            end if

            if (.not.periodic(2)) then
                if (xstart(2).eq.1) then
                    do kc=xstart(3),xend(3)
                        do ic=xstart(1),xend(1)
                            vy(ic,1 ,kc) = vy(ic,1 ,kc) + (gdiv/gcar)*glpcys(ic,kc)
                        enddo
                    enddo
                end if
                if (xend(2).eq.nym) then
                    do kc=xstart(3),xend(3)
                        do ic=xstart(1),xend(1)
                            vy(ic,ny,kc) = vy(ic,ny,kc) - (gdiv/gcar)*glpcye(ic,kc)
                        enddo
                    enddo
                endif
            end if

            if (.not.periodic(3)) then
                if (xstart(3).eq.1) then
                    do jc=xstart(2),xend(2)
                        do ic=xstart(1),xend(1)
                            vz(ic,jc,1 ) = vz(ic,jc,1 ) + (gdiv/gcar)*glpczs(ic,jc)
                        enddo
                    enddo
                end if
                if (xend(3).eq.nzm) then
                    do jc=xstart(2),xend(2)
                        do ic=xstart(1),xend(1)
                            vz(ic,jc,nz) = vz(ic,jc,nz) - (gdiv/gcar)*glpcze(ic,jc)
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

    end if

    return

end subroutine CorrectGlobalDivergence