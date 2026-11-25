!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: BoundaryRoutinesPC.F90                         !
!    CONTAINS: subroutine CalcBoundaryConditionPC         !
!              subroutine AddBoundaryRHSTermsPC           !
!              subroutine ImposeBoundaryConditionPC       !
!                                                         !
!    PURPOSE: Compute boundary conditions for the         !
!    pressure correction before the pressure solve, add   !
!    the corresponding RHS terms to the pressure Poisson  !
!    equation, and impose them on the pressure correction !
!    after the pressure solve.                            !
!                                                         !
!    The ghost value of the pressure correction is set as !
!    dph(ghost) = s*dph(inner) + pcor, with s = -1 for    !
!    DIRICHLET and s = +1 for NEUMANN. The offsets "pcor" !
!    are chosen such that after the pressure update, the  !
!    face pressure (DIRICHLET) or the face-normal         !
!    pressure gradient (NEUMANN) equals the value set by  !
!    the PrBc(X/Y/Z)(s/e) subroutines.                    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcBoundaryConditionPC

    use decomp_2d
    use param
    use local_arrays, only: pr
    use boundary_arrays
    use pressure_solver, only: pcbcxs, pcbcxe, pcbcys, pcbcye, pcbczs, pcbcze

    implicit none

    integer :: ic,jc,kc
    real    :: pval

    ! DIRICHLET: twice the deficit between the required and current face pressure
    ! NEUMANN:   (current - required) face-normal pressure gradient, times the spacing
    ! The DIRICHLET offsets are completed in AddBoundaryRHSTermsPC

    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                !! Boundary condition at x-axis start
                call PrBcXs(jc,kc,pval)
                if (pcbcxs.eq.DIRICHLET) pcorxs(jc,kc) = 2.0*pval - (pr(1  ,jc,kc) + pr(0  ,jc,kc))
                if (pcbcxs.eq.NEUMANN)   pcorxs(jc,kc) = (pr(1  ,jc,kc) - pr(0  ,jc,kc)) - pval*dxm(1 )
                !! Boundary condition at x-axis end
                call PrBcXe(jc,kc,pval)
                if (pcbcxe.eq.DIRICHLET) pcorxe(jc,kc) = 2.0*pval - (pr(nxm,jc,kc) + pr(nx ,jc,kc))
                if (pcbcxe.eq.NEUMANN)   pcorxe(jc,kc) = pval*dxm(nx) - (pr(nx ,jc,kc) - pr(nxm,jc,kc))
            end do
        end do
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis start
                    call PrBcYs(ic,kc,pval)
                    if (pcbcys.eq.DIRICHLET) pcorys(ic,kc) = 2.0*pval - (pr(ic,1  ,kc) + pr(ic,0  ,kc))
                    if (pcbcys.eq.NEUMANN)   pcorys(ic,kc) = (pr(ic,1  ,kc) - pr(ic,0  ,kc)) - pval*dym(1 )
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis end
                    call PrBcYe(ic,kc,pval)
                    if (pcbcye.eq.DIRICHLET) pcorye(ic,kc) = 2.0*pval - (pr(ic,nym,kc) + pr(ic,ny ,kc))
                    if (pcbcye.eq.NEUMANN)   pcorye(ic,kc) = pval*dym(ny) - (pr(ic,ny ,kc) - pr(ic,nym,kc))
                end do
            end do
        end if
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis start
                    call PrBcZs(ic,jc,pval)
                    if (pcbczs.eq.DIRICHLET) pcorzs(ic,jc) = 2.0*pval - (pr(ic,jc,1  ) + pr(ic,jc,0  ))
                    if (pcbczs.eq.NEUMANN)   pcorzs(ic,jc) = (pr(ic,jc,1  ) - pr(ic,jc,0  )) - pval*dzm(1 )
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis end
                    call PrBcZe(ic,jc,pval)
                    if (pcbcze.eq.DIRICHLET) pcorze(ic,jc) = 2.0*pval - (pr(ic,jc,nzm) + pr(ic,jc,nz ))
                    if (pcbcze.eq.NEUMANN)   pcorze(ic,jc) = pval*dzm(nz) - (pr(ic,jc,nz ) - pr(ic,jc,nzm))
                end do
            end do
        end if
    end if

end subroutine CalcBoundaryConditionPC

subroutine AddBoundaryRHSTermsPC

    use decomp_2d
    use param
    use local_arrays, only: rhx
    use boundary_arrays
    use pressure_solver, only: pcbcxs, pcbcxe, pcbcys, pcbcye, pcbczs, pcbcze

    implicit none

    integer :: ic,jc,kc
    real    :: be

    ! Match the premultiplier with the one used in CorrectPressure
    be = 0.5*al*dt/rey

    ! Complete the DIRICHLET offsets with the rotational pressure update, so that the
    ! face pressure attains its boundary value exactly after CorrectPressure. At the
    ! boundary cells, the Laplacian of the pressure correction equals the unmodified
    ! RHS, so this must be done before any RHS terms are added (edges and corners).

    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                if (pcbcxs.eq.DIRICHLET) pcorxs(jc,kc) = pcorxs(jc,kc) + 2.0*be*rhx(1  ,jc,kc)
                if (pcbcxe.eq.DIRICHLET) pcorxe(jc,kc) = pcorxe(jc,kc) + 2.0*be*rhx(nxm,jc,kc)
            end do
        end do
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    if (pcbcys.eq.DIRICHLET) pcorys(ic,kc) = pcorys(ic,kc) + 2.0*be*rhx(ic,1  ,kc)
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    if (pcbcye.eq.DIRICHLET) pcorye(ic,kc) = pcorye(ic,kc) + 2.0*be*rhx(ic,nym,kc)
                end do
            end do
        end if
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    if (pcbczs.eq.DIRICHLET) pcorzs(ic,jc) = pcorzs(ic,jc) + 2.0*be*rhx(ic,jc,1  )
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    if (pcbcze.eq.DIRICHLET) pcorze(ic,jc) = pcorze(ic,jc) + 2.0*be*rhx(ic,jc,nzm)
                end do
            end do
        end if
    end if

    ! Add RHS terms from the ghost value offsets (same form for DIRICHLET and NEUMANN)

    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                rhx(1  ,jc,kc) = rhx(1  ,jc,kc) - amphi(1  )*pcorxs(jc,kc)
                rhx(nxm,jc,kc) = rhx(nxm,jc,kc) - apphi(nxm)*pcorxe(jc,kc)
            end do
        end do
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    rhx(ic,1  ,kc) = rhx(ic,1  ,kc) - amphj(1  )*pcorys(ic,kc)
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    rhx(ic,nym,kc) = rhx(ic,nym,kc) - apphj(nym)*pcorye(ic,kc)
                end do
            end do
        end if
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    rhx(ic,jc,1  ) = rhx(ic,jc,1  ) - amphk(1  )*pcorzs(ic,jc)
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    rhx(ic,jc,nzm) = rhx(ic,jc,nzm) - apphk(nzm)*pcorze(ic,jc)
                end do
            end do
        end if
    end if

end subroutine AddBoundaryRHSTermsPC

subroutine ImposeBoundaryConditionPC

    use decomp_2d
    use param
    use local_arrays, only: dphhalo
    use boundary_arrays
    use pressure_solver, only: pcbcxs, pcbcxe, pcbcys, pcbcye, pcbczs, pcbcze

    implicit none

    integer :: ic,jc,kc

    ! Halo exchange doesn't handle the pencils in x-direction, manually corrected.
    if (periodic(1)) then
        dphhalo(0 ,:,:) = dphhalo(nxm,:,:)
        dphhalo(nx,:,:) = dphhalo(1  ,:,:)
    else
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                !! Boundary condition at x-axis start
                if (pcbcxs.eq.DIRICHLET) dphhalo(0 ,jc,kc) = -dphhalo(1  ,jc,kc) + pcorxs(jc,kc)
                if (pcbcxs.eq.NEUMANN)   dphhalo(0 ,jc,kc) =  dphhalo(1  ,jc,kc) + pcorxs(jc,kc)
                !! Boundary condition at x-axis end
                if (pcbcxe.eq.DIRICHLET) dphhalo(nx,jc,kc) = -dphhalo(nxm,jc,kc) + pcorxe(jc,kc)
                if (pcbcxe.eq.NEUMANN)   dphhalo(nx,jc,kc) =  dphhalo(nxm,jc,kc) + pcorxe(jc,kc)
            end do
        end do
    end if

    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis start
                    if (pcbcys.eq.DIRICHLET) dphhalo(ic,0 ,kc) = -dphhalo(ic,1  ,kc) + pcorys(ic,kc)
                    if (pcbcys.eq.NEUMANN)   dphhalo(ic,0 ,kc) =  dphhalo(ic,1  ,kc) + pcorys(ic,kc)
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis end
                    if (pcbcye.eq.DIRICHLET) dphhalo(ic,ny,kc) = -dphhalo(ic,nym,kc) + pcorye(ic,kc)
                    if (pcbcye.eq.NEUMANN)   dphhalo(ic,ny,kc) =  dphhalo(ic,nym,kc) + pcorye(ic,kc)
                end do
            end do
        end if
    end if

    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis start
                    if (pcbczs.eq.DIRICHLET) dphhalo(ic,jc,0 ) = -dphhalo(ic,jc,1  ) + pcorzs(ic,jc)
                    if (pcbczs.eq.NEUMANN)   dphhalo(ic,jc,0 ) =  dphhalo(ic,jc,1  ) + pcorzs(ic,jc)
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis end
                    if (pcbcze.eq.DIRICHLET) dphhalo(ic,jc,nz) = -dphhalo(ic,jc,nzm) + pcorze(ic,jc)
                    if (pcbcze.eq.NEUMANN)   dphhalo(ic,jc,nz) =  dphhalo(ic,jc,nzm) + pcorze(ic,jc)
                end do
            end do
        end if
    end if

end subroutine ImposeBoundaryConditionPC
