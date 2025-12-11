!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: BoundaryRoutinesVY.F90                         !
!    CONTAINS: subroutine CalcBoundaryConditionVY         !
!              subroutine AddBoundaryRHSTermsVY           !
!              subroutine ImposeInternalBoundaryVY        !
!              subroutine ImposeExternalBoundaryVY        !
!                                                         !
!    PURPOSE: Compute boundary conditions for vy before   !
!    implicit solve step, and impose boundary conditions  !
!    after the implicit solve and pressure correction.    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcBoundaryConditionVY

    use decomp_2d
    use param
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc
    integer :: btyp
    real    :: bval

    ! Set the coefficients to be used for tridiagonal solvers
    !! For Dirichlet boundary condition
    coefxs(DIRICHLET,2)     = -1.0
    coefxe(DIRICHLET,2)     = -1.0
    coefys(DIRICHLET,2)     =  0.0
    coefye(DIRICHLET,2)     =  0.0
    coefzs(DIRICHLET,2)     = -1.0
    coefze(DIRICHLET,2)     = -1.0
    !! For Neumann boundary condition
    coefxs(NEUMANN,2)       =  1.0
    coefxe(NEUMANN,2)       =  1.0
    coefys(NEUMANN,2)       =  1.0
    coefye(NEUMANN,2)       =  1.0
    coefzs(NEUMANN,2)       =  1.0
    coefze(NEUMANN,2)       =  1.0
    !! For non-reflecting/advective/Sommerfeld boundary condition
    coefxs(SOMMERFELD,2)    =  (cvel*al*dt)/((cvel*al*dt) + dxm(1  ))
    coefxe(SOMMERFELD,2)    =  (cvel*al*dt)/((cvel*al*dt) + dxm(nx ))
    coefys(SOMMERFELD,2)    =  (cvel*al*dt)/((cvel*al*dt) + dyc(1  ))
    coefye(SOMMERFELD,2)    =  (cvel*al*dt)/((cvel*al*dt) + dyc(nym))
    coefzs(SOMMERFELD,2)    =  (cvel*al*dt)/((cvel*al*dt) + dzm(1  ))
    coefze(SOMMERFELD,2)    =  (cvel*al*dt)/((cvel*al*dt) + dzm(nz ))

    ! Store coefficient arrays for implicit solver 
    !! For x-pencils
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                call VyBcXs(jc,kc,btyp,bval)
                cfbcxs(jc,kc,2) = coefxs(btyp,2)
                call VyBcXe(jc,kc,btyp,bval)
                cfbcxe(jc,kc,2) = coefxe(btyp,2)
            end do
        end do
    end if
    !! For y-pencils
    if (.not.periodic(2)) then
        do kc = ystart(3),yend(3)
            do ic = ystart(1),yend(1)
                call VyBcYs(ic,kc,btyp,bval)
                cfbcys(ic,kc,2) = coefys(btyp,2)
                call VyBcYe(ic,kc,btyp,bval)
                cfbcye(ic,kc,2) = coefye(btyp,2)
            end do
        end do
    end if
    !! For z-pencils
    if (.not.periodic(3)) then
        do jc = zstart(2),zend(2)
            do ic = zstart(1),zend(1)
                call VyBcZs(ic,jc,btyp,bval)
                cfbczs(ic,jc,2) = coefzs(btyp,2)
                call VyBcZe(ic,jc,btyp,bval)
                cfbcze(ic,jc,2) = coefze(btyp,2)
            end do
        end do
    end if

    ! Store boundary type and boundary value
    !! For x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
            do jc = xstart(2),xend(2)+lvlhalo
                call VyBcXs(jc,kc,btypxs(jc,kc,2), bvalxs(jc,kc,2))
                !!! Store this value needed for imposing B.C.s later
                if (btypxs(jc,kc,2).eq.SOMMERFELD) bvalxs(jc,kc,2) = cvel*al*dt/dxm(1  )
                call VyBcXe(jc,kc,btypxe(jc,kc,2), bvalxe(jc,kc,2))
                !!! Store this value needed for imposing B.C.s later
                if (btypxe(jc,kc,2).eq.SOMMERFELD) bvalxe(jc,kc,2) = cvel*al*dt/dxm(nx )
            end do
        end do
    end if
    !! For y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VyBcYs(ic,kc,btypys(ic,kc,2), bvalys(ic,kc,2))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypys(ic,kc,2).eq.SOMMERFELD) bvalys(ic,kc,2) = cvel*al*dt/dyc(1  )
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VyBcYe(ic,kc,btypye(ic,kc,2), bvalye(ic,kc,2))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypye(ic,kc,2).eq.SOMMERFELD) bvalye(ic,kc,2) = cvel*al*dt/dyc(nym)
                end do
            end do
        end if
    end if
    !! For z-boundaries
    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VyBcZs(ic,jc,btypzs(ic,jc,2), bvalzs(ic,jc,2))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypzs(ic,jc,2).eq.SOMMERFELD) bvalzs(ic,jc,2) = cvel*al*dt/dzm(1  )
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VyBcZe(ic,jc,btypze(ic,jc,2), bvalze(ic,jc,2))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypze(ic,jc,2).eq.SOMMERFELD) bvalze(ic,jc,2) = cvel*al*dt/dzm(nz )
                end do
            end do
        end if
    end if

end subroutine CalcBoundaryConditionVY

subroutine AddBoundaryRHSTermsVY

    use decomp_2d
    use param
    use local_arrays, only: vy,rhx
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc
    real    :: beta

    ! Match the premultiplier with the one used in ImplicitAndUpdateVy
    beta = 0.5*al*dt/rey

    ! Add RHS terms to the x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                !! Boundary condition at x-axis start
                if (btypxs(jc,kc,2).eq.DIRICHLET)  rhx(1  ,jc,kc) = rhx(1  ,jc,kc) + beta*am1ci(1  )*(2.0*bvalxs(jc,kc,2) - (vy(1,jc,kc) + vy(0,jc,kc)))
                if (btypxs(jc,kc,2).eq.NEUMANN)    rhx(1  ,jc,kc) = rhx(1  ,jc,kc) - beta*am1ci(1  )*(bvalxs(jc,kc,2)*dxm(1) - (vy(1,jc,kc) - vy(0,jc,kc)))
                if (btypxs(jc,kc,2).eq.SOMMERFELD) rhx(1  ,jc,kc) = rhx(1  ,jc,kc) + beta*am1ci(1  )*(coefxs(SOMMERFELD,2)*(vy(1,jc,kc) - vy(0,jc,kc)))
                !! Boundary condition at x-axis end
                if (btypxe(jc,kc,2).eq.DIRICHLET)  rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1ci(nxm)*(2.0*bvalxe(jc,kc,2) - (vy(nxm,jc,kc) + vy(nx,jc,kc)))
                if (btypxe(jc,kc,2).eq.NEUMANN)    rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1ci(nxm)*(bvalxe(jc,kc,2)*dxm(nx) + (vy(nxm,jc,kc) - vy(nx,jc,kc)))
                if (btypxe(jc,kc,2).eq.SOMMERFELD) rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1ci(nxm)*(coefxe(SOMMERFELD,2)*(vy(nxm,jc,kc) - vy(nx,jc,kc)))
            end do
        end do
    end if

    ! Add RHS terms to the y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis start
                    if (btypys(ic,kc,2).eq.DIRICHLET)  rhx(ic,2  ,kc) = rhx(ic,2  ,kc) + beta*am2sj(2  )*(bvalys(ic,kc,2) - vy(ic,1,kc))
                    if (btypys(ic,kc,2).eq.NEUMANN)    rhx(ic,2  ,kc) = rhx(ic,2  ,kc) - beta*am2sj(2  )*(bvalys(ic,kc,2)*dyc(1) - (vy(ic,2,kc) - vy(ic,1,kc)))
                    if (btypys(ic,kc,2).eq.SOMMERFELD) rhx(ic,2  ,kc) = rhx(ic,2  ,kc) + beta*am2sj(2  )*(coefys(SOMMERFELD,2)*(vy(ic,2,kc) - vy(ic,1,kc)))
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis end
                    if (btypye(ic,kc,2).eq.DIRICHLET)  rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2sj(nym)*(bvalye(ic,kc,2) - vy(ic,ny,kc))
                    if (btypye(ic,kc,2).eq.NEUMANN)    rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2sj(nym)*(bvalye(ic,kc,2)*dyc(nym) + (vy(ic,nym,kc) - vy(ic,ny,kc)))
                    if (btypye(ic,kc,2).eq.SOMMERFELD) rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2sj(nym)*(coefye(SOMMERFELD,2)*(vy(ic,nym,kc) - vy(ic,ny,kc)))
                end do
            end do
        end if
    end if

    ! Add RHS terms to the z-boundaries
    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis start
                    if (btypzs(ic,jc,2).eq.DIRICHLET)  rhx(ic,jc,1  ) = rhx(ic,jc,1  ) + beta*am3ck(1  )*(2.0*bvalzs(ic,jc,2) - (vy(ic,jc,1) + vy(ic,jc,0)))
                    if (btypzs(ic,jc,2).eq.NEUMANN)    rhx(ic,jc,1  ) = rhx(ic,jc,1  ) - beta*am3ck(1  )*(bvalzs(ic,jc,2)*dzm(1) - (vy(ic,jc,1) - vy(ic,jc,0)))
                    if (btypzs(ic,jc,2).eq.SOMMERFELD) rhx(ic,jc,1  ) = rhx(ic,jc,1  ) + beta*am3ck(1  )*(coefzs(SOMMERFELD,2)*(vy(ic,jc,1) - vy(ic,jc,0)))
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis end
                    if (btypze(ic,jc,2).eq.DIRICHLET)  rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3ck(nzm)*(2.0*bvalze(ic,jc,2) - (vy(ic,jc,nzm) + vy(ic,jc,nz)))
                    if (btypze(ic,jc,2).eq.NEUMANN)    rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3ck(nzm)*(bvalze(ic,jc,2)*dzm(nz) + (vy(ic,jc,nzm) - vy(ic,jc,nz)))
                    if (btypze(ic,jc,2).eq.SOMMERFELD) rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3ck(nzm)*(coefze(SOMMERFELD,2)*(vy(ic,jc,nzm) - vy(ic,jc,nz)))
                end do
            end do
        end if
    end if
        
end subroutine AddBoundaryRHSTermsVY

subroutine ImposeInternalBoundaryVY

    use decomp_2d
    use param
    use local_arrays, only: vy
    use boundary_arrays

    implicit none

    integer :: ic,kc

    ! Impose boundary values on x-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis start
                    if (btypys(ic,kc,2).eq.DIRICHLET) then 
                        vy(ic,1  ,kc) = bvalys(ic,kc,2)
                        glpcys(ic,kc) = .false.
                    else if (btypys(ic,kc,2).eq.NEUMANN) then 
                        vy(ic,1  ,kc) = (vy(ic,2,kc) - bvalys(ic,kc,2)*dyc(1))
                        glpcys(ic,kc) = .true.
                    else if (btypys(ic,kc,2).eq.SOMMERFELD) then 
                        vy(ic,1  ,kc) = (vy(ic,2,kc)*bvalys(ic,kc,2) + vy(ic,1,kc))/(1.0 + bvalys(ic,kc,2))
                        glpcys(ic,kc) = .true.
                    end if
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis start
                    if (btypye(ic,kc,2).eq.DIRICHLET) then 
                        vy(ic,ny ,kc) = bvalye(ic,kc,2)
                        glpcye(ic,kc) = .false.
                    else if (btypye(ic,kc,2).eq.NEUMANN) then 
                        vy(ic,ny ,kc) = (vy(ic,nym,kc) + bvalye(ic,kc,2)*dyc(nym))
                        glpcye(ic,kc) = .true.
                    else if (btypye(ic,kc,2).eq.SOMMERFELD) then 
                        vy(ic,ny ,kc) = (vy(ic,nym,kc)*bvalye(ic,kc,2) + vy(ic,ny,kc))/(1.0 + bvalye(ic,kc,2))
                        glpcye(ic,kc) = .true.
                    end if
                end do
            end do
        end if
    end if

end subroutine ImposeInternalBoundaryVY

subroutine ImposeExternalBoundaryVY

    use decomp_2d
    use param
    use local_arrays, only: vy
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc

    ! Impose boundary values on x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
            do jc = xstart(2),xend(2)+lvlhalo
                !! Boundary condition at x-axis start
                if (btypxs(jc,kc,2).eq.DIRICHLET)  vy(0 ,jc,kc) = (2.0*bvalxs(jc,kc,2) - vy(1,jc,kc))
                if (btypxs(jc,kc,2).eq.NEUMANN)    vy(0 ,jc,kc) = (vy(1,jc,kc) - bvalxs(jc,kc,2)*dxm(1))
                if (btypxs(jc,kc,2).eq.SOMMERFELD) vy(0 ,jc,kc) = (vy(1,jc,kc)*bvalxs(jc,kc,2) + vy(0,jc,kc))/(1.0 + bvalxs(jc,kc,2))
                !! Boundary condition at x-axis end
                if (btypxe(jc,kc,2).eq.DIRICHLET)  vy(nx,jc,kc) = (2.0*bvalxe(jc,kc,2) - vy(nxm,jc,kc))
                if (btypxe(jc,kc,2).eq.NEUMANN)    vy(nx,jc,kc) = (vy(nxm,jc,kc) + bvalxe(jc,kc,2)*dxm(nx))
                if (btypxe(jc,kc,2).eq.SOMMERFELD) vy(nx,jc,kc) = (vy(nxm,jc,kc)*bvalxe(jc,kc,2) + vy(nx,jc,kc))/(1.0 + bvalxe(jc,kc,2))
            end do
        end do
    end if

    ! Impose boundary values on z-boundaries
    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    !! Boundary condition at z-axis start
                    if (btypzs(ic,jc,2).eq.DIRICHLET)  vy(ic,jc,0 ) = (2.0*bvalzs(ic,jc,2) - vy(ic,jc,1))
                    if (btypzs(ic,jc,2).eq.NEUMANN)    vy(ic,jc,0 ) = (vy(ic,jc,1) - bvalzs(ic,jc,2)*dzm(1))
                    if (btypzs(ic,jc,2).eq.SOMMERFELD) vy(ic,jc,0 ) = (vy(ic,jc,1)*bvalzs(ic,jc,2) + vy(ic,jc,0))/(1.0 + bvalzs(ic,jc,2))
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    !! Boundary condition at z-axis end
                    if (btypze(ic,jc,2).eq.DIRICHLET)  vy(ic,jc,nz) = (2.0*bvalze(ic,jc,2) - vy(ic,jc,nzm))
                    if (btypze(ic,jc,2).eq.NEUMANN)    vy(ic,jc,nz) = (vy(ic,jc,nzm) + bvalze(ic,jc,2)*dzm(nz))
                    if (btypze(ic,jc,2).eq.SOMMERFELD) vy(ic,jc,nz) = (vy(ic,jc,nzm)*bvalze(ic,jc,2) + vy(ic,jc,nz))/(1.0 + bvalze(ic,jc,2))
                end do
            end do
        end if
    end if

end subroutine ImposeExternalBoundaryVY