!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: BoundaryRoutinesVX.F90                         !
!    CONTAINS: subroutine CalcBoundaryConditionVX         !
!              subroutine AddBoundaryRHSTermsVX           !
!              subroutine ImposeInternalBoundaryVX        !
!              subroutine ImposeExternalBoundaryVX        !
!                                                         !
!    PURPOSE: Compute boundary conditions for vx before   !
!    implicit solve step, and impose boundary conditions  !
!    after the implicit solve and pressure correction.    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcBoundaryConditionVX

    use decomp_2d
    use param
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc
    integer :: btyp
    real    :: bval

    ! Set the coefficients to be used for tridiagonal solvers
    !! For Dirichlet boundary condition
    coefxs(DIRICHLET,1)     =  0.0
    coefxe(DIRICHLET,1)     =  0.0
    coefys(DIRICHLET,1)     = -1.0
    coefye(DIRICHLET,1)     = -1.0
    coefzs(DIRICHLET,1)     = -1.0
    coefze(DIRICHLET,1)     = -1.0
    !! For Neumann boundary condition
    coefxs(NEUMANN,1)       =  1.0
    coefxe(NEUMANN,1)       =  1.0
    coefys(NEUMANN,1)       =  1.0
    coefye(NEUMANN,1)       =  1.0
    coefzs(NEUMANN,1)       =  1.0
    coefze(NEUMANN,1)       =  1.0
    !! For non-reflecting/advective/Sommerfeld boundary condition
    coefxs(SOMMERFELD,1)    =  (cvel*al*dt)/((cvel*al*dt) + dxc(1  ))
    coefxe(SOMMERFELD,1)    =  (cvel*al*dt)/((cvel*al*dt) + dxc(nxm))
    coefys(SOMMERFELD,1)    =  (cvel*al*dt)/((cvel*al*dt) + dym(1  ))
    coefye(SOMMERFELD,1)    =  (cvel*al*dt)/((cvel*al*dt) + dym(ny ))
    coefzs(SOMMERFELD,1)    =  (cvel*al*dt)/((cvel*al*dt) + dzm(1  ))
    coefze(SOMMERFELD,1)    =  (cvel*al*dt)/((cvel*al*dt) + dzm(nz ))

    ! Store coefficient arrays for implicit solver 
    !! For x-pencils
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                call VxBcXs(jc,kc,btyp,bval)
                cfbcxs(jc,kc,1) = coefxs(btyp,1)
                call VxBcXe(jc,kc,btyp,bval)
                cfbcxe(jc,kc,1) = coefxe(btyp,1)
            end do
        end do
    end if
    !! For y-pencils
    if (.not.periodic(2)) then
        do kc = ystart(3),yend(3)
            do ic = ystart(1),yend(1)
                call VxBcYs(ic,kc,btyp,bval)
                cfbcys(ic,kc,1) = coefys(btyp,1)
                call VxBcYe(ic,kc,btyp,bval)
                cfbcye(ic,kc,1) = coefye(btyp,1)
            end do
        end do
    end if
    !! For z-pencils
    if (.not.periodic(3)) then
        do jc = zstart(2),zend(2)
            do ic = zstart(1),zend(1)
                call VxBcZs(ic,jc,btyp,bval)
                cfbczs(ic,jc,1) = coefzs(btyp,1)
                call VxBcZe(ic,jc,btyp,bval)
                cfbcze(ic,jc,1) = coefze(btyp,1)
            end do
        end do
    end if

    ! Store boundary type and boundary value
    !! For x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                call VxBcXs(jc,kc,btypxs(jc,kc,1), bvalxs(jc,kc,1))
                !!! Store this value needed for imposing B.C.s later
                if (btypxs(jc,kc,1).eq.SOMMERFELD) bvalxs(jc,kc,1) = cvel*al*dt/dxc(1  )
                call VxBcXe(jc,kc,btypxe(jc,kc,1), bvalxe(jc,kc,1))
                !!! Store this value needed for imposing B.C.s later
                if (btypxe(jc,kc,1).eq.SOMMERFELD) bvalxe(jc,kc,1) = cvel*al*dt/dxc(nxm)
            end do
        end do
    end if
    !! For y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    call VxBcYs(ic,kc,btypys(ic,kc,1), bvalys(ic,kc,1))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypys(ic,kc,1).eq.SOMMERFELD) bvalys(ic,kc,1) = cvel*al*dt/dym(1  )
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    call VxBcYe(ic,kc,btypye(ic,kc,1), bvalye(ic,kc,1))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypye(ic,kc,1).eq.SOMMERFELD) bvalye(ic,kc,1) = cvel*al*dt/dym(ny )
                end do
            end do
        end if
    end if
    !! For z-boundaries
    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    call VxBcZs(ic,jc,btypzs(ic,jc,1), bvalzs(ic,jc,1))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypzs(ic,jc,1).eq.SOMMERFELD) bvalzs(ic,jc,1) = cvel*al*dt/dzm(1  )
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    call VxBcZe(ic,jc,btypze(ic,jc,1),bvalze(ic,jc,1))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypze(ic,jc,1).eq.SOMMERFELD) bvalze(ic,jc,1) = cvel*al*dt/dzm(nz )
                end do
            end do
        end if
    end if

end subroutine CalcBoundaryConditionVX

subroutine AddBoundaryRHSTermsVX

    use decomp_2d
    use param
    use local_arrays, only: vx,rhx
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc
    real    :: beta

    ! Match the premultiplier with the one used in ImplicitAndUpdateVx
    beta = 0.5*al*dt/rey

    ! Add RHS terms to the x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                !! Boundary condition at x-axis start
                if (btypxs(jc,kc,1).eq.DIRICHLET)  rhx(2  ,jc,kc) = rhx(2  ,jc,kc) + beta*am1si(2  )*(bvalxs(jc,kc,1) - vx(1,jc,kc))
                if (btypxs(jc,kc,1).eq.NEUMANN)    rhx(2  ,jc,kc) = rhx(2  ,jc,kc) - beta*am1si(2  )*(bvalxs(jc,kc,1)*dxc(1) - (vx(2,jc,kc) - vx(1,jc,kc)))
                if (btypxs(jc,kc,1).eq.SOMMERFELD) rhx(2  ,jc,kc) = rhx(2  ,jc,kc) + beta*am1si(2  )*(coefxs(SOMMERFELD,1)*(vx(2,jc,kc) - vx(1,jc,kc)))
                !! Boundary condition at x-axis end
                if (btypxe(jc,kc,1).eq.DIRICHLET)  rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1si(nxm)*(bvalxe(jc,kc,1) - vx(nx,jc,kc))
                if (btypxe(jc,kc,1).eq.NEUMANN)    rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1si(nxm)*(bvalxe(jc,kc,1)*dxc(nxm) + (vx(nxm,jc,kc) - vx(nx,jc,kc)))
                if (btypxe(jc,kc,1).eq.SOMMERFELD) rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1si(nxm)*(coefxe(SOMMERFELD,1)*(vx(nxm,jc,kc) - vx(nx,jc,kc)))
            end do
        end do
    end if

    ! Add RHS terms to the y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis start
                    if (btypys(ic,kc,1).eq.DIRICHLET)  rhx(ic,1  ,kc) = rhx(ic,1  ,kc) + beta*am2cj(1  )*(2.0*bvalys(ic,kc,1) - (vx(ic,1,kc) + vx(ic,0,kc)))
                    if (btypys(ic,kc,1).eq.NEUMANN)    rhx(ic,1  ,kc) = rhx(ic,1  ,kc) - beta*am2cj(1  )*(bvalys(ic,kc,1)*dym(1) - (vx(ic,1,kc) - vx(ic,0,kc)))
                    if (btypys(ic,kc,1).eq.SOMMERFELD) rhx(ic,1  ,kc) = rhx(ic,1  ,kc) + beta*am2cj(1  )*(coefys(SOMMERFELD,1)*(vx(ic,1,kc) - vx(ic,0,kc)))
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis end
                    if (btypye(ic,kc,1).eq.DIRICHLET)  rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2cj(nym)*(2.0*bvalye(ic,kc,1) - (vx(ic,nym,kc) + vx(ic,ny,kc)))
                    if (btypye(ic,kc,1).eq.NEUMANN)    rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2cj(nym)*(bvalye(ic,kc,1)*dym(ny) + (vx(ic,nym,kc) - vx(ic,ny,kc)))
                    if (btypye(ic,kc,1).eq.SOMMERFELD) rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2cj(nym)*(coefye(SOMMERFELD,1)*(vx(ic,nym,kc) - vx(ic,ny,kc)))
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
                    if (btypzs(ic,jc,1).eq.DIRICHLET)  rhx(ic,jc,1  ) = rhx(ic,jc,1  ) + beta*am3ck(1  )*(2.0*bvalzs(ic,jc,1) - (vx(ic,jc,1) + vx(ic,jc,0)))
                    if (btypzs(ic,jc,1).eq.NEUMANN)    rhx(ic,jc,1  ) = rhx(ic,jc,1  ) - beta*am3ck(1  )*(bvalzs(ic,jc,1)*dzm(1) - (vx(ic,jc,1) - vx(ic,jc,0)))
                    if (btypzs(ic,jc,1).eq.SOMMERFELD) rhx(ic,jc,1  ) = rhx(ic,jc,1  ) + beta*am3ck(1  )*(coefzs(SOMMERFELD,1)*(vx(ic,jc,1) - vx(ic,jc,0)))
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis end
                    if (btypze(ic,jc,1).eq.DIRICHLET)  rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3ck(nzm)*(2.0*bvalze(ic,jc,1) - (vx(ic,jc,nzm) + vx(ic,jc,nz)))
                    if (btypze(ic,jc,1).eq.NEUMANN)    rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3ck(nzm)*(bvalze(ic,jc,1)*dzm(nz) + (vx(ic,jc,nzm) - vx(ic,jc,nz)))
                    if (btypze(ic,jc,1).eq.SOMMERFELD) rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3ck(nzm)*(coefze(SOMMERFELD,1)*(vx(ic,jc,nzm) - vx(ic,jc,nz)))
                end do
            end do
        end if
    end if
        
end subroutine AddBoundaryRHSTermsVX

subroutine ImposeInternalBoundaryVX

    use decomp_2d
    use param
    use local_arrays, only: vx
    use boundary_arrays

    implicit none

    integer :: jc,kc

    ! Impose boundary values on x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                !! Boundary condition at x-axis start
                if (btypxs(jc,kc,1).eq.DIRICHLET) then 
                    vx(1  ,jc,kc) = bvalxs(jc,kc,1)
                    glpcxs(jc,kc) = .false.
                else if (btypxs(jc,kc,1).eq.NEUMANN) then 
                    vx(1  ,jc,kc) = (vx(2,jc,kc) - bvalxs(jc,kc,1)*dxc(1))
                    glpcxs(jc,kc) = .true.
                else if (btypxs(jc,kc,1).eq.SOMMERFELD) then 
                    vx(1  ,jc,kc) = (vx(2,jc,kc)*bvalxs(jc,kc,1) + vx(1,jc,kc))/(1.0 + bvalxs(jc,kc,1))
                    glpcxs(jc,kc) = .true.
                end if
                !! Boundary condition at x-axis end
                if (btypxe(jc,kc,1).eq.DIRICHLET) then 
                    vx(nx ,jc,kc) = bvalxe(jc,kc,1)
                    glpcxe(jc,kc) = .false.
                else if (btypxe(jc,kc,1).eq.NEUMANN) then 
                    vx(nx ,jc,kc) = (vx(nxm,jc,kc) + bvalxe(jc,kc,1)*dxc(nxm))
                    glpcxe(jc,kc) = .true.
                else if (btypxe(jc,kc,1).eq.SOMMERFELD) then 
                    vx(nx ,jc,kc) = (vx(nxm,jc,kc)*bvalxe(jc,kc,1) + vx(nx,jc,kc))/(1.0 + bvalxe(jc,kc,1))
                    glpcxe(jc,kc) = .true.
                end if
            end do
        end do
    end if

end subroutine ImposeInternalBoundaryVX

subroutine ImposeExternalBoundaryVX

    use decomp_2d
    use param
    use local_arrays, only: vx
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc

    ! Impose boundary values on y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    !! Boundary condition at y-axis start
                    if (btypys(ic,kc,1).eq.DIRICHLET)  vx(ic,0 ,kc) = (2.0*bvalys(ic,kc,1) - vx(ic,1,kc))
                    if (btypys(ic,kc,1).eq.NEUMANN)    vx(ic,0 ,kc) = (vx(ic,1,kc) - bvalys(ic,kc,1)*dym(1))
                    if (btypys(ic,kc,1).eq.SOMMERFELD) vx(ic,0 ,kc) = (vx(ic,1,kc)*bvalys(ic,kc,1) + vx(ic,0,kc))/(1.0 + bvalys(ic,kc,1))
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3)-lvlhalo,xend(3)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    !! Boundary condition at y-axis end
                    if (btypye(ic,kc,1).eq.DIRICHLET)  vx(ic,ny,kc) = (2.0*bvalye(ic,kc,1) - vx(ic,nym,kc))
                    if (btypye(ic,kc,1).eq.NEUMANN)    vx(ic,ny,kc) = (vx(ic,nym,kc) + bvalye(ic,kc,1)*dym(ny))
                    if (btypye(ic,kc,1).eq.SOMMERFELD) vx(ic,ny,kc) = (vx(ic,nym,kc)*bvalye(ic,kc,1) + vx(ic,ny,kc))/(1.0 + bvalye(ic,kc,1))
                end do
            end do
        end if
    end if

    ! Impose boundary values on z-boundaries
    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    !! Boundary condition at z-axis start
                    if (btypzs(ic,jc,1).eq.DIRICHLET)  vx(ic,jc,0 ) = (2.0*bvalzs(ic,jc,1) - vx(ic,jc,1))
                    if (btypzs(ic,jc,1).eq.NEUMANN)    vx(ic,jc,0 ) = (vx(ic,jc,1) - bvalzs(ic,jc,1)*dzm(1))
                    if (btypzs(ic,jc,1).eq.SOMMERFELD) vx(ic,jc,0 ) = (vx(ic,jc,1)*bvalzs(ic,jc,1) + vx(ic,jc,0))/(1.0 + bvalzs(ic,jc,1))
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                do ic = xstart(1),xend(1)+lvlhalo
                    !! Boundary condition at z-axis end
                    if (btypze(ic,jc,1).eq.DIRICHLET)  vx(ic,jc,nz) = (2.0*bvalze(ic,jc,1) - vx(ic,jc,nzm))
                    if (btypze(ic,jc,1).eq.NEUMANN)    vx(ic,jc,nz) = (vx(ic,jc,nzm) + bvalze(ic,jc,1)*dzm(nz))
                    if (btypze(ic,jc,1).eq.SOMMERFELD) vx(ic,jc,nz) = (vx(ic,jc,nzm)*bvalze(ic,jc,1) + vx(ic,jc,nz))/(1.0 + bvalze(ic,jc,1))
                end do
            end do
        end if
    end if

end subroutine ImposeExternalBoundaryVX