!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: BoundaryRoutinesVZ.F90                         !
!    CONTAINS: subroutine CalcBoundaryConditionVZ         !
!              subroutine AddBoundaryRHSTermsVZ           !
!              subroutine ImposeInternalBoundaryVZ        !
!              subroutine ImposeExternalBoundaryVZ        !
!                                                         !
!    PURPOSE: Compute boundary conditions for vz before   !
!    implicit solve step, and impose boundary conditions  !
!    after the implicit solve and pressure correction.    !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CalcBoundaryConditionVZ  

    use decomp_2d
    use param
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc
    integer :: btyp
    real    :: bval

    ! Set the coefficients to be used for tridiagonal solvers
    !! For Dirichlet boundary condition
    coefxs(DIRICHLET,3)     = -1.0
    coefxe(DIRICHLET,3)     = -1.0
    coefys(DIRICHLET,3)     = -1.0
    coefye(DIRICHLET,3)     = -1.0
    coefzs(DIRICHLET,3)     =  0.0
    coefze(DIRICHLET,3)     =  0.0
    !! For Neumann boundary condition
    coefxs(NEUMANN,3)       =  1.0
    coefxe(NEUMANN,3)       =  1.0
    coefys(NEUMANN,3)       =  1.0
    coefye(NEUMANN,3)       =  1.0
    coefzs(NEUMANN,3)       =  1.0
    coefze(NEUMANN,3)       =  1.0
    !! For non-reflecting/advective/Sommerfeld boundary condition
    coefxs(SOMMERFELD,3)    =  (cvel*al*dt)/((cvel*al*dt) + dxm(1  ))
    coefxe(SOMMERFELD,3)    =  (cvel*al*dt)/((cvel*al*dt) + dxm(nx ))
    coefys(SOMMERFELD,3)    =  (cvel*al*dt)/((cvel*al*dt) + dym(1  ))
    coefye(SOMMERFELD,3)    =  (cvel*al*dt)/((cvel*al*dt) + dym(ny ))
    coefzs(SOMMERFELD,3)    =  (cvel*al*dt)/((cvel*al*dt) + dzc(1  ))
    coefze(SOMMERFELD,3)    =  (cvel*al*dt)/((cvel*al*dt) + dzc(nzm))

    ! Store coefficient arrays for implicit solver 
    !! For x-pencils
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                call VzBcXs(jc,kc,btyp,bval)
                cfbcxs(jc,kc,3) = coefxs(btyp,3)
                call VzBcXe(jc,kc,btyp,bval)
                cfbcxe(jc,kc,3) = coefxe(btyp,3)
            end do
        end do
    end if
    !! For y-pencils
    if (.not.periodic(2)) then
        do kc = ystart(3),yend(3)
            do ic = ystart(1),yend(1)
                call VzBcYs(ic,kc,btyp,bval)
                cfbcys(ic,kc,3) = coefys(btyp,3)
                call VzBcYe(ic,kc,btyp,bval)
                cfbcye(ic,kc,3) = coefye(btyp,3)
            end do
        end do
    end if
    !! For z-pencils
    if (.not.periodic(3)) then
        do jc = zstart(2),zend(2)
            do ic = zstart(1),zend(1)
                call VzBcZs(ic,jc,btyp,bval)
                cfbczs(ic,jc,3) = coefzs(btyp,3)
                call VzBcZe(ic,jc,btyp,bval)
                cfbcze(ic,jc,3) = coefze(btyp,3)
            end do
        end do
    end if

    ! Store boundary type and boundary value
    !! For x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)+lvlhalo
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                call VzBcXs(jc,kc,btypxs(jc,kc,3), bvalxs(jc,kc,3))
                !!! Store this value needed for imposing B.C.s later
                if (btypxs(jc,kc,3).eq.SOMMERFELD) bvalxs(jc,kc,3) = cvel*al*dt/dxm(1  )
                call VzBcXe(jc,kc,btypxe(jc,kc,3), bvalxe(jc,kc,3))
                !!! Store this value needed for imposing B.C.s later
                if (btypxe(jc,kc,3).eq.SOMMERFELD) bvalxe(jc,kc,3) = cvel*al*dt/dxm(nx )
            end do
        end do
    end if
    !! For y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VzBcYs(ic,kc,btypys(ic,kc,3), bvalys(ic,kc,3))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypys(ic,kc,3).eq.SOMMERFELD) bvalys(ic,kc,3) = cvel*al*dt/dym(1  )
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VzBcYe(ic,kc,btypye(ic,kc,3), bvalye(ic,kc,3))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypye(ic,kc,3).eq.SOMMERFELD) bvalye(ic,kc,3) = cvel*al*dt/dym(ny )
                end do
            end do
        end if
    end if
    !! For z-boundaries
    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VzBcZs(ic,jc,btypzs(ic,jc,3), bvalzs(ic,jc,3))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypzs(ic,jc,3).eq.SOMMERFELD) bvalzs(ic,jc,3) = cvel*al*dt/dzc(1  )
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    call VzBcZe(ic,jc,btypze(ic,jc,3), bvalze(ic,jc,3))
                    !!! Store this value needed for imposing B.C.s later
                    if (btypze(ic,jc,3).eq.SOMMERFELD) bvalze(ic,jc,3) = cvel*al*dt/dzc(nzm)
                end do
            end do
        end if
    end if

end subroutine CalcBoundaryConditionVZ

subroutine AddBoundaryRHSTermsVZ

    use decomp_2d
    use param
    use local_arrays, only: vz,rhx
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc
    real    :: beta

    ! Match the premultiplier with the one used in ImplicitAndUpdateVz
    beta = 0.5*al*dt/rey

    ! Add RHS terms to the x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                !! Boundary condition at x-axis start
                if (btypxs(jc,kc,3).eq.DIRICHLET)  rhx(1  ,jc,kc) = rhx(1  ,jc,kc) + beta*am1ci(1  )*(2.0*bvalxs(jc,kc,3) - (vz(1,jc,kc) + vz(0,jc,kc)))
                if (btypxs(jc,kc,3).eq.NEUMANN)    rhx(1  ,jc,kc) = rhx(1  ,jc,kc) - beta*am1ci(1  )*(bvalxs(jc,kc,3)*dxm(1) - (vz(1,jc,kc) - vz(0,jc,kc)))
                if (btypxs(jc,kc,3).eq.SOMMERFELD) rhx(1  ,jc,kc) = rhx(1  ,jc,kc) + beta*am1ci(1  )*(coefxs(SOMMERFELD,3)*(vz(1,jc,kc) - vz(0,jc,kc)))
                !! Boundary condition at x-axis end
                if (btypxe(jc,kc,3).eq.DIRICHLET)  rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1ci(nxm)*(2.0*bvalxe(jc,kc,3) - (vz(nxm,jc,kc) + vz(nx,jc,kc)))
                if (btypxe(jc,kc,3).eq.NEUMANN)    rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1ci(nxm)*(bvalxe(jc,kc,3)*dxm(nx) + (vz(nxm,jc,kc) - vz(nx,jc,kc)))
                if (btypxe(jc,kc,3).eq.SOMMERFELD) rhx(nxm,jc,kc) = rhx(nxm,jc,kc) + beta*ap1ci(nxm)*(coefxe(SOMMERFELD,3)*(vz(nxm,jc,kc) - vz(nx,jc,kc)))
            end do
        end do
    end if

    ! Add RHS terms to the y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis start
                    if (btypys(ic,kc,3).eq.DIRICHLET)  rhx(ic,1  ,kc) = rhx(ic,1  ,kc) + beta*am2cj(1  )*(2.0*bvalys(ic,kc,3) - (vz(ic,1,kc) + vz(ic,0,kc)))
                    if (btypys(ic,kc,3).eq.NEUMANN)    rhx(ic,1  ,kc) = rhx(ic,1  ,kc) - beta*am2cj(1  )*(bvalys(ic,kc,3)*dym(1) - (vz(ic,1,kc) - vz(ic,0,kc)))
                    if (btypys(ic,kc,3).eq.SOMMERFELD) rhx(ic,1  ,kc) = rhx(ic,1  ,kc) + beta*am2cj(1  )*(coefys(SOMMERFELD,3)*(vz(ic,1,kc) - vz(ic,0,kc)))
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at y-axis end
                    if (btypye(ic,kc,3).eq.DIRICHLET)  rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2cj(nym)*(2.0*bvalye(ic,kc,3) - (vz(ic,nym,kc) + vz(ic,ny,kc)))
                    if (btypye(ic,kc,3).eq.NEUMANN)    rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2cj(nym)*(bvalye(ic,kc,3)*dym(ny) + (vz(ic,nym,kc) - vz(ic,ny,kc)))
                    if (btypye(ic,kc,3).eq.SOMMERFELD) rhx(ic,nym,kc) = rhx(ic,nym,kc) + beta*ap2cj(nym)*(coefye(SOMMERFELD,3)*(vz(ic,nym,kc) - vz(ic,ny,kc)))
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
                    if (btypzs(ic,jc,3).eq.DIRICHLET)  rhx(ic,jc,2  ) = rhx(ic,jc,2  ) + beta*am3sk(2  )*(bvalzs(ic,jc,3) - vz(ic,jc,1))
                    if (btypzs(ic,jc,3).eq.NEUMANN)    rhx(ic,jc,2  ) = rhx(ic,jc,2  ) - beta*am3sk(2  )*(bvalzs(ic,jc,3)*dzc(1) - (vz(ic,jc,2) - vz(ic,jc,1)))
                    if (btypzs(ic,jc,3).eq.SOMMERFELD) rhx(ic,jc,2  ) = rhx(ic,jc,2  ) + beta*am3sk(2  )*(coefzs(SOMMERFELD,3)*(vz(ic,jc,2) - vz(ic,jc,1)))
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis end
                    if (btypze(ic,jc,3).eq.DIRICHLET)  rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3sk(nzm)*(bvalze(ic,jc,3) - vz(ic,jc,nz))
                    if (btypze(ic,jc,3).eq.NEUMANN)    rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3sk(nzm)*(bvalze(ic,jc,3)*dzc(nzm) + (vz(ic,jc,nzm) - vz(ic,jc,nz)))
                    if (btypze(ic,jc,3).eq.SOMMERFELD) rhx(ic,jc,nzm) = rhx(ic,jc,nzm) + beta*ap3sk(nzm)*(coefze(SOMMERFELD,3)*(vz(ic,jc,nzm) - vz(ic,jc,nz)))
                end do
            end do
        end if
    end if
        
end subroutine AddBoundaryRHSTermsVZ

subroutine ImposeInternalBoundaryVZ

    use decomp_2d
    use param
    use local_arrays, only: vz
    use boundary_arrays

    implicit none

    integer :: ic,jc

    ! Impose boundary values on z-boundaries
    if (.not.periodic(3)) then
        if (xstart(3).eq.1) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis start
                    if (btypzs(ic,jc,3).eq.DIRICHLET)  then 
                        vz(ic,jc,1  ) = bvalzs(ic,jc,3)
                        glpczs(ic,jc) = .false.
                    else if (btypzs(ic,jc,3).eq.NEUMANN) then 
                        vz(ic,jc,1  ) = (vz(ic,jc,2) - bvalzs(ic,jc,3)*dzc(1))
                        glpczs(ic,jc) = .true.
                    else if (btypzs(ic,jc,3).eq.SOMMERFELD) then 
                        vz(ic,jc,1  ) = (vz(ic,jc,2)*bvalzs(ic,jc,3) + vz(ic,jc,1))/(1.0 + bvalzs(ic,jc,3))
                        glpczs(ic,jc) = .true.
                    end if
                end do
            end do
        end if
        if (xend(3).eq.nzm) then
            do jc = xstart(2),xend(2)
                do ic = xstart(1),xend(1)
                    !! Boundary condition at z-axis end
                    if (btypze(ic,jc,3).eq.DIRICHLET) then 
                        vz(ic,jc,nz ) = bvalze(ic,jc,3) 
                        glpcze(ic,jc) = .false.
                    else if (btypze(ic,jc,3).eq.NEUMANN) then 
                        vz(ic,jc,nz ) = (vz(ic,jc,nzm) + bvalze(ic,jc,3)*dzc(nzm))
                        glpcze(ic,jc) = .true.
                    else if (btypze(ic,jc,3).eq.SOMMERFELD) then 
                        vz(ic,jc,nz ) = (vz(ic,jc,nzm)*bvalze(ic,jc,3) + vz(ic,jc,nz))/(1.0 + bvalze(ic,jc,3))
                        glpcze(ic,jc) = .true.
                    end if
                end do
            end do
        end if
    end if

end subroutine ImposeInternalBoundaryVZ

subroutine ImposeExternalBoundaryVZ

    use decomp_2d
    use param
    use local_arrays, only: vz
    use boundary_arrays

    implicit none

    integer :: ic,jc,kc

    ! Impose boundary values on x-boundaries
    if (.not.periodic(1)) then
        do kc = xstart(3),xend(3)+lvlhalo
            do jc = xstart(2)-lvlhalo,xend(2)+lvlhalo
                !! Boundary condition at x-axis start
                if (btypxs(jc,kc,3).eq.DIRICHLET)  vz(0 ,jc,kc) = (2.0*bvalxs(jc,kc,3) - vz(1,jc,kc))
                if (btypxs(jc,kc,3).eq.NEUMANN)    vz(0 ,jc,kc) = (vz(1,jc,kc) - bvalxs(jc,kc,3)*dxm(1))
                if (btypxs(jc,kc,3).eq.SOMMERFELD) vz(0 ,jc,kc) = (vz(1,jc,kc)*bvalxs(jc,kc,3) + vz(0,jc,kc))/(1.0 + bvalxs(jc,kc,3))
                !! Boundary condition at x-axis end
                if (btypxe(jc,kc,3).eq.DIRICHLET)  vz(nx,jc,kc) = (2.0*bvalxe(jc,kc,3) - vz(nxm,jc,kc))
                if (btypxe(jc,kc,3).eq.NEUMANN)    vz(nx,jc,kc) = (vz(nxm,jc,kc) + bvalxe(jc,kc,3)*dxm(nx))
                if (btypxe(jc,kc,3).eq.SOMMERFELD) vz(nx,jc,kc) = (vz(nxm,jc,kc)*bvalxe(jc,kc,3) + vz(nx,jc,kc))/(1.0 + bvalxe(jc,kc,3))
            end do
        end do
    end if

    ! Impose boundary values on y-boundaries
    if (.not.periodic(2)) then
        if (xstart(2).eq.1) then
            do kc = xstart(3),xend(3)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    !! Boundary condition at y-axis start
                    if (btypys(ic,kc,3).eq.DIRICHLET)  vz(ic,0 ,kc) = (2.0*bvalys(ic,kc,3) - vz(ic,1,kc))
                    if (btypys(ic,kc,3).eq.NEUMANN)    vz(ic,0 ,kc) = (vz(ic,1,kc) - bvalys(ic,kc,3)*dym(1))
                    if (btypys(ic,kc,3).eq.SOMMERFELD) vz(ic,0 ,kc) = (vz(ic,1,kc)*bvalys(ic,kc,3) + vz(ic,0,kc))/(1.0 + bvalys(ic,kc,3))
                end do
            end do
        end if
        if (xend(2).eq.nym) then
            do kc = xstart(3),xend(3)+lvlhalo
                do ic = xstart(1)-lvlhalo,xend(1)+lvlhalo
                    !! Boundary condition at y-axis end
                    if (btypye(ic,kc,3).eq.DIRICHLET)  vz(ic,ny,kc) = (2.0*bvalye(ic,kc,3) - vz(ic,nym,kc))
                    if (btypye(ic,kc,3).eq.NEUMANN)    vz(ic,ny,kc) = (vz(ic,nym,kc) + bvalye(ic,kc,3)*dym(ny))
                    if (btypye(ic,kc,3).eq.SOMMERFELD) vz(ic,ny,kc) = (vz(ic,nym,kc)*bvalye(ic,kc,3) + vz(ic,ny,kc))/(1.0 + bvalye(ic,kc,3))
                end do
            end do
        end if
    end if

end subroutine ImposeExternalBoundaryVZ