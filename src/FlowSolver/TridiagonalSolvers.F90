!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
!    FILE:     TridiagonalSolvers.F90                          !
!                                                              !
!    CONTAINS: subroutine TridiagonalSolverXM                  !
!              subroutine TridiagonalSolverXC                  !
!              subroutine TridiagonalSolverYM                  !
!              subroutine TridiagonalSolverYC                  !
!              subroutine TridiagonalSolverZM                  !
!              subroutine TridiagonalSolverZC                  !
!                                                              !
!              subroutine TridiagonalSolverXR                  !
!              subroutine TridiagonalSolverYR                  !
!              subroutine TridiagonalSolverZR                  !
!                                                              !
!              subroutine ThomasBatched                        !
!              subroutine ThomasBatchedPeriodic                !
!                                                              !
!    PURPOSE: subroutines that are used for                    !
!             solving tridiagonal systems for diffusive        !
!             terms of the Navier-Stokes or Pressure-Poisson   !
!             equations implicitly                             !
!                                                              !
!    All solvers call the batched Thomas algorithm, which      !
!    solves many independent lines at once. Arrays are passed  !
!    to it by their first element, e.g. rhy(ystart(1),1,kc),   !
!    together with the memory distance between two rows of a   !
!    line (ld). In X a single line is solved per call, in Y    !
!    all lines of a kc plane, and in Z all lines of a jc row.  !
!                                                              !
!    NAMING                                                    !
!    ------                                                    !
!    Suffix M : variable on cell centres in the solve          !
!               direction (e.g. vy, vz, tp in X), rows 1:n     !
!    Suffix C : variable on cell faces in the solve direction  !
!               (e.g. vx in X). Without periodicity, face 1    !
!               lies on the boundary and is not solved for,    !
!               rows 2:n are solved                            !
!    Suffix R : pressure correction (Helmholtz) equation       !
!                                                              !
!    am, ac, ap  : sub-, main and super-diagonal coefficients  !
!    sh          : shift of the main diagonal per line         !
!    dcs, dce    : boundary correction of the first/last main  !
!                  diagonal entry per line                     !
!    suffix il/jl/kl : line along X/Y/Z (as in the old code)   !
!                                                              !
!    WHY BATCHED                                               !
!    -----------                                               !
!    Along a line, the Thomas algorithm is a recurrence (row k !
!    needs row k-1) and cannot be vectorised. Different lines  !
!    are independent, so the kernels loop over the lines in    !
!    the innermost loop. These are contiguous in memory for    !
!    Y and Z pencils, so the compiler vectorises them.         !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------- FOR MOMENTUM EQUATIONS -------------------!
!                                                              !
! Each routine solves, along one direction, the factor of the  !
! approximately factorised implicit diffusion operator:        !
!                                                              !
!   (1 - betadx*D2) dq = rhs,   betadx = 0.5*al*dt*difc        !
!                                                              !
! where D2 is the second derivative with stencil coefficients  !
! am/ac/ap (e.g. am1ci, ac1ci, ap1ci from param) and difc is   !
! the diffusion coefficient (1/Re for velocity, 1/(Re*Pr) for  !
! temperature). The rows are used unscaled:                    !
!                                                              !
!   -betadx*am(k)*dq(k-1) + (1 - betadx*ac(k))*dq(k)           !
!   -betadx*ap(k)*dq(k+1) = rhs(k)                             !
!                                                              !
! Boundary conditions: the ghost value outside the domain is   !
! written as dq(0) = cf*dq(1), with cf = cfbc*s/e computed in  !
! CalcBoundaryCondition* (bcid selects vx/vy/vz/tp). The ghost !
! term am(1)*dq(0) then becomes cf*am(1)*dq(1), which is added !
! to the first main diagonal entry (dcs). The same holds at    !
! the end of the line (dce). For periodic directions the       !
! periodic kernel is used and no boundary terms are needed.    !
!                                                              !
! The momentum equations have no diagonal shift (sh = 0).      !
!                                                              !
!--------------------------------------------------------------!

subroutine TridiagonalSolverXM(bcid,difc)

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: rhx
    use boundary_arrays, only: cfbcxs,cfbcxe
    use tridiagonal_buffer

    implicit none

    integer, intent(in) :: bcid     ! Boundary condition index
    real, intent(in)    :: difc     ! Diffusion coefficient (1/Re or 1/(Re*Pr))

    real                :: amil(nxm),acil(nxm),apil(nxm)
    real                :: shil(1),dcis(1),dcie(1)
    real                :: betadx
    integer             :: ic,jc,kc

    betadx = 0.5*al*dt*difc

    ! Coefficients of (1 - betadx*D2) along X for a cell-centred variable
    do ic = 1,nxm
        amil(ic) = -am1ci(ic)*betadx
        acil(ic) =  1.0 - ac1ci(ic)*betadx
        apil(ic) = -ap1ci(ic)*betadx
    end do

    ! No diagonal shift for the momentum equations
    shil(1) = 0.0

    if (periodic(1)) then

        ! One line per call: rhx(1:nxm,jc,kc) is contiguous (nl = 1, ld = 1)
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                call ThomasBatchedPeriodic(1,nxm,1,amil,acil,apil,shil,rhx(1,jc,kc),twk1(1),twk2(1))
            end do
        end do

    else

        ! One line per call, with the boundary conditions of line (jc,kc)
        ! folded into its first and last main diagonal entry
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                dcis(1) = cfbcxs(jc,kc,bcid)*amil(1)
                dcie(1) = cfbcxe(jc,kc,bcid)*apil(nxm)
                call ThomasBatched(1,nxm,1,amil,acil,apil,shil,dcis,dcie,rhx(1,jc,kc),twk1(1))
            end do
        end do

    end if

    return

end subroutine TridiagonalSolverXM

subroutine TridiagonalSolverXC(bcid,difc)

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: rhx
    use boundary_arrays, only: cfbcxs,cfbcxe
    use tridiagonal_buffer

    implicit none

    integer, intent(in) :: bcid     ! Boundary condition index
    real, intent(in)    :: difc     ! Diffusion coefficient (1/Re or 1/(Re*Pr))

    real                :: amil(nxm),acil(nxm),apil(nxm)
    real                :: shil(1),dcis(1),dcie(1)
    real                :: betadx
    integer             :: ic,jc,kc

    betadx = 0.5*al*dt*difc

    ! Coefficients of (1 - betadx*D2) along X for a face-centred variable
    do ic = 1,nxm
        amil(ic) = -am1si(ic)*betadx
        acil(ic) =  1.0 - ac1si(ic)*betadx
        apil(ic) = -ap1si(ic)*betadx
    end do

    ! No diagonal shift for the momentum equations
    shil(1) = 0.0

    if (periodic(1)) then

        ! One line per call: rhx(1:nxm,jc,kc) is contiguous (nl = 1, ld = 1)
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                call ThomasBatchedPeriodic(1,nxm,1,amil,acil,apil,shil,rhx(1,jc,kc),twk1(1),twk2(1))
            end do
        end do

    else

        ! Face ic = 1 lies on the boundary and is not solved for: its value is
        ! imposed by the boundary routines, so its increment is zero. Rows 2:nxm
        ! are solved by passing all arrays from their 2nd element on, and the
        ! boundary condition is folded into row 2 through amil(2).
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)
                dcis(1) = cfbcxs(jc,kc,bcid)*amil(2)
                dcie(1) = cfbcxe(jc,kc,bcid)*apil(nxm)
                call ThomasBatched(1,nxm-1,1,amil(2),acil(2),apil(2),shil,dcis,dcie,rhx(2,jc,kc),twk1(1))
                rhx(1,jc,kc) = 0.0
            end do
        end do

    end if

    return

end subroutine TridiagonalSolverXC

subroutine TridiagonalSolverYM(bcid,difc)

    use param
    use decomp_2d, only: ystart,yend,ysize
    use local_arrays, only: rhy
    use boundary_arrays, only: cfbcys,cfbcye
    use tridiagonal_buffer

    implicit none

    integer, intent(in) :: bcid     ! Boundary condition index
    real, intent(in)    :: difc     ! Diffusion coefficient (1/Re or 1/(Re*Pr))

    real                :: amjl(nym),acjl(nym),apjl(nym)
    real                :: shjl(ysize(1)),dcjs(ysize(1)),dcje(ysize(1))
    real                :: betadx
    integer             :: jc,kc

    betadx = 0.5*al*dt*difc

    ! Coefficients of (1 - betadx*D2) along Y for a cell-centred variable
    do jc = 1,nym
        amjl(jc) = -am2cj(jc)*betadx
        acjl(jc) =  1.0 - ac2cj(jc)*betadx
        apjl(jc) = -ap2cj(jc)*betadx
    end do

    ! No diagonal shift for the momentum equations
    shjl(:) = 0.0

    if (periodic(2)) then

        ! All ysize(1) lines of a kc plane per call: row jc of line ic is
        ! rhy(ic,jc,kc), so lines are contiguous and rows are ysize(1) apart
        do kc = ystart(3),yend(3)
            call ThomasBatchedPeriodic(ysize(1),nym,ysize(1),amjl,acjl,apjl,shjl,rhy(ystart(1),1,kc),twk1(1),twk2(1))
        end do

    else

        ! All lines of a kc plane per call, with the boundary conditions of
        ! each line ic folded into its first and last main diagonal entry
        do kc = ystart(3),yend(3)
            dcjs(:) = cfbcys(ystart(1):yend(1),kc,bcid)*amjl(1)
            dcje(:) = cfbcye(ystart(1):yend(1),kc,bcid)*apjl(nym)
            call ThomasBatched(ysize(1),nym,ysize(1),amjl,acjl,apjl,shjl,dcjs,dcje,rhy(ystart(1),1,kc),twk1(1))
        end do

    end if

    return

end subroutine TridiagonalSolverYM

subroutine TridiagonalSolverYC(bcid,difc)

    use param
    use decomp_2d, only: ystart,yend,ysize
    use local_arrays, only: rhy
    use boundary_arrays, only: cfbcys,cfbcye
    use tridiagonal_buffer

    implicit none

    integer, intent(in) :: bcid     ! Boundary condition index
    real, intent(in)    :: difc     ! Diffusion coefficient (1/Re or 1/(Re*Pr))

    real                :: amjl(nym),acjl(nym),apjl(nym)
    real                :: shjl(ysize(1)),dcjs(ysize(1)),dcje(ysize(1))
    real                :: betadx
    integer             :: jc,kc

    betadx = 0.5*al*dt*difc

    ! Coefficients of (1 - betadx*D2) along Y for a face-centred variable
    do jc = 1,nym
        amjl(jc) = -am2sj(jc)*betadx
        acjl(jc) =  1.0 - ac2sj(jc)*betadx
        apjl(jc) = -ap2sj(jc)*betadx
    end do

    ! No diagonal shift for the momentum equations
    shjl(:) = 0.0

    if (periodic(2)) then

        ! All ysize(1) lines of a kc plane per call: row jc of line ic is
        ! rhy(ic,jc,kc), so lines are contiguous and rows are ysize(1) apart
        do kc = ystart(3),yend(3)
            call ThomasBatchedPeriodic(ysize(1),nym,ysize(1),amjl,acjl,apjl,shjl,rhy(ystart(1),1,kc),twk1(1),twk2(1))
        end do

    else

        ! Face jc = 1 lies on the boundary and is not solved for: its value is
        ! imposed by the boundary routines, so its increment is zero. Rows 2:nym
        ! are solved by passing all arrays from their 2nd element on, and the
        ! boundary condition is folded into row 2 through amjl(2).
        do kc = ystart(3),yend(3)
            dcjs(:) = cfbcys(ystart(1):yend(1),kc,bcid)*amjl(2)
            dcje(:) = cfbcye(ystart(1):yend(1),kc,bcid)*apjl(nym)
            call ThomasBatched(ysize(1),nym-1,ysize(1),amjl(2),acjl(2),apjl(2),shjl,dcjs,dcje,rhy(ystart(1),2,kc),twk1(1))
            rhy(ystart(1):yend(1),1,kc) = 0.0
        end do

    end if

    return

end subroutine TridiagonalSolverYC

subroutine TridiagonalSolverZM(bcid,difc)

    use param
    use decomp_2d, only: zstart,zend,zsize
    use local_arrays, only: rhz
    use boundary_arrays, only: cfbczs,cfbcze
    use tridiagonal_buffer

    implicit none

    integer, intent(in) :: bcid     ! Boundary condition index
    real, intent(in)    :: difc     ! Diffusion coefficient (1/Re or 1/(Re*Pr))

    real                :: amkl(nzm),ackl(nzm),apkl(nzm)
    real                :: shkl(zsize(1)),dcks(zsize(1)),dcke(zsize(1))
    real                :: betadx
    integer             :: jc,kc,ldz

    betadx = 0.5*al*dt*difc

    ! Coefficients of (1 - betadx*D2) along Z for a cell-centred variable
    do kc = 1,nzm
        amkl(kc) = -am3ck(kc)*betadx
        ackl(kc) =  1.0 - ac3ck(kc)*betadx
        apkl(kc) = -ap3ck(kc)*betadx
    end do

    ! No diagonal shift for the momentum equations
    shkl(:) = 0.0

    ! All zsize(1) lines of a jc row per call: row kc of line ic is
    ! rhz(ic,jc,kc), so lines are contiguous and rows are ldz apart
    ldz     = zsize(1)*zsize(2)

    if (periodic(3)) then

        do jc = zstart(2),zend(2)
            call ThomasBatchedPeriodic(zsize(1),nzm,ldz,amkl,ackl,apkl,shkl,rhz(zstart(1),jc,1),twk1(1),twk2(1))
        end do

    else

        ! All lines of a jc row per call, with the boundary conditions of
        ! each line ic folded into its first and last main diagonal entry
        do jc = zstart(2),zend(2)
            dcks(:) = cfbczs(zstart(1):zend(1),jc,bcid)*amkl(1)
            dcke(:) = cfbcze(zstart(1):zend(1),jc,bcid)*apkl(nzm)
            call ThomasBatched(zsize(1),nzm,ldz,amkl,ackl,apkl,shkl,dcks,dcke,rhz(zstart(1),jc,1),twk1(1))
        end do

    end if

    return

end subroutine TridiagonalSolverZM

subroutine TridiagonalSolverZC(bcid,difc)

    use param
    use decomp_2d, only: zstart,zend,zsize
    use local_arrays, only: rhz
    use boundary_arrays, only: cfbczs,cfbcze
    use tridiagonal_buffer

    implicit none

    integer, intent(in) :: bcid     ! Boundary condition index
    real, intent(in)    :: difc     ! Diffusion coefficient (1/Re or 1/(Re*Pr))

    real                :: amkl(nzm),ackl(nzm),apkl(nzm)
    real                :: shkl(zsize(1)),dcks(zsize(1)),dcke(zsize(1))
    real                :: betadx
    integer             :: jc,kc,ldz

    betadx = 0.5*al*dt*difc

    ! Coefficients of (1 - betadx*D2) along Z for a face-centred variable
    do kc = 1,nzm
        amkl(kc) = -am3sk(kc)*betadx
        ackl(kc) =  1.0 - ac3sk(kc)*betadx
        apkl(kc) = -ap3sk(kc)*betadx
    end do

    ! No diagonal shift for the momentum equations
    shkl(:) = 0.0

    ! All zsize(1) lines of a jc row per call: row kc of line ic is
    ! rhz(ic,jc,kc), so lines are contiguous and rows are ldz apart
    ldz     = zsize(1)*zsize(2)

    if (periodic(3)) then

        do jc = zstart(2),zend(2)
            call ThomasBatchedPeriodic(zsize(1),nzm,ldz,amkl,ackl,apkl,shkl,rhz(zstart(1),jc,1),twk1(1),twk2(1))
        end do

    else

        ! Face kc = 1 lies on the boundary and is not solved for: its value is
        ! imposed by the boundary routines, so its increment is zero. Rows 2:nzm
        ! are solved by passing all arrays from their 2nd element on, and the
        ! boundary condition is folded into row 2 through amkl(2).
        do jc = zstart(2),zend(2)
            dcks(:) = cfbczs(zstart(1):zend(1),jc,bcid)*amkl(2)
            dcke(:) = cfbcze(zstart(1):zend(1),jc,bcid)*apkl(nzm)
            call ThomasBatched(zsize(1),nzm-1,ldz,amkl(2),ackl(2),apkl(2),shkl,dcks,dcke,rhz(zstart(1),jc,2),twk1(1))
            rhz(zstart(1):zend(1),jc,1) = 0.0
        end do

    end if

    return

end subroutine TridiagonalSolverZC

!--------------- FOR PRESSURE POISSON EQUATION  ---------------!
!                                                              !
! Before these routines are called, the right-hand side has    !
! been transformed to the eigenvector basis of the discrete    !
! second derivative in the two other directions (see           !
! MultiplyEigenvectorsBW_* in SolvePressureCorrection.F90).    !
! Each line then satisfies a 1D Helmholtz equation along the   !
! remaining direction:                                         !
!                                                              !
!   (D2 + eig_a + eig_b) p = rhs                               !
!                                                              !
! where D2 has the stencil coefficients amph*/acph*/apph* and  !
! eig_a + eig_b is the sum of the eigenvalues of the two other !
! directions belonging to this line (diagonal shift sh).       !
!                                                              !
! Boundary conditions (cell-centred, ghost cell p(0)):         !
!   NEUMANN   : p(0) =  p(1)  ->  dcs = +amph(1)               !
!   DIRICHLET : p(0) = -p(1)  ->  dcs = -amph(1)               !
! and likewise at the end of the line with apph(n) (dce).      !
!                                                              !
! If all directions are periodic or Neumann-Neumann, the line  !
! carrying the zero eigenvalue of the other two directions is  !
! singular (pressure is defined up to a constant). That line   !
! is solved with one value pinned to zero (izmx/izmy/izmz).    !
! Because the right-hand side is compatible (the global        !
! divergence is corrected beforehand), the dropped equation is !
! satisfied automatically and the solution is exact:           !
!   non-periodic : last value pinned, rows 1:n-1 solved        !
!   periodic     : first value pinned, rows 2:n solved (both   !
!                  periodic couplings multiply the first value !
!                  and vanish, leaving an ordinary system)     !
!                                                              !
!--------------------------------------------------------------!

subroutine TridiagonalSolverXR

    use param
    use decomp_2d, only: xstart,xend
    use local_arrays, only: rhx
    use boundary_arrays, only: DIRICHLET,NEUMANN
    use pressure_solver, only: eigy,eigz,pcbcxs,pcbcxe,izmx,izmy,izmz
    use tridiagonal_buffer

    implicit none

    real                :: shil(1),dcis(1),dcie(1),dzro(1)
    real                :: prst,pren
    integer             :: jc,kc
    logical             :: zmod

    dzro(1) = 0.0

    if (periodic(1)) then

        ! One line per call: rhx(1:nxm,jc,kc) is contiguous (nl = 1, ld = 1)
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)

                ! Shift: eigenvalues of the transformed Y and Z directions
                shil(1) = eigy(jc) + eigz(kc)
                ! Singular line: zero mode in Y and Z, and X itself has a zero mode
                zmod    = (izmx.gt.0).and.(jc.eq.izmy).and.(kc.eq.izmz)

                if (zmod) then
                    ! Singular line: pin ic = 1 to zero and solve rows 2:nxm
                    ! (the periodic couplings multiply rhx(1) and drop out)
                    call ThomasBatched(1,nxm-1,1,amphi(2),acphi(2),apphi(2),shil,dzro,dzro,rhx(2,jc,kc),twk1(1))
                    rhx(1,jc,kc) = 0.0
                else
                    call ThomasBatchedPeriodic(1,nxm,1,amphi,acphi,apphi,shil,rhx(1,jc,kc),twk1(1),twk2(1))
                end if

            end do
        end do

    else

        ! Ghost cell factor: +1 for Neumann (p(0) = p(1)), -1 for Dirichlet (p(0) = -p(1))
        if (pcbcxs.eq.DIRICHLET) prst = -1.0
        if (pcbcxs.eq.NEUMANN)   prst =  1.0
        if (pcbcxe.eq.DIRICHLET) pren = -1.0
        if (pcbcxe.eq.NEUMANN)   pren =  1.0

        ! Boundary corrections are the same for every line
        dcis(1) = prst*amphi(1)
        dcie(1) = pren*apphi(nxm)

        ! One line per call: rhx(1:nxm,jc,kc) is contiguous (nl = 1, ld = 1)
        do kc = xstart(3),xend(3)
            do jc = xstart(2),xend(2)

                ! Shift: eigenvalues of the transformed Y and Z directions
                shil(1) = eigy(jc) + eigz(kc)
                ! Singular line: zero mode in Y and Z, and X itself has a zero mode
                zmod    = (izmx.gt.0).and.(jc.eq.izmy).and.(kc.eq.izmz)

                if (zmod) then
                    ! Singular line: pin ic = nxm to zero and solve rows 1:nxm-1
                    ! (no end boundary term, the coupling to rhx(nxm) drops out)
                    call ThomasBatched(1,nxm-1,1,amphi,acphi,apphi,shil,dcis,dzro,rhx(1,jc,kc),twk1(1))
                    rhx(nxm,jc,kc) = 0.0
                else
                    call ThomasBatched(1,nxm,1,amphi,acphi,apphi,shil,dcis,dcie,rhx(1,jc,kc),twk1(1))
                end if

            end do
        end do

    end if

    return

end subroutine TridiagonalSolverXR

subroutine TridiagonalSolverYR

    use param
    use decomp_2d, only: ystart,yend,ysize
    use local_arrays, only: rhy
    use boundary_arrays, only: DIRICHLET,NEUMANN
    use pressure_solver, only: eigx,eigz,pcbcys,pcbcye,izmx,izmy,izmz
    use tridiagonal_buffer

    implicit none

    real                :: shjl(ysize(1)),dcjs(ysize(1)),dcje(ysize(1)),dzro(1)
    real                :: rsav(nym)
    real                :: prst,pren
    integer             :: ic,kc,lz
    logical             :: zmod

    dzro(1) = 0.0
    lz      = izmx - ystart(1) + 1      ! Local line index of the zero mode in X

    if (.not.periodic(2)) then
        ! Ghost cell factor: +1 for Neumann (p(0) = p(1)), -1 for Dirichlet (p(0) = -p(1))
        if (pcbcys.eq.DIRICHLET) prst = -1.0
        if (pcbcys.eq.NEUMANN)   prst =  1.0
        if (pcbcye.eq.DIRICHLET) pren = -1.0
        if (pcbcye.eq.NEUMANN)   pren =  1.0
        ! Boundary corrections are the same for every line
        dcjs(:) = prst*amphj(1)
        dcje(:) = pren*apphj(nym)
    end if

    do kc = ystart(3),yend(3)

        ! Shift of line ic: eigenvalues of the transformed X and Z directions
        do ic = ystart(1),yend(1)
            shjl(ic-ystart(1)+1) = eigx(ic) + eigz(kc)
        end do

        ! The singular line ic = izmx is in this plane if kc = izmz and izmx is on this rank.
        ! Save its right-hand side: the batched call below also processes it, but gives a
        ! meaningless result for that (singular) line, which is overwritten afterwards.
        zmod = (izmy.gt.0).and.(kc.eq.izmz).and.(izmx.ge.ystart(1)).and.(izmx.le.yend(1))
        if (zmod) rsav(1:nym) = rhy(izmx,1:nym,kc)

        ! All ysize(1) lines of this kc plane per call (rows are ysize(1) apart)
        if (periodic(2)) then
            call ThomasBatchedPeriodic(ysize(1),nym,ysize(1),amphj,acphj,apphj,shjl,rhy(ystart(1),1,kc),twk1(1),twk2(1))
        else
            call ThomasBatched(ysize(1),nym,ysize(1),amphj,acphj,apphj,shjl,dcjs,dcje,rhy(ystart(1),1,kc),twk1(1))
        end if

        if (zmod) then
            ! Singular line: restore its right-hand side and solve it with one value pinned to zero
            ! (a single line, nl = 1, still with rows ysize(1) apart in memory)
            rhy(izmx,1:nym,kc) = rsav(1:nym)
            if (periodic(2)) then
                ! Pin jc = 1, solve rows 2:nym without boundary terms
                call ThomasBatched(1,nym-1,ysize(1),amphj(2),acphj(2),apphj(2),shjl(lz),dzro,dzro,rhy(izmx,2,kc),twk1(1))
                rhy(izmx,1,kc) = 0.0
            else
                ! Pin jc = nym, solve rows 1:nym-1 with the start boundary term only
                call ThomasBatched(1,nym-1,ysize(1),amphj,acphj,apphj,shjl(lz),dcjs(lz),dzro,rhy(izmx,1,kc),twk1(1))
                rhy(izmx,nym,kc) = 0.0
            end if
        end if

    end do

    return

end subroutine TridiagonalSolverYR

subroutine TridiagonalSolverZR

    use param
    use decomp_2d, only: zstart,zend,zsize
    use local_arrays, only: rhz
    use boundary_arrays, only: DIRICHLET,NEUMANN
    use pressure_solver, only: eigx,eigy,pcbczs,pcbcze,izmx,izmy,izmz
    use tridiagonal_buffer

    implicit none

    real                :: shkl(zsize(1)),dcks(zsize(1)),dcke(zsize(1)),dzro(1)
    real                :: rsav(nzm)
    real                :: prst,pren
    integer             :: ic,jc,lz,ldz
    logical             :: zmod

    dzro(1) = 0.0
    lz      = izmx - zstart(1) + 1      ! Local line index of the zero mode in X
    ldz     = zsize(1)*zsize(2)

    if (.not.periodic(3)) then
        ! Ghost cell factor: +1 for Neumann (p(0) = p(1)), -1 for Dirichlet (p(0) = -p(1))
        if (pcbczs.eq.DIRICHLET) prst = -1.0
        if (pcbczs.eq.NEUMANN)   prst =  1.0
        if (pcbcze.eq.DIRICHLET) pren = -1.0
        if (pcbcze.eq.NEUMANN)   pren =  1.0
        ! Boundary corrections are the same for every line
        dcks(:) = prst*amphk(1)
        dcke(:) = pren*apphk(nzm)
    end if

    do jc = zstart(2),zend(2)

        ! Shift of line ic: eigenvalues of the transformed X and Y directions
        do ic = zstart(1),zend(1)
            shkl(ic-zstart(1)+1) = eigx(ic) + eigy(jc)
        end do

        ! The singular line ic = izmx is in this row if jc = izmy and izmx is on this rank.
        ! Save its right-hand side: the batched call below also processes it, but gives a
        ! meaningless result for that (singular) line, which is overwritten afterwards.
        zmod = (izmz.gt.0).and.(jc.eq.izmy).and.(izmx.ge.zstart(1)).and.(izmx.le.zend(1))
        if (zmod) rsav(1:nzm) = rhz(izmx,jc,1:nzm)

        ! All zsize(1) lines of this jc row per call (rows are ldz apart)
        if (periodic(3)) then
            call ThomasBatchedPeriodic(zsize(1),nzm,ldz,amphk,acphk,apphk,shkl,rhz(zstart(1),jc,1),twk1(1),twk2(1))
        else
            call ThomasBatched(zsize(1),nzm,ldz,amphk,acphk,apphk,shkl,dcks,dcke,rhz(zstart(1),jc,1),twk1(1))
        end if

        if (zmod) then
            ! Singular line: restore its right-hand side and solve it with one value pinned to zero
            ! (a single line, nl = 1, still with rows ldz apart in memory)
            rhz(izmx,jc,1:nzm) = rsav(1:nzm)
            if (periodic(3)) then
                ! Pin kc = 1, solve rows 2:nzm without boundary terms
                call ThomasBatched(1,nzm-1,ldz,amphk(2),acphk(2),apphk(2),shkl(lz),dzro,dzro,rhz(izmx,jc,2),twk1(1))
                rhz(izmx,jc,1) = 0.0
            else
                ! Pin kc = nzm, solve rows 1:nzm-1 with the start boundary term only
                call ThomasBatched(1,nzm-1,ldz,amphk,acphk,apphk,shkl(lz),dcks(lz),dzro,rhz(izmx,jc,1),twk1(1))
                rhz(izmx,jc,nzm) = 0.0
            end if
        end if

    end do

    return

end subroutine TridiagonalSolverZR

!------------------ BATCHED THOMAS ALGORITHM ------------------!
!                                                              !
! Solves nl independent tridiagonal systems of ng rows each.   !
! Row k of line l is stored at rr(l,k), rows are ld apart in   !
! memory. All lines share the off-diagonals am(k), ap(k). The  !
! main diagonal of line l is ac(k) + sh(l), with dcs(l) added  !
! to the first row and dce(l) added to the last row:           !
!                                                              !
!   am(k)*x(k-1) + (ac(k)+sh(l))*x(k) + ap(k)*x(k+1) = rr(l,k) !
!                                                              !
! No pivoting: the systems must be diagonally dominant, which  !
! holds for all momentum and pressure Helmholtz equations.     !
! ib(nl,ng) is a work buffer for the inverse pivots.           !
!                                                              !
!--------------------------------------------------------------!

subroutine ThomasBatched(nl,ng,ld,am,ac,ap,sh,dcs,dce,rr,ib)

    implicit none

    integer, intent(in) :: nl,ng,ld                 ! Number of lines, rows per line, memory distance between rows
    real, intent(in)    :: am(ng),ac(ng),ap(ng)     ! Sub-, main and super-diagonal, shared by all lines
    real, intent(in)    :: sh(nl),dcs(nl),dce(nl)   ! Per line: diagonal shift, first/last row correction
    real, intent(inout) :: rr(ld,*)                 ! In: right-hand side, out: solution (row k of line l at rr(l,k))
    real, intent(out)   :: ib(nl,ng)                ! Work buffer: inverse of the pivots

    integer             :: l,k
    real                :: w                        ! Elimination multiplier

    ! Single row (e.g. a face-centred solve with only 2 cells, or a pinned zero-mode
    ! line with 2 rows): both boundary corrections act on the same diagonal entry
    if (ng.eq.1) then
        do l = 1,nl
            rr(l,1) = rr(l,1)/(ac(1) + sh(l) + dcs(l) + dce(l))
        end do
        return
    end if

    ! Forward elimination: row k-1 times w = am(k)/pivot(k-1) is subtracted from
    ! row k, which removes the sub-diagonal entry am(k). The new pivot of row k is
    ! diag(k) - w*ap(k-1). The inverse pivots are stored in ib, so the backward
    ! substitution needs only multiplications. In all loops below, the line index
    ! l is innermost: the lines are independent and contiguous in memory, so the
    ! compiler vectorises over them.

    ! First row, including the start boundary correction
    do l = 1,nl
        ib(l,1) = 1.0/(ac(1) + sh(l) + dcs(l))
    end do

    ! Interior rows
    do k = 2,ng-1
        do l = 1,nl
            w       = am(k)*ib(l,k-1)
            ib(l,k) = 1.0/(ac(k) + sh(l) - w*ap(k-1))
            rr(l,k) = rr(l,k) - w*rr(l,k-1)
        end do
    end do

    ! Last row, including the end boundary correction. Dividing by its pivot
    ! already gives the solution in the last row.
    do l = 1,nl
        w        = am(ng)*ib(l,ng-1)
        ib(l,ng) = 1.0/(ac(ng) + sh(l) + dce(l) - w*ap(ng-1))
        rr(l,ng) = (rr(l,ng) - w*rr(l,ng-1))*ib(l,ng)
    end do

    ! Backward substitution: x(k) = (rhs(k) - ap(k)*x(k+1))/pivot(k)
    do k = ng-1,1,-1
        do l = 1,nl
            rr(l,k) = (rr(l,k) - ap(k)*rr(l,k+1))*ib(l,k)
        end do
    end do

    return

end subroutine ThomasBatched

!-------------- BATCHED PERIODIC THOMAS ALGORITHM -------------!
!                                                              !
! Same as ThomasBatched, but the systems are periodic: row 1   !
! couples to x(ng) through am(1) and row ng couples to x(1)    !
! through ap(ng). The corners are written as a rank-one        !
! update u*v^T of a tridiagonal matrix and removed with the    !
! Sherman-Morrison formula:                                    !
!                                                              !
!   u = (gm, 0, ..., 0, ap(ng)),  v = (1, 0, ..., 0, am(1)/gm) !
!                                                              !
! with gm = -(ac(1)+sh(l)). Both the right-hand side and u are !
! eliminated in the same sweep. zz(nl,ng) is a work buffer     !
! for u. Requires ng >= 3.                                     !
!                                                              !
!--------------------------------------------------------------!

subroutine ThomasBatchedPeriodic(nl,ng,ld,am,ac,ap,sh,rr,ib,zz)

    implicit none

    integer, intent(in) :: nl,ng,ld                 ! Number of lines, rows per line, memory distance between rows
    real, intent(in)    :: am(ng),ac(ng),ap(ng)     ! Sub-, main and super-diagonal, shared by all lines
    real, intent(in)    :: sh(nl)                   ! Per line: diagonal shift
    real, intent(inout) :: rr(ld,*)                 ! In: right-hand side, out: solution (row k of line l at rr(l,k))
    real, intent(out)   :: ib(nl,ng),zz(nl,ng)      ! Work buffers: inverse pivots, and z = T^-1 u

    integer             :: l,k
    real                :: w,vn                     ! Elimination multiplier, last entry of v
    real                :: gm(nl),qf(nl)            ! Per line: gm (see above), Sherman-Morrison factor

    ! The periodic matrix A is split as A = T + u*v^T, where T is tridiagonal
    ! (A with its first and last diagonal entries modified and the corners
    ! removed). With y = T^-1 rhs and z = T^-1 u, the solution is
    !
    !   x = y - (v.y)/(1 + v.z) * z
    !
    ! y and z are computed together in one forward/backward sweep over T.

    ! First row: T(1,1) = A(1,1) - u(1)*v(1) = A(1,1) - gm, and u(1) = gm
    do l = 1,nl
        gm(l)   = -(ac(1) + sh(l))
        ib(l,1) = 1.0/(ac(1) + sh(l) - gm(l))
        zz(l,1) = gm(l)
    end do

    ! Interior rows: u is zero here, only the fill-in from row 1 propagates into z
    do k = 2,ng-1
        do l = 1,nl
            w       = am(k)*ib(l,k-1)
            ib(l,k) = 1.0/(ac(k) + sh(l) - w*ap(k-1))
            rr(l,k) = rr(l,k) - w*rr(l,k-1)
            zz(l,k) = -w*zz(l,k-1)
        end do
    end do

    ! Last row: T(ng,ng) = A(ng,ng) - u(ng)*v(ng) = A(ng,ng) - ap(ng)*am(1)/gm,
    ! and u(ng) = ap(ng)
    do l = 1,nl
        w        = am(ng)*ib(l,ng-1)
        ib(l,ng) = 1.0/(ac(ng) + sh(l) - ap(ng)*am(1)/gm(l) - w*ap(ng-1))
        rr(l,ng) = (rr(l,ng) - w*rr(l,ng-1))*ib(l,ng)
        zz(l,ng) = (ap(ng)   - w*zz(l,ng-1))*ib(l,ng)
    end do

    ! Backward substitution, for y (in rr) and z (in zz)
    do k = ng-1,1,-1
        do l = 1,nl
            rr(l,k) = (rr(l,k) - ap(k)*rr(l,k+1))*ib(l,k)
            zz(l,k) = (zz(l,k) - ap(k)*zz(l,k+1))*ib(l,k)
        end do
    end do

    ! Sherman-Morrison factor (v.y)/(1 + v.z), with v = (1, 0, ..., 0, am(1)/gm)
    do l = 1,nl
        vn    = am(1)/gm(l)
        qf(l) = (rr(l,1) + vn*rr(l,ng))/(1.0 + zz(l,1) + vn*zz(l,ng))
    end do

    ! Correction x = y - qf*z
    do k = 1,ng
        do l = 1,nl
            rr(l,k) = rr(l,k) - qf(l)*zz(l,k)
        end do
    end do

    return

end subroutine ThomasBatchedPeriodic
