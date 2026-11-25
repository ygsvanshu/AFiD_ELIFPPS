!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                              !
!    FILE: BoundaryConditions.F90                              !
!    CONTAINS: subroutines V(x/y/z)Bc(X/Y/Z)(s/e),             !
!              subroutines GlFl(X/Y/Z)(s/e),                   !
!              subroutine PrBcTp and                           !
!              subroutines PrBc(X/Y/Z)(s/e)                    !
!                                                              !
!    PURPOSE: Sets the boundary conditions for each velocity   !
!    component Vx, Vy, Vz in three directions X, Y, Z at the   !
!    start (s) or end (e) of the axis, the flux at the start   !
!    (s) or end (e) of the axis, and the boundary conditions   !
!    for the pressure correction.                              !
!                                                              !
!    All B.C. functions take two spatial coordinate indices as !
!    inputs, to set boundary type "btyp" which can be          !
!    "DIRICHLET", "NEUMANN", or "SOMMERFELD". The value at the !
!    boundary is set by "bval".                                !
!                                                              !
!    All flux functions take two spatial coordinate indices as !
!    inputs, to set the flux correction magnitude "gval" at    !
!    the boundary.                                             !
!                                                              !
!    The pressure correction boundary type of each face is set !
!    in "PrBcTp", and can be "DIRICHLET" or "NEUMANN"          !
!    (default). It is uniform over each face. If any face is   !
!    DIRICHLET, the global flux correction is not applied. All !
!    pressure B.C. functions take two spatial coordinate       !
!    indices as inputs, to set "pval", which is the pressure   !
!    at the boundary for DIRICHLET faces, or the pressure      !
!    gradient along the positive axis direction at the         !
!    boundary for NEUMANN faces. For NEUMANN faces, a non-zero !
!    "pval" changes the boundary-normal velocity in the        !
!    pressure correction step, and the global flux correction  !
!    accounts for the resulting net flux.                      !
!                                                              !
!    Additional information such as coordinates, grid spacing, !
!    time, and time-step information can also be accessed in   !
!    the subroutines. The variable "time" gives the flow time  !
!    at the end of substep, and the product of the variables   !
!    "al*dt" gives the substep time duration.                  !
!                                                              !
!    For Vx grid, the coordinates are xc(i), ym(j), zm(k)      !
!    For Vy grid, the coordinates are xm(i), yc(j), zm(k)      !
!    For Vz grid, the coordinates are xm(i), ym(j), zc(k)      !
!    For Vx grid, the spacings are dxc(i), dym(j), dzm(k)      !
!    For Vy grid, the spacings are dxm(i), dyc(j), dzm(k)      !
!    For Vz grid, the spacings are dxm(i), dym(j), dzc(k)      !
!    For pressure, the coordinates are xm(i), ym(j), zm(k)     !
!    For pressure, the spacings are dxm(i), dym(j), dzm(k)     !
!                                                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!================ VX ================!

subroutine VxBcXs(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VxBcXs

subroutine VxBcXe(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VxBcXe

subroutine VxBcYs(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VxBcYs

subroutine VxBcYe(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VxBcYe

subroutine VxBcZs(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VxBcZs

subroutine VxBcZe(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VxBcZe

!================ VY ================!

subroutine VyBcXs(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VyBcXs

subroutine VyBcXe(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 1.0

end subroutine VyBcXe

subroutine VyBcYs(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VyBcYs

subroutine VyBcYe(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VyBcYe

subroutine VyBcZs(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VyBcZs

subroutine VyBcZe(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VyBcZe

!================ VZ ================!

subroutine VzBcXs(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VzBcXs

subroutine VzBcXe(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VzBcXe

subroutine VzBcYs(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VzBcYs

subroutine VzBcYe(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VzBcYe

subroutine VzBcZs(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VzBcZs

subroutine VzBcZe(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = DIRICHLET
    bval = 0.0

end subroutine VzBcZe

!=============== Flux ===============!

subroutine GlFlXs(j,k,gval)

    implicit none

    integer, intent(in)     :: j,k
    real,    intent(out)    :: gval

    gval = 0.0

end subroutine GlFlXs

subroutine GlFlXe(j,k,gval)

    implicit none

    integer, intent(in)     :: j,k
    real,    intent(out)    :: gval

    gval = 0.0

end subroutine GlFlXe   

subroutine GlFlYs(i,k,gval)

    implicit none

    integer, intent(in)     :: i,k
    real,    intent(out)    :: gval

    gval = 0.0

end subroutine GlFlYs

subroutine GlFlYe(i,k,gval)

    implicit none

    integer, intent(in)     :: i,k
    real,    intent(out)    :: gval

    gval = 0.0

end subroutine GlFlYe

subroutine GlFlZs(i,j,gval)

    implicit none

    integer, intent(in)     :: i,j
    real,    intent(out)    :: gval

    gval = 0.0

end subroutine GlFlZs

subroutine GlFlZe(i,j,gval)

    implicit none

    integer, intent(in)     :: i,j
    real,    intent(out)    :: gval

    gval = 0.0

end subroutine GlFlZe

!============= Pressure =============!

subroutine PrBcTp(xs,xe,ys,ye,zs,ze)

    use boundary_arrays, only: DIRICHLET, NEUMANN

    implicit none

    integer, intent(inout)  :: xs,xe,ys,ye,zs,ze

    ! By default everything is set to NEUMANN, changes are made below if needed

end subroutine PrBcTp

subroutine PrBcXs(j,k,pval)

    implicit none

    integer, intent(in)     :: j,k
    real,    intent(out)    :: pval

    pval = 0.0

end subroutine PrBcXs

subroutine PrBcXe(j,k,pval)

    implicit none

    integer, intent(in)     :: j,k
    real,    intent(out)    :: pval

    pval = 0.0

end subroutine PrBcXe

subroutine PrBcYs(i,k,pval)

    implicit none

    integer, intent(in)     :: i,k
    real,    intent(out)    :: pval

    pval = 0.0

end subroutine PrBcYs

subroutine PrBcYe(i,k,pval)

    implicit none

    integer, intent(in)     :: i,k
    real,    intent(out)    :: pval

    pval = 0.0

end subroutine PrBcYe

subroutine PrBcZs(i,j,pval)

    implicit none

    integer, intent(in)     :: i,j
    real,    intent(out)    :: pval

    pval = 0.0

end subroutine PrBcZs

subroutine PrBcZe(i,j,pval)

    implicit none

    integer, intent(in)     :: i,j
    real,    intent(out)    :: pval

    pval = 0.0

end subroutine PrBcZe