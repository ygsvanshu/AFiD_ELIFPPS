!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: BoundaryConditions.F90                         !
!    CONTAINS: subroutines V(x/y/z)Bc(X/Y/Z)(s/e)         !
!                                                         !
!    PURPOSE: Sets the boundary conditions for each       !
!    velocity component Vx, Vy, Vz in three directions    !
!    X, Y, Z at the start (s) or end (e) of the axis.     !
!                                                         !
!    All B.C. functions take two spatial coordinate       !
!    indices as inputs, to set boundary type "btyp" which !
!    can be "DIRICHLET", "NEUMANN", OR "SOMMERFELD".      !
!    The value at the boundary is set by bval. Additional !
!    information such as coordinates, grid spacing, time, !
!    and time-step infomration can also be accessed in    !
!    the subroutines. The variable "time" gives the flow  ! 
!    time at the end of substep, and the product of the   !
!    variables "al*dt" gives the substep time duration    !
!    For Vx grid, the coordinates are xc(i), ym(j), zm(k) !
!    For Vy grid, the coordinates are xm(i), yc(j), zm(k) !
!    For Vz grid, the coordinates are xm(i), ym(j), zc(k) !
!    For Vx grid, the spacings are dxm(i), dyc(j), dzc(k) !
!    For Vy grid, the spacings are dxc(i), dym(j), dzc(k) !
!    For Vz grid, the spacings are dxc(i), dyc(j), dzm(k) !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!================ VX ================!

subroutine VxBcXs(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VxBcXs

subroutine VxBcXe(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VxBcXe

subroutine VxBcYs(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VxBcYs

subroutine VxBcYe(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VxBcYe

subroutine VxBcZs(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VxBcZs

subroutine VxBcZe(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
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

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VyBcXs

subroutine VyBcXe(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VyBcXe

subroutine VyBcYs(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VyBcYs

subroutine VyBcYe(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VyBcYe

subroutine VyBcZs(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VyBcZs

subroutine VyBcZe(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
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

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VzBcXs

subroutine VzBcXe(j,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: j,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VzBcXe

subroutine VzBcYs(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VzBcYs

subroutine VzBcYe(i,k,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,k
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VzBcYe

subroutine VzBcZs(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VzBcZs

subroutine VzBcZe(i,j,btyp,bval)

    use param
    use boundary_arrays, only: DIRICHLET, NEUMANN, SOMMERFELD

    implicit none

    integer, intent(in)     :: i,j
    integer, intent(out)    :: btyp
    real,    intent(out)    :: bval

    btyp = SOMMERFELD
    bval = 0.0

end subroutine VzBcZe