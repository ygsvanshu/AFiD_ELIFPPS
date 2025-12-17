!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: CreateInitialConditions.F90                    !
!    CONTAINS: subroutine CreateInitialConditions         !
!                                                         !
!    PURPOSE: Initialization routine. Sets initial        !
!    conditions for pressure and velocity components.     !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine CreateInitialConditions

    use decomp_2d, only: xstart,xend
    use param
    use local_arrays, only: vx,vy,vz,pr
    
    implicit none

    integer :: i,j,k,l
    real    :: rval
    real    :: xxx,yyy,zzz

    pr(:,:,:) = 0.0

    call random_seed()

    do j = xstart(2),xend(2)
        do k = xstart(3),xend(3)

            do i = 1,nx

                xxx = xc(i)/xlen
                yyy = ym(j)/ylen
                zzz = zm(k)/zlen

                vy(i,j,k) = 0.0

                do l = 1,epsnum
                    call random_number(rval)
                    rval      = 2.0*eps*(rval-0.5)
                    rval      = rval*sin(2.0*pi*l*xxx)*abs(sin(2.0*pi*l*xxx))
                    rval      = rval*sin(2.0*pi*l*yyy)*abs(sin(2.0*pi*l*yyy))
                    rval      = rval*sin(2.0*pi*l*zzz)*abs(sin(2.0*pi*l*zzz))
                    vx(i,j,k) = vx(i,j,k) + rval
                    call random_number(rval)
                end do

            end do

            do i = 0,nx

                xxx = xm(i)/xlen
                yyy = yc(j)/ylen
                zzz = zm(k)/zlen

                vy(i,j,k) = 0.0

                do l = 1,epsnum
                    call random_number(rval)
                    rval      = 2.0*eps*(rval-0.5)
                    rval      = rval*sin(2.0*pi*l*xxx)*abs(sin(2.0*pi*l*xxx))
                    rval      = rval*sin(2.0*pi*l*yyy)*abs(sin(2.0*pi*l*yyy))
                    rval      = rval*sin(2.0*pi*l*zzz)*abs(sin(2.0*pi*l*zzz))
                    vy(i,j,k) = vy(i,j,k) + rval
                end do

                xxx = xm(i)/xlen
                yyy = ym(j)/ylen
                zzz = zc(k)/zlen

                vz(i,j,k) = 0.0

                do l = 1,epsnum
                    call random_number(rval)
                    rval      = 2.0*eps*(rval-0.5)
                    rval      = rval*sin(2.0*pi*l*xxx)*abs(sin(2.0*pi*l*xxx))
                    rval      = rval*sin(2.0*pi*l*yyy)*abs(sin(2.0*pi*l*yyy))
                    rval      = rval*sin(2.0*pi*l*zzz)*abs(sin(2.0*pi*l*zzz))
                    vz(i,j,k) = vz(i,j,k) + rval
                    call random_number(rval)
                end do

            end do

        end do
    end do

    return

end subroutine CreateInitialConditions