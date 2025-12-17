!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: InitVariables.F90                              !
!    CONTAINS: subroutine InitVariables                   !
!                                                         !
!    PURPOSE: Initialization routine. Sets to zero all    !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine InitVariables

    use decomp_2d
    use param
    use local_arrays
    use boundary_arrays
    use stat_arrays
    use movie_arrays
    use AuxiliaryRoutines

    implicit none

    !-------------------------------------------------
    ! Arrays for grid making
    !-------------------------------------------------

    call AllocateReal1DArray(xc,1,nx)
    call AllocateReal1DArray(xm,0,nx)

    call AllocateReal1DArray(dxc,0,nx)
    call AllocateReal1DArray(dxm,1,nx)

    call AllocateReal1DArray(ixcm,1,nx)
    call AllocateReal1DArray(ixcp,1,nx)

    call AllocateReal1DArray(ap1ci,1,nxm)
    call AllocateReal1DArray(ac1ci,1,nxm)
    call AllocateReal1DArray(am1ci,1,nxm)

    call AllocateReal1DArray(ap1si,1,nx)
    call AllocateReal1DArray(ac1si,1,nx)
    call AllocateReal1DArray(am1si,1,nx)

    call AllocateReal1DArray(amphi,1,nxm)
    call AllocateReal1DArray(acphi,1,nxm)
    call AllocateReal1DArray(apphi,1,nxm)

    call AllocateReal1DArray(ak1,1,nx)
    call AllocateReal1DArray(ao,1,nx)

    call AllocateReal1DArray(yc,1,ny)
    call AllocateReal1DArray(ym,0,ny)

    call AllocateReal1DArray(dyc,0,ny)
    call AllocateReal1DArray(dym,1,ny)

    call AllocateReal1DArray(iycm,1,ny)
    call AllocateReal1DArray(iycp,1,ny)

    call AllocateReal1DArray(ap2cj,1,nym)
    call AllocateReal1DArray(ac2cj,1,nym)
    call AllocateReal1DArray(am2cj,1,nym)

    call AllocateReal1DArray(ap2sj,1,ny)
    call AllocateReal1DArray(ac2sj,1,ny)
    call AllocateReal1DArray(am2sj,1,ny)

    call AllocateReal1DArray(amphj,1,nym)
    call AllocateReal1DArray(acphj,1,nym)
    call AllocateReal1DArray(apphj,1,nym)

    call AllocateReal1DArray(ak2,1,ny)
    call AllocateReal1DArray(ap,1,ny)

    call AllocateReal1DArray(zc,1,nz)
    call AllocateReal1DArray(zm,0,nz)

    call AllocateReal1DArray(dzc,0,nz)
    call AllocateReal1DArray(dzm,1,nz)

    call AllocateReal1DArray(izcm,1,nz)
    call AllocateReal1DArray(izcp,1,nz)

    call AllocateReal1DArray(ap3ck,1,nzm)
    call AllocateReal1DArray(ac3ck,1,nzm)
    call AllocateReal1DArray(am3ck,1,nzm)

    call AllocateReal1DArray(ap3sk,1,nz)
    call AllocateReal1DArray(ac3sk,1,nz)
    call AllocateReal1DArray(am3sk,1,nz)

    call AllocateReal1DArray(amphk,1,nzm)
    call AllocateReal1DArray(acphk,1,nzm)
    call AllocateReal1DArray(apphk,1,nzm)

    call AllocateReal1DArray(ak3,1,nz)
    call AllocateReal1DArray(aq,1,nz)

    !-------------------------------------------------
    ! Arrays for statistical profiles
    !-------------------------------------------------

    if (save1d) then

        call AllocateReal1DArray(vx_m1_xp,1,nx)
        call AllocateReal1DArray(vx_m2_xp,1,nx)
        call AllocateReal1DArray(vy_m1_xp,0,nx)
        call AllocateReal1DArray(vy_m2_xp,0,nx)
        call AllocateReal1DArray(vz_m1_xp,0,nx)
        call AllocateReal1DArray(vz_m2_xp,0,nx)

        call AllocateReal1DArray(vx_m1_yp,0,ny)
        call AllocateReal1DArray(vx_m2_yp,0,ny)
        call AllocateReal1DArray(vy_m1_yp,1,ny)
        call AllocateReal1DArray(vy_m2_yp,1,ny)
        call AllocateReal1DArray(vz_m1_yp,0,ny)
        call AllocateReal1DArray(vz_m2_yp,0,ny)

        call AllocateReal1DArray(vx_m1_zp,0,nz)
        call AllocateReal1DArray(vx_m2_zp,0,nz)
        call AllocateReal1DArray(vy_m1_zp,0,nz)
        call AllocateReal1DArray(vy_m2_zp,0,nz)
        call AllocateReal1DArray(vz_m1_zp,1,nz)
        call AllocateReal1DArray(vz_m2_zp,1,nz)

    end if

    !-------------------------------------------------
    ! Arrays for boundaries
    !-------------------------------------------------

    ! Arrays to store boundary types (for default pencil)
    call AllocateInteger3DArray(btypxs,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    call AllocateInteger3DArray(btypxe,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    if (xstart(2).eq.1) call AllocateInteger3DArray(btypys,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    if (xend(2).eq.nym) call AllocateInteger3DArray(btypye,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    if (xstart(3).eq.1) call AllocateInteger3DArray(btypzs,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,1,3)
    if (xend(3).eq.nzm) call AllocateInteger3DArray(btypze,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,1,3)
    ! Arrays to store boundary values (for default pencil)
    call AllocateReal3DArray(bvalxs,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    call AllocateReal3DArray(bvalxe,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    if (xstart(2).eq.1) call AllocateReal3DArray(bvalys,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    if (xend(2).eq.nym) call AllocateReal3DArray(bvalye,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo,1,3)
    if (xstart(3).eq.1) call AllocateReal3DArray(bvalzs,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,1,3)
    if (xend(3).eq.nzm) call AllocateReal3DArray(bvalze,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,1,3)
    ! Arrays to store coefficients for implicit solver
    call AllocateReal3DArray(cfbcxs,xstart(2),xend(2),xstart(3),xend(3),1,3)
    call AllocateReal3DArray(cfbcxe,xstart(2),xend(2),xstart(3),xend(3),1,3)
    call AllocateReal3DArray(cfbcys,ystart(1),yend(1),ystart(3),yend(3),1,3)
    call AllocateReal3DArray(cfbcye,ystart(1),yend(1),ystart(3),yend(3),1,3)
    call AllocateReal3DArray(cfbczs,zstart(1),zend(1),zstart(2),zend(2),1,3)
    call AllocateReal3DArray(cfbcze,zstart(1),zend(1),zstart(2),zend(2),1,3)
    ! Arrays to store whether global flux balance is allowed
    call AllocateLogical2DArray(glpcxs,xstart(2),xend(2),xstart(3),xend(3))
    call AllocateLogical2DArray(glpcxe,xstart(2),xend(2),xstart(3),xend(3))
    if (xstart(2).eq.1) call AllocateLogical2DArray(glpcys,xstart(1),xend(1),xstart(3),xend(3))
    if (xend(2).eq.nym) call AllocateLogical2DArray(glpcye,xstart(1),xend(1),xstart(3),xend(3))
    if (xstart(3).eq.1) call AllocateLogical2DArray(glpczs,xstart(1),xend(1),xstart(2),xend(2))
    if (xend(3).eq.nzm) call AllocateLogical2DArray(glpcze,xstart(1),xend(1),xstart(2),xend(2))

    !-------------------------------------------------
    ! Arrays for movie slices
    !-------------------------------------------------

    if (save2d) then
        call AllocateReal2DArray(mslx,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
        call AllocateReal2DArray(msly,0,nx,xstart(3)-lvlhalo,xend(3)+lvlhalo)
        call AllocateReal2DArray(mslz,0,nx,xstart(2)-lvlhalo,xend(2)+lvlhalo)
    end if

    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------

    call AllocateReal3DArray(vx,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(vy,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(vz,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(pr,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)
    call AllocateReal3DArray(dphhalo,xstart(1)-lvlhalo,xend(1)+lvlhalo,xstart(2)-lvlhalo,xend(2)+lvlhalo,xstart(3)-lvlhalo,xend(3)+lvlhalo)

    !-----------------------------------------------
    ! Arrays without ghost cells
    !-----------------------------------------------

    call AllocateReal3DArray(dph,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

    call AllocateReal3DArray(adx,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(ady,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(adz,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    
    call AllocateReal3DArray(rux,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(ruy,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(ruz,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))

    if (pre_diff) then
        call AllocateReal3DArray(dfx,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal3DArray(dfy,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
        call AllocateReal3DArray(dfz,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    end if

    call AllocateReal3DArray(rhx,xstart(1),xend(1),xstart(2),xend(2),xstart(3),xend(3))
    call AllocateReal3DArray(rhy,ystart(1),yend(1),ystart(2),yend(2),ystart(3),yend(3))
    call AllocateReal3DArray(rhz,zstart(1),zend(1),zstart(2),zend(2),zstart(3),zend(3))

    return

end subroutine InitVariables

