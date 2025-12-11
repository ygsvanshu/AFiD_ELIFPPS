!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: DeallocateVariables.F90                        !
!    CONTAINS: subroutine DeallocateVariables             !
!                                                         ! 
!    PURPOSE: Finalization routine. Deallocates all       !
!     variables used in the code                          !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DeallocateVariables

    use decomp_2d
    use param
    use local_arrays
    use boundary_arrays
    use stat_arrays
    use movie_arrays
    use AuxiliaryRoutines

    implicit none
    
    call DestroyReal1DArray(xc)
    call DestroyReal1DArray(xm)

    call DestroyReal1DArray(dxc)
    call DestroyReal1DArray(dxm)

    call DestroyReal1DArray(ixcm)
    call DestroyReal1DArray(ixcp)

    call DestroyReal1DArray(ap1ci)
    call DestroyReal1DArray(ac1ci)
    call DestroyReal1DArray(am1ci)

    call DestroyReal1DArray(ap1si)
    call DestroyReal1DArray(ac1si)
    call DestroyReal1DArray(am1si)

    call DestroyReal1DArray(amphi)
    call DestroyReal1DArray(acphi)
    call DestroyReal1DArray(apphi)

    call DestroyReal1DArray(ak1)
    call DestroyReal1DArray(ao)


    call DestroyReal1DArray(yc)
    call DestroyReal1DArray(ym)

    call DestroyReal1DArray(dyc)
    call DestroyReal1DArray(dym)

    call DestroyReal1DArray(iycm)
    call DestroyReal1DArray(iycp)

    call DestroyReal1DArray(ap2cj)
    call DestroyReal1DArray(ac2cj)
    call DestroyReal1DArray(am2cj)

    call DestroyReal1DArray(ap2sj)
    call DestroyReal1DArray(ac2sj)
    call DestroyReal1DArray(am2sj)

    call DestroyReal1DArray(amphj)
    call DestroyReal1DArray(acphj)
    call DestroyReal1DArray(apphj)

    call DestroyReal1DArray(ak2)
    call DestroyReal1DArray(ap)


    call DestroyReal1DArray(zc)
    call DestroyReal1DArray(zm)

    call DestroyReal1DArray(dzc)
    call DestroyReal1DArray(dzm)

    call DestroyReal1DArray(izcm)
    call DestroyReal1DArray(izcp)

    call DestroyReal1DArray(ap3ck)
    call DestroyReal1DArray(ac3ck)
    call DestroyReal1DArray(am3ck)

    call DestroyReal1DArray(ap3sk)
    call DestroyReal1DArray(ac3sk)
    call DestroyReal1DArray(am3sk)

    call DestroyReal1DArray(amphk)
    call DestroyReal1DArray(acphk)
    call DestroyReal1DArray(apphk)

    call DestroyReal1DArray(ak3)
    call DestroyReal1DArray(aq)

    !-------------------------------------------------
    ! Arrays for statistical profiles
    !-------------------------------------------------

    if (save1d) then

        call DestroyReal1DArray(vx_m1_xp)
        call DestroyReal1DArray(vx_m2_xp)
        call DestroyReal1DArray(vy_m1_xp)
        call DestroyReal1DArray(vy_m2_xp)
        call DestroyReal1DArray(vz_m1_xp)
        call DestroyReal1DArray(vz_m2_xp)

        call DestroyReal1DArray(vx_m1_yp)
        call DestroyReal1DArray(vx_m2_yp)
        call DestroyReal1DArray(vy_m1_yp)
        call DestroyReal1DArray(vy_m2_yp)
        call DestroyReal1DArray(vz_m1_yp)
        call DestroyReal1DArray(vz_m2_yp)

        call DestroyReal1DArray(vx_m1_zp)
        call DestroyReal1DArray(vx_m2_zp)
        call DestroyReal1DArray(vy_m1_zp)
        call DestroyReal1DArray(vy_m2_zp)
        call DestroyReal1DArray(vz_m1_zp)
        call DestroyReal1DArray(vz_m2_zp)

    end if

    !-------------------------------------------------
    ! Arrays for boundaries
    !-------------------------------------------------

    ! Arrays to store boundary types (for default pencil)
    call DestroyInteger3DArray(btypxs)
    call DestroyInteger3DArray(btypxe)
    if (xstart(2).eq.1) call  DestroyInteger3DArray(btypys)
    if (xend(2).eq.nym) call  DestroyInteger3DArray(btypye)
    if (xstart(3).eq.1) call  DestroyInteger3DArray(btypzs)
    if (xend(3).eq.nzm) call  DestroyInteger3DArray(btypze)
    ! Arrays to store boundary values (for default pencil)
    call DestroyReal3DArray(bvalxs)
    call DestroyReal3DArray(bvalxe)
    if (xstart(2).eq.1) call DestroyReal3DArray(bvalys)
    if (xend(2).eq.nym) call DestroyReal3DArray(bvalye)
    if (xstart(3).eq.1) call DestroyReal3DArray(bvalzs)
    if (xend(3).eq.nzm) call DestroyReal3DArray(bvalze)
    ! Arrays to store coefficients for implicit solver
    call DestroyReal3DArray(cfbcxs)
    call DestroyReal3DArray(cfbcxe)
    call DestroyReal3DArray(cfbcys)
    call DestroyReal3DArray(cfbcye)
    call DestroyReal3DArray(cfbczs)
    call DestroyReal3DArray(cfbcze)
    ! Arrays to store whether global flux balance is allowed
    call DestroyLogical2DArray(glpcxs)
    call DestroyLogical2DArray(glpcxe)
    if (xstart(2).eq.1) call DestroyLogical2DArray(glpcys)
    if (xend(2).eq.nym) call DestroyLogical2DArray(glpcye)
    if (xstart(3).eq.1) call DestroyLogical2DArray(glpczs)
    if (xend(3).eq.nzm) call DestroyLogical2DArray(glpcze)

    !-------------------------------------------------
    ! Arrays for movie slices
    !-------------------------------------------------

    if (save2d) then

        call DestroyReal1DArray(slcx)
        call DestroyReal1DArray(slcy)
        call DestroyReal1DArray(slcz)

        call DestroyInteger2DArray(sidx)
        call DestroyInteger2DArray(sidy)
        call DestroyInteger2DArray(sidz)

        call DestroyReal2DArray(scfx)
        call DestroyReal2DArray(scfy)
        call DestroyReal2DArray(scfz)

        call DestroyReal2DArray(mslx)
        call DestroyReal2DArray(msly)
        call DestroyReal2DArray(mslz)

    end if

    !-------------------------------------------------
    ! Arrays with ghost cells
    !-------------------------------------------------
    
    call DestroyReal3DArray(vx)
    call DestroyReal3DArray(vy)
    call DestroyReal3DArray(vz)
    call DestroyReal3DArray(pr)

    call DestroyReal3DArray(dphhalo)

    !-----------------------------------------------
    ! Arrays without ghost cells
    !-----------------------------------------------

    call DestroyReal3DArray(dph)

    call DestroyReal3DArray(adx)
    call DestroyReal3DArray(ady)
    call DestroyReal3DArray(adz)

    call DestroyReal3DArray(rux)
    call DestroyReal3DArray(ruy)
    call DestroyReal3DArray(ruz)

    if (pre_diff) then
        call DestroyReal3DArray(dfx)
        call DestroyReal3DArray(dfy)
        call DestroyReal3DArray(dfz)
    end if

    call DestroyReal3DArray(rhx)
    call DestroyReal3DArray(rhy)
    call DestroyReal3DArray(rhz)

    return 

end subroutine DeallocateVariables