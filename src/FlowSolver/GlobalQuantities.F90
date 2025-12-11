!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: GlobalQuantities.F90                           !
!    CONTAINS: subroutine GlobalQuantities                !
!                                                         ! 
!    PURPOSE: Calculate maximum, mean, rms velocity       !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GlobalQuantities

    use param
    use local_arrays, only: vy,vx,vz
    use decomp_2d, only: xstart,xend
    use mpih

    implicit none

    integer             :: kc,kp
    integer             :: jc,jp
    integer             :: ic,ip
    real                :: cvol,fvol
    real                :: res_dummy
    logical             :: fexist

    vmax(1) = 0.0
    vmax(2) = 0.0
    vmax(3) = 0.0

    vavg(1) = 0.0
    vavg(2) = 0.0
    vavg(3) = 0.0

    vrms(1) = 0.0
    vrms(2) = 0.0
    vrms(3) = 0.0

    fvol = xlen*zlen*ylen

    do kc=xstart(3),xend(3)
        kp = kc+1
        do jc=xstart(2),xend(2)
            jp = jc+1
            do ic=1,nxm
                ip = ic+1

                cvol = dxc(ic)*dyc(jc)*dzc(kc)

                vmax(1) = max(vmax(1),abs(vx(ic,jc,kc)),abs(vx(ip,jc,kc)))
                vmax(2) = max(vmax(2),abs(vy(ic,jc,kc)),abs(vy(ic,jp,kc)))
                vmax(3) = max(vmax(3),abs(vz(ic,jc,kc)),abs(vz(ic,jc,kp)))

                vavg(1) = vavg(1) + 0.5*(cvol/fvol)*(vx(ic,jc,kc) + vx(ip,jc,kc))
                vavg(2) = vavg(2) + 0.5*(cvol/fvol)*(vy(ic,jc,kc) + vy(ic,jp,kc))
                vavg(3) = vavg(3) + 0.5*(cvol/fvol)*(vz(ic,jc,kc) + vz(ic,jc,kp))

                vrms(1) = vrms(1) + 0.5*(cvol/fvol)*(vx(ic,jc,kc)**2 + vx(ip,jc,kc)**2)
                vrms(2) = vrms(2) + 0.5*(cvol/fvol)*(vy(ic,jc,kc)**2 + vy(ic,jp,kc)**2)
                vrms(3) = vrms(3) + 0.5*(cvol/fvol)*(vz(ic,jc,kc)**2 + vz(ic,jc,kp)**2)
                
            enddo
        enddo
    enddo

    call MpiMaxRealScalar(vmax(1),res_dummy)
    vmax(1) = res_dummy
    call MpiMaxRealScalar(vmax(2),res_dummy)
    vmax(2) = res_dummy
    call MpiMaxRealScalar(vmax(3),res_dummy)
    vmax(3) = res_dummy

    call MpiSumRealScalar(vavg(1),res_dummy)
    vavg(1) = res_dummy
    call MpiSumRealScalar(vavg(2),res_dummy)
    vavg(2) = res_dummy
    call MpiSumRealScalar(vavg(3),res_dummy)
    vavg(3) = res_dummy

    call MpiSumRealScalar(vrms(1),res_dummy)
    vrms(1) = res_dummy
    call MpiSumRealScalar(vrms(2),res_dummy)
    vrms(2) = res_dummy
    call MpiSumRealScalar(vrms(3),res_dummy)
    vrms(3) = res_dummy

    vmag    = dsqrt(sum(vrms))
    vrms(1) = dsqrt(vrms(1))
    vrms(2) = dsqrt(vrms(2))
    vrms(3) = dsqrt(vrms(3))

    if(ismaster) then

        inquire(file='Results/global.out',exist=fexist)
        open(94,file='Results/global.out',status='unknown',position='append',access='sequential')
        if (.not.fexist) write(94,'(11(A16,X))') 'Time', 'Max_Vx', 'Max_Vy', 'Max_Vz', 'Avg_Vx', 'Avg_Vy', 'Avg_Vz', 'RMS_Vx', 'RMS_Vy', 'RMS_Vz', 'RMS_VxVyVz'
        write(94,'(11(ES16.8,X))') time, vmax(1), vmax(2), vmax(3), vavg(1), vavg(2), vavg(3), vrms(1), vrms(2), vrms(3), vmag
        close(94)

    endif

    return

end subroutine GlobalQuantities
