!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: GlobalParticleStatistics.F90                   !
!    CONTAINS: subroutine GlobalParticleStatistics        !
!                                                         !
!    PURPOSE: To update the global statistics of particle !
!    velocities and counts                                !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine GlobalParticleStatistics

    use mpih
    use param
    use decomp_2d, only: xstart,xend
    use lagrangian_point_particle

    implicit none

    character*200   :: filename
    logical         :: gbound,exists
    integer         :: kc,kp
    integer         :: jc,jp
    integer         :: ic,ip
    integer         :: ndir,nlpp
    integer         :: int_dummy
    real            :: cvol,fvol
    real            :: res_dummy

    ! Initialize particle velocity statistics to zero
    lpp_vmax(:) = 0.0
    lpp_vavg(:) = 0.0
    lpp_vrms(:) = 0.0
    ! Initialize body force statistics to zero
    lpp_bmax(:) = 0.0
    lpp_bavg(:) = 0.0
    lpp_brms(:) = 0.0

    ! Get the total counts of spawned and exited particles
    !! Call MPI_Reduce to get the total spawn and exit counts
    call MPI_REDUCE(lpp_spwn,int_dummy,1,MPI_INTEGER,MPI_SUM,0,mpi_comm,mpi_ierr)
    tot_spwn = tot_spwn + int_dummy
    call MPI_REDUCE(lpp_exit,int_dummy,1,MPI_INTEGER,MPI_SUM,0,mpi_comm,mpi_ierr)
    tot_exit = tot_exit + int_dummy
    !! Call MPI_Allreduce to get the total active particles
    call MPI_ALLREDUCE((lpp_actv-sub_exit),int_dummy,1,MPI_INTEGER,MPI_SUM,mpi_comm,mpi_ierr)
    tot_actv = int_dummy
    !! Reset the count of spawned and exited particles
    lpp_spwn = 0
    lpp_exit = 0

    ! Check if there are any globally active particles to write statistical data for
    if (tot_actv.gt.0) then

        !! Compute particle maximum, mean, and RMS velocities
        do ndir = 1,3
            
            do nlpp = 1,lpp_actv
                if (abs(lpp_list(nlpp)%lpp_vel(ndir)) .gt. lpp_vmax(ndir)) lpp_vmax(ndir) = abs(lpp_list(nlpp)%lpp_vel(ndir))
                lpp_vavg(ndir) = lpp_vavg(ndir) + (lpp_list(nlpp)%lpp_vel(ndir))/tot_actv
                lpp_vrms(ndir) = lpp_vrms(ndir) + (lpp_list(nlpp)%lpp_vel(ndir)**2.0)/tot_actv
            end do

            call MPI_REDUCE(lpp_vmax(ndir),res_dummy,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpi_comm,mpi_ierr)
            lpp_vmax(ndir) = res_dummy
            call MPI_REDUCE(lpp_vavg(ndir),res_dummy,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm,mpi_ierr)
            lpp_vavg(ndir) = res_dummy
            call MPI_REDUCE(lpp_vrms(ndir),res_dummy,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm,mpi_ierr)
            lpp_vrms(ndir) = sqrt(res_dummy)

        end do

        fvol = xlen*zlen*ylen

        !! Compute the maximum, mean, and RMS body force values
        do kc=xstart(3),xend(3)
            kp = kc+1
            do jc=xstart(2),xend(2)
                jp = jc+1
                do ic=1,nxm
                    ip = ic+1

                    cvol = dxc(ic)*dyc(jc)*dzc(kc)

                    lpp_bmax(1) = max(lpp_bmax(1),abs(lpp_bdfx(ic,jc,kc)),abs(lpp_bdfx(ip,jc,kc)))
                    lpp_bmax(2) = max(lpp_bmax(2),abs(lpp_bdfy(ic,jc,kc)),abs(lpp_bdfy(ic,jp,kc)))
                    lpp_bmax(3) = max(lpp_bmax(3),abs(lpp_bdfz(ic,jc,kc)),abs(lpp_bdfz(ic,jc,kp)))

                    lpp_bavg(1) = lpp_bavg(1) + 0.5*(cvol/fvol)*(lpp_bdfx(ic,jc,kc) + lpp_bdfx(ip,jc,kc))
                    lpp_bavg(2) = lpp_bavg(2) + 0.5*(cvol/fvol)*(lpp_bdfy(ic,jc,kc) + lpp_bdfy(ic,jp,kc))
                    lpp_bavg(3) = lpp_bavg(3) + 0.5*(cvol/fvol)*(lpp_bdfz(ic,jc,kc) + lpp_bdfz(ic,jc,kp))

                    lpp_brms(1) = lpp_brms(1) + 0.5*(cvol/fvol)*(lpp_bdfx(ic,jc,kc)**2.0 + lpp_bdfx(ip,jc,kc)**2.0)
                    lpp_brms(2) = lpp_brms(2) + 0.5*(cvol/fvol)*(lpp_bdfy(ic,jc,kc)**2.0 + lpp_bdfy(ic,jp,kc)**2.0)
                    lpp_brms(3) = lpp_brms(3) + 0.5*(cvol/fvol)*(lpp_bdfz(ic,jc,kc)**2.0 + lpp_bdfz(ic,jc,kp)**2.0)
                    
                enddo
            enddo
        enddo

        do ndir = 1,3

            call MPI_REDUCE(lpp_bmax(ndir),res_dummy,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,mpi_comm,mpi_ierr)
            lpp_bmax(ndir) = res_dummy
            call MPI_REDUCE(lpp_bavg(ndir),res_dummy,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm,mpi_ierr)
            lpp_bavg(ndir) = res_dummy
            call MPI_REDUCE(lpp_brms(ndir),res_dummy,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,mpi_comm,mpi_ierr)
            lpp_brms(ndir) = sqrt(res_dummy)

        end do

    end if

    ! Write the global data to file
    if (ismaster) then

        filename = trim("Results/particle_stats.out")
        inquire(file=filename,exist=exists)
        open(unit=140,file=filename,access="sequential",status="unknown",position="append")
        if (.not.exists) write (140,"((A16,X),3(A16,X),18(A16,X))") "Time","Injected","Active","Exited","Re_Max_Vx","Re_Max_Vy","Re_Max_Vz","Re_Avg_Vx","Re_Avg_Vy","Re_Avg_Vz","Re_RMS_Vx","Re_RMS_Vy","Re_RMS_Vz","BF_Max_Vx","BF_Max_Vy","BF_Max_Vz","BF_Avg_Vx","BF_Avg_Vy","BF_Avg_Vz","BF_RMS_Vx","BF_RMS_Vy","BF_RMS_Vz"
        write(140,"((ES16.8,X),3(I16,X),18(ES16.8,X))") time,tot_spwn,tot_actv,tot_exit,lpp_vmax(1),lpp_vmax(2),lpp_vmax(3),lpp_vavg(1),lpp_vavg(2),lpp_vavg(3),lpp_vrms(1),lpp_vrms(2),lpp_vrms(3),lpp_bmax(1),lpp_bmax(2),lpp_bmax(3),lpp_bavg(1),lpp_bavg(2),lpp_bavg(3),lpp_brms(1),lpp_brms(2),lpp_brms(3)
        close(140)

    end if

    return

end subroutine GlobalParticleStatistics