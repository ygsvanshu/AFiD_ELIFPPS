!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ReadParticleInputs.F90                         !
!    CONTAINS: subroutine ReadParticleInputs              !
!                                                         !
!    PURPOSE: Subroutine to read inputs related to the    !
!    particle solver
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine ReadParticleInputs

    use mpih
    use param
    use lagrangian_point_particle

    implicit none
    
    logical             :: exists
    integer             :: narg,ierror
    integer,parameter   :: nimp = 100
    character*200       :: filename,line
    character*100       :: ss(nimp)
    character*1         :: stringdummy1

    filename = "Inputs/particle.in"

    inquire(file=filename,exist=exists)
    if (exists) then 

        particle = .true.
    
        open (unit=15,file=filename,status='old')

        ierror = 0
        do while(ierror.eq.0)

            read (15,'(a200)',iostat=ierror) line

            if (line(1:3) == '101') then
                !###### DRAG MODEL ######
                call scan_string(line,1,ss,narg)
                read (ss(1),*) lpp_dmod

            else if (line(1:3) == '102') then
                !###### SLIP CORRECTION ######
                call scan_string(line,1,ss,narg)
                stringdummy1 = ss(1)
                if ('y' == stringdummy1) then
                    lpp_scor = .true.
                else if ('n' == stringdummy1) then
                    lpp_scor = .false.
                else
                    write (*,*) "ERROR: Input value of parameter SCOR not valid"
                    write (*,*) "       ===> valid values 'n' or 'y'           "
                    call stop_config
                end if

            else if (line(1:3) == '201') then
                !###### EULERIAN - LAGRANGIAN COUPLING ######
                call scan_string(line,2,ss,narg)
                read (ss(1),*) e2l_mult
                read (ss(2),*) l2e_mult

            else if (line(1:3) == '301') then
                !###### GRAVITY ######
                call scan_string(line,3,ss,narg)
                read (ss(1),*) lpp_grav(1)
                read (ss(2),*) lpp_grav(2)
                read (ss(3),*) lpp_grav(3)

            else if (line(1:3) == '401') then
                !###### TIME STEP LIMITS ######
                call scan_string(line,1,ss,narg)
                read (ss(1),*) lpp_stol

            else if (line(1:3) == '402') then
                !###### TIME STEP LIMITS ######
                call scan_string(line,1,ss,narg)
                read (ss(1),*) lpp_clim

            else if (line(1:3) == '403') then
                !###### TIME STEP LIMITS ######
                call scan_string(line,1,ss,narg)
                read (ss(1),*) lpp_tlim

            else if (line(1:3) == '501') then
                !###### SAVING PARTICLE EXIT EVENTS ######
                call scan_string(line,3,ss,narg)
                stringdummy1 = ss(1)
                if ('y' == stringdummy1) then
                    pex_save = .true.
                else if ('n' == stringdummy1) then
                    pex_save = .false.
                else
                    write (*,*) "ERROR: Input value of parameter PSAV not valid"
                    write (*,*) "       ===> valid values 'n' or 'y'           "
                    call stop_config
                end if
                read (ss(2),*) pex_ssta
                read (ss(3),*) pex_sfrq
            
            else if (line(1:3) == '502') then
                !###### SAVING PARTICLE SNAPSHOTS ######
                call scan_string(line,3,ss,narg)
                stringdummy1 = ss(1)
                if ('y' == stringdummy1) then
                    lpp_save = .true.
                else if ('n' == stringdummy1) then
                    lpp_save = .false.
                else
                    write (*,*) "ERROR: Input value of parameter PSAV not valid"
                    write (*,*) "       ===> valid values 'n' or 'y'           "
                    call stop_config
                end if
                read (ss(2),*) lpp_ssta
                read (ss(3),*) lpp_sfrq

            end if

        end do

        pre_diff = lpp_scor

        close (15)

    end if

end subroutine ReadParticleInputs