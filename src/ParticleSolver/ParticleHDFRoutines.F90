!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: ParticleHDFRoutines.F90                        !
!    CONTAINS: subroutine HdfParallelWriteParticle1D      !
!    CONTAINS: subroutine HdfParallelReadParticle1D       !
!    CONTAINS: subroutine HdfParallelAppendExit1D         !
!                                                         !
!    PURPOSE: Routines to read/write particle data or     !
!    exit events to HDF5 files.                           !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HdfParallelWriteParticle1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5
    use lagrangian_point_particle, only: particle_data

    implicit none

    character*200, intent(in)           :: filename,dsetname
    type(particle_data), intent(in)     :: var(1:en1-st1+1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,link_id,plist_id,dset_id(0:9)
    integer(HID_T)                      :: filespace(0:9)
    integer(HID_T)                      :: memspace1,memspace2
    integer(HSIZE_T)                    :: dims0(1),dims1(1),dims2(2)
    integer(HSIZE_T)                    :: data_count1(1),data_count2(2)  
    integer(HSIZE_T)                    :: data_offset1(1),data_offset2(2)
    integer                             :: i,j
    integer                             :: hdf_error
    logical                             :: exists
    integer                             :: indx(1:en1-st1+1)
    real                                :: sclr(1:en1-st1+1)
    real                                :: vctr(1:en1-st1+1,1:3)

    ! For count
    dims0(1) = 1

    ! For scalar data
    dims1(1) = n1e-n1s+1
    data_count1(1) = en1-st1+1
    data_offset1(1) = st1-n1s

    ! For vector data
    dims2(1) = n1e-n1s+1
    dims2(2) = 3

    data_count2(1) = en1-st1+1
    data_count2(2) = 3

    data_offset2(1) = st1-n1s
    data_offset2(2) = 0
    
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    inquire(file=filename,exist=exists)
    call MPI_ALLREDUCE(MPI_IN_PLACE,exists,1,MPI_LOGICAL,MPI_LOR,comm,mpi_ierr)
    if (exists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
    end if
    call h5pclose_f(plist_id,hdf_error)

    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/lpp_num")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_idx")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/lpp_lft")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/lpp_dia")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/lpp_den")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/lpp_rey")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/lpp_pos")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/lpp_vel")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/acc_old")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/acc_now")

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)

    ! Write the number of particles
    call h5screate_simple_f(1,dims0,filespace(0),hdf_error)
    call h5dcreate_f(file_id,trim(dsetname)//"/lpp_num",H5T_NATIVE_INTEGER,filespace(0),dset_id(0),hdf_error,lcpl_id=link_id)
    call h5dwrite_f(dset_id(0),H5T_NATIVE_INTEGER,(n1e-n1s+1),dims0,hdf_error)
    call h5dclose_f(dset_id(0),hdf_error)

    ! Write only if there is atleast one particle (in total/global sense)
    if (n1e.ge.n1s) then

        call h5screate_simple_f(1,dims1,filespace(1),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(2),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(3),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(4),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(5),hdf_error)
        call h5screate_simple_f(2,dims2,filespace(6),hdf_error)
        call h5screate_simple_f(2,dims2,filespace(7),hdf_error)
        call h5screate_simple_f(2,dims2,filespace(8),hdf_error)
        call h5screate_simple_f(2,dims2,filespace(9),hdf_error)

        call h5screate_simple_f(1,data_count1,memspace1,hdf_error)
        call h5screate_simple_f(2,data_count2,memspace2,hdf_error)

        call h5dcreate_f(file_id,trim(dsetname)//"/src_idx",H5T_NATIVE_INTEGER,filespace(1),dset_id(1),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/lpp_lft",H5T_NATIVE_DOUBLE, filespace(2),dset_id(2),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/lpp_dia",H5T_NATIVE_DOUBLE, filespace(3),dset_id(3),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/lpp_den",H5T_NATIVE_DOUBLE, filespace(4),dset_id(4),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/lpp_rey",H5T_NATIVE_DOUBLE, filespace(5),dset_id(5),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/lpp_pos",H5T_NATIVE_DOUBLE, filespace(6),dset_id(6),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/lpp_vel",H5T_NATIVE_DOUBLE, filespace(7),dset_id(7),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/acc_old",H5T_NATIVE_DOUBLE, filespace(8),dset_id(8),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/acc_now",H5T_NATIVE_DOUBLE, filespace(9),dset_id(9),hdf_error,lcpl_id=link_id)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)

        ! Write source indices
        do i = 1,en1-st1+1
            indx(i) = var(i)%src_idx
        end do
        call h5sselect_hyperslab_f(filespace(1),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(1),H5T_NATIVE_INTEGER,indx,data_count1,hdf_error,file_space_id=filespace(1),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write particle lifetimes
        do i = 1,en1-st1+1
            sclr(i) = var(i)%lpp_lft
        end do
        call h5sselect_hyperslab_f(filespace(2),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(2),H5T_NATIVE_DOUBLE, sclr,data_count1,hdf_error,file_space_id=filespace(2),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write particle diameters
        do i = 1,en1-st1+1
            sclr(i) = var(i)%lpp_dia
        end do
        call h5sselect_hyperslab_f(filespace(3),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(3),H5T_NATIVE_DOUBLE, sclr,data_count1,hdf_error,file_space_id=filespace(3),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write particle densities
        do i = 1,en1-st1+1
            sclr(i) = var(i)%lpp_den
        end do
        call h5sselect_hyperslab_f(filespace(4),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(4),H5T_NATIVE_DOUBLE, sclr,data_count1,hdf_error,file_space_id=filespace(4),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write particle Reynolds number
        do i = 1,en1-st1+1
            sclr(i) = var(i)%lpp_rey
        end do
        call h5sselect_hyperslab_f(filespace(5),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(5),H5T_NATIVE_DOUBLE, sclr,data_count2,hdf_error,file_space_id=filespace(5),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write particle positions
        do i = 1,en1-st1+1
            do j = 1,3
                vctr(i,j) = var(i)%lpp_pos(j)
            end do
        end do
        call h5sselect_hyperslab_f(filespace(6),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dwrite_f(dset_id(6),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(6),mem_space_id=memspace2,xfer_prp=plist_id)
        ! Write particle velocities
        do i = 1,en1-st1+1
            do j = 1,3
                vctr(i,j) = var(i)%lpp_vel(j)
            end do
        end do
        call h5sselect_hyperslab_f(filespace(7),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dwrite_f(dset_id(7),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(7),mem_space_id=memspace2,xfer_prp=plist_id)
        ! Write particle accelerations (previous)
        do i = 1,en1-st1+1
            do j = 1,3
                vctr(i,j) = var(i)%acc_old(i)
            end do
        end do
        call h5sselect_hyperslab_f(filespace(8),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dwrite_f(dset_id(8),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(8),mem_space_id=memspace2,xfer_prp=plist_id)
        ! Write particle accelerations (current)
        do i = 1,en1-st1+1
            do j = 1,3
                vctr(i,j) = var(i)%acc_now(i)
            end do
        end do
        call h5sselect_hyperslab_f(filespace(9),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dwrite_f(dset_id(9),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(9),mem_space_id=memspace2,xfer_prp=plist_id)
        
        call h5pclose_f(plist_id,hdf_error)

        do i = 1,9 
            call h5dclose_f(dset_id(i),hdf_error)
            call h5sclose_f(filespace(i),hdf_error)
        end do

        call h5sclose_f(memspace2,hdf_error)
        call h5sclose_f(memspace1,hdf_error)

    end if

    call h5pclose_f(link_id,hdf_error)

    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelWriteParticle1D

subroutine HdfParallelReadParticle1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5
    use lagrangian_point_particle, only: particle_data

    implicit none

    character*200, intent(in)           :: filename,dsetname
    type(particle_data), intent(out)    :: var(1:en1-st1+1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,plist_id,dset_id(9)
    integer(HID_T)                      :: filespace(9)
    integer(HID_T)                      :: memspace1,memspace2
    integer(HSIZE_T)                    :: dims1(1),dims2(2)
    integer(HSIZE_T)                    :: data_count1(1),data_count2(2)  
    integer(HSIZE_T)                    :: data_offset1(1),data_offset2(2)
    integer                             :: i,j
    integer                             :: hdf_error
    logical                             :: exists
    integer                             :: indx(1:en1-st1+1)
    real                                :: sclr(1:en1-st1+1)
    real                                :: vctr(1:en1-st1+1,1:3)

    ! For scalar data
    dims1(1) = n1e-n1s+1
    data_count1(1) = en1-st1+1
    data_offset1(1) = st1-n1s

    ! For vector data
    dims2(1) = n1e-n1s+1
    dims2(2) = 3

    data_count2(1) = en1-st1+1
    data_count2(2) = 3

    data_offset2(1) = st1-n1s
    data_offset2(2) = 0
    
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    inquire(file=filename,exist=exists)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)

    ! Read only if there is atleast one particle (in total/global sense)
    if (n1e.ge.n1s) then

        call h5screate_simple_f(1,data_count1,memspace1,hdf_error)
        call h5screate_simple_f(2,data_count2,memspace2,hdf_error)

        call h5dopen_f(file_id,trim(dsetname)//"/src_idx",dset_id(1),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/lpp_lft",dset_id(2),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/lpp_dia",dset_id(3),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/lpp_den",dset_id(4),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/lpp_rey",dset_id(5),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/lpp_pos",dset_id(6),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/lpp_vel",dset_id(7),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/acc_old",dset_id(8),hdf_error)
        call h5dopen_f(file_id,trim(dsetname)//"/acc_now",dset_id(9),hdf_error)
        
        call h5dget_space_f(dset_id(1),filespace(1),hdf_error)
        call h5dget_space_f(dset_id(2),filespace(2),hdf_error)
        call h5dget_space_f(dset_id(3),filespace(3),hdf_error)
        call h5dget_space_f(dset_id(4),filespace(4),hdf_error)
        call h5dget_space_f(dset_id(5),filespace(5),hdf_error)
        call h5dget_space_f(dset_id(6),filespace(6),hdf_error)
        call h5dget_space_f(dset_id(7),filespace(7),hdf_error)
        call h5dget_space_f(dset_id(8),filespace(8),hdf_error)
        call h5dget_space_f(dset_id(9),filespace(9),hdf_error)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)

        ! Read source indices
        call h5sselect_hyperslab_f(filespace(1),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dread_f(dset_id(1),H5T_NATIVE_INTEGER,indx,dims1,hdf_error,file_space_id=filespace(1),mem_space_id=memspace1,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            var(i)%src_idx = indx(i)
        end do
        ! Read particle lifetimes
        call h5sselect_hyperslab_f(filespace(2),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dread_f(dset_id(2),H5T_NATIVE_DOUBLE, sclr,dims1,hdf_error,file_space_id=filespace(2),mem_space_id=memspace1,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            var(i)%lpp_lft = sclr(i)
        end do
        ! Read particle diameters
        call h5sselect_hyperslab_f(filespace(3),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dread_f(dset_id(3),H5T_NATIVE_DOUBLE, sclr,dims1,hdf_error,file_space_id=filespace(3),mem_space_id=memspace1,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            var(i)%lpp_dia = sclr(i)
        end do
        ! Read particle densities
        call h5sselect_hyperslab_f(filespace(4),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dread_f(dset_id(4),H5T_NATIVE_DOUBLE, sclr,dims1,hdf_error,file_space_id=filespace(4),mem_space_id=memspace1,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            var(i)%lpp_den = sclr(i)
        end do
        ! Read particle Reynolds number
        call h5sselect_hyperslab_f(filespace(5),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dread_f(dset_id(5),H5T_NATIVE_DOUBLE, sclr,dims1,hdf_error,file_space_id=filespace(5),mem_space_id=memspace1,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            var(i)%lpp_rey = sclr(i)
        end do
        ! Read particle positions
        call h5sselect_hyperslab_f(filespace(6),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dread_f(dset_id(6),H5T_NATIVE_DOUBLE, vctr,dims2,hdf_error,file_space_id=filespace(6),mem_space_id=memspace2,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            do j = 1,3
                var(i)%lpp_pos(j) = vctr(i,j)
            end do
        end do
        ! Read particle velocities
        call h5sselect_hyperslab_f(filespace(7),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dread_f(dset_id(7),H5T_NATIVE_DOUBLE, vctr,dims2,hdf_error,file_space_id=filespace(7),mem_space_id=memspace2,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            do j = 1,3
                var(i)%lpp_vel(j) = vctr(i,j)
            end do
        end do
        ! Read particle accelerations (previous)
        call h5sselect_hyperslab_f(filespace(8),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dread_f(dset_id(8),H5T_NATIVE_DOUBLE, vctr,dims2,hdf_error,file_space_id=filespace(8),mem_space_id=memspace2,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            do j = 1,3
                var(i)%acc_old(i) = vctr(i,j)
            end do
        end do
        ! Read particle accelerations (current)
        call h5sselect_hyperslab_f(filespace(9),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dread_f(dset_id(9),H5T_NATIVE_DOUBLE, vctr,dims2,hdf_error,file_space_id=filespace(9),mem_space_id=memspace2,xfer_prp=plist_id)
        do i = 1,en1-st1+1
            do j = 1,3
                var(i)%acc_now(i)  = vctr(i,j)
            end do
        end do
        
        call h5pclose_f(plist_id,hdf_error)

        do i = 1,9 
            call h5dclose_f(dset_id(i),hdf_error)
            call h5sclose_f(filespace(i),hdf_error)
        end do

        call h5sclose_f(memspace2,hdf_error)
        call h5sclose_f(memspace1,hdf_error)

    end if
    
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelReadParticle1D

subroutine HdfParallelWriteSource1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5
    use lagrangian_point_particle, only: particle_source

    implicit none

    character*200, intent(in)           :: filename,dsetname
    type(particle_source), intent(in)   :: var(1:en1-st1+1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,link_id,plist_id,dset_id(0:8)
    integer(HID_T)                      :: filespace(0:8)
    integer(HID_T)                      :: memspace1,memspace2
    integer(HSIZE_T)                    :: dims0(1),dims1(1),dims2(2)
    integer(HSIZE_T)                    :: data_count1(1),data_count2(2)  
    integer(HSIZE_T)                    :: data_offset1(1),data_offset2(2)
    integer                             :: i,j
    integer                             :: hdf_error
    logical                             :: exists
    integer                             :: indx(1:en1-st1+1)
    real                                :: sclr(1:en1-st1+1)
    real                                :: vctr(1:en1-st1+1,1:3)

    ! For count
    dims0(1) = 1

    ! For scalar data
    dims1(1) = n1e-n1s+1
    data_count1(1) = en1-st1+1
    data_offset1(1) = st1-n1s

    ! For vector data
    dims2(1) = n1e-n1s+1
    dims2(2) = 3

    data_count2(1) = en1-st1+1
    data_count2(2) = 3

    data_offset2(1) = st1-n1s
    data_offset2(2) = 0
    
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    inquire(file=filename,exist=exists)
    call MPI_ALLREDUCE(MPI_IN_PLACE,exists,1,MPI_LOGICAL,MPI_LOR,comm,mpi_ierr)
    if (exists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
    end if
    call h5pclose_f(plist_id,hdf_error)

    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_idx")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_sta")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_end")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_frq")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_dia")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_den")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_pos")
    call HdfDeleteExistingDataset(file_id,trim(dsetname)//"/src_vel")

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)

    ! Write the number of sources
    call h5screate_simple_f(1,dims0,filespace(0),hdf_error)
    call h5dcreate_f(file_id,trim(dsetname)//"/src_num",H5T_NATIVE_INTEGER,filespace(0),dset_id(0),hdf_error,lcpl_id=link_id)
    call h5dwrite_f(dset_id(0),H5T_NATIVE_INTEGER,(n1e-n1s+1),dims0,hdf_error)
    call h5dclose_f(dset_id(0),hdf_error)

    ! Write only if there is atleast one source (in total/global sense)
    if (n1e.ge.n1s) then

        call h5screate_simple_f(1,dims1,filespace(1),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(2),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(3),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(4),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(5),hdf_error)
        call h5screate_simple_f(1,dims1,filespace(6),hdf_error)
        call h5screate_simple_f(2,dims2,filespace(7),hdf_error)
        call h5screate_simple_f(2,dims2,filespace(8),hdf_error)

        call h5screate_simple_f(1,data_count1,memspace1,hdf_error)
        call h5screate_simple_f(2,data_count2,memspace2,hdf_error)

        call h5dcreate_f(file_id,trim(dsetname)//"/src_idx",H5T_NATIVE_INTEGER,filespace(1),dset_id(1),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/src_sta",H5T_NATIVE_DOUBLE, filespace(2),dset_id(2),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/src_end",H5T_NATIVE_DOUBLE, filespace(3),dset_id(3),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/src_frq",H5T_NATIVE_DOUBLE, filespace(4),dset_id(4),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/src_dia",H5T_NATIVE_DOUBLE, filespace(5),dset_id(5),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/src_den",H5T_NATIVE_DOUBLE, filespace(6),dset_id(6),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/src_pos",H5T_NATIVE_DOUBLE, filespace(7),dset_id(7),hdf_error,lcpl_id=link_id)
        call h5dcreate_f(file_id,trim(dsetname)//"/src_vel",H5T_NATIVE_DOUBLE, filespace(8),dset_id(8),hdf_error,lcpl_id=link_id)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)

        ! Write source indices
        do i = 1,en1-st1+1
            indx(i) = var(i)%src_idx
        end do
        call h5sselect_hyperslab_f(filespace(1),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(1),H5T_NATIVE_INTEGER,indx,data_count1,hdf_error,file_space_id=filespace(1),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write source injection start times
        do i = 1,en1-st1+1
            sclr(i) = var(i)%src_sta
        end do
        call h5sselect_hyperslab_f(filespace(2),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(2),H5T_NATIVE_DOUBLE, sclr,data_count1,hdf_error,file_space_id=filespace(2),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write source injection end times
        do i = 1,en1-st1+1
            sclr(i) = var(i)%src_end
        end do
        call h5sselect_hyperslab_f(filespace(3),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(3),H5T_NATIVE_DOUBLE, sclr,data_count1,hdf_error,file_space_id=filespace(3),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write source injection frequency
        do i = 1,en1-st1+1
            sclr(i) = var(i)%src_frq
        end do
        call h5sselect_hyperslab_f(filespace(4),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(4),H5T_NATIVE_DOUBLE, sclr,data_count1,hdf_error,file_space_id=filespace(4),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write source diameters
        do i = 1,en1-st1+1
            sclr(i) = var(i)%src_dia
        end do
        call h5sselect_hyperslab_f(filespace(5),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(5),H5T_NATIVE_DOUBLE, sclr,data_count1,hdf_error,file_space_id=filespace(5),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write source densities
        do i = 1,en1-st1+1
            sclr(i) = var(i)%src_den
        end do
        call h5sselect_hyperslab_f(filespace(6),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
        call h5dwrite_f(dset_id(6),H5T_NATIVE_DOUBLE, sclr,data_count2,hdf_error,file_space_id=filespace(6),mem_space_id=memspace1,xfer_prp=plist_id)
        ! Write source positions
        do i = 1,en1-st1+1
            do j = 1,3
                vctr(i,j) = var(i)%src_pos(j)
            end do
        end do
        call h5sselect_hyperslab_f(filespace(7),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dwrite_f(dset_id(7),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(7),mem_space_id=memspace2,xfer_prp=plist_id)
        ! Write source velocities
        do i = 1,en1-st1+1
            do j = 1,3
                vctr(i,j) = var(i)%src_vel(j)
            end do
        end do
        call h5sselect_hyperslab_f(filespace(8),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
        call h5dwrite_f(dset_id(8),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(8),mem_space_id=memspace2,xfer_prp=plist_id)
        
        call h5pclose_f(plist_id,hdf_error)

        do i = 1,8
            call h5dclose_f(dset_id(i),hdf_error)
            call h5sclose_f(filespace(i),hdf_error)
        end do

        call h5sclose_f(memspace2,hdf_error)
        call h5sclose_f(memspace1,hdf_error)

    end if

    call h5pclose_f(link_id,hdf_error)

    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelWriteSource1D

subroutine HdfParallelCreateParticleExit1D(filename,dsetname,comm)

    use mpih
    use hdf5
    use lagrangian_point_particle, only: particle_exit

    implicit none

    character*200, intent(in)           :: filename,dsetname
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,plist_id,link_id,chnk_id1,chnk_id2,dset_id(0:6)
    integer(HID_T)                      :: filespace(0:6)
    integer(HID_T)                      :: H5T_DERIVED_STRING
    integer(HSIZE_T)                    :: strlength
    integer(HSIZE_T)                    :: dims0(1),dims1(1),dims2(2)
    integer(HSIZE_T)                    :: maxdims1(1),maxdims2(2)
    integer(HSIZE_T)                    :: chunkdims1(1),chunkdims2(2)
    integer                             :: datasize
    integer                             :: i
    integer                             :: hdf_error
    logical                             :: exists

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    inquire(file=filename,exist=exists)
    call MPI_ALLREDUCE(MPI_IN_PLACE,exists,1,MPI_LOGICAL,MPI_LOR,comm,mpi_ierr)
    if (exists) then
        call h5pclose_f(plist_id,hdf_error)
    else

        ! Create the derived datatype for character/string
        strlength = 1
        call h5tcopy_f(H5T_FORTRAN_S1, H5T_DERIVED_STRING, mpi_ierr)
        call h5tset_size_f(H5T_DERIVED_STRING, strlength, mpi_ierr)

        ! For initialization of the write start index
        datasize = 0

        ! Dimension of exit count
        dims0(1) = 1

        ! For initialization of scalar data
        dims1(1) = 0
        maxdims1(1) = H5S_UNLIMITED_F
        chunkdims1(1) = 1024

        ! For initialization of vector data
        dims2(1) = 0
        dims2(2) = 3

        maxdims2(1) = H5S_UNLIMITED_F
        maxdims2(2) = 3
        
        chunkdims2(1) = 1024
        chunkdims2(2) = 3
    
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error,access_prp=plist_id)
        
        call h5pclose_f(plist_id,hdf_error)
        
        call h5screate_simple_f(1,dims0,filespace(0),hdf_error)
        
        call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
        call h5pset_create_inter_group_f(link_id,1,hdf_error)
        call h5dcreate_f(file_id,trim(dsetname)//"/pex_num",H5T_NATIVE_INTEGER,filespace(0),dset_id(0),hdf_error,lcpl_id=link_id)
        call h5pclose_f(link_id,hdf_error)
        
        call h5dwrite_f(dset_id(0),H5T_NATIVE_INTEGER,datasize,dims0,hdf_error)

        call h5screate_simple_f(1,dims1,filespace(1),hdf_error,maxdims1)
        call h5screate_simple_f(1,dims1,filespace(2),hdf_error,maxdims1)
        call h5screate_simple_f(1,dims1,filespace(3),hdf_error,maxdims1)
        call h5screate_simple_f(1,dims1,filespace(4),hdf_error,maxdims1)
        call h5screate_simple_f(2,dims2,filespace(5),hdf_error,maxdims2)
        call h5screate_simple_f(2,dims2,filespace(6),hdf_error,maxdims2)

        call h5pcreate_f(H5P_DATASET_CREATE_F,chnk_id1,hdf_error)
        call h5pset_chunk_f(chnk_id1,1,chunkdims1,hdf_error)

        call h5pcreate_f(H5P_DATASET_CREATE_F,chnk_id2,hdf_error)
        call h5pset_chunk_f(chnk_id2,2,chunkdims2,hdf_error)

        call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
        call h5pset_create_inter_group_f(link_id,1,hdf_error)

        call h5dcreate_f(file_id,trim(dsetname)//"/src_idx",H5T_NATIVE_INTEGER,filespace(1),dset_id(1),hdf_error,lcpl_id=link_id,dcpl_id=chnk_id1)
        call h5dcreate_f(file_id,trim(dsetname)//"/pex_pln",H5T_DERIVED_STRING,filespace(2),dset_id(2),hdf_error,lcpl_id=link_id,dcpl_id=chnk_id1)
        call h5dcreate_f(file_id,trim(dsetname)//"/pex_lft",H5T_NATIVE_DOUBLE ,filespace(3),dset_id(3),hdf_error,lcpl_id=link_id,dcpl_id=chnk_id1)
        call h5dcreate_f(file_id,trim(dsetname)//"/pex_eft",H5T_NATIVE_DOUBLE ,filespace(4),dset_id(4),hdf_error,lcpl_id=link_id,dcpl_id=chnk_id1)
        call h5dcreate_f(file_id,trim(dsetname)//"/pex_pos",H5T_NATIVE_DOUBLE ,filespace(5),dset_id(5),hdf_error,lcpl_id=link_id,dcpl_id=chnk_id2)
        call h5dcreate_f(file_id,trim(dsetname)//"/pex_vel",H5T_NATIVE_DOUBLE ,filespace(6),dset_id(6),hdf_error,lcpl_id=link_id,dcpl_id=chnk_id2)

        call h5pclose_f(link_id,hdf_error)

        call h5pclose_f(chnk_id1,hdf_error)
        call h5pclose_f(chnk_id2,hdf_error)

        do i = 0,6
            call h5dclose_f(dset_id(i),hdf_error)
            call h5sclose_f(filespace(i),hdf_error)
        end do

        call h5fclose_f(file_id,hdf_error)

    end if

end subroutine HdfParallelCreateParticleExit1D

subroutine HdfParallelAppendParticleExit1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5
    use lagrangian_point_particle, only: particle_exit

    implicit none

    character*200, intent(in)           :: filename,dsetname
    type(particle_exit), intent(in)     :: var(1:en1-st1+1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,plist_id,dset_id(0:6)
    integer(HID_T)                      :: filespace(0:6)
    integer(HID_T)                      :: memspace1,memspace2
    integer(HID_T)                      :: H5T_DERIVED_STRING
    integer(HSIZE_T)                    :: strlength
    integer(HSIZE_T)                    :: dims0(1),dims1(1),dims2(2)
    integer(HSIZE_T)                    :: data_count1(1),data_count2(2)  
    integer(HSIZE_T)                    :: data_offset1(1),data_offset2(2)
    integer                             :: datasize
    integer                             :: i,j
    integer                             :: hdf_error
    logical                             :: exists
    character                           :: chrr(1:en1-st1+1)
    integer                             :: indx(1:en1-st1+1)
    real                                :: sclr(1:en1-st1+1)
    real                                :: vctr(1:en1-st1+1,1:3)

    ! Create the derived datatype for character/string
    strlength = 1
    call h5tcopy_f(H5T_FORTRAN_S1, H5T_DERIVED_STRING, mpi_ierr)
    call h5tset_size_f(H5T_DERIVED_STRING, strlength, mpi_ierr)

    ! Dimension of exit count
    dims0(1) = 1

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    
    call h5dopen_f(file_id,trim(dsetname)//"/pex_num",dset_id(0),hdf_error)
    call h5dget_space_f(dset_id(0),filespace(0),hdf_error)
    call h5dread_f(dset_id(0),H5T_NATIVE_INTEGER,datasize,dims0,hdf_error)

    call h5dopen_f(file_id,trim(dsetname)//"/src_idx",dset_id(1),hdf_error)
    call h5dopen_f(file_id,trim(dsetname)//"/pex_pln",dset_id(2),hdf_error)
    call h5dopen_f(file_id,trim(dsetname)//"/pex_lft",dset_id(3),hdf_error)
    call h5dopen_f(file_id,trim(dsetname)//"/pex_eft",dset_id(4),hdf_error)
    call h5dopen_f(file_id,trim(dsetname)//"/pex_pos",dset_id(5),hdf_error)
    call h5dopen_f(file_id,trim(dsetname)//"/pex_vel",dset_id(6),hdf_error)

    ! For scalar data
    dims1(1) = n1e-n1s+1 + datasize
    data_count1(1) = en1-st1+1
    data_offset1(1) = st1-n1s + datasize

    ! For vector data
    dims2(1) = n1e-n1s+1 + datasize
    dims2(2) = 3

    data_count2(1) = en1-st1+1
    data_count2(2) = 3

    data_offset2(1) = st1-n1s + datasize
    data_offset2(2) = 0

    ! Update data start for next append 
    datasize = n1e-n1s+1 + datasize
    
    call h5dwrite_f(dset_id(0),H5T_NATIVE_INTEGER,datasize,dims0,hdf_error)

    ! Extend datasets to accommodate new exit events to be added 
    call h5dset_extent_f(dset_id(1),dims1,hdf_error)
    call h5dset_extent_f(dset_id(2),dims1,hdf_error)
    call h5dset_extent_f(dset_id(3),dims1,hdf_error)
    call h5dset_extent_f(dset_id(4),dims1,hdf_error)
    call h5dset_extent_f(dset_id(5),dims2,hdf_error)
    call h5dset_extent_f(dset_id(6),dims2,hdf_error)

    ! Get new filespace
    do i = 1,6
        call h5dget_space_f(dset_id(i),filespace(i),hdf_error)
    end do

    call h5screate_simple_f(1,data_count1,memspace1,hdf_error)
    call h5screate_simple_f(2,data_count2,memspace2,hdf_error)

    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)

    ! Append source indices corresponding to exit events
    do i = 1,en1-st1+1
        indx(i) = var(i)%src_idx
    end do
    call h5sselect_hyperslab_f(filespace(1),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
    call h5dwrite_f(dset_id(1),H5T_NATIVE_INTEGER,indx,data_count1,hdf_error,file_space_id=filespace(1),mem_space_id=memspace1,xfer_prp=plist_id)
    ! Append exit plane
    do i = 1,en1-st1+1
        chrr(i) = var(i)%pex_pln
    end do
    call h5sselect_hyperslab_f(filespace(2),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
    call h5dwrite_f(dset_id(2),H5T_DERIVED_STRING,chrr,data_count1,hdf_error,file_space_id=filespace(2),mem_space_id=memspace1,xfer_prp=plist_id)
    ! Append particle lifetime at exit
    do i = 1,en1-st1+1
        sclr(i) = var(i)%pex_lft
    end do
    call h5sselect_hyperslab_f(filespace(3),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
    call h5dwrite_f(dset_id(3),H5T_NATIVE_DOUBLE ,sclr,data_count1,hdf_error,file_space_id=filespace(3),mem_space_id=memspace1,xfer_prp=plist_id)
    ! Append flow time at exit
    do i = 1,en1-st1+1
        sclr(i) = var(i)%pex_eft
    end do
    call h5sselect_hyperslab_f(filespace(4),H5S_SELECT_SET_F,data_offset1,data_count1,hdf_error)
    call h5dwrite_f(dset_id(4),H5T_NATIVE_DOUBLE ,sclr,data_count1,hdf_error,file_space_id=filespace(4),mem_space_id=memspace1,xfer_prp=plist_id)
    ! Append particle positions at exit
    do i = 1,en1-st1+1
        do j = 1,3
            vctr(i,j) = var(i)%pex_pos(j)
        end do
    end do
    call h5sselect_hyperslab_f(filespace(5),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
    call h5dwrite_f(dset_id(5),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(5),mem_space_id=memspace2,xfer_prp=plist_id)
    ! Append particle velocities at exit
    do i = 1,en1-st1+1
        do j = 1,3
            vctr(i,j) = var(i)%pex_vel(j)
        end do
    end do
    call h5sselect_hyperslab_f(filespace(6),H5S_SELECT_SET_F,data_offset2,data_count2,hdf_error)
    call h5dwrite_f(dset_id(6),H5T_NATIVE_DOUBLE, vctr,data_count2,hdf_error,file_space_id=filespace(6),mem_space_id=memspace2,xfer_prp=plist_id)
    
    do i = 0,6
        call h5dclose_f(dset_id(i),hdf_error)
        call h5sclose_f(filespace(i),hdf_error)
    end do

    call h5sclose_f(memspace2,hdf_error)
    call h5sclose_f(memspace1,hdf_error)

    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelAppendParticleExit1D