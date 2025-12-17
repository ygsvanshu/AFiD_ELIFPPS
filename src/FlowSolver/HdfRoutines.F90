!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         ! 
!    FILE: HdfRoutines.F90                                !
!    CONTAINS: subroutine hdf_read_serial_1d              !
!                                                         ! 
!    PURPOSE: I/O routines.                               !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine HdfStart
    use hdf5
    implicit none
    integer :: hdf_error
    call h5open_f(hdf_error)
end subroutine HdfStart
  
subroutine HdfClose
    use hdf5
    implicit none
    integer :: hdf_error
    call h5close_f(hdf_error)
end subroutine HdfClose

subroutine HdfClean(filename)
    use hdf5
    implicit none
    integer         :: hdf_error
    character*200   :: filename
    integer(HID_T)  :: file_id
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    call h5fclose_f(file_id,hdf_error)
end subroutine HdfClean

subroutine HdfDeleteExistingDataset(file_id,dsetname)

    use hdf5

    implicit none
    
    character*200, intent(in)   :: dsetname
    integer(HID_T), intent(in)  :: file_id
    integer                     :: slashpos,checkpos,dsetnlen
    integer                     :: hdf_error
    logical                     :: exists

    dsetnlen = len(trim(dsetname))
    exists   = .true.
    slashpos = 1
    checkpos = 1
    if (index(dsetname,"/").eq.1) checkpos = 2
    do while ((checkpos.lt.dsetnlen).and.(exists))
        slashpos = index(dsetname(checkpos:),"/")
        if (slashpos.eq.0) slashpos = dsetnlen - checkpos + 2
        call h5lexists_f(file_id,dsetname(:slashpos+checkpos-2),exists,hdf_error)
        checkpos = checkpos + slashpos
    end do
    if (exists) call h5ldelete_f(file_id,dsetname,hdf_error)

end subroutine HdfDeleteExistingDataset

!==================================================================================================

subroutine HdfSerialWriteRealScalar(filename,dsetname,var)

    use hdf5

    implicit none

    character*200, intent(in)   :: filename,dsetname
    real, intent(in)            :: var

    integer(HID_T)              :: file_id,link_id,dset_id
    integer(HID_T)              :: filespace
    integer(HSIZE_T)            :: dims(1)
    integer                     :: hdf_error
    logical                     :: exists

    dims(1)  = 1

    inquire(file=filename,exist=exists)
    if (exists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    end if

    call HdfDeleteExistingDataset(file_id,dsetname)

    call h5screate_simple_f(1,dims,filespace,hdf_error)

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,hdf_error,lcpl_id=link_id)
    call h5pclose_f(link_id,hdf_error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var,dims,hdf_error)
    call h5dclose_f(dset_id,hdf_error)

    call h5sclose_f(filespace,hdf_error)

    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialWriteRealScalar

subroutine HdfSerialWriteIntScalar(filename,dsetname,var)

    use hdf5

    implicit none

    character*200, intent(in)   :: filename,dsetname
    integer, intent(in)         :: var

    integer(HID_T)              :: file_id,link_id,dset_id
    integer(HID_T)              :: filespace 
    integer(HSIZE_T)            :: dims(1)
    integer                     :: hdf_error
    logical                     :: exists

    dims(1)  = 1

    inquire(file=filename,exist=exists)
    if (exists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    end if

    call HdfDeleteExistingDataset(file_id,dsetname)

    call h5screate_simple_f(1,dims,filespace,hdf_error)

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_INTEGER,filespace,dset_id,hdf_error,lcpl_id=link_id)
    call h5pclose_f(link_id,hdf_error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,var,dims,hdf_error)
    call h5dclose_f(dset_id,hdf_error)

    call h5sclose_f(filespace,hdf_error)

    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialWriteIntScalar

subroutine HdfSerialWriteReal1D(filename,dsetname,var,st1,en1)

    use hdf5

    implicit none

    character*200, intent(in)   :: filename,dsetname
    real, intent(in)            :: var(st1:en1)
    integer, intent(in)         :: st1,en1

    integer(HID_T)              :: file_id,link_id,dset_id
    integer(HID_T)              :: filespace
    integer(HSIZE_T)            :: dims(1)
    integer                     :: hdf_error
    logical                     :: exists

    dims(1)  = en1-st1+1

    inquire(file=filename,exist=exists)
    if (exists) then
        call h5fopen_f(filename,H5F_ACC_RDWR_F,file_id,hdf_error)
    else
        call h5fcreate_f(filename,H5F_ACC_TRUNC_F,file_id,hdf_error)
    end if

    call HdfDeleteExistingDataset(file_id,dsetname)

    call h5screate_simple_f(1,dims,filespace,hdf_error)

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,hdf_error,lcpl_id=link_id)
    call h5pclose_f(link_id,hdf_error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var(st1:en1),dims,hdf_error)
    call h5dclose_f(dset_id,hdf_error)

    call h5sclose_f(filespace,hdf_error)

    call h5fclose_f(file_id,hdf_error)
    
end subroutine HdfSerialWriteReal1D

subroutine HdfParallelWriteInt1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5

    implicit none

    character*200, intent(in)           :: filename,dsetname
    integer, intent(in)                 :: var(st1:en1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,link_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(1)
    integer(HSIZE_T)                    :: data_count(1)  
    integer(HSIZE_T)                    :: data_offset(1)
    integer                             :: hdf_error
    logical                             :: exists

    dims(1) = n1e-n1s+1

    data_count(1) = en1-st1+1

    data_offset(1) = st1-n1s

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

    call HdfDeleteExistingDataset(file_id,dsetname)

    call h5screate_simple_f(1,dims,filespace,hdf_error)
    call h5screate_simple_f(1,data_count,memspace,hdf_error)

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_INTEGER,filespace,dset_id,hdf_error,lcpl_id=link_id)
    call h5pclose_f(link_id,hdf_error)

    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dwrite_f(dset_id,H5T_NATIVE_INTEGER,var(st1:en1),data_count,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    
    call h5dclose_f(dset_id,hdf_error)
    
    call h5sclose_f(memspace,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfParallelWriteInt1D

subroutine HdfParallelWriteReal1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5

    implicit none

    character*200, intent(in)           :: filename,dsetname
    real, intent(in)                    :: var(st1:en1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,link_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(1)
    integer(HSIZE_T)                    :: data_count(1)  
    integer(HSIZE_T)                    :: data_offset(1)
    integer                             :: hdf_error
    logical                             :: exists

    dims(1) = n1e-n1s+1

    data_count(1) = en1-st1+1

    data_offset(1) = st1-n1s

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

    call HdfDeleteExistingDataset(file_id,dsetname)

    call h5screate_simple_f(1,dims,filespace,hdf_error)
    call h5screate_simple_f(1,data_count,memspace,hdf_error)

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,hdf_error,lcpl_id=link_id)
    call h5pclose_f(link_id,hdf_error)

    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var(st1:en1),data_count,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    
    call h5dclose_f(dset_id,hdf_error)
    
    call h5sclose_f(memspace,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfParallelWriteReal1D

subroutine HdfParallelWriteReal2D(filename,dsetname,var,st1,en1,st2,en2,n1s,n1e,n2s,n2e,comm)

    use mpih
    use hdf5

    implicit none

    character*200, intent(in)           :: filename,dsetname
    real, intent(in)                    :: var(st1:en1,st2:en2)
    integer, intent(in)                 :: st1,en1,st2,en2
    integer, intent(in)                 :: n1s,n1e,n2s,n2e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,link_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(2)
    integer(HSIZE_T)                    :: data_count(2)  
    integer(HSIZE_T)                    :: data_offset(2)
    integer                             :: hdf_error
    logical                             :: exists

    dims(1) = n1e-n1s+1
    dims(2) = n2e-n2s+1

    data_count(1) = en1-st1+1
    data_count(2) = en2-st2+1

    data_offset(1) = st1-n1s
    data_offset(2) = st2-n2s

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

    call HdfDeleteExistingDataset(file_id,dsetname)

    call h5screate_simple_f(2,dims,filespace,hdf_error)
    call h5screate_simple_f(2,data_count,memspace,hdf_error)

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,hdf_error,lcpl_id=link_id)
    call h5pclose_f(link_id,hdf_error)

    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var(st1:en1,st2:en2),data_count,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)
    
    call h5dclose_f(dset_id,hdf_error)
    
    call h5sclose_f(memspace,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfParallelWriteReal2D

subroutine HdfParallelWriteReal3D(filename,dsetname,var,st1,en1,st2,en2,st3,en3,n1s,n1e,n2s,n2e,n3s,n3e,comm)

    use mpih
    use hdf5
      
    implicit none

    character*200, intent(in)           :: filename,dsetname
    real, intent(in)                    :: var(st1:en1,st2:en2,st3:en3)
    integer, intent(in)                 :: st1,en1,st2,en2,st3,en3
    integer, intent(in)                 :: n1s,n1e,n2s,n2e,n3s,n3e
    integer, intent(in)                 :: comm
    
    integer(HID_T)                      :: file_id,link_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(3)
    integer(HSIZE_T)                    :: data_count(3)
    integer(HSIZE_T)                    :: data_offset(3)
    integer                             :: hdf_error
    logical                             :: exists

    dims(1) = n1e-n1s+1
    dims(2) = n2e-n2s+1
    dims(3) = n3e-n3s+1

    data_count(1) = en1-st1+1
    data_count(2) = en2-st2+1
    data_count(3) = en3-st3+1

    data_offset(1) = st1-n1s
    data_offset(2) = st2-n2s
    data_offset(3) = st3-n3s

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

    call HdfDeleteExistingDataset(file_id,dsetname)

    call h5screate_simple_f(3,dims,filespace,hdf_error)
    call h5screate_simple_f(3,data_count,memspace,hdf_error) 

    call h5pcreate_f(H5P_LINK_CREATE_F,link_id,hdf_error)
    call h5pset_create_inter_group_f(link_id,1,hdf_error)
    call h5dcreate_f(file_id,dsetname,H5T_NATIVE_DOUBLE,filespace,dset_id,hdf_error,lcpl_id=link_id)
    call h5pclose_f(link_id,hdf_error)
    
    call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

    call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
    call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
    
    call h5dwrite_f(dset_id,H5T_NATIVE_DOUBLE,var(st1:en1,st2:en2,st3:en3),dims,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
    
    call h5pclose_f(plist_id,hdf_error)

    call h5dclose_f(dset_id,hdf_error)

    call h5sclose_f(memspace,hdf_error)
    call h5sclose_f(filespace,hdf_error)
    
    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelWriteReal3D

!==================================================================================================

subroutine HdfSerialReadRealScalar(filename,dsetname,var)

    use hdf5

    implicit none

    character*200, intent(in)   :: filename,dsetname
    real, intent(out)           :: var

    integer(HID_T)              :: file_id,dset_id
    integer(HSIZE_T)            :: dims(1)
    integer                     :: hdf_error
    logical                     :: exists

    var = 0.0

    dims(1)  = 1

    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error)
    call h5lexists_f(file_id,dsetname,exists,hdf_error)
    if (exists) then
        call h5dopen_f(file_id,dsetname,dset_id,hdf_error)
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,var,dims,hdf_error)
        call h5dclose_f(dset_id,hdf_error)
    end if
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialReadRealScalar

subroutine HdfSerialReadIntScalar(filename,dsetname,var)

    use hdf5

    implicit none

    character*200, intent(in)   :: filename,dsetname
    integer, intent(out)        :: var
    
    integer(HID_T)              :: file_id,dset_id
    integer(HSIZE_T)            :: dims(1)
    integer                     :: hdf_error
    logical                     :: exists
    
    var = 0

    dims(1)  = 1

    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error)
    call h5lexists_f(file_id,dsetname,exists,hdf_error)
    if (exists) then
        call h5dopen_f(file_id,dsetname,dset_id,hdf_error)
        call h5dread_f(dset_id,H5T_NATIVE_INTEGER,var,dims,hdf_error)
        call h5dclose_f(dset_id,hdf_error)
    end if
    call h5fclose_f(file_id,hdf_error)
    
end subroutine HdfSerialReadIntScalar

subroutine HdfSerialReadReal1D(filename,dsetname,var,st1,en1)

    use hdf5

    implicit none

    character*200, intent(in)   :: filename,dsetname
    real, intent(out)           :: var(st1:en1)
    integer, intent(in)         :: st1,en1
    
    integer(HID_T)              :: file_id,dset_id
    integer(HSIZE_T)            :: dims(1)
    integer                     :: hdf_error
    logical                     :: exists

    var(:) = 0.0
    
    dims(1)  = en1-st1+1

    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error)
    call h5lexists_f(file_id,dsetname,exists,hdf_error)
    if (exists) then
        call h5dopen_f(file_id,dsetname,dset_id,hdf_error)
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,var,dims,hdf_error)
        call h5dclose_f(dset_id,hdf_error)
    end if
    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfSerialReadReal1D

subroutine HdfParallelReadInt1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5

    implicit none

    character*200, intent(in)           :: filename,dsetname
    integer, intent(out)                :: var(st1:en1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(1)
    integer(HSIZE_T)                    :: data_count(1)
    integer(HSIZE_T)                    :: data_offset(1) 
    integer                             :: hdf_error
    logical                             :: exists

    var(st1:en1) = 0

    !RO   Sort out MPI definitions

    dims(1) = n1e-n1s+1

    data_count(1) = en1-st1+1

    data_offset(1) = st1-n1s

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)

    call h5lexists_f(file_id,dsetname,exists,hdf_error)
    if (exists) then

        call h5dopen_f(file_id,dsetname,dset_id,hdf_error)
        call h5dget_space_f(dset_id,filespace,hdf_error)

        call h5screate_simple_f(1,data_count,memspace,hdf_error) 
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
        call h5dread_f(dset_id,H5T_NATIVE_INTEGER,var(st1:en1),dims,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
        call h5pclose_f(plist_id,hdf_error)
        call h5dclose_f(dset_id,hdf_error)
        
        call h5sclose_f(memspace,hdf_error)
        call h5sclose_f(filespace,hdf_error)

    end if

    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfParallelReadInt1D

subroutine HdfParallelReadReal1D(filename,dsetname,var,st1,en1,n1s,n1e,comm)

    use mpih
    use hdf5

    implicit none

    character*200, intent(in)           :: filename,dsetname
    real, intent(out)                   :: var(st1:en1)
    integer, intent(in)                 :: st1,en1
    integer, intent(in)                 :: n1s,n1e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(1)
    integer(HSIZE_T)                    :: data_count(1)
    integer(HSIZE_T)                    :: data_offset(1) 
    integer                             :: hdf_error
    logical                             :: exists

    var(st1:en1) = 0.0

    !RO   Sort out MPI definitions

    dims(1) = n1e-n1s+1

    data_count(1) = en1-st1+1

    data_offset(1) = st1-n1s

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)

    call h5lexists_f(file_id,dsetname,exists,hdf_error)
    if (exists) then

        call h5dopen_f(file_id,dsetname,dset_id,hdf_error)
        call h5dget_space_f(dset_id,filespace,hdf_error)

        call h5screate_simple_f(1,data_count,memspace,hdf_error) 
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,var(st1:en1),dims,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
        call h5pclose_f(plist_id,hdf_error)
        call h5dclose_f(dset_id,hdf_error)
        
        call h5sclose_f(memspace,hdf_error)
        call h5sclose_f(filespace,hdf_error)

    end if

    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfParallelReadReal1D

subroutine HdfParallelReadReal2D(filename,dsetname,var,st1,en1,st2,en2,n1s,n1e,n2s,n2e,comm)

    use mpih
    use hdf5

    implicit none

    character*200, intent(in)           :: filename,dsetname
    real, intent(out)                   :: var(st1:en1,st2:en2)
    integer, intent(in)                 :: st1,en1,st2,en2
    integer, intent(in)                 :: n1s,n1e,n2s,n2e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(2)
    integer(HSIZE_T)                    :: data_count(2)
    integer(HSIZE_T)                    :: data_offset(2) 
    integer                             :: hdf_error
    logical                             :: exists

    var(st1:en1,st2:en2) = 0.0

    !RO   Sort out MPI definitions

    dims(1) = n1e-n1s+1
    dims(2) = n2e-n2s+1

    data_count(1) = en1-st1+1
    data_count(2) = en2-st2+1

    data_offset(1) = st1-n1s
    data_offset(2) = st2-n2s

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)

    call h5lexists_f(file_id,dsetname,exists,hdf_error)
    if (exists) then

        call h5dopen_f(file_id,dsetname,dset_id,hdf_error)
        call h5dget_space_f(dset_id,filespace,hdf_error)

        call h5screate_simple_f(2,data_count,memspace,hdf_error) 
        call h5sselect_hyperslab_f(filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,var(st1:en1,st2:en2),dims,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
        call h5pclose_f(plist_id,hdf_error)
        call h5dclose_f(dset_id,hdf_error)
        
        call h5sclose_f(memspace,hdf_error)
        call h5sclose_f(filespace,hdf_error)

    end if

    call h5fclose_f(file_id,hdf_error)
      
end subroutine HdfParallelReadReal2D

subroutine HdfParallelReadReal3D(filename,dsetname,var,st1,en1,st2,en2,st3,en3,n1s,n1e,n2s,n2e,n3s,n3e,comm)
    
    use mpih
    use hdf5
    
    implicit none

    character*200, intent(in)           :: filename,dsetname
    real, intent(out)                   :: var(st1:en1,st2:en2,st3:en3)
    integer, intent(in)                 :: st1,en1,st2,en2,st3,en3
    integer, intent(in)                 :: n1s,n1e,n2s,n2e,n3s,n3e
    integer, intent(in)                 :: comm

    integer(HID_T)                      :: file_id,dset_id,plist_id
    integer(HID_T)                      :: memspace,filespace
    integer(HSIZE_T)                    :: dims(3)
    integer(HSIZE_T)                    :: data_count(3)
    integer(HSIZE_T)                    :: data_offset(3)
    integer                             :: hdf_error
    logical                             :: exists

    var(st1:en1,st2:en2,st3:en3) = 0.0
    
    dims(1) = n1e-n1s+1
    dims(2) = n2e-n2s+1
    dims(3) = n3e-n3s+1

    data_count(1) = en1-st1+1
    data_count(2) = en2-st2+1
    data_count(3) = en3-st3+1

    data_offset(1) = st1-n1s
    data_offset(2) = st2-n2s
    data_offset(3) = st3-n3s

    call h5pcreate_f(H5P_FILE_ACCESS_F,plist_id,hdf_error)
    call h5pset_fapl_mpio_f(plist_id,comm,MPI_INFO_NULL,hdf_error)
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdf_error,access_prp=plist_id)
    call h5pclose_f(plist_id,hdf_error)

    call h5lexists_f(file_id,dsetname,exists,hdf_error)
    if (exists) then
    
        call h5dopen_f(file_id,dsetname,dset_id,hdf_error)
        call h5dget_space_f(dset_id,filespace,hdf_error)

        call h5screate_simple_f(3,data_count,memspace,hdf_error) 
        call h5sselect_hyperslab_f (filespace,H5S_SELECT_SET_F,data_offset,data_count,hdf_error)

        call h5pcreate_f(H5P_DATASET_XFER_F,plist_id,hdf_error) 
        call h5pset_dxpl_mpio_f(plist_id,H5FD_MPIO_COLLECTIVE_F,hdf_error)
        call h5dread_f(dset_id,H5T_NATIVE_DOUBLE,var(st1:en1,st2:en2,st3:en3),dims,hdf_error,file_space_id=filespace,mem_space_id=memspace,xfer_prp=plist_id)
        call h5pclose_f(plist_id,hdf_error)
        call h5dclose_f(dset_id,hdf_error)
        
        call h5sclose_f(memspace,hdf_error)
        call h5sclose_f(filespace,hdf_error)

    end if

    call h5fclose_f(file_id,hdf_error)

end subroutine HdfParallelReadReal3D