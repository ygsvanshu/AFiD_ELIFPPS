!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                         !
!    FILE: AuxiliaryRoutines.F90                          !
!    CONTAINS: subroutines Allocate*,Destroy*             !
!                                                         !
!    PURPOSE: Auxiliary routines used for memory allocs   !
!    and memory freeing                                   !
!                                                         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module AuxiliaryRoutines

contains

    subroutine AllocateInteger1DArray(var, st1, en1)

        use decomp_2d

        implicit none

        integer, intent(in)                                 :: st1, en1
        integer, allocatable, dimension(:), intent(inout)   :: var
        integer                                             :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate (var(st1:en1), stat=alloc_stat)
        if (alloc_stat /= 0) then
            errorcode = 8
            call decomp_2d_abort(errorcode, 'Memory allocation failed when creating new arrays')
        end if
        var = 0

        return

    end subroutine AllocateInteger1DArray

    subroutine AllocateReal1DArray(var, st1, en1)

        use decomp_2d

        implicit none

        integer, intent(in)                             :: st1, en1
        real, allocatable, dimension(:), intent(inout)  :: var
        integer                                         :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate (var(st1:en1), stat=alloc_stat)
        if (alloc_stat /= 0) then
            errorcode = 8
            call decomp_2d_abort(errorcode, 'Memory allocation failed when creating new arrays')
        end if
        var = 0.0

        return

    end subroutine AllocateReal1DArray

    subroutine AllocateLogical2DArray(var, st1, en1, st2, en2)

        use decomp_2d

        implicit none

        integer, intent(in)                                     :: st1, en1, st2, en2
        logical, allocatable, dimension(:, :), intent(inout)    :: var
        integer                                                 :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate (var(st1:en1, st2:en2), stat=alloc_stat)
        if (alloc_stat /= 0) then
            errorcode = 8
            call decomp_2d_abort(errorcode, 'Memory allocation failed when creating new arrays')
        end if
        var = .false.

        return

    end subroutine AllocateLogical2DArray

    subroutine AllocateInteger2DArray(var, st1, en1, st2, en2)

        use decomp_2d

        implicit none

        integer, intent(in)                                     :: st1, en1, st2, en2
        integer, allocatable, dimension(:, :), intent(inout)    :: var
        integer                                                 :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate (var(st1:en1, st2:en2), stat=alloc_stat)
        if (alloc_stat /= 0) then
            errorcode = 8
            call decomp_2d_abort(errorcode, 'Memory allocation failed when creating new arrays')
        end if
        var = 0

        return

    end subroutine AllocateInteger2DArray

    subroutine AllocateReal2DArray(var, st1, en1, st2, en2)

        use decomp_2d

        implicit none

        integer, intent(in)                                 :: st1, en1, st2, en2
        real, allocatable, dimension(:, :), intent(inout)   :: var
        integer                                             :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate (var(st1:en1, st2:en2), stat=alloc_stat)
        if (alloc_stat /= 0) then
            errorcode = 8
            call decomp_2d_abort(errorcode, 'Memory allocation failed when creating new arrays')
        end if
        var = 0.0

        return

    end subroutine AllocateReal2DArray

    subroutine AllocateInteger3DArray(var, st1, en1, st2, en2, st3, en3)

        use decomp_2d

        implicit none

        integer, intent(in)                                     :: st1, en1, st2, en2, st3, en3
        integer, allocatable, dimension(:, :, :), intent(inout) :: var
        integer                                                 :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate (var(st1:en1, st2:en2, st3:en3), stat=alloc_stat)
        if (alloc_stat /= 0) then
            errorcode = 8
            call decomp_2d_abort(errorcode, 'Memory allocation failed when creating new arrays')
        end if
        var = 0

        return

    end subroutine AllocateInteger3DArray

    subroutine AllocateReal3DArray(var, st1, en1, st2, en2, st3, en3)

        use decomp_2d

        implicit none

        integer, intent(in)                                     :: st1, en1, st2, en2, st3, en3
        real, allocatable, dimension(:, :, :), intent(inout)    :: var
        integer                                                 :: alloc_stat, errorcode

        if (.not. allocated(var)) allocate (var(st1:en1, st2:en2, st3:en3), stat=alloc_stat)
        if (alloc_stat /= 0) then
            errorcode = 8
            call decomp_2d_abort(errorcode, 'Memory allocation failed when creating new arrays')
        end if
        var = 0.0

        return

    end subroutine AllocateReal3DArray

    subroutine AllocateRealFFTArray(var, dcmp, pdir)

        use decomp_2d

        implicit none

        real, allocatable, intent(inout)    :: var(:,:,:)
        type(decomp_info), intent(in)       :: dcmp
        character, intent(in)               :: pdir
        integer                             :: ists

        if (.not.allocated(var)) then

            if (pdir.eq.'x') then
                allocate(var(dcmp%xst(1):dcmp%xen(1),dcmp%xst(2):dcmp%xen(2),dcmp%xst(3):dcmp%xen(3)),stat=ists)
            else if (pdir.eq.'y') then
                allocate(var(dcmp%yst(1):dcmp%yen(1),dcmp%yst(2):dcmp%yen(2),dcmp%yst(3):dcmp%yen(3)),stat=ists)
            else if (pdir.eq.'z') then
                allocate(var(dcmp%zst(1):dcmp%zen(1),dcmp%zst(2):dcmp%zen(2),dcmp%zst(3):dcmp%zen(3)),stat=ists)
            end if

            if (ists.ne.0) call decomp_2d_abort(8, 'Memory allocation failed when creating new arrays')

            var = 0.0

        end if

    end subroutine AllocateRealFFTArray

    subroutine AllocateComplexFFTArray(var, dcmp, pdir)

        use decomp_2d

        implicit none

        complex, allocatable, intent(inout) :: var(:,:,:)
        type(decomp_info), intent(in)       :: dcmp
        character, intent(in)               :: pdir
        integer                             :: ists

        if (.not.allocated(var)) then

            if (pdir.eq.'x') then
                allocate(var(dcmp%xst(1):dcmp%xen(1),dcmp%xst(2):dcmp%xen(2),dcmp%xst(3):dcmp%xen(3)),stat=ists)
            else if (pdir.eq.'y') then
                allocate(var(dcmp%yst(1):dcmp%yen(1),dcmp%yst(2):dcmp%yen(2),dcmp%yst(3):dcmp%yen(3)),stat=ists)
            else if (pdir.eq.'z') then
                allocate(var(dcmp%zst(1):dcmp%zen(1),dcmp%zst(2):dcmp%zen(2),dcmp%zst(3):dcmp%zen(3)),stat=ists)
            end if

            if (ists.ne.0) call decomp_2d_abort(8, 'Memory allocation failed when creating new arrays')

            var = complex(0.0,0.0)

        end if

    end subroutine AllocateComplexFFTArray

    !----------------------------------------------------------------------------------------------------------------------!

    subroutine DestroyInteger1DArray(var)

        use decomp_2d

        implicit none

        integer, allocatable, dimension(:), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyInteger1DArray

    subroutine DestroyReal1DArray(var)

        use decomp_2d

        implicit none

        real, allocatable, dimension(:), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyReal1DArray

    subroutine DestroyLogical2DArray(var)

        use decomp_2d

        implicit none

        logical, allocatable, dimension(:, :), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyLogical2DArray

    subroutine DestroyInteger2DArray(var)

        use decomp_2d

        implicit none

        integer, allocatable, dimension(:, :), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyInteger2DArray

    subroutine DestroyReal2DArray(var)

        use decomp_2d

        implicit none

        real, allocatable, dimension(:, :), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyReal2DArray

    subroutine DestroyInteger3DArray(var)

        use decomp_2d

        implicit none

        integer, allocatable, dimension(:, :, :), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyInteger3DArray
    
    subroutine DestroyReal3DArray(var)

        use decomp_2d

        implicit none

        real, allocatable, dimension(:, :, :), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyReal3DArray

    subroutine DestroyRealFFTArray(var)

        use decomp_2d

        implicit none

        real, allocatable, dimension(:, :, :), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyRealFFTArray

    subroutine DestroyComplexFFTArray(var)

        use decomp_2d

        implicit none

        complex, allocatable, dimension(:, :, :), intent(inout) :: var

        if (allocated(var)) deallocate (var)

        return

    end subroutine DestroyComplexFFTArray

end module AuxiliaryRoutines