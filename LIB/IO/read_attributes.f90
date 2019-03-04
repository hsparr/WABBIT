!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name read_attributes.f90
!> \version 0.5
!> \author sm
!
!> \brief read attributes saved in a hdf5-file
! = log ======================================================================================
!
!> \date 02/02/18 - create
!

subroutine read_attributes(fname, lgt_n, time, iteration, domain, bs, tc_length, dim)

    implicit none
    !> file name
    character(len=*), intent(in)                  :: fname
    !> number of active blocks (required to allocate light data, prior to reading)
    integer(kind=ik), intent(out)                 :: lgt_n
    !> time (to be read from file)
    real(kind=rk), intent(out)                    :: time
    !> iteration (to be read from file)
    integer(kind=ik), intent(out)                 :: iteration
    !> blocksize in the file (required to allocate light data, prior to reading)
    integer(kind=ik), dimension(3), intent(out)   :: Bs
    !> length of treecodes in the file (required to allocate light data, prior to reading)
    integer(kind=ik), intent(out)                 :: tc_length
    !> data dimensionality (2 or 3)
    integer(kind=ik), intent(out)                 :: dim
    !> domain size
    real(kind=rk), dimension(3), intent(out)      :: domain

    integer(kind=ik), dimension(1)                :: iiteration, number_blocks
    real(kind=rk), dimension(1)                   :: ttime
    integer(hid_t)                                :: file_id
    integer(kind=ik)                              :: datarank, Nb
    integer(kind=hsize_t)                         :: size_field(1:4)
    integer(hsize_t), dimension(2)                :: dims_treecode

    call check_file_exists(fname)

    ! open the file
    call open_file_hdf5( trim(adjustl(fname)), file_id, .false.)
    ! read attributes
    call read_attribute(file_id, "blocks", "domain-size", domain)
    call read_attribute(file_id, "blocks", "time", ttime)
    time = ttime(1)
    call read_attribute(file_id, "blocks", "iteration", iiteration)
    iteration = iiteration(1)
    call read_attribute(file_id, "blocks", "total_number_blocks", number_blocks)
    lgt_n = number_blocks(1)

    !---------------------------------------------------------------------------
    ! Number of blocks and blocksize
    !---------------------------------------------------------------------------
    ! check if we deal with 2D or 3D data
    call get_rank_datafield(file_id, "blocks", datarank)
    if (datarank == 3) then
        ! 2D data
        call get_size_datafield( datarank, file_id, "blocks", size_field(1:datarank))
        Bs(1) = int( size_field(1), kind=ik)
        Bs(2) = int( size_field(2), kind=ik)
        Bs(3) = 1
        Nb = int( size_field(3), kind=ik)
        domain(3) = 0.0_rk
        dim = 2

    elseif (datarank == 4) then
        ! 3D data
        call get_size_datafield( datarank, file_id, "blocks", size_field(1:datarank))
        Bs(1) = int( size_field(1), kind=ik)
        Bs(2) = int( size_field(2), kind=ik)
        Bs(3) = int( size_field(3), kind=ik)
        Nb = int( size_field(4), kind=ik)
        dim = 3

    else
        ! crazy data
        call abort(33321, "Datarank neither 2d nor 3d..that is unusual.")

    endif

    if ( Nb /= lgt_n ) then
        ! the number of blocks stored in metadata and the dimensionality of the
        ! array do not match.
        write(*,*) "Nb= ", Nb, "lgt_n= ", lgt_n
        call abort(333139, "the number of blocks stored in metadata and the dimensionality of the array do not match.")
    endif

    !---------------------------------------------------------------------------
    ! length of treecodes in file
    !---------------------------------------------------------------------------
    ! NOTE: we do store only the treecode, not the level or refinement status
    ! so the length of this array is indeed the treecode length, and not treecode_length+2
    call get_size_datafield(2, file_id, "block_treecode", dims_treecode)
    tc_length = int(dims_treecode(1), kind=ik)



    ! close file and HDF5 library
    call close_file_hdf5(file_id)

end subroutine read_attributes
