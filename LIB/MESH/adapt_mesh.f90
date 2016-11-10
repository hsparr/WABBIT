! ********************************************************************************************
! WABBIT
! ============================================================================================
! name: adapt_mesh.f90
! version: 0.4
! author: msr
!
! mesh adapting main function
!
! input:    - params, light and heavy data
! output:   - light and heavy data arrays
!
! = log ======================================================================================
!
! 10/11/16 - switch to v0.4
! ********************************************************************************************

subroutine adapt_mesh( params, block_list, block_data, neighbor_list )

!---------------------------------------------------------------------------------------------
! modules

    ! global parameters
    use module_params

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    ! user defined parameter structure
    type (type_params), intent(in)      :: params
    ! light data array
    integer(kind=ik), intent(inout)     :: block_list(:, :)
    ! heavy data array - block data
    real(kind=rk), intent(inout)        :: block_data(:, :, :, :)
    ! neighbor list
    integer(kind=ik), intent(in)        :: neighbor_list(:)

    ! loop variables
    integer(kind=ik)                    :: k

!---------------------------------------------------------------------------------------------
! interfaces

    interface
        subroutine threshold_block( params, block_list, block_data, neighbor_list )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
        end subroutine threshold_block

        subroutine ensure_gradedness( block_list, neighbor_list, N, max_treelevel )
            use module_params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            integer(kind=ik), intent(in)                :: neighbor_list(:)
            integer(kind=ik), intent(in)                :: N, max_treelevel
        end subroutine ensure_gradedness

        subroutine ensure_completeness( block_list, max_treelevel )
            use module_params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            integer(kind=ik), intent(in)                :: max_treelevel
        end subroutine ensure_completeness

        subroutine coarse_mesh( params, block_list, block_data )
            use module_params
            type (type_params), intent(in)              :: params
            integer(kind=ik), intent(inout)             :: block_list(:, :)
            real(kind=rk), intent(inout)                :: block_data(:, :, :, :)
        end subroutine coarse_mesh

    end interface

!---------------------------------------------------------------------------------------------
! variables initialization

!---------------------------------------------------------------------------------------------
! main body

    ! maximal number of loops to coarsen the mesh == one block go down from max_treelevel to min_treelevel
    do k = 1, (params%max_treelevel - params%min_treelevel)

        ! check where to coarsen (refinement done with safety zone)
        call threshold_block( params, block_list, block_data, neighbor_list )

        ! unmark blocks that cannot be coarsened due to gradedness
        call ensure_gradedness( block_list, neighbor_list, params%number_blocks, params%max_treelevel )

        ! ensure completeness
        call ensure_completeness( block_list, params%max_treelevel )

        ! adapt the mesh
        !call interpolate_mesh()

    end do

end subroutine adapt_mesh
