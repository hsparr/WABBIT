
! ============================================================================================
!> \name coarsening_indicator.f90
!> \version 0.5
!> \author engels
!> \brief Set coarsening status for all active blocks, different methods possible
!
!> \details This routine sets the coarsening flag for all blocks. We allow for different
!! mathematical methods (everywhere / random) currently not very complex, but expected to grow
!! in the future.
!! \n
!! ------------------ \n
!! Refinement status: \n
!! ------------------ \n
!! +1 refine \n
!! 0 do nothing \n
!! -1 block wants to coarsen (ignoring other constraints, such as gradedness) \n
!! -2 block will coarsen and be merged with her sisters \n
!! ------------------ \n
!! \n
!! = log ======================================================================================
!! \n
!! 29/05/2018 create
! ********************************************************************************************
subroutine grid_coarsening_indicator( time, params, lgt_block, hvy_block, hvy_tmp, lgt_active, lgt_n, &
  hvy_active, hvy_n, indicator, iteration, hvy_neighbor)
  !---------------------------------------------------------------------------------------------
  ! modules
    use module_indicators

    implicit none
    real(kind=rk), intent(in)           :: time
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data.
    real(kind=rk), intent(inout)        :: hvy_tmp(:, :, :, :, :)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> how to choose blocks for refinement
    character(len=*), intent(in)        :: indicator
    !> coarsening iteration index. coarsening is done until the grid has reached
    !! the steady state; therefore, this routine is called several times during the
    !! mesh adaptation. Random coarsening (used for testing) is done only in the first call.
    integer(kind=ik), intent(in)        :: iteration
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)


    ! local variables
    integer(kind=ik) :: k, Jmax, neq, lgt_id, Bs, g, mpierr
    ! local block spacing and origin
    real(kind=rk) :: dx(1:3), x0(1:3), tmp(1:params%n_eqn)
    real(kind=rk) :: norm(1:params%n_eqn)


    Jmax = params%max_treelevel
    neq = params%n_eqn
    Bs = params%Bs
    g = params%n_ghosts

    !> reset refinement status to "stay" on all blocks
    do k = 1, lgt_n
        lgt_block( lgt_active(k), Jmax + idx_refine_sts ) = 0
    enddo

    !---------------------------------------------------------------------------
    !> Compute vorticity (if required)
    !---------------------------------------------------------------------------
    !! it is a little unfortunate that we have to compute the vorticity here and not in
    !! block_coarsening_indicator, where it should be. but after computing it, we have to synch
    !! its vorticity ghost nodes in order to apply the detail operator to the entire
    !! vorticity field (incl gost nodes)
    if (params%coarsening_indicator=="threshold-vorticity") then
        ! loop over my active hvy data
        do k = 1, hvy_n
            ! get lgt id of block
            call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

            ! some indicators may depend on the grid (e.g. to compute the vorticity), hence
            ! we pass the spacing and origin of the block
            call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

            ! actual computation of vorticity on the block
            call compute_vorticity( hvy_block(:,:,:,1,hvy_active(k)), hvy_block(:,:,:,2,hvy_active(k)), hvy_block(:,:,:,3,hvy_active(k)), &
            dx, Bs, g, params%order_discretization, hvy_tmp(:,:,:,1:3,hvy_active(k)) )
        enddo

        ! note here we synch hvy_tmp (=vorticity) and not hvy_block
        call sync_ghosts( params, lgt_block, hvy_tmp(:,:,:,1:3,:), hvy_neighbor, hvy_active, hvy_n )
    endif

    !---------------------------------------------------------------------------
    !> Compute normalization for eps, if desired.
    !---------------------------------------------------------------------------
    !! versions <30.05.2018 used fixed eps for all qtys of the state vector, but that is not very smart
    !! as each qty can have different mangitudes. If the switch eps_normalized is on, we compute here
    !! the vector of normalization factors for each qty that adaptivity will bebased on (state vector
    !! or vorticity). We currently use the L_infty norm. I have made bad experience with L_2 norm
    !! (to be checked...)

    !! default norm (e.g. for compressible navier-stokes) is 1 so in this case eps
    !! is an absolute value.
    norm = 1.0_rk

    if (params%eps_normalized .and. params%physics_type=="ACM-new") then
        norm = 0.0_rk
        if (params%coarsening_indicator=="threshold-vorticity") then
            ! loop over my active hvy data
            do k = 1, hvy_n
                norm(1) = max(norm(1),  maxval(abs(hvy_tmp(:, :, :, 1, hvy_active(k)))) )
            enddo
        else
            do k = 1, hvy_n
                ! call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )
                ! call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )
                ! norm(1) = norm(1) + sum(hvy_block(:,:,:,1,hvy_active(k))**2)*dx(1)*dx(2)

                ! max over velocities
                norm(1) = max( norm(1), maxval(abs(hvy_block(:,:,:,1:neq-1,hvy_active(k)))) )
                ! pressure
                norm(neq) = max( norm(neq), maxval(abs(hvy_block(:,:,:,neq,hvy_active(k)))) )
            enddo
            ! isotropy: uz=uy=ux
            norm(1:neq-1) = norm(1)
            ! norm(2) = sqrt( norm(1) )
        endif

        ! note that a check norm>1.0e-10 is in threshold-block
        tmp = norm
        call MPI_ALLREDUCE(tmp, norm, neq, MPI_DOUBLE_PRECISION, MPI_MAX, WABBIT_COMM, mpierr)

        ! during dev is it useful to know what the normalization is, if that is active
        if (params%rank == 0) then
            open(14,file='eps_norm.t',status='unknown',position='append')
            write (14,'(10(g15.8,1x))') time, norm, params%eps
            close(14)
        endif

    endif


    !---------------------------------------------------------------------------
    !> evaluate coarsing criterion on all blocks
    !---------------------------------------------------------------------------
    ! loop over all my blocks
    do k = 1, hvy_n
        ! get lgt id of block
        call hvy_id_to_lgt_id( lgt_id, hvy_active(k), params%rank, params%number_blocks )

        ! some indicators may depend on the grid, hence
        ! we pass the spacing and origin of the block (as we have to compute vorticity
        ! here, this can actually be omitted.)
        call get_block_spacing_origin( params, lgt_id, lgt_block, x0, dx )

        ! evaluate the criterion on this block.
        call block_coarsening_indicator( params, hvy_block(:,:,:,1:neq,hvy_active(k)), &
        hvy_tmp(:,:,:,1:neq,hvy_active(k)), dx, x0, indicator, iteration, &
        lgt_block(lgt_id, Jmax + idx_refine_sts), norm )
    enddo

    !---------------------------------------------------------------------------
    !> force blocks on maximum refinement level to coarsen, if parameter is set
    !---------------------------------------------------------------------------
    if (params%force_maxlevel_dealiasing) then
        do k = 1, lgt_n
            if (lgt_block(lgt_active(k), Jmax + idx_mesh_lvl) == params%max_treelevel) then
                ! force blocks on maxlevel to coarsen
                lgt_block(lgt_active(k), Jmax + idx_refine_sts) = -1
            endif
        enddo
    endif


    !> after modifying all refinement statusses, we need to synchronize light data
    call synchronize_lgt_data( params, lgt_block, refinement_status_only=.true. )

end subroutine grid_coarsening_indicator
