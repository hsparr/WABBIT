!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name final_stage_RK.f90
!> \version 0.5
!> \author sm
!
!> \brief final stage of Runge-Kutta time step. Gives back data field  at t+dt
!
!>
!! for the RK4 the final stage looks like this: \n
!! data_field(t+dt) = data_field(t) + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4)
!!
!! input:
!!           - params
!!           - heavy data
!!           - time step dt
!!           - coefficients for Runge Kutta
!!
!! output:
!!           - hvy_work
!!
!!
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!!
!!
!! = log ======================================================================================
!! \n
!! 23/05/17 - create
!
!**********************************************************************************************

subroutine final_stage_RK(params, dt, hvy_work, hvy_block, hvy_active, hvy_n, rk_coeffs)

!---------------------------------------------------------------------------------------------
! modules

!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> dt
    real(kind=rk), intent(in)           :: dt
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(inout)        :: hvy_work(:, :, :, :, :, :)

    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n


    ! array containing Runge-Kutta coefficients
    real(kind=rk), intent(in)           :: rk_coeffs(:,:)

    ! loop variables
    integer(kind=ik)                    :: dF, k, j, Neqn, Bs, g, z1, z2

    Neqn  = params%n_eqn
    Bs    = params%Bs
    g     = params%n_ghosts

    if (params%dim==2) then
        z1 = 1
        z2 = 1
    else
        z1 = g+1
        z2 = Bs+g
    endif


    ! loop over all active heavy data blocks
    do k = 1, hvy_n
        !u_n = u_n +...
        hvy_block( g+1:Bs+g, g+1:Bs+g, z1:z2, 1:Neqn, hvy_active(k)) = hvy_work( g+1:Bs+g, g+1:Bs+g, z1:z2, 1:Neqn,hvy_active(k), 1)

        do j = 2, size(rk_coeffs, 2)
            if ( abs(rk_coeffs(size(rk_coeffs, 1),j)) < 1e-8_rk) then
            else
                ! ... dt*(b1*k1 + b2*k2+ ..)
                ! rk_coeffs(size(rk_coeffs,1)) , since we want to access last line,
                ! e.g. b1 = butcher(last line,2)
                hvy_block( g+1:Bs+g, g+1:Bs+g, z1:z2, 1:Neqn, hvy_active(k)) = hvy_block( g+1:Bs+g, g+1:Bs+g, z1:z2, 1:Neqn, hvy_active(k)) &
                       + dt*rk_coeffs(size(rk_coeffs,1),j) * hvy_work( g+1:Bs+g, g+1:Bs+g, z1:z2, 1:Neqn, hvy_active(k), j)
            end if
        end do
    end do


end subroutine final_stage_RK
