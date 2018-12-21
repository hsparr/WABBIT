!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name set_RK_input.f90
!> \version 0.5
!> \author sm
!> \brief set input for Runge Kutta time stepper
!>
!! gives back the input for the RHS (from which in the final stage the next
!! time step is computed).\n
!!
!! k_j = RHS(t+dt*c_j,  datafield(t) + dt*sum(a_jl*k_l))
!! (e.g. k3 = RHS(t+dt*c_3, data_field(t) + dt*(a31*k1+a32*k2)) ) \n
!!
!! This routine is in charge of setting the input and saving it in the hvy_work and hvy_block array \n
!! butcher table, e.g.
!!
!! |   |    |    |   |
!! |---|----|----|---|
!! | 0 | 0  | 0  |  0|
!! |c2 | a21| 0  |  0|
!! |c3 | a31| a32|  0|
!! | 0 | b1 | b2 | b3|
!!
! = log ======================================================================================
!
!> \date  22/05/17 - create
! ********************************************************************************************

subroutine set_RK_input(dt, params, rk_coeffs, j, hvy_block, hvy_work, hvy_active, hvy_n)
!---------------------------------------------------------------------------------------------
! variables

    implicit none

    !> dt
    real(kind=rk), intent(in)           :: dt
    !> array containing Runge-Kutta coefficients
    real(kind=rk), intent(in)           :: rk_coeffs(:)
    !> loop variable
    integer(kind=ik), intent(in)        :: j
    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> heavy data array - block data
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy work data array - block data
    real(kind=rk), intent(in)           :: hvy_work(:, :, :, :, :, :)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(in)        :: hvy_n

    ! loop variables
    integer(kind=ik)                    :: l, Neqn, k, Bs, g, z1, z2

!---------------------------------------------------------------------------------------------
! variables initialization
    Neqn  = params%n_eqn
    Bs    = params%Bs
    g     = params%n_ghosts

!---------------------------------------------------------------------------------------------
! main body

    if (params%dim==2) then
        z1 = 1
        z2 = 1
    else
        z1 = g+1
        z2 = Bs+g
    endif

    ! first: k_j = RHS(data_field(t) + ...
    ! loop over all active heavy data blocks
    do k = 1, hvy_n
        ! first slot in hvy_work is previous time step
        hvy_block(g+1:Bs+g,g+1:Bs+g,z1:z2,:,hvy_active(k)) = hvy_work(g+1:Bs+g,g+1:Bs+g,z1:z2,:,hvy_active(k),1)
    end do

    do l = 2, j
        ! check if coefficient is zero - if so, avoid loop over all data fields and active blocks
        if (abs(rk_coeffs(l)) > 1e-8_rk) then
            ! loop over all active heavy data blocks
            do k = 1, hvy_n
                ! new input for computation of k-coefficients
                ! k_j = RHS((t+dt*c_j, data_field(t) + sum(a_jl*k_l))
                hvy_block(g+1:Bs+g, g+1:Bs+g, z1:z2, :, hvy_active(k)) = hvy_block(g+1:Bs+g, g+1:Bs+g, z1:z2, :, hvy_active(k)) &
                + dt * rk_coeffs(l) * hvy_work(g+1:Bs+g, g+1:Bs+g, z1:z2, :, hvy_active(k), l)

            end do
        end if
    end do


end subroutine set_RK_input
