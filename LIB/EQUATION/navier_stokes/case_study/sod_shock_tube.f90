

subroutine add_sod_shock_tube(penalization, x0, dx, Bs, g ,phi)


    implicit none

    ! grid
    integer(kind=ik), intent(in)                     :: Bs, g
    !> penalization term including mask
    real(kind=rk), dimension(:,:,:), intent(inout)   :: penalization
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)          :: x0, dx
    !> statevector
    real(kind=rk), dimension(:,:,:), intent(in)      :: phi
      !> statevector
    real(kind=rk)                                    :: mask

    ! auxiliary variables
    real(kind=rk)                                    :: x, y, h, domain_size(3)
    ! preasure,density velocities
    real(kind=rk)                                    :: p(Bs+2*g,Bs+2*g),rho(Bs+2*g,Bs+2*g), &
                                                        u(Bs+2*g,Bs+2*g),v(Bs+2*g,Bs+2*g)
    ! loop variables
    integer(kind=ik)                                 :: ix, iy,n
    ! left and right boundary
    real(kind=rk)                                    :: rho_R,rho_L,p_L,p_R,x_L,x_R,rho_ref,p_ref

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(penalization,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask         = 0.0_rk
    penalization = 0.0_rk

    ! parameter for smoothing function (width)
    h       = 1.5_rk*max(dx(1), dx(2))

    rho         = phi(:,:,1)**2
    u           = phi(:,:,2)/phi(:,:,1)
    v           = phi(:,:,3)/phi(:,:,1)
    p           = phi(:,:,4)
    domain_size = params_ns%domain_size

    ! left boundary

    x_L      =0.1_rk*domain_size(1)
    p_L      = 1.0_rk
    rho_L    = 1.0_rk

    ! right boundary

    x_R      =domain_size(1)*0.9_rk
    rho_R    = 0.125_rk
    p_R      = 0.1_rk

    ! parameter for smoothing function (width)

    do ix=1, Bs+2*g
      x = dble(ix-(g+1)) * dx(1) + x0(1)
      call continue_periodic(x,domain_size(1))
      if (x<domain_size(1)*0.5_rk) then
        mask      = smoothstep(x-x_L,h)
        rho_ref   = rho_L
        p_ref     = p_L
      else
        mask      = smoothstep(x_R-x,h)
        rho_ref   = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),rho_R,rho_L)
        p_ref     = transition(x,0.925_rk*domain_size(1),0.05_rk*domain_size(1),p_R  ,p_L    )
      endif

      ! density
      penalization(ix,:,1)= penalization(ix,:,1) + C_sp_inv*mask * ( rho(ix,:) - rho_ref )
      ! x-velocity
      penalization(ix,:,2)= penalization(ix,:,2) + C_sp_inv*mask * ( rho(ix,:)*u(ix,:) )
      ! y-velocity
      penalization(ix,:,3)= penalization(ix,:,3) + C_sp_inv*mask * ( rho(ix,:)*v(ix,:) )
      ! preasure
      penalization(ix,:,4)= penalization(ix,:,4) + C_sp_inv*mask *( p(ix,:) - p_ref )
    end do


end subroutine add_sod_shock_tube
