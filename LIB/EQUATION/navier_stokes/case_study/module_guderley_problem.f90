!----------------------------------------------------------------
!> Implementation of Guderley Problem
!> \details
!> \version 05.02.2019
!> \ author H. Sparr
!----------------------------------------------------------------

module module_guderley_problem

  use module_navier_stokes_params
  use module_precision
  use module_ini_files_parser_mpi
  use module_ns_penalization
  use module_helpers


  implicit none

  !**********************************************************************************************
  ! only this functions are visible outside this module
  PUBLIC :: read_params_guderley_problem, set_inicond_guderley_problem, draw_guderley_sponges
  !**********************************************************************************************
  ! make everything private if not explicitly marked public
  PRIVATE
  !**********************************************************************************************
  type :: type_Guderley_Problem_params
      character(len=80) :: name
      real(kind=rk)  ::    rho_0, u_0(3), p_0 ,  rho_2a, p_2a, u_2a 
  end type type_Guderley_Problem_params

  !------------------------------------------------------------
  type(type_Guderley_Problem_params), save :: Guderley_Problem
  !------------------------------------------------------------ 

contains


!=========================================================================================
! INITIALIZATIONs
!=========================================================================================

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> \brief reads parameters for mask function from file
subroutine read_params_guderley_problem( params,FILE )
  use module_navier_stokes_params
    implicit none
    !--------------------------------------------
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params
    !---------------------------------------------
    ! add sponge parameters for the in and outflow of the domain
    call init_simple_sponge(FILE)

    if (params%mpirank==0) then
     write(*,*)
     write(*,*)
     write(*,*) "PARAMS: Guderley Problem!"
     write(*,'(" ------------------------")')
   endif
!   call read_param_mpi(FILE, 'Guderley_Problem', 'name', Guderley_Problem%name, 'no_Guderley' )

   ! geometry of the solid obstacle
!    call read_param_mpi(FILE, 'Guderley_Problem','density_0',Guderley_Problem%rho_0, 1.0_rk )
!    call read_param_mpi(FILE, 'Guderley_Problem','velocity_0', Guderley_Problem%u_0, 3.0_rk )
!    call read_param_mpi(FILE, 'Guderley_Problem','pressure_0', Guderley_Problem%p_0, 1e-15 )
!    call read_param_mpi(FILE, 'Guderley_Problem','densitry_2a', Guderley_Problem%rho_2a, 1.0_rk )
!    call read_param_mpi(FILE, 'Guderley_Problem','velocity_2a', Guderley_Problem%u_2a, 2.0_rk) 
!    call read_param_mpi(FILE, 'Guderley_Problem', 'pressure_2a',Guderley_Problem%p_2a , 2.0_rk)
    ! bottom half part Ly/2>r
!    Guderley_Problem%rho_0 = 1
!    Guderley_Problem%u_0   = 0
!    Guderley_Problem%p_0   = 1e-15
    !
    !top half part Ly/2< r
!    rho_2a = 2.1044
!    u_2a   = -1.7743
!    p_2a   = 1.9725

end subroutine read_params_guderley_problem
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!==========================================================================
subroutine set_inicond_guderley_problem (x_0, delta_x, Bs, g, u)
     implicit none
     !--------------------------------------------------------------
     integer(kind=ik), intent(in)      :: Bs, g                 !< grid parameter
     real(kind=rk),  intent(in)        :: x_0(3), delta_x(3)    !< cylindrical coordinates
     real(kind=rk)                     :: dr, dx, x0, r0
     real(kind=rk), intent(inout)      :: u(:,:,:,:)            !< Statevector for t=0
     !---------------------------------------------------------------
     real(kind=rk), allocatable        :: mask(:,:,:)
     integer(kind=rk)                  :: ir, ix 
     real(kind=rk)                     :: r, rho_0, rho_2a, p_0, p_2a, u_0(3), u_2a(3)
     real(kind=rk)                     :: u_bottom(size(u,4)), u_top(size(u,4))    
     
     rho_0 = params_ns%initial_density
     p_0   = params_ns%initial_pressure
     u_0   = params_ns%initial_velocity   


      do ix = g+1, Bs + g
        do ir = g+1, Bs + g
          ! this if is necessary to not overwrite the privious values
          r  = abs(dble(ix-(g+1)) * dr + r0)
         ! if (r<=params_ns%domain_size(2)) then
             u_bottom(:)     = u_0  
             u_bottom(rhoF)  = rho_0    
             u_bottom(pF)    = p_0      
!          else
!             u_top(:)     = u_2a  
!             u_top(rhoF)  = rho_2a   
!             u_top(pF)    = p_2a       
!          end if  
        end do
      end do

end subroutine set_inicond_guderley_problem
!==========================================================================



!==========================================================================
subroutine draw_guderley_sponges(mask, x_0, delta_x, Bs, g )
    implicit none
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk),  intent(out)                               :: mask
    !> spacing and origin of block
    real(kind=rk), intent(in)                                 :: x_0(3), delta_x(3)
    real(kind=rk)                                             :: r0, dr, x0, dx

    call cartesian2cylinder(x_0(1:2),delta_x(1:2),dx,dr,x0,r0)
  !  call sponge_2D(mask(:,:,1),r0, dr, x0, dx, Bs, g)


end subroutine draw_guderley_sponges
!==========================================================================




!==========================================================================
!> This function adds a penalization term
!> to navier stokes equations
!subroutine guderley_penalization2D(Bs, g, x_0, delta_x, mask, phi_ref)
!      implicit none
      ! -----------------------------------------------------------------
!      integer(kind=ik), intent(in)  :: Bs, g          !< grid parameter
!      real(kind=rk), intent(in)     :: x_0(2), delta_x(2)   !< coordinates of block and block spacinf
!      real(kind=rk), intent(inout)  :: phi_ref(:,:,:) !< reference values of penalized volume
!      real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
!      real(kind=rk)                 ::  dr, dx, x0, r0
      ! -----------------------------------------------------------------
!      real(kind=rk) ::  r, x
!      real(kind=rk),save,allocatable :: tmp_mask(:,:)    !< mask function
!      integer(kind=ik):: ix, ir
!      call cartesian2cylinder(x_0, delta_x, dx, dr, x0, r0)
 
    
!      if (.not. allocated(tmp_mask)) allocate(tmp_mask(Bs+2*g,Bs+2*g))
!      if (size(mask,1) /= Bs+2*g) call abort(7109,"wrong array size, there's pirates, captain!")
!      if (size(phi_ref,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")
      ! reset mask array
!      mask    = 0.0_rk
!      phi_ref = 0.0_rk
!      u_0     = params_ns%initial_velocity
!      gamma_  = params_ns%gamma_   

      ! sponge for in and outflow
      !---------------------------
  !    call sponge_2D(tmp_mask,r0,dr, x0, dx, Bs, g)

!      do ix = g+1, Bs + g
!        do ir = g+1, Bs + g
          ! this if is necessary to not overwrite the privious values
!          r  = abs(dble(ix-(g+1)) * dr + r0)
!          if (tmp_mask(ix,ir) > 0.0_rk) then
            ! mask of the inlet and outlet sponge
          !  mask(iz,ir,rhoF) = C_sp_inv*tmp_mask(iz,ir)
          !  mask(iz,ir,UxF ) = C_sp_inv*tmp_mask(iz,ir)
          !  mask(iz,ir,UyF ) = C_sp_inv*tmp_mask(iz,ir)
!            mask(ix,ir,pF  ) = C_sp_inv*tmp_mask(ix,ir)
            ! values of the sponge inlet
            !phi_ref(iz,ir,rhoF)= rho0
          !  phi_ref(iz,ir,UxF) = 0.0_rk
          !  phi_ref(iz,ir,UyF) = 0.0_rk
!            if ( r== params_ns%domain_size(1)*0.5_rk ) then
!              phi_ref(ix,ir,rhoF) = rho_0
!              phi_ref(ix,ir,UyF)  = u_0(1)
!              phi_ref(ix,ir,UyF)  = u_0(2)
!              phi_ref(ix,ir,pF)   = p_0
!              else
!              if( x<  params_ns%domain_size(1)*0.5_rk ) then
!                phi_ref(ix,ir,rhoF) = rho_0*((gamma_+1)/(gamma_-1)) 
               ! phi_ref(ix,ir,pF)   =  
!              end if       
!           end if
!          end if
!        end do
!      end do

!end subroutine guderley_penalization2D

!########################################################################
!> Transform computational coordinates to cylindrical coords
subroutine cartesian2cylinder(x_0, delta_x, dx, dr, x0, r0)

  implicit none
  !--------------------------------------------------------
  real(kind=rk), intent(in) :: x_0(2), delta_x(2)     !<cartesian coords
  real(kind=rk), intent(out) :: dx, dr, r0, x0  !<cylindrical coords
  !--------------------------------------------------------
  ! Set the lattice spacing and coordinate origin:

  x0 = x_0(1)
  dx = delta_x(1)
  dr = delta_x(2)
  ! The coordinate origin is dependent on the coordinate system
  if ( params_ns%coordinates=="cylindrical" ) then
    ! The total grid is shifted by R_min, which accounts for the
    ! infinitesimal cylinder centered arround the symmetrie axis.
    ! The origin is therefore shifted to (x,r) = (0, R_min)
    r0 = x_0(2) + params_ns%R_min
  else
    ! The coordinate system is centered at (x,r) = (0, L_r/2) and -L_r/2<r<L_r/2
    r0 = x_0(2) - params_ns%domain_size(2)*0.5_rk
  endif

end subroutine cartesian2cylinder
!########################################################################





end module module_guderley_problem


