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
  PUBLIC :: read_params_guderley_problem,  draw_guderley_sponges, &
            guderley_penalization2D
  !**********************************************************************************************

 


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
!    call init_simple_sponge(FILE)

    if (params%mpirank==0) then
     write(*,*)
     write(*,*)
     write(*,*) "PARAMS: Guderley Problem!"
     write(*,'(" ------------------------")')
   endif
!   call read_param_mpi(FILE, 'Guderley_Problem', 'name', Guderley_Problem%name, 'no_Guderley' )

   ! geometry of the solid obstacle

end subroutine read_params_guderley_problem
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Allocate and compute mask for 2D/3D funnel. The different parts of the mask can be
!> colored if the boolean mask_is_colored is true. If the boolean is false then mask returns
!> only the mask of the solid objects (like walls and plates)
!subroutine  draw_skimmer(x_0, delta_x, Bs, g, mask, mask_is_colored)
!end subroutine draw_guderley
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++








!==========================================================================
subroutine draw_guderley_sponges(mask, x_0, delta_x, Bs, g )
  implicit none
  ! -----------------------------------------------------------------
  integer(kind=ik), intent(in)  ::  g        !< grid parameter
  integer(kind=ik),dimension(3),intent(in) :: Bs
  real(kind=rk), intent(in)     :: x_0(3), delta_x(3) !< coordinates of block and block spacinf
  real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
  real(kind=rk) :: r0, dr, dx, x0      !< cylindrical coordinates
  real(kind=rk) :: chi
  ! -----------------------------------------------------------------

    call cartesian2cylinder(x_0(1:2),delta_x(1:2),dx,dr,x0,r0)
!    call sponge_2D(mask(:,:,1),r0, dr, x0, dx, Bs, g)

!    mask(ix,ir,:)     = 0.0_rk

    ! Guderley wall at the top
    !  -------------------------

!    chi = draw_guderley_wall(x,r,guderley,h)
!    if (chi>0.0_rk) then
!       mask(ix,ir,UxF)   = mask(ix,ir,UxF) + chi
!       mask(ix,ir,UyF)   = mask(ix,ir,UyF) + chi
!       mask(ix,ir,pF)   = mask(ix,ir,pF) + chi
!
!    end if



end subroutine draw_guderley_sponges
!==========================================================================







!==========================================================================
!> This function adds a penalization term
!> to navier stokes equations
subroutine guderley_penalization2D(Bs, g, x_0, delta_x, mask, phi_ref)
      implicit none
      ! -----------------------------------------------------------------
      integer(kind=ik), intent(in)  ::  g          !< grid parameter
      integer(kind=ik), dimension(3), intent(in) :: Bs
      real(kind=rk), intent(in)     :: x_0(2), delta_x(2)   !< coordinates of block and block spacinf
      real(kind=rk), intent(inout)  :: phi_ref(:,:,:) !< reference values of penalized volume
!      real(kind=rk), intent(inout)  :: phi(:,:,:) 
      real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
      real(kind=rk)                 ::  dr, dx, x0, r0
      ! -----------------------------------------------------------------
      real(kind=rk) ::  r, x
      real(kind=rk),save,allocatable :: tmp_mask(:,:)    !< mask function
      integer(kind=ik):: ix, ir
      real(kind=rk)    :: rho
      call cartesian2cylinder(x_0, delta_x, dx, dr, x0, r0)
      
      if (.not. allocated(tmp_mask)) allocate(tmp_mask(Bs(1)+2*g,Bs(2)+2*g))
      if (size(mask,1) /= Bs(1)+2*g) call abort(7109,"wrong array size, there's pirates, captain!")
      if (size(phi_ref,1) /= Bs(1)+2*g) call abort(777109,"wrong array size, there's pirates, captain!") 
     
  !     call sponge_2D(tmp_mask, x0, dx, Bs, g)
      
       
 
      do ir = g+1, Bs(2) + g
        do ix = g+1, Bs(1) + g
          ! this if is necessary to not overwrite the privious values
          r  = dble(ir-(g+1)) * dr + r0

!          rho = phi(ix,ir,rhoF)    

          if (tmp_mask(ix,ir) > 0.0_rk) then
            ! mask of the inlet and outlet sponge
          !  mask(ix,ir,rhoF) = C_sp_inv*tmp_mask(ix,ir)
          !  mask(ix,ir,UxF ) = C_sp_inv*tmp_mask(ix,ir)
          !  mask(ix,ir,UyF ) = C_sp_inv*tmp_mask(ix,ir)
            mask(ix,ir,pF  ) = C_sp_inv*tmp_mask(ix,ir)
            ! values of the sponge inlet
            !phi_ref(ix,ir,rhoF)= rho0
          !  phi_ref(ix,ir,UxF) = 0.0_rk
          !  phi_ref(ix,ir,UyF) = 0.0_rk
!            if (r > params_ns%domain_size(2)*0.05_rk) then
!              phi_ref(ix,ir,pF) = rho*Rs*guderley%temperature 
!            end if
          end if
        end do
      end do

end subroutine guderley_penalization2D
!=========================================================================


    
!=============================================================================

!function draw_walls(x,r,guderley,h)

!    real(kind=rk), intent(in)       ::x,r,h
!    type(type_guderley), intent(in) ::guderley

!    real(kind=rk)                   ::mask, draw_walls

!    if (r> domain_size(2)*0.05) then
!       mask = mask+hardstep(x)
    
!    end if

!    draw_walls= mask

!end function draw_walls


!=========================================================================

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

