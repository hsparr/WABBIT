
!-----------------------------------------------------------------
!> \file
!> \brief
!! Module of 2D/3D ion skimmer
!> \version 17.10.2018
!> \author H. Sparr
!-----------------------------------------------------------------

module module_skimmer

  !-------------------------------------------------------
  ! modules
  use module_navier_stokes_params
  use module_precision
  use module_ns_penalization
  use module_ini_files_parser_mpi, only : read_param_mpi
  use mpi
  !--------------------------------------------------------
  ! variables

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: integrate_over_pump_area,read_params_skimmer,mean_quantity,draw_skimmer, &
            set_inicond_skimmer,skimmer_penalization2D
  !**********************************************************************************************

!  real(kind=rk),    allocatable,     save        :: mask(:,:,:)
  character(len=80),save    :: mask_geometry!273.15_rk
  logical      ,save        :: smooth_mask, use_sponge
  !real(kind=rk),save        :: C_eta_inv, C_sp_inv, L_sponge
  real(kind=rk),save        :: domain_size(3)=0.0_rk
  ! radius of domain (Ly/2)
  real(kind=rk),save        :: R_domain
  real(kind=rk),save        :: Rs,gamma_


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! identifyers of the different parts of the skimmer
! they are used in the array mask_color
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer(kind=2),parameter :: color_capillary  =6
  integer(kind=2),parameter :: color_outlet     =5
!  integer(kind=2),parameter :: color_plates     =4
  integer(kind=2),parameter :: color_walls      =3
  integer(kind=2),parameter :: color_pumps      =2
  integer(kind=2),parameter :: color_pumps_sink =1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++



  type :: type_skimmer
      real(kind=rk)       ::outer_diameter         ! outer diameter
      real(kind=rk)       ::max_inner_diameter     ! maximal inner diameter
      real(kind=rk)       ::min_inner_diameter    =-1.0_rk ! minimal inner diameter
      real(kind=rk)       ::temperatur            =-1.0_rk ! temperatur of plates

      real(kind=rk)       ::length                =-1.0_rk ! total length of funnel
      real(kind=rk)       ::slope                 =-1.0_rk ! slope of funnel
      real(kind=rk)       ::offset(2)             =-1.0_rk ! offset of funnel in x and y

      ! parameters of flow inlet outlet
      real(kind=rk)       ::pump_diameter  =-1.0_rk
      real(kind=rk)       ::pump_x_center  =-1.0_rk
      real(kind=rk)       ::jet_radius     =-1.0_rk        ! cappilary inner Radius
      real(kind=rk)       ::r_out_cappilary=-1.0_rk         ! cappilary outer Radus
      real(kind=rk)       ::wall_thickness =-1.0_rk           !

      real(kind=rk)       ::inlet_velocity(3)       !
      real(kind=rk)       ::inlet_density       !
      real(kind=rk)       ::inlet_pressure       !
      real(kind=rk)       ::outlet_pressure       !
      real(kind=rk)       ::outlet_density
      real(kind=rk)       ::pump_speed       !
      real(kind=rk)       ::pump_density      !
      real(kind=rk)       ::pump_pressure     !
  end type type_skimmer


  !------------------------------------------------
  type(type_skimmer)   , save :: skimmer
  !------------------------------------------------

contains

  include 'skimmer2D.f90'
  



    !> \brief reads parameters for mask function from file
    subroutine read_params_skimmer(params,FILE)

      implicit none

      ! character(len=*), intent(in) :: filename
      type(inifile) , intent(inout) :: FILE
      !> params structure of navier stokes
      type(type_params_ns),intent(inout)  :: params

      real(kind=rk)                 :: dmax,dmin
      ! inifile structure
      !type(inifile) :: FILE
      !call read_ini_file_mpi(FILE, filename, .true.)

      ! READ IN geometry
      ! ----------------
      call read_param_mpi(FILE, 'skimmer', 'outer_diameter'        , skimmer%outer_diameter, R_domain*0.5_rk )
      call read_param_mpi(FILE, 'skimmer', 'maximal_inner_diameter', dmax, domain_size(2)/3.0_rk )
      call read_param_mpi(FILE, 'skimmer', 'minimal_inner_diameter', dmin, domain_size(2)/4.0_rk )
      call read_param_mpi(FILE, 'skimmer', 'Temperatur_of_plates'  , skimmer%temperatur, 300.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'jet_diameter'          , skimmer%jet_radius, R_domain*0.5_rk)
      !call read_param_mpi(FILE, 'skimmer', 'pump_diameter'         , skimmer%pump_diameter, domain_size(1)/5.0_rk)
      !optional values
       ! this parameters are global in skimmer module!
      Rs         =params%Rs
      gamma_     =params%gamma_
      domain_size=params%domain_size
      R_domain   =params%domain_size(2)*0.5_rk
      C_sp_inv   =1.0_rk/params%C_sp
      C_eta_inv   =1.0_rk/params%C_eta
      skimmer%max_inner_diameter   = dmax
      skimmer%min_inner_diameter   = dmin
      params%inicond_width        = dmax
      skimmer%wall_thickness       = 0.05*domain_size(1)
      skimmer%length               = domain_size(1)*0.95_rk-skimmer%wall_thickness*2.0_rk
      skimmer%slope                = (dmax - dmin)
 !   else....-> umschreiben für den skimmer
      ! convert diameter slope to slope in y=slope*x
 !     skimmer%slope  = skimmer%slope*0.5_rk
 !    skimmer%length = funnel%first_plate_thickness &
 !     + funnel%plates_thickness*(funnel%nr_plates-1)+funnel%plates_distance*(funnel%nr_plates-1)
 !     if ( funnel%length+2*funnel%wall_thickness>domain_size(1)) then
!        write(*,*) "funnel length + 2*walls=",funnel%length
!        call abort(7543,'your funnel does not fit in the vacuum chamber! you are a bad experimentalist! try again')
!      end if

!    end if
    skimmer%jet_radius           = skimmer%jet_radius/2.0_rk !inner radius of cappilary
    skimmer%r_out_cappilary      = skimmer%jet_radius*3.0_rk  !outer radius of cappilary
    skimmer%pump_density         = 0
    ! we pump the full right half of the chamber:
    skimmer%pump_diameter        = (domain_size(1)-2.0_rk*skimmer%wall_thickness)*0.5_rk
    skimmer%pump_x_center        = domain_size(1)-skimmer%wall_thickness-skimmer%pump_diameter*0.5_rk   ! the pump is located at the right half of the skimmer chamber
    ! READ IN Capillary inlet flow
    ! ----------------------------
    skimmer%inlet_velocity=(/ 0.0_rk, 0.0_rk, 0.0_rk /)
    call read_param_mpi(FILE, 'skimmer', 'inlet_velocity'  , skimmer%inlet_velocity(1:params%dim) &
                                                          , skimmer%inlet_velocity(1:params%dim))
    call read_param_mpi(FILE, 'skimmer', 'inlet_density'   , skimmer%inlet_density  , 1.0_rk )
    call read_param_mpi(FILE, 'skimmer', 'inlet_pressure'  , skimmer%inlet_pressure, 1.0_rk )
    call read_param_mpi(FILE, 'skimmer', 'pump_speed'      , skimmer%pump_speed, 30.0_rk )
    call read_param_mpi(FILE, 'skimmer', 'outlet_pressure' , skimmer%outlet_pressure, 1.0_rk)
    skimmer%outlet_density=skimmer%outlet_pressure/(params_ns%Rs*skimmer%temperatur)
    if (skimmer%length         >domain_size(1)-2.0_rk*skimmer%wall_thickness .or. &
    skimmer%outer_diameter >domain_size(2)-2.0_rk*skimmer%wall_thickness) then
      call abort(5032,"ERROR [skimmer.f90]:funnel is larger then simulation domain!")
    endif
    if ( skimmer%pump_diameter>domain_size(2)-2*skimmer%wall_thickness ) then
        call abort(3464,"ERROR [module_skimmer]: your pump diameter is larger then the vacuum chamber!!")
    end if
  
end subroutine read_params_skimmer


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Allocate and compute mask for 2D/3D funnel. The different parts of the mask can be
!> colored if the boolean mask_is_colored is true. If the boolean is false then mask returns
!> only the mask of the solid objects (like walls and plates)
subroutine  draw_skimmer(x0, dx, Bs, g, mask, mask_is_colored)
  implicit none
  ! -----------------------------------------------------------------
  integer(kind=ik), intent(in)  :: Bs, g        !< grid parameter
  real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
  real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
  logical, optional, intent(in) :: mask_is_colored
  integer(kind=2),allocatable   :: mask_color(:,:,:)!< identifyers of mask parts 
  logical, save :: is_colored =.false.
  real(kind=rk), allocatable  :: mask_tmp(:,:,:,:)    !< mask function for the statevector
  ! -----------------------------------------------------------------
  if (size(mask,1) /= Bs+2*g) call abort(127109,"wrong array size!")
  ! if variable is present the default (false) is overwritten by the input
  if( present(mask_is_colored)) is_colored=mask_is_colored
  ! allocate and compute mask and colored mask
  if (params_ns%dim==3) then
    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
    if (.not. allocated(mask_tmp))        allocate(mask_tmp(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g,5))
    mask_tmp    = 0.0_rk
    mask_color  = 0
!    call  draw_skimmer3D(x0, dx, Bs, g, mask_tmp, mask_color)
!    call  draw_sponge3D(x0,dx,Bs,g,mask_tmp,mask_color)
  else
    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs+2*g, 1:Bs+2*g, 1))
    if (.not. allocated(mask_tmp))        allocate(mask_tmp(1:Bs+2*g, 1:Bs+2*g, 1,4))
    mask_tmp    = 0.0_rk
    mask_color  = 0
    ! call  draw_sponge2D(x0,dx,Bs,g,mask_tmp(:,:,1,:), mask_color(:,:,1))
    call  draw_skimmer2D(x0, dx, Bs, g, mask_tmp(:,:,1,:), mask_color(:,:,1))
  endif

  ! mask coloring is optional, which is mainly used for plotting the different parts
  ! of the skimmer in paraview
  if (is_colored) then
    mask(:,:,:)= real(mask_color(:,:,:),kind=rk)
  else
    ! if the mask is not colored we use the mask of the solid obstacles
    mask(:,:,:) = mask_tmp(:,:,:,UxF)
  endif

end subroutine draw_skimmer
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> \brief Set the initial condition of a specific case
  !> \details
   subroutine set_inicond_skimmer(x0, dx, Bs, g, u)
       implicit none
       ! -----------------------------------------------------------------
       integer(kind=ik), intent(in)  :: Bs, g        !< grid parameter
       real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
       real(kind=rk), intent(inout)  :: u(:,:,:,:)    !< Statevector for t=0
       ! -----------------------------------------------------------------
       real(kind=rk),allocatable:: mask(:,:,:)
       real(kind=rk)            :: y_rel,p_init, rho_init,u_init(3),T_init,b
       integer(kind=ik)         :: iy

       p_init    =params_ns%initial_pressure
       rho_init  =params_ns%initial_density
       u_init    =params_ns%initial_velocity
       T_init    =params_ns%initial_temp

       allocate(mask(size(u,1), size(u,2), size(u,3)))
      ! set velocity field u(x)=1 for x in mask
      ! u(x)=(1-mask(x))*u0 to make sure that flow is zero at mask values
      call draw_skimmer(x0, dx, Bs, g, mask)
      u( :, :, :, pF) = p_init
      u( :, :, :, rhoF) = sqrt(rho_init)
      u( :, :, :, UxF) = ( 1 - mask ) * u_init(1)*sqrt(rho_init) !flow in x
      u( :, :, :, UyF) = (1-mask)*u_init(2)*sqrt(rho_init) !flow in y
      if (params_ns%dim==3) then
        u( :, :, :, UzF) = (1-mask)*u_init(2)*sqrt(rho_init) !flow in z
      endif
      ! if ( params_ns%geometry=="funnel" ) then
      !   do iy=g+1, Bs+g
      !       !initial y-velocity negative in lower half and positive in upper half
      !       y_rel = dble(iy-(g+1)) * dx(2) + x0(2) - params_ns%domain_size(2)*0.5_rk
      !       b=tanh(y_rel*2.0_rk/(params_ns%inicond_width))
      !       u( :, iy, 1, UyF) = (1-mask(:,iy,1))*b*u_init(2)*sqrt(rho_init)
      !   enddo
      ! endif
    end subroutine set_inicond_skimmer
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  subroutine integrate_over_pump_area(u,g,Bs,x0,dx,integral,area)
      implicit none
      !---------------------------------------------------------------
      integer(kind=ik), intent(in):: Bs, g            !< grid parameter (g ghostnotes,Bs Bulk)
      real(kind=rk), intent(in)   :: u(:,:,:,:)       !< statevector in PURE VARIABLES \f$ (rho,u,v,w,p) \f$
      real(kind=rk),  intent(in)  :: x0(3), dx(3)     !< spacing and origin of block
      real(kind=rk),intent(out)   :: integral(5), area!< mean values
      !---------------------------------------------------------------

      if ( params_ns%dim==2 ) then
        call integrate_over_pump_area2D(u(:,:,1,:),g,Bs,x0(1:2),dx(1:2),integral(1:4),area)
      end if

  end subroutine integrate_over_pump_area





  subroutine mean_quantity(integral,area)
      !> area of taking the mean
      real(kind=rk),intent(in)    :: area
      !> integral over the area
      real(kind=rk),intent(inout) :: integral(1:)

      ! temporary values
      real(kind=rk),allocatable,save :: tmp(:)
      real(kind=rk)                  :: A
      integer(kind=ik)               :: mpierr,Nq


      Nq = size(integral,1)
      if ( .not. allocated(tmp) ) allocate(tmp(Nq))

      tmp=integral

      ! integrate over all procs
      call MPI_ALLREDUCE(tmp  ,integral, Nq , MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      call MPI_ALLREDUCE(area ,A       , 1  , MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

      if ( .not. abs(A) > 0) then
        call abort(24636,"Error [skimmer.f90]: only chuck norris can devide by zero!!")
      endif

      !devide by the area of the region
      integral = integral / A
      skimmer%pump_density = integral(rhoF)
      skimmer%pump_pressure = integral(pF)
  end subroutine mean_quantity









end module module_skimmer
