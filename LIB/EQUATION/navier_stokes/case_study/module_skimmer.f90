

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
  PUBLIC :: integrate_over_pump_area_skimmer,read_params_skimmer,mean_quantity_skimmer,draw_skimmer, &
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
  integer(kind=2),parameter :: color_capillary    =6
  integer(kind=2),parameter :: color_outlet       =5
  integer(kind=2),parameter :: color_plates       =2
  integer(kind=2),parameter :: color_walls        =3
  integer(kind=2),parameter :: color_pumps_2      =8
  integer(kind=2),parameter :: color_pumps_1      =4
  integer(kind=2),parameter :: color_pumps_sink_1 =1
  integer(kind=2),parameter :: color_pumps_sink_2 =7


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++

  type :: type_skimmer_plate
    real(kind=rk) :: x0(1:2)
    real(kind=rk) :: width
    real(kind=rk) :: r_in
    real(kind=rk) :: r_out
    real(kind=rk) :: height  
  end type type_skimmer_plate
 

  type :: type_skimmer
      real(kind=rk)       ::outer_diameter         ! outer diameter
      real(kind=rk)       ::max_inner_diameter     ! maximal inner diameter
      real(kind=rk)       ::max_inner_diameter_2     ! maximal inner diameter skimmer 2
      real(kind=rk)       ::min_inner_diameter    =-1.0_rk ! minimal inner diameter
      integer(kind=ik)    ::nr_plates             =0_ik ! Number of plates
      real(kind=rk)       ::plates_distance       =-1.0_rk ! distance between origin of plates
      real(kind=rk)       ::plates_distance_2       =-1.0_rk ! distance between origin of plates
      real(kind=rk)       ::plates_thickness      =-1.0_rk !
      real(kind=rk)       ::first_plate_thickness =-1.0_rk
      real(kind=rk)       ::temperatur            =-1.0_rk ! temperatur of plates
      real(kind=rk)       ::alpha_1               =-1.0_rk !
      real(kind=rk)       ::alpha_2               =-1.0_rk !
      real(kind=rk)       ::l_sk1_in              =-1.0_rk !
      real(kind=rk)       ::l_sk2_in             =-1.0_rk !
 
      real(kind=rk)       ::length                =-1.0_rk ! total length of funnel
      real(kind=rk)       ::slope                 =-1.0_rk ! slope of funnel
      real(kind=rk)       ::offset(2)             =-1.0_rk ! offset of funnel in x and y

      ! parameters of flow inlet outlet
      real(kind=rk)       ::jet_radius     =-1.0_rk        ! cappilary inner Radius
      real(kind=rk)       ::r_out_cappilary=-1.0_rk         ! cappilary outer Radus
      real(kind=rk)       ::wall_thickness_x =-1.0_rk           !Wall_thicknes in x-direct
      real(kind=rk)       ::wall_thickness_y =-1.0_rk       ! Wall_thicknes in y-direct

      real(kind=rk)       ::inlet_velocity(3)       !
      real(kind=rk)       ::inlet_density       !
      real(kind=rk)       ::inlet_pressure       !
      real(kind=rk)       ::outlet_pressure       !
      real(kind=rk)       ::outlet_density
      real(kind=rk)       ::pump_speed_1       !
      real(kind=rk)       ::pump_speed_2       !
      real(kind=rk)       ::pump_density_1      !
      real(kind=rk)       ::pump_density_2      !
      real(kind=rk)       ::pump_pressure_1     !
      real(kind=rk)       ::pump_pressure_2     !
      type(type_skimmer_plate), allocatable:: plate(:)
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
      !> params structure of navier stokekimm
      type(type_params_ns),intent(inout)  :: params

      integer(kind=ik)              :: nr_focus_plates
      ! inifile structure
      !type(inifile) :: FILE
      !call read_ini_file_mpi(FILE, filename, .true.)

      ! READ IN geometry
      ! ----------------
      call read_param_mpi(FILE, 'skimmer', 'outer_diameter'        , skimmer%outer_diameter, R_domain*0.5_rk )
      call read_param_mpi(FILE, 'skimmer', 'maximal_inner_diameter', skimmer%max_inner_diameter, domain_size(2)/3.0_rk )
      call read_param_mpi(FILE, 'skimmer', 'maximal_inner_diameter_skimmer_2',skimmer%max_inner_diameter_2, domain_size(2)/3.0_rk )
      call read_param_mpi(FILE, 'skimmer', 'minimal_inner_diameter', skimmer%min_inner_diameter, domain_size(2)/4.0_rk )
      call read_param_mpi(FILE, 'skimmer', 'angle_alpha_1'      ,skimmer%alpha_1 , 30.0_rk)            
      call read_param_mpi(FILE, 'skimmer', 'angle_alpha_2'         ,skimmer%alpha_2 , 30.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'skimmer_1_funnel_length_in' , skimmer%l_sk1_in, 1.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'skimmer_2_funnel_length_in', skimmer%l_sk2_in, 1.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'Number_of_plates'      , skimmer%nr_plates, 30 )
      call read_param_mpi(FILE, 'skimmer', 'Number_of_focus_plates', nr_focus_plates, 15)
      call read_param_mpi(FILE, 'skimmer', 'Temperatur_of_plates'  , skimmer%temperatur, 300.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'jet_diameter'          , skimmer%jet_radius, R_domain*0.5_rk)
      !optional values
      call read_param_mpi(FILE, 'skimmer', 'plates_thickness'         , skimmer%plates_thickness, 1.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'first_plate_thickness'    , skimmer%first_plate_thickness, 1.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'distance_first_skimmer_Wallwest'  , skimmer%plates_distance, 1.0_rk)
      call read_param_mpi(FILE, 'skimmer', 'distance_first_skimmer_to2nd_skimmer'  , skimmer%plates_distance_2, 1.0_rk)
       ! this parameters are global in skimmer module!
      Rs         =params%Rs
      gamma_     =params%gamma_
      domain_size=params%domain_size
      R_domain   =params%domain_size(2)*0.5_rk
      C_sp_inv   =1.0_rk/params%C_sp
      C_eta_inv   =1.0_rk/params%C_eta
      params%inicond_width        = skimmer%max_inner_diameter
      skimmer%wall_thickness_x       = 0.05*domain_size(1)
      skimmer%wall_thickness_y       = 0.05*domain_size(2)
      skimmer%length               = domain_size(1)*0.95_rk-skimmer%wall_thickness_x*2.0_rk

     if ( skimmer%plates_thickness  < 0.0_rk .or. &
      skimmer%plates_distance   <0.0_rk .or. &
      skimmer%plates_thickness  <0.0_rk ) then !default values
      skimmer%length               = domain_size(1)*0.95_rk-skimmer%wall_thickness_x*2.0_rk
      skimmer%plates_thickness     = skimmer%length/(2.0_rk*skimmer%nr_plates)
      skimmer%first_plate_thickness= skimmer%plates_thickness
      skimmer%plates_distance      = (skimmer%length-skimmer%nr_plates*skimmer%plates_thickness)/(skimmer%nr_plates-1)
      else
      ! convert diameter slope to slope in y=slope*x
      skimmer%length = skimmer%first_plate_thickness &
      + skimmer%plates_thickness*(skimmer%nr_plates-1)+skimmer%plates_distance*(skimmer%nr_plates-1)
      if ( skimmer%length+2*skimmer%wall_thickness_x>domain_size(1)) then
        write(*,*) "skimmer length + 2*walls=",skimmer%length
        call abort(7543,'your skimmer does not fit in the vacuum chamber! you are a bad experimentalist! try again')
      end if

    end if
    ! convert degrees to radians
    if ( skimmer%alpha_1>0.0_rk .and. skimmer%alpha_1<90.0_rk .or. & 
         skimmer%alpha_2>0.0_rk .and. skimmer%alpha_2<90.0_rk) then
     skimmer%alpha_1=skimmer%alpha_1*PI/180.0_rk
     skimmer%alpha_2=skimmer%alpha_2*PI/180.0_rk 
    else
     call abort(45756,"somebody has to go back to preeshool! 0< angle <90")
    end if   
    skimmer%jet_radius           = skimmer%jet_radius/2.0_rk !inner radius of cappilary
    skimmer%r_out_cappilary      = skimmer%jet_radius*1.5_rk  !outer radius of cappilary
    skimmer%pump_density_1         = 0
    skimmer%pump_density_2         = 0
    ! READ IN Capillary inlet flow
    ! ----------------------------
    skimmer%inlet_velocity=(/ 0.0_rk, 0.0_rk, 0.0_rk /)
    call read_param_mpi(FILE, 'skimmer', 'inlet_velocity'  , skimmer%inlet_velocity(1:params%dim) &
                                                          , skimmer%inlet_velocity(1:params%dim))
    call read_param_mpi(FILE, 'skimmer', 'inlet_density'   , skimmer%inlet_density  , 1.0_rk )
    call read_param_mpi(FILE, 'skimmer', 'inlet_pressure'  , skimmer%inlet_pressure, 1.0_rk )
    call read_param_mpi(FILE, 'skimmer', 'pump_speed_1'      , skimmer%pump_speed_1, 30.0_rk )
    call read_param_mpi(FILE, 'skimmer', 'pump_speed_2'      , skimmer%pump_speed_2, 30.0_rk )
    call read_param_mpi(FILE, 'skimmer', 'outlet_pressure' , skimmer%outlet_pressure, 1.0_rk)
    skimmer%outlet_density=skimmer%outlet_pressure/(params_ns%Rs*skimmer%temperatur)
    if (skimmer%length         >domain_size(1)-2.0_rk*skimmer%wall_thickness_x .or. &
    skimmer%outer_diameter >domain_size(2)-2.0_rk*skimmer%wall_thickness_y) then
      call abort(5032,"ERROR [skimmer.f90]:skimmer is larger then simulation domain!")
    endif
    !initialice geometry of ion skimmer plates
    call init_plates(skimmer) 
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
    call  draw_skimmer2D(x0, dx, Bs, g, mask_tmp(:,:,1,:), mask_color(:,:,1))
    call  draw_sponge2D(x0,dx,Bs,g,mask_tmp(:,:,1,:), mask_color(:,:,1))
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
      ! if ( params_ns%geometry=="skimmer" ) then
      !   do iy=g+1, Bs+g
      !       !initial y-velocity negative in lower half and positive in upper half
      !       y_rel = dble(iy-(g+1)) * dx(2) + x0(2) - params_ns%domain_size(2)*0.5_rk
      !       b=tanh(y_rel*2.0_rk/(params_ns%inicond_width))
      !       u( :, iy, 1, UyF) = (1-mask(:,iy,1))*b*u_init(2)*sqrt(rho_init)
      !   enddo
      ! endif
    end subroutine set_inicond_skimmer
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



  subroutine integrate_over_pump_area_skimmer(u,g,Bs,x0,dx,integral)
      implicit none
      !---------------------------------------------------------------
      integer(kind=ik), intent(in):: Bs, g            !< grid parameter (g ghostnotes,Bs Bulk)
      real(kind=rk), intent(in)   :: u(:,:,:,:)       !< statevector in PURE VARIABLES \f$ (rho,u,v,w,p) \f$
      real(kind=rk),  intent(in)  :: x0(3), dx(3)     !< spacing and origin of block
      real(kind=rk),intent(out)   :: integral(10)      !< mean values
      !---------------------------------------------------------------

      if ( params_ns%dim==2 ) then
        call integrate_over_pump_area_skimmer2D(u(:,:,1,:),g,Bs,x0(1:2),dx(1:2),integral(1:10))
      end if
  end subroutine integrate_over_pump_area_skimmer





  subroutine mean_quantity_skimmer(integral)
      !> integral over the area
      real(kind=rk),intent(inout) :: integral(10)

      ! temporary values
      real(kind=rk),allocatable,save :: tmp(:)
      real(kind=rk)                  :: A_1, A_2
      integer(kind=ik)               :: mpierr,Nq
     
               

      Nq = size(integral,1)
      if ( .not. allocated(tmp) ) allocate(tmp(Nq))

      tmp=integral

      ! integrate over all procs
      call MPI_ALLREDUCE(tmp  ,integral, Nq , MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

      A_1 = integral(5)
      A_2 = integral(10)
!      write (*,*) "A_1=", A_1  
!      write (*,*) "A_2=", A_2
      if ( .not.  A_1>  0 .or. .not. A_2> 0 ) then
        call abort(24636,"Error [skimmer.f90]: only chuck norris can devide by zero!!")
      endif
      !devide by the area of the region
      skimmer%pump_density_1 = integral(1)/ A_1
      skimmer%pump_pressure_1 = integral(4)/A_1
      skimmer%pump_density_2 = integral(6)/A_2
      skimmer%pump_pressure_2 = integral(9)/A_2
      write (*,*) "pump_pressure_1=", skimmer%pump_pressure_1  
      write (*,*) "pump_pressure_2=", skimmer%pump_pressure_2
  
  end subroutine mean_quantity_skimmer









end module module_skimmer
