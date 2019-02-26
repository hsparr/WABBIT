!#########################################################################################
!                     2D SKIMMER IMPLEMENTATION
!#########################################################################################
subroutine  skimmer_penalization2D(Bs, g, x_0, delta_x, phi, mask, phi_ref)
  use module_helpers
  implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)     :: x_0(2), delta_x(2)   !< coordinates of block and block spacinf
    real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:) !< reference values of penalized volume
    real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
    integer(kind=2), allocatable,save:: mask_color(:,:)!< identifyers of mask parts (plates etc)
    logical                       :: mesh_was_adapted=.true.
    real(kind=rk)                 :: x0, dx, r0, dr   
    ! -----------------------------------------------------------------
    call cartesian2cylinder(x_0, delta_x, dx, dr, x0, r0)

    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs+2*g, 1:Bs+2*g))
    !!> todo implement function check_if_mesh_adapted (true/false) in adapt mesh
    if ( mesh_was_adapted .eqv. .true. ) then
      ! dont switch the order of draw_skimmer3D and draw_sponge3D,
      ! because mask and color are reset in the draw_skimmer
      call draw_skimmer2D(r0, dr, x0, dx, Bs, g, mask, mask_color)
      call draw_sponge2D(r0, dr, x0, dx, Bs, g, mask, mask_color)
    end if

    call compute_penal2D(mask_color,mask,phi, r0, dr, x0, dx, Bs, g ,phi_ref)

end subroutine  skimmer_penalization2D



subroutine draw_skimmer2D(r0, dr, x0, dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)                :: r0, dr, x0, dx   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, r, h
    real(kind=rk)     :: chi
    integer(kind=ik)  :: ix, ir, n ! loop variables
  ! -----------------------------------------------------------------

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx, dr)

    ! smooth width in x and y direction
    do ir=g+1, Bs+g
      ! y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(dble(ir-(g+1)) *dr +r0)
       
       do ix=g+1, Bs+g
            x = dble(ix-(g+1)) * dx + x0

            !=============================
!     /\     reset the mask function!
!    /  \    caution! Reseting should only be done once
!   /    \   for optimal performance
!  / stop \
! +--------+
          !  mask_color(ix,iy)=0
          !  mask(ix,iy,:)=0.0_rk
            mask_color(ix,ir)=0
            mask(ix,ir,:)=0.0_rk
            !=============================


            ! Skimmer
            ! ------
            chi = draw_skimmer_plates(x,r,skimmer,h)
            if (chi>0.0_rk) then
              mask_color(ix,ir)  = color_plates
              mask(ix,ir,2:4)    = mask(ix,ir,2:4) + chi
            endif                                           ! the temperature of the funnel

            ! Walls
            ! -----
            chi = draw_walls(x,r,skimmer,h)
            if (chi>0.0_rk) then                       ! default values on the walls
              mask_color(ix,ir)  = color_walls
              mask(ix,ir,2:4)    = mask(ix,ir,2:4) + chi
            endif                                           ! the temperature of the funnel

       end do
    end do

end subroutine  draw_skimmer2D




subroutine draw_sponge2D(r0, dr, x0 , dx, Bs, g, mask, mask_color)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)                :: r0, dr, x0, dx   
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, r, h
    real(kind=rk)     :: chi
    integer(kind=ik)  :: ix, ir,n ! loop variables
  ! -----------------------------------------------------------------

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx, dr)

    ! smooth width in x and y direction
    do ir=g+1, Bs+g
      ! y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(dble(ir-(g+1)) * dr +r0)
       do ix=g+1, Bs+g
            x = dble(ix-(g+1)) * dx + x0 
                                      ! the temperature of the skimmer

            ! Outlet flow: PUMPS
            ! ------------------
            ! pump volume flow
!            chi=  draw_pumps_volume_flow(x,r,skimmer,h)
!            if (chi>0.0_rk) then
!              if (x> skimmer%wall_thickness_x .and. x< skimmer%plate(1)%x0(1)) then
!              mask(ix,ir,2:3)   = mask(ix,ir,2:3)+chi
!              mask_color(ix,ir) = color_pumps_1
!              endif 
!              if( x>skimmer%plate(1)%x0(1)+skimmer%plates_thickness .and. x< skimmer%plate(2)%x0(1)) then 
!              mask(ix,ir,2:3)   = mask(ix,ir,2:3)+chi
!              mask_color(ix,ir) = color_pumps_2
!              endif
!            endif

            ! mass and energy sink
!            chi=  draw_pumps_sink(x,r,skimmer,h)
!            if (chi>0.0_rk) then
!              if (x> skimmer%wall_thickness_x .and. x< skimmer%plate(1)%x0(1)) then
!              mask_color(ix,ir) = color_pumps_sink_1
!              mask(ix,ir,1) = mask(ix,ir,1)+chi
!              mask(ix,ir,4) = mask(ix,ir,4)+chi
!              endif
!              if( x>skimmer%plate(1)%x0(1)+skimmer%plates_thickness .and. x< skimmer%plate(2)%x0(1)) then 
!              mask_color(ix,ir) = color_pumps_sink_2
!              mask(ix,ir,1) = mask(ix,ir,1)+chi
!              mask(ix,ir,4) = mask(ix,ir,4)+chi
!              endif
!             endif

            ! Inlet flow: Capillary
            ! ---------------------
            chi=  draw_jet(x,r,skimmer,h)
            if (chi>0.0_rk) then
              mask_color(ix,ir) = color_capillary
              mask(ix,ir,1:4)   =  mask(ix,ir,1:4)+chi
            endif

            ! Outlet flow: Transition to 2pump
            ! ---------------------
!            chi=  draw_outlet(x,r,skimmer,h)
            !   chi=  draw_sink(x,y,skimmer,h)
!              if (chi>0.0_rk) then
!                mask_color(ix,ir) = color_outlet
!                mask(ix,ir,1)     = mask(ix,ir,1)+chi
!                mask(ix,ir,3:4)   = mask(ix,ir,3:4)+chi
!              endif

       end do
    end do

end subroutine  draw_sponge2D





!> Computes the 2D skimmer mask with reference values of the penalized system
subroutine  compute_penal2D(mask_color,mask,phi,r0, dr, x0, dx, Bs, g ,phi_ref)
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: Bs, g          !< grid parameter
    real(kind=rk), intent(in)     :: r0, dr, x0, dx   !< coordinates of block and block spacing
    integer(kind=2), intent(inout):: mask_color(:,:)!< identifyers of mask parts (plates etc)
    real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: mask(:,:,:)     !< mask
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:)  !< funnel penalization term
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h,velocity
    real(kind=rk)     :: rho,chi,v_ref_1,v_ref_2,dq,u,v,p,C_inv
    integer(kind=ik)  :: ix, ir,n                                    ! loop variables
    real(kind=rk)     :: velocity_pump_1,rho_pump_1,pressure_pump_1,velocity_pump_2,rho_pump_2,pressure_pump_2, &    ! outlets and inlets
                      rho_capillary,u_capillary,v_capillary,p_capillary, &
                      p_2nd_pump_stage,rho_2nd_pump_stage
    real(kind=rk)     ::jet_smooth_width,pump_smooth_width_1 ,pump_smooth_width_2 ! smooth width of jet and pumpsponges (boundarylayerthickness)
    ! -----------------------------------------------------------------

    ! initialice parameters
    phi_ref     = 0.0_rk
    u_capillary       =skimmer%inlet_velocity(1)
    v_capillary       =skimmer%inlet_velocity(2)
    rho_capillary     =skimmer%inlet_density
    rho_pump_1        =skimmer%pump_density_1
    rho_pump_2        =skimmer%pump_density_2
    p_capillary       =skimmer%inlet_pressure
    velocity_pump_1   =skimmer%pump_speed_1
    velocity_pump_2   =skimmer%pump_speed_2
    pressure_pump_1   =skimmer%pump_pressure_1
    pressure_pump_2   =skimmer%pump_pressure_2
    p_2nd_pump_stage  =skimmer%outlet_pressure
    rho_2nd_pump_stage=skimmer%outlet_density

    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx, dr)
    if (3*dr<=0.05_rk*skimmer%jet_radius) then
      jet_smooth_width = 0.05_rk*skimmer%jet_radius
    else
      jet_smooth_width = 3*dr
      !call abort('ERROR [skimmer.f90]: discretication constant dy to large')
    endif

  
    ! smooth width in x and y direction
    do ir=1, Bs+2*g
      ! y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(dble(ir-(g+1))* dr + r0)
       do ix=1, Bs+2*g
            x = dble(ix-(g+1)) * dx + x0
            rho = phi(ix,ir,rhoF)
            u   = phi(ix,ir,UxF)
            v   = phi(ix,ir,UyF)
            p   = phi(ix,ir,pF)

            C_inv=C_sp_inv

            !solid obstacles: walls and plates
            ! ------
            if  (mask_color(ix,ir) == color_plates &
            .or. mask_color(ix,ir) == color_walls ) then
              Phi_ref(ix,ir,rhoF) = params_ns%initial_density
              Phi_ref(ix,ir,UxF) = 0.0_rk                     ! no velocity in x
              Phi_ref(ix,ir,UyF) = 0.0_rk                     ! no velocity in y
              Phi_ref(ix,ir,pF) = rho*Rs*skimmer%temperatur   ! pressure set according to
              C_inv=C_eta_inv
            endif                                           ! the temperature of the skimmer 

            ! Outlet flow: PUMPS
            ! ------------------
            
!            if (mask_color(ix,ir) == color_pumps_1)then            
!                v_ref_1=velocity_pump_1
!               Phi_ref(ix,ir,UxF) = 0
!               C_inv=C_eta_inv
!              if (r0>0.0_rk) then
!                Phi_ref(ix,ir,UyF) = rho*v_ref_1
!              else
!                Phi_ref(ix,ir,UyF) = -rho*v_ref_1
!              endif
!            endif 
          
!            if( mask_color(ix,ir)== color_pumps_2) then
!             v_ref_2=velocity_pump_2
!               Phi_ref(ix,ir,UxF) = 0
!               C_inv=C_eta_inv
!               if (r0> 0.0_rk) then
!                Phi_ref(ix,ir,UyF) = rho*v_ref_2
!              else
!                Phi_ref(ix,ir,UyF) = -rho*v_ref_2
!              endif
!            endif
            ! mass and energy sink
!            if (mask_color(ix,ir)==color_pumps_sink_1 ) then
!              Phi_ref(ix,ir,rhoF) = rho_pump_1 
!             Phi_ref(ix,ir,pF) = pressure_pump_1
!              C_inv=C_eta_inv
!            endif
!            if (mask_color(ix,ir)==color_pumps_sink_2 ) then
!              Phi_ref(ix,ir,rhoF) = rho_pump_2 
!              Phi_ref(ix,ir,pF) = pressure_pump_2
!              C_inv=C_eta_inv
!           endif

            ! Inlet flow: Capillary
            ! ---------------------
            if (mask_color(ix,ir)==color_capillary) then
              dq               =jet_stream(r,skimmer%jet_radius,jet_smooth_width)
              C_inv=C_sp_inv
              Phi_ref(ix,ir,rhoF) =  rho_capillary
              Phi_ref(ix,ir,UxF) =  rho_capillary*u_capillary*dq
              Phi_ref(ix,ir,UyF) =  rho_capillary*v_capillary
              Phi_ref(ix,ir,pF) =  p_capillary  !rho*Rs*skimmer%temperatur * (1 - dq) + p_capillary * dq
            endif

            ! Outlet flow: Transition to 2pump
            ! ---------------------
!              if (mask_color(ix,ir)==color_outlet) then
!                Phi_ref(ix,ir,rhoF) = rho_2nd_pump_stage
                !Phi_ref(ix,ir,UxF) = 0
!                Phi_ref(ix,ir,UyF) = 0
!                Phi_ref(ix,ir,pF) = p_2nd_pump_stage
!                C_inv=C_sp_inv!*1e2
!              endif
              ! add penalization strength to mask
              mask(ix,ir,:)=C_inv*mask(ix,ir,:)
       end do
    end do
end subroutine  compute_penal2D







!> \brief Compute mask term for stacked rings of the ion skimmer
function draw_skimmer_plates(x,r,skimmer,h)

  real(kind=rk)     , intent(in)   :: x, r, h
  type(type_skimmer) , intent(in)  :: skimmer
  real(kind=rk)                    :: draw_skimmer_plates
  integer(kind=ik)                 :: n
  ! distance frommin_inner_diameterter of cylinder
  draw_skimmer_plates=0.0_rk

  ! loop over all plates
  do n=1,skimmer%nr_plates

      draw_skimmer_plates=draw_skimmer_plates+draw_plate(x,r,skimmer%plate(n),h)

  enddo

end function draw_skimmer_plates


!> \brief Compute mask term for single ring-plate
function draw_plate(x,r,plate,h)

  real(kind=rk), intent(in)          :: x, r, h
  type(type_skimmer_plate),intent(in) ::plate

  real(kind=rk)                      ::draw_plate,delta_r1,delta_r2,r2,r1,r3,r4,x_sk_1 ,height, mask,x_sk_2
  
  mask=0.0_rk  


  height    = skimmer%min_inner_diameter/2
  x_sk_1      = (x-(skimmer%plate(1)%x0(1)-skimmer%l_sk1_in+plate%width))
 if(x< skimmer%plate(1)%x0(1)) then
 
   r1        = tan(skimmer%alpha_1)*x_sk_1+height 
else
   r1 = skimmer%plate(1)%x0(1) 
 endif

   r2        = tan(skimmer%alpha_2)*x_sk_1+height 
  delta_r1     = skimmer%plate(1)%r_out-skimmer%plate(1)%r_in

 if (x > plate%width+skimmer%plate(1)%x0(1)-skimmer%l_sk1_in  .and. x <  skimmer%plate(1)%x0(1)+plate%width.and. &
    r <r1  .and. r > r2 ) then 
     mask=1
 endif 
if (skimmer%nr_plates==2) then
  x_sk_2    = (x-(skimmer%plate(2)%x0(1)-skimmer%l_sk2_in+plate%width))
 if(x< skimmer%plate(2)%x0(1)) then
 
   r3        = tan(skimmer%alpha_1)*x_sk_2+height 
 else
   r3 = skimmer%plate(2)%x0(1) 
 endif

   r4        = tan(skimmer%alpha_2)*x_sk_2+height 
  delta_r2     = skimmer%plate(2)%r_out-skimmer%plate(2)%r_in

 if (x > plate%width+skimmer%plate(2)%x0(1)-skimmer%l_sk2_in  .and. x <  skimmer%plate(2)%x0(1)+plate%width.and. &
    r <r3  .and. r > r4 ) then 
     mask=1
 endif
endif
 if (mask>1)then
   mask=1
 endif

!  draw_plate  = soft_bump(x,plate%x0(1),plate%width,h)*soft_bump(r,plate%r_in,delta_r,h)
!  if (skimmer%nr_plates==2) then
  draw_plate = hard_bump(x,skimmer%plate(1)%x0(1),plate%width)*hard_bump(r,skimmer%plate(1)%r_out,delta_r1)+ & 
               hard_bump(x,skimmer%plate(2)%x0(1),plate%width)*hard_bump(r,skimmer%plate(2)%r_out,delta_r2)+mask 
              ! hardstep(y1)*smoothstep(y2,h)
!  else 
!  draw_plate = hard_bump(x,plate%x0(1),plate%width)*hard_bump(r,plate%r_in,delta_r1)+mask
!  endif
end function draw_plate


function draw_walls(x,r,skimmer,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_skimmer),intent(in)          ::skimmer

  real(kind=rk)                         ::  mask, draw_walls

  mask=0.0_rk

  ! wall in south and north (i.e. in radial direction)
  if ( x>skimmer%wall_thickness_x .and. x<skimmer%plate(1)%x0(1) .or. &
       x>skimmer%plate(1)%x0(1)+skimmer%plates_thickness .and. x<skimmer%plate(2)%x0(1) ) then 
        mask=mask+hardstep(R_domain-skimmer%wall_thickness_y-r)
  else
         mask=mask+hardstep(R_domain-0.333_rk*skimmer%wall_thickness_y-r)
  endif

  ! wall in east
  !mask=mask+smoothstep(x-skimmer%wall_thickness,h)
  if ( r > skimmer%jet_radius) then
      mask=mask+hardstep(x-skimmer%wall_thickness_x)
  else
      mask=mask+hardstep(x-skimmer%wall_thickness_x*0.5_rk)
  endif
  ! attach cappilary to wall in EAST
  if (  r > skimmer%jet_radius  ) then
     mask=mask+hardstep(x+skimmer%wall_thickness_x-skimmer%plate(1)%x0(1)+skimmer%l_sk1_in-skimmer%plates_thickness+0.5_rk*h)*hardstep(r-skimmer%r_out_cappilary*1.5_rk)

         !mask=mask+smoothstep(x,h)*smoothstep(r-skimmer%r_out_cappilary,h)
  endif


  ! wall in WEST
!  if ( r > skimmer%min_inner_diameter*0.5_rk) then
        !  mask=mask+smoothstep(domain_size(1)-x-skimmer%wall_thickness,h)

         mask=mask+hardstep(skimmer%plate(2)%x0(1)+skimmer%plates_thickness-x)
!  else
!        mask=mask+hardstep(domain_size(1)-x-skimmer%wall_thickness_x*0.5_rk)
        ! mask=mask+smoothstep(domain_size(1)-x-skimmer%wall_thickness*0.5_rk,h)
!  endif

   ! is needed because mask off walls overlap
  if (mask>1) then
         mask=1
  endif

  draw_walls=mask
end function draw_walls

!function draw_pumps_volume_flow(x,r,skimmer,h)

!  real(kind=rk),    intent(in)          :: x, r, h
!  type(type_skimmer),intent(in)         ::skimmer

! real(kind=rk)                         ::  mask, draw_pumps_volume_flow,r0,width

!  mask  =0
!  width =skimmer%wall_thickness_x*0.333_rk
!  r0    =(R_domain-skimmer%wall_thickness_y)
!  if ( x>skimmer%wall_thickness_x .and. x<skimmer%plate(1)%x0(1) .or. &
!       x>skimmer%plate(1)%x0(1)+skimmer%plates_thickness .and. x< skimmer%plate(2)%x0(1)) then 
!       mask = soft_bump2(r,r0,width,h)
!  endif



!  draw_pumps_volume_flow=mask
!end function draw_pumps_volume_flow

!function draw_pumps_sink(x,r,skimmer,h)

!  real(kind=rk),    intent(in)          :: x, r, h
!  type(type_skimmer),intent(in)         ::skimmer

!  real(kind=rk)                         ::r0,mask, draw_pumps_sink,width,depth
 
!  draw_pumps_sink  = 0.0_rk
!  r0    =(R_domain-skimmer%wall_thickness_y*0.666_rk)
!  depth =skimmer%wall_thickness_y*0.5_rk

!  if ( x>skimmer%wall_thickness_x .and. x<skimmer%plate(1)%x0(1) .or. &
!     x>skimmer%plate(1)%x0(1)+skimmer%plates_thickness .and. x< skimmer%plate(2)%x0(1)) then
!      mask =soft_bump2(r,r0,depth,h)
!  endif
!  draw_pumps_sink=mask 

!end function draw_pumps_sink

function draw_jet(x,r,skimmer,h)

  real(kind=rk),    intent(in)          :: x, r, h
  type(type_skimmer),intent(in)         ::skimmer

  real(kind=rk)                         ::draw_jet,length_of_jet

  length_of_jet= -skimmer%wall_thickness_x*0.5_rk+h+skimmer%plate(1)%x0(1)-skimmer%l_sk1_in+skimmer%plates_thickness-0.006_rk

  !if (r< skimmer%jet_radius) then
    ! wall in EAST
    draw_jet=soft_bump2(x,skimmer%wall_thickness_x*0.5_rk,length_of_jet,h)*smoothstep(r-skimmer%jet_radius+h,h)

    !draw_jet=smoothstep(x-skimmer%wall_thickness,h)-smoothstep(x-skimmer%wall_thickness/2.0_rk,h)
  !else
  !endif

end function draw_jet


!function draw_outlet(x,r,skimmer,h)

!  real(kind=rk),    intent(in)          :: x, r, h
!  type(type_skimmer),intent(in)         ::skimmer

!  real(kind=rk)                         ::draw_outlet

         ! wall in WEST
   !!!draw_outlet=smoothstep(domain_size(1)-x-skimmer%wall_thickness,h)
!    draw_outlet=soft_bump2(x,skimmer%plate(2)%x0(1)+skimmer%plates_thickness,domain_size(1)-skimmer%plate(2)%x0(1)-skimmer%plates_thickness-skimmer%wall_thickness_x,h) &
!            *smoothstep(r-skimmer%max_inner_diameter_2*0.5_rk+skimmer%wall_thickness_y*0.5_rk,h)
          !  *smoothstep(r-(skimmer%max_inner_diameter_2-skimmer%wall_thickness_x-h-skimmer%min_inner_diameter),h)
       !  x_sk_2    = (x-(skimmer%plate(2)%x0(1)-skimmer%l_sk2_in+skimmer%plate%width))
     ! if(x< skimmer%plate(2)%x0(1)) then
 
  !   r3        = tan(skimmer%alpha_1)*x_sk_2+height 
  ! else
  !   r3 = skimmer%plate(2)%x0(1) 
  ! endif

  !   r4        = tan(skimmer%alpha_2)*x_sk_2+height 
  !  delta_r2     = skimmer%plate(2)%r_out-skimmer%plate(2)%r_in

  !   if (x > plate%width+skimmer%plate(2)%x0(1)-skimmer%l_sk2_in  .and. x <  skimmer%plate(2)%x0(1)+plate%width.and. &
  !    r <r3  .and. r > r4 ) then 
  !     mask =1
  ! endif

!end function draw_outlet


!function draw_sink(x,r,skimmer,h)

!  real(kind=rk),    intent(in)          :: x, r, h
!  type(type_skimmer),intent(in)         ::skimmer

!  real(kind=rk)                         ::draw_sink,radius

!  radius     = sqrt((x-domain_size(1)+skimmer%wall_thickness_x*0.6_rk)**2+r**2)
!  draw_sink  = smoothstep(r-skimmer%min_inner_diameter*0.4_rk,h)


!end function draw_sink




!==========================================================================
    !> \brief initialization of all plates in the skimmer
    !> \todo insert picture
    subroutine init_plates(skimmer)

      implicit none
      !> geometric parameters of ion funnel
      type(type_skimmer), intent(inout) :: skimmer

      real(kind=rk)                   :: distance,length,length_focus,width
      integer(kind=ik)                :: n
      type(type_skimmer_plate)         :: plate

      allocate(skimmer%plate(skimmer%nr_plates))

      distance    =skimmer%plates_distance
      length      =skimmer%length
      width       =skimmer%plates_thickness
      ! length of focus area in skimmer (+distance +wdith because the last plate is in the orifice)
      length_focus           = distance+width+length - (skimmer%max_inner_diameter-skimmer%min_inner_diameter)
      ! origin of skimmer
      skimmer%offset=(/skimmer%wall_thickness_x+distance, &
                       R_domain/)
      if(skimmer%offset(1)<skimmer%wall_thickness_x) then
       call abort(13457,'Error [module_mask.f90]: your skimmer is to long')
      endif

      ! initialicd all plates
      ! ---------------------
      ! first plate is often fat and ugly (:
      skimmer%plate(1)%x0    =skimmer%offset                 ! coordinates of the midpoint of the skimmer
      skimmer%plate(1)%width =skimmer%first_plate_thickness   ! as the name says its a fat plate
      skimmer%plate(1)%r_out =skimmer%outer_diameter/2        ! outer diameter of the ring
      skimmer%plate(1)%r_in  =skimmer%max_inner_diameter/2    ! inner diameter of the ring
          
      ! all the other plates are similar
      do n=2,skimmer%nr_plates
        plate%x0(1)     = skimmer%plate(n-1)%x0(n-1)+skimmer%plate(1)%width+skimmer%plates_distance_2
        plate%x0(2)     = skimmer%offset(2)
        plate%width     = width
        plate%r_in      =skimmer%max_inner_diameter_2/2    ! inner diameter of the ring
        plate%r_out     =skimmer%outer_diameter/2
        skimmer%plate(n)=plate
     enddo

    end subroutine init_plates
!==========================================================================



!==========================================================================
  !> Integrates the flow field close to the pump outlet,
  !> and computes the area of the intagration region.
  subroutine integrate_over_pump_area_skimmer2D(u,g,Bs,x0,dx,integral)
      implicit none

      !> grid parameter (g ghostnotes,Bs Bulk)
      integer(kind=ik), intent(in)                     :: Bs, g
      !> density,pressure
      real(kind=rk), dimension(1:,1:,1:), intent(in)   :: u
      !> spacing and origin of block
      real(kind=rk), dimension(2), intent(in)          :: x0, dx
      !> mean density
      real(kind=rk),intent(out)                      :: integral(10) 

      real(kind=rk)                                    :: h,r,y,x,r0,width
      !temporal data field
      real(kind=rk),dimension(10)                       :: tmp

      integer(kind=ik)                                 :: ix,ir

!        write (*,*) "integral=" ,integral 
       h  = 1.5_rk*max(dx(1), dx(2))
      ! calculate mean density close to the pump
      width =skimmer%wall_thickness_y
      tmp   =  0.0_rk
      r0    =(R_domain-2*skimmer%wall_thickness_y)
       do ir=g+1, Bs+g
!         if (params%coordinate == "cylydrical") then
         ! y = dble(iy-(g+1)) * dx(2) + x0(2)
!          r = abs(dble(ir-(g+1)) * dr +r0
!          if  ( r>r0 .and. r<r0+width) then
!           do ix=g+1, Bs+g
!              x = dble(ix-(g+1)) * dx(1) + x0(1)
!                if ( x>skimmer%wall_thickness_x .and. x<skimmer%plate(1)%x0(1)) then 
!                tmp(1:4)  = (tmp(1:4)+ u(ix,iy,:))*r
!                tmp(5)    = (tmp(5)  + 1.0_rk)*r
!              else if (x>skimmer%plate(1)%x0(1)+skimmer%plates_thickness .and. x< skimmer%plate(2)%x0(1)) then
!                tmp(6:9)  = (tmp(6:9)+ u(ix,iy,:)) *r
!                tmp(10)    = (tmp(10)  + 1.0_rk) *r
!               endif
!            enddo    
!         else
          y = dble(ir-(g+1)) * dx(2) + x0(2)
          r = abs(y-R_domain)
         if  ( r>r0 .and. r<r0+width) then
          do ix=g+1, Bs+g
              x = dble(ix-(g+1)) * dx(1) + x0(1)
                if ( x>skimmer%wall_thickness_x .and. x<skimmer%plate(1)%x0(1)) then 
                tmp(1:4)  = tmp(1:4)+ u(ix,ir,:)
                tmp(5)    = tmp(5)  + 1.0_rk
              else if (x>skimmer%plate(1)%x0(1)+skimmer%plates_thickness .and. x< skimmer%plate(2)%x0(1)) then
                tmp(6:9)  = tmp(6:9)+ u(ix,ir,:)
                tmp(10)    = tmp(10)  + 1.0_rk
               endif
         
           enddo
          endif
!         endif
        enddo
        integral  = integral + tmp *dx(1)*dx(2)
      !  write (*,*) "integral=" , integral
  end subroutine integrate_over_pump_area_skimmer2D
  !==========================================================================

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
