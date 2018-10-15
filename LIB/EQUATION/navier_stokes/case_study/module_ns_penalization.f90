
!-----------------------------------------------------------------
!> \file
!> \brief
!! Module for volume penaliza \f$\chi(x,t)\f$
!> \details
!> This module implements the mask function on each Block
!!          \f[
!!                           \chi(x,t)\quad \forall\;  x\in\mathcal{–}^l
!!            \f]
!!                        for Volume penalization
!> To increase peformance there are certain tasks:
!> \todo
!>       + include update flag: block_updated
!>                - if true the mask needs to be computed again
!>                - if false the grid has not changed on the proc rank and the mask function stays
!!                  the same
!!       + maybe it is better to have a global mask on the finest grid level and coarsen it for
!!         lower levels on the specific Blocks
!!
!> \version 23.2.2018
!> \author P.Krah
!-----------------------------------------------------------------

module module_ns_penalization

  !---------------------------------------------------------------------------------------------
  ! modules
  use module_navier_stokes_params
  use module_precision
  use module_ini_files_parser_mpi
  use module_helpers, only:smoothstep
  use mpi
  !---------------------------------------------------------------------------------------------
  ! variables

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: init_penalization,smoothstep,hardstep,soft_bump, &
            soft_bump2,hard_bump,jet_stream,add_penalization_term, &
            transition,draw_free_outlet_wall, &
            init_simple_sponge,draw_simple_sponge
  !**********************************************************************************************

!  real(kind=rk),    allocatable,     save        :: mask(:,:,:)
  character(len=80),                 save        :: mask_geometry!273.15_rk
  logical           ,                save        :: smooth_mask, use_sponge
  real(kind=rk),public             , save        :: C_eta_inv,C_sp_inv,L_sponge
  real(kind=rk),                     save        :: domain_size(3)=0.0_rk
  ! radius of domain (Ly/2)
  real(kind=rk),                     save        :: R_domain
  real(kind=rk),                     save        :: Rs,gamma_

contains


!> \brief reads parameters for mask function from file
subroutine init_penalization( params,FILE )
    use module_navier_stokes_params
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params

     if (params%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: penalization!"
      write(*,'(" ----------------------")')
    endif
    ! =============================================================================
    ! parameters needed for ns_physics module
    ! -----------------------------------------------------------------------------
    call read_param_mpi(FILE, 'VPM', 'penalization', params%penalization, .false.)
    if (.not.params%penalization)  return

    call read_param_mpi(FILE, 'Sponge', 'C_sponge',  params%C_sp, 0.01_rk )
    call read_param_mpi(FILE, 'VPM', 'C_eta', params%C_eta, 0.01_rk )
    ! =============================================================================
    ! parameters needed for penalization only
    ! -----------------------------------------------------------------------------
    call read_param_mpi(FILE, 'VPM', 'smooth_mask', smooth_mask, .true.)
    domain_size=params%domain_size
    gamma_=params%gamma_
    Rs    =params%Rs
    C_eta_inv=1.0_rk/params%C_eta
    C_sp_inv =1.0_rk/params%C_sp
    R_domain =domain_size(2)*0.5_rk

end subroutine init_penalization







!==========================================================================
!> This function adds a penalization term
!> to navier stokes equations
subroutine add_penalization_term(rhs,penalization,phi)

    !-------------------------------------------------------
    !> discrete spatial operator
    real(kind=rk), intent(inout)    :: rhs(:,:,:,:)
    !> state variables and penalizationtem
    real(kind=rk), intent(in)       :: phi(:,:,:,:),penalization(:,:,:,:)
    !--------------------------------------------------------

    ! preasure
    rhs(:,:,:,pF)  =rhs(:,:,:,pF) -                         penalization(:,:,:,pF)
   ! density
   rhs(:,:,:,rhoF)=rhs(:,:,:,rhoF)- 0.5_rk/ phi(:,:,:,rhoF)*penalization(:,:,:,rhoF)
   ! x-velocity
   rhs(:,:,:,UxF) =rhs(:,:,:,UxF) - 1.0_rk/phi(:,:,:,rhoF)* penalization(:,:,:,UxF)
   ! y-velocity
   rhs(:,:,:,UyF) =rhs(:,:,:,UyF) - 1.0_rk/phi(:,:,:,rhoF)* penalization(:,:,:,UyF)
   ! z-velocity
   if ( params_ns%dim==3 ) then
     rhs(:,:,:,UzF) =rhs(:,:,:,UzF) -1.0_rk/phi(:,:,:,rhoF)*penalization(:,:,:,UzF)
   end if

end subroutine add_penalization_term
!==========================================================================




!==========================================================================
  !> \brief This subroutine returns the value f of a step function \n
  !> The sharp step function is 1 if delta<=0 and 0 if delta>0 \n
  !> h is the semi-size of the smoothing area, so \n
  !> f is 1 if delta<=0-h \n
  !> f is 0 if delta>0+h \n
    function hardstep(delta)

      implicit none
      real(kind=rk), intent(in)  :: delta

      real(kind=rk)              :: hardstep
      !-------------------------------------------------
      ! cos shaped smoothing (compact in phys.space)
      !-------------------------------------------------
      if (delta<=0) then
        hardstep= 1.0_rk
      else
        hardstep = 0.0_rk
      endif

    end function hardstep
!==========================================================================



!> This function computes a normalized boxcar function
!   f(x)
!   |
! 1 |       ________
!   |      |        |
! 0 |______|_ _ _ _ |____________________x
!         x0        x0+width

function hard_bump(x,x0,width)

  real(kind=rk), intent(in)      :: x, x0, width
  real(kind=rk)                  :: hard_bump

    if (x>=x0 .and. x<=x0+width) then
        hard_bump = 1.0
    else
        hard_bump = 0.0
    endif

end function hard_bump


!> \brief
!> Thif function is the smooth version of hard_bump
!   f(x)
!   |
! 1 |       ________
!   |      |        |
! 0 |______|_ _ _ _ |____________________x
!         x0        x0+width
!> \details
!> NOTE: the smoothing layer is centered at x0 and x0+width (compare to soft_bump2)
function soft_bump(x,x0,width,h)

  real(kind=rk), intent(in)      :: x, x0, h, width
  real(kind=rk)                  :: soft_bump,d

    d=x-x0

    if (d>=h .and. d<=width-h) then
        soft_bump = 1.0
    elseif ((-h<d) .and.(d<+h)) then
        soft_bump = 0.5 * (1.0 - dcos((d+h) * pi / (2.0*h)) )
    elseif (d>width-h .and. d<width+h) then
        soft_bump = 0.5 * (1.0 - dcos((d-h-width) * pi / (-2.0*(h))) )
    else
        soft_bump = 0.0
    endif

end function soft_bump



!> \brief
!> Thif function is the smooth version of hard_bump
!   f(x)
!   |
! 1 |       ________
!   |      |        |
! 0 |______|_ _ _ _ |____________________x
!         x0        x0+width
!> \details
!> NOTE: the smoothing layer at x0 is starting x0 and ends at x0+h ,
!> respectively the smoothing arround x0+width starts at x0+width-h and ends at x0+width
!> (compare to soft_bump2)
function soft_bump2(x,x0,width,h)

  real(kind=rk), intent(in)      :: x, x0, h, width
  real(kind=rk)                  :: soft_bump2,max_R,smooth_width,radius

  soft_bump2=soft_bump(x,x0+h,width-2.0_rk*h,h)

end function soft_bump2

function jet_stream(radius,max_R,smooth_width)

  real(kind=rk), intent(in)      :: radius, smooth_width,max_R
  real(kind=rk)                  :: jet_stream

    jet_stream=0.5_rk*(1-tanh((radius-(max_R-0.5_rk*smooth_width))*2*PI/smooth_width ))
end function jet_stream


!> \brief computes a smooth transition between val_left and val_right
function transition(x,x0,trans_width,val_left,val_rigth)
  !> position x
  real(kind=rk), intent(in)     :: x
  !> beginning point of transition
  real(kind=rk), intent(in)     :: x0
  !> length of transition region
  real(kind=rk), intent(in)     :: trans_width
  !> values at the left and right of the transition region
  real(kind=rk), intent(in)     ::val_left,val_rigth

  !--------------------------------------------
  real(kind=rk)                  :: s,units,transition


    ! convert transition range from 2pi to trans_width
    units  = trans_width/(2*PI)
    ! transition function 0<s<1
    s         = 0.5_rk+0.5_rk * tanh((x-x0-0.5_rk*trans_width)/units)

    transition= s * val_rigth + (1-s) * val_left
end function transition





!> \brief reads parameters for sponge initializing a sipmle sponge form file
!> \details The simple sponge starts at x=0 and ends at x=L_sponge
subroutine init_simple_sponge( params,FILE )
    use module_navier_stokes_params
    implicit none
    !> pointer to inifile
    type(inifile) ,intent(inout)       :: FILE
   !> params structure of navier stokes
    type(type_params_ns),intent(inout)  :: params

     if (params%mpirank==0) then
      write(*,*)
      write(*,*)
      write(*,*) "PARAMS: sponge!"
      write(*,'(" --------------")')
    endif

    call read_param_mpi(FILE, 'Sponge', 'use_sponge', use_sponge, .false. )

    if (use_sponge) then
      call read_param_mpi(FILE, 'Sponge', 'L_sponge', L_sponge, 0.1_rk )
      if ( L_sponge<1e-10 ) then
        L_sponge=0.1*params%domain_size(1)
      end if
      call read_param_mpi(FILE, 'Sponge', 'C_sponge', C_sp_inv, 1.0e-2_rk )
      C_sp_inv=1.0_rk/C_sp_inv
    endif

end subroutine init_simple_sponge


!> \brief creates the mask of a simple sponge
subroutine draw_simple_sponge(sponge, x0, dx, Bs, g)

    implicit none
    !-----------------------------------------------------------------
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> sponge term for every grid point of this block
    real(kind=rk), dimension(2*g+Bs, 2*g+Bs), intent(out)     :: sponge
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx
    !-------------------------------------------------------------------
    ! auxiliary variables
    real(kind=rk)                                             :: x
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

    ! reset sponge array
    sponge = 0.0_rk
    if ( use_sponge .eqv. .false. ) then
      return
    end if


        do ix=1, Bs+2*g
          x = dble(ix-(g+1)) * dx(1) + x0(1)
          do iy=1, Bs+2*g
            sponge(ix,iy)=soft_bump2(x,(domain_size(1)-L_sponge),L_sponge,1.5*min(dx(1),dx(2)))
          end do
        end do

end subroutine draw_simple_sponge






subroutine draw_free_outlet_wall(mask, x0, dx, Bs, g )

    implicit none
    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)              :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, y,h,r,Delta_r
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

!---------------------------------------------------------------------------------------------
! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(745109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask  = 0.0_rk


!---------------------------------------------------------------------------------------------
! main body
    Delta_r =Domain_Size(2)*0.475_rk

    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do ix=1, Bs+2*g
        x = dble(ix-(g+1)) * dx(1) + x0(1)
        do iy=1, Bs+2*g
           y = dble(iy-(g+1)) * dx(2) + x0(2)
           r = abs(y-Domain_Size(2)*0.5_rk)

           mask(ix,iy) = smoothstep(Delta_r-r,h)

       end do
    end do
end subroutine draw_free_outlet_wall
!==========================================================================







end module module_ns_penalization
