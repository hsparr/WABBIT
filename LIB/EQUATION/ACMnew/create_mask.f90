!-------------------------------------------------------------------------------
! 2 possibilities:
!   - grid_qty is given: we possibly use a part of grid_qty as time-independent mask function
!   - grid_qty is not given: all masks are created, incl non-moving parts
!-------------------------------------------------------------------------------
subroutine create_mask_3D( time, x0, dx, Bs, g, mask, us, grid_qty )
    implicit none

    ! grid
    integer(kind=ik), intent(in) :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:,:), intent(inout) :: mask
    real(kind=rk), dimension(:,:,:,:), intent(inout) :: us
    real(kind=rk), dimension(:,:,:,:), optional, intent(in) :: grid_qty
    !> spacing and origin of block
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3), time

    integer(kind=2), allocatable, save :: mask_color(:,:,:)

    ! usually, the routine should not be called with no penalization, but if it still
    ! happens, do nothing.
    if ( params_acm%penalization .eqv. .false.) return


    if (size(mask,1) /= Bs+2*g .or. size(mask,2) /= Bs+2*g .or. size(mask,3) /= Bs+2*g ) then
        call abort(777107, "mask: wrong array size, there's pirates, captain!")
    endif

    if (size(us,4) /= 3 ) then
        call abort(777108, "us: wrong array size, there's pirates, captain!")
    endif


    mask = 0.0_rk
    us = 0.0_rk


    select case (params_acm%geometry)
    case ('Insect')
        !-----------------------------------------------------------------------
        ! INSECT MODULE
        !-----------------------------------------------------------------------
        if (.not. allocated(mask_color)) allocate(mask_color(1:Bs+2*g,1:Bs+2*g,1:Bs+2*g))

        ! case I: the body is fixed and alread created on "grid_qty"
        if ( Insect%body_moves == "no" .and. present(grid_qty) ) then
            ! the grid_qty already contains the body mask + fixed obstacles. Note that the insect
            ! module needs to add the wings to this existing data, not the other way around (i.e. not
            ! adding the body afterwards)
            mask = grid_qty(:,:,:,1)
            mask_color = int( grid_qty(:,:,:,5), kind=2 )

            ! add the wings to the existing mask. note: the "old" wings from the previous
            ! time step are already deleted in grid_qty (in update_grid_qtys_ACM, the mask is deleted)
            ! Hence, here, you do not have to delete again.
            call Draw_Insect( time, Insect, x0, dx, mask(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g), &
                mask_color(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g), us(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,1:3), &
                with_body = .false., with_wings = .true., delete_before_drawing = .false. )

        ! case II: the body moves OR the grid qty is not available: create everything
        else
            ! note the shift in origin: we pass the coordinates of point (1,1,1) since the insect module cannot
            ! know that the first g points are in fact ghost nodes...
            call Draw_Insect( time, Insect, x0, dx, mask(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g), &
                mask_color(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g), us(g+1:Bs+g,g+1:Bs+g,g+1:Bs+g,1:3) )
        endif

    case ('none')
        mask = 0.0_rk
        us = 0.0_rk

    case default
        call abort(120001,"ERROR: geometry for 3d VPM is unknown "//params_acm%geometry)

    end select

end subroutine create_mask_3D

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine create_mask_2D( time, x0, dx, Bs, g, mask, us )
    implicit none

    ! grid
    integer(kind=ik), intent(in) :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(inout) :: mask
    real(kind=rk), dimension(:,:,:), intent(inout) :: us
    !> spacing and origin of block
    real(kind=rk), intent(in) :: x0(1:2), dx(1:2), time

    ! some cheap checks
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")
    if (size(us,3) /= 2) call abort(777209,"wrong array size, there's pirates, captain!")

    mask = 0.0_rk
    us = 0.0_rk

    ! usually, the routine should not be called with no penalization, but if it still
    ! happens, do nothing.
    if ( params_acm%penalization .eqv. .false.) return


    select case (params_acm%geometry)
    case ('cylinder')
        call draw_cylinder( mask, x0, dx, Bs, g )

    case ('two-cylinders')
        call draw_two_cylinders( mask, x0, dx, Bs, g )

    case ('none')
        mask = 0.0_rk

    case default
        call abort(120002,"ERROR: geometry for 2d VPM is unknown"//params_acm%geometry)

    end select

end subroutine create_mask_2D

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine draw_cylinder(mask, x0, dx, Bs, g )

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in)                              :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)                   :: x0, dx

    ! auxiliary variables
    real(kind=rk)                                             :: x, y, r, h
    ! loop variables
    integer(kind=ik)                                          :: ix, iy

    !---------------------------------------------------------------------------------------------
    ! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk

    !---------------------------------------------------------------------------------------------
    ! main body


    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
        y = dble(iy-(g+1)) * dx(2) + x0(2) - params_acm%x_cntr(2)
        do ix=1, Bs+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1) - params_acm%x_cntr(1)
            ! distance from center of cylinder
            r = dsqrt(x*x + y*y)
            if (params_acm%smooth_mask) then
                mask(ix,iy) = smoothstep(r, params_acm%R_cyl, h)
            else
                ! if point is inside the cylinder, set mask to 1
                if (r <= params_acm%R_cyl) then
                    mask(ix,iy) = 1.0_rk
                else
                    mask(ix,iy) = 0.0_rk
                end if
            end if
        end do
    end do

end subroutine draw_cylinder

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------

subroutine draw_two_cylinders( mask, x0, dx, Bs, g)

    use module_params
    use module_precision

    implicit none

    ! grid
    integer(kind=ik), intent(in)                   :: Bs, g
    !> mask term for every grid point of this block
    real(kind=rk), dimension(:,:), intent(out)     :: mask
    !> spacing and origin of block
    real(kind=rk), dimension(2), intent(in)        :: x0, dx

    ! auxiliary variables
    real(kind=rk)         :: x1, x2, y1, y2, R, cx1, cx2, cy1,&
    cy2, r_1, r_2, h, mask1, mask2
    ! loop variables
    integer(kind=ik)      :: ix, iy

    !---------------------------------------------------------------------------------------------
    ! variables initialization
    if (size(mask,1) /= Bs+2*g) call abort(777109,"wrong array size, there's pirates, captain!")

    ! reset mask array
    mask = 0.0_rk
    mask1 = 0.0_rk
    mask2 = 0.0_rk

    !---------------------------------------------------------------------------------------------
    ! main body

    ! center of the first cylinder
    cx1 = 0.5884_rk*params_acm%domain_size(1)
    cy1 = 0.4116_rk*params_acm%domain_size(2)

    ! center of the second cylinder
    cx2 = 0.4116_rk*params_acm%domain_size(1)
    cy2 = 0.5884_rk*params_acm%domain_size(2)

    ! radius of the cylinders
    R = params_acm%R_cyl
    ! parameter for smoothing function (width)
    h = 1.5_rk*max(dx(1), dx(2))

    do iy=1, Bs+2*g
        y1 = dble(iy-(g+1)) * dx(2) + x0(2) - cy1
        y2 = dble(iy-(g+1)) * dx(2) + x0(2) - cy2
        do ix=1, Bs+2*g
            x1 = dble(ix-(g+1)) * dx(1) + x0(1) - cx1
            x2 = dble(ix-(g+1)) * dx(1) + x0(1) - cx2
            ! distance from center of cylinder 1
            r_1 = dsqrt(x1*x1 + y1*y1)
            ! distance from center of cylinder 2
            r_2 = dsqrt(x2*x2 + y2*y2)
            if (params_acm%smooth_mask) then
                mask1 = smoothstep( r_1, R, h)
                mask2 = smoothstep( r_2, R, h)
                mask(ix,iy) = mask1 + mask2
            else
                ! if point is inside one of the cylinders, set mask to 1
                if (r_1 <= R) then
                    mask(ix,iy) = 1.0_rk
                elseif ( r_2 <= R) then
                    mask(ix,iy) = 1.0_rk
                else
                    mask(ix,iy) = 0.0_rk
                end if
            end if
        end do
    end do


end subroutine draw_two_cylinders
