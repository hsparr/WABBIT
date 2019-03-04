!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name init_physics_modules.f90
!> \version 0.5
!> \author engels
!
!> \brief Wrapper to call the parameter-initialization routines for the physics module in use.
!!
!! NOTE: this is *NOT* a part of module_params.f90 in order to avoid circulars in the makefile
!!
!! input:
!!           - params struct to know what module is used.
!!           - filename of the ini file to be read
!! output:
!!           - parameters are stored in the modules themselves.
!!
!! = log ======================================================================================
!! \n
!! 17/12/17 - create
!
!**********************************************************************************************

subroutine init_physics_modules( params, filename, n_gridQ)
  ! of course, we have to load all available physics modules here.
  use module_physics_metamodule
  ! NOTE: this is *NOT* a part of module_params.f90 in order to avoid circulars in the makefile
  ! therefore load also params module.
  use module_params

  implicit none

  type (type_params), intent(in) :: params
  ! number of grid-dependent (and not time-dependend qtys) is decided by the physics modules
  integer(kind=ik), intent(out) :: n_gridQ
  character(len=*), intent(in) :: filename

  if (params%rank==0) then
    write(*,'(80("-"))')
    write(*,*) "Initializing physics modules"
    write(*,'(80("-"))')
  endif

  n_gridQ = 0

  ! call the initialization routines for the physics module that is in use
  call READ_PARAMETERS_meta( params%physics_type, filename, n_gridQ )

end subroutine
