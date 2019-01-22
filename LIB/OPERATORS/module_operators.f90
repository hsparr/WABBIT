!> \file
!> \callgraph
! ********************************************************************************************
! WABBIT
! ============================================================================================
!> \name module_operators.f90
!> \version 0.5
!> \author sm
!
!> \brief module for all operator routines
!
!>
!! = log ======================================================================================
!! \n
!! 28/7/17 - create
! *********************************************************************************************

module module_operators

!---------------------------------------------------------------------------------------------
! modules

use mpi
! global parameters
use module_params
! timing module
use module_timing
! use mesh module, since we want to compute origin/spacing of blocks
! use module_mesh, only : get_block_spacing_origin
implicit none


PRIVATE
!**********************************************************************************************
! These are the important routines that are visible to WABBIT:
!**********************************************************************************************
PUBLIC :: compute_vorticity, divergence


contains

    ! include "volume_integral.f90"
    include "compute_vorticity.f90"
    include "divergence.f90"



end module module_operators
