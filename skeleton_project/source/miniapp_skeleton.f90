!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page miniapp skeleton Program

!> @brief Main program used to illustrate how to write LFRic miniapps.

!> @details Calls init, run and finalise routines from a driver module

program miniapp_skeleton

  use cli_mod,                     only : get_initial_filename
  use skeleton_miniapp_driver_mod, only : initialise, run, finalise

  implicit none

  character(:), allocatable :: filename

  call get_initial_filename( filename )
  call initialise( filename )
  deallocate( filename )

  call run()

  call finalise()

end program miniapp_skeleton
