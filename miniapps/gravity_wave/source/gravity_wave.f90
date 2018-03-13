!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page gravity_wave Gravity Wave miniapp
!> Test program for the automatic generation of boundary condition enforcement
!> by PSyclone.
!>
!> @brief Main program used to simulate the linear gravity waves equations.

program gravity_wave

  use cli_mod,                 only: get_initial_filename
  use gravity_wave_driver_mod, only: initialise, run, finalise

  implicit none

  character(:), allocatable :: filename

  call get_initial_filename( filename )
  call initialise( filename )
  deallocate( filename )

  call run()

  call finalise()

end program gravity_wave
