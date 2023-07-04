!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page io_dev io_dev Miniapp
!> Test program for the XIOS IO implementation.

!> @brief Main program used to test XIOS setup and output of a field.

program io_dev

  use cli_mod,                only: get_initial_filename
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_log_mod,         only: init_logger, final_logger
  use driver_timer_mod,       only: init_timers, final_timers
  use io_dev_mod,             only: io_dev_required_namelists
  use io_dev_driver_mod,      only: initialise, run, finalise
  use io_dev_data_mod,        only: io_dev_data_type
  use mpi_mod,                only: global_mpi

  implicit none

  character(*), parameter :: program_name = "io_dev"

  type (io_dev_data_type)   :: model_data

  character(:), allocatable :: filename

  call init_comm( "io_dev" )
  call get_initial_filename( filename )
  call init_config( filename, io_dev_required_namelists )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_timers( program_name )
  call init_collections()
  deallocate( filename )


  call initialise( model_data, global_mpi, program_name )
  call run( model_data, program_name )
  call finalise( model_data, program_name )


  call final_collections()
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm()

end program io_dev
