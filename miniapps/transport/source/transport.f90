!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page transport Transport miniapp
!> Program file for running transport miniapp. Subroutine calls include initialise_transport(),
!> run_transport() and finalise_transport().
program transport

  use cli_mod,                only: get_initial_filename
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_log_mod,         only: init_logger, final_logger
  use driver_timer_mod,       only: init_timers, final_timers
  use log_mod,                only: log_event, log_level_trace
  use mpi_mod,                only: global_mpi
  use transport_mod,          only: transport_required_namelists
  use transport_driver_mod,   only: initialise_transport, &
                                    run_transport,        &
                                    finalise_transport

  implicit none

  character(*), parameter :: program_name = "transport"

  character(:), allocatable :: filename

  call init_comm( program_name )
  call get_initial_filename( filename )
  call init_config( filename, transport_required_namelists )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_timers( program_name )
  call init_collections()
  deallocate( filename )


  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call initialise_transport( global_mpi, program_name )
  call log_event( 'Running ' // program_name // ' ...', log_level_trace )
  call run_transport()
  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise_transport( program_name )


  call final_collections()
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm()

end program transport
