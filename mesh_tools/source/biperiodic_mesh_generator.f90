!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Biperiodic mesh generator
!>
!> @brief   Utility to generate a biperiodic surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          biperiodic_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program biperiodic_mesh_generator

  use biperiodic_mesh_generator_config_mod,                              &
                         only : read_biperiodic_mesh_generator_namelist, &
                                mesh_name, cells_in_x, cells_in_y,       &
                                cell_width, cell_height,                 &
                                mesh_filename
  use cli_mod,           only : get_initial_filename
  use constants_mod,     only : i_def, str_def
  use ESMF
  use genbiperiodic_mod, only : genbiperiodic_type
  use io_utility_mod,    only : open_file, close_file
  use log_mod,           only : log_scratch_space, log_event, &
                                LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use ncdf_quad_mod,     only : ncdf_quad_type
  use ugrid_2d_mod,      only : ugrid_2d_type
  use ugrid_file_mod,    only : ugrid_file_type

  implicit none

  type(ESMF_VM)  :: vm
  integer(i_def) :: rc

  character(:), allocatable :: filename
  integer(i_def)            :: namelist_unit

  type(genbiperiodic_type)            :: bpgen
  type(ugrid_2d_type)                 :: ugrid_2d
  class(ugrid_file_type), allocatable :: ugrid_file
  integer(i_def)                      :: fsize

  character(str_def) :: rchar
  character(str_def) :: fmt_str

  ! Start up ESMF
  call ESMF_Initialize( vm=vm, defaultlogfilename="biperiodic.log", &
                        logkindflag=ESMF_LOGKIND_SINGLE, rc=rc )
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', &
                                          LOG_LEVEL_ERROR )

  ! Read mesh generation namelist from file
  call get_initial_filename( filename )
  namelist_unit = open_file( filename )
  call read_biperiodic_mesh_generator_namelist( namelist_unit, vm, 0 )
  call close_file( namelist_unit )
  deallocate( filename )

  ! Create object to manipulate UGRID conforming NetCDF file
  allocate(ncdf_quad_type::ugrid_file)
  call ugrid_2d%set_file_handler(ugrid_file)

  ! Create object which can generate the biperiodic mesh from
  ! specified inputs.
  bpgen = genbiperiodic_type( mesh_name,              &
                              cells_in_x, cells_in_y, &
                              cell_width, cell_height )

  ! Mesh for the UGRID conforming NetCDF file is
  ! set by the generator passed to it
  call log_event( 'Generating biperiodic mesh, "' // trim(mesh_name)// &
                  '", with...', LOG_LEVEL_INFO )

  fmt_str = '(A,T16,I0)'

  write(log_scratch_space, fmt_str) "  Cells in x:", cells_in_x
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  write(log_scratch_space, fmt_str) "  Cells in y:", cells_in_y
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  fmt_str = '(A,T16,A)'

  write(rchar, '(F6.1)') cell_width
  write(log_scratch_space, fmt_str) &
      "  Cell width:",trim(adjustl(rchar))
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  write(rchar, '(F6.1)') cell_height
  write(log_scratch_space, fmt_str) &
      "  Cell height:",trim(adjustl(rchar))
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  ! Mesh for the UGRID conforming NetCDF file is
  ! set by the generator passed to it
  call ugrid_2d%set_by_generator( bpgen )

  ! Now the write out mesh to the NetCDF file
  call ugrid_2d%write_to_file( trim(mesh_filename) )

  call log_event( "...generation complete.", LOG_LEVEL_INFO )
  inquire(file=mesh_filename, size=fsize)
  write( log_scratch_space, '(A,I0,A)' )                       &
      'Writing ugrid mesh to '//trim(adjustl(mesh_filename))// &
      ' - ', fsize, ' bytes written.'
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  call ESMF_Finalize(rc=rc)

end program biperiodic_mesh_generator
