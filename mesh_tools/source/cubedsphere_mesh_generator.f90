!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @mainpage Cubedsphere mesh generator
!> @brief   Utility to generate a cubedsphere surface mesh and write to a file
!>          which conforms to the UGRID format convention.
!> @details Usage:
!>
!>          cubedsphere_mesh_generator <filename>
!>          filename - Controlling namelist file
!>
!-----------------------------------------------------------------------------
program cubedsphere_mesh_generator

  use cli_mod,         only : get_initial_filename
  use constants_mod,   only : i_def, str_def
  use cubedsphere_mesh_generator_config_mod,                            &
                       only : read_cubedsphere_mesh_generator_namelist, &
                              mesh_name, edge_cells, smooth_passes,     &
                              mesh_filename
  use ESMF
  use gencube_ps_mod,  only : gencube_ps_type
  use io_utility_mod,  only : open_file, close_file
  use log_mod,         only : log_scratch_space, log_event, &
                              LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use ncdf_quad_mod,   only : ncdf_quad_type
  use ugrid_file_mod,  only : ugrid_file_type
  use ugrid_2d_mod,    only : ugrid_2d_type

  implicit none

  type(ESMF_VM)  :: vm
  integer(i_def) :: rc

  character(:), allocatable :: filename
  integer(i_def)            :: namelist_unit

  type(gencube_ps_type)               :: csgen
  type(ugrid_2d_type)                 :: ugrid_2d
  class(ugrid_file_type), allocatable :: ugrid_file
  integer(i_def)                      :: fsize

  ! Start up ESMF
  call ESMF_Initialize(vm=vm, defaultlogfilename="cubedsphere.log", &
                       logkindflag=ESMF_LOGKIND_SINGLE, rc=rc)
  if (rc /= ESMF_SUCCESS) call log_event( 'Failed to initialise ESMF.', &
                                          LOG_LEVEL_ERROR )

  ! Read mesh generation namelist from file
  call get_initial_filename( filename )
  namelist_unit = open_file( filename )
  call read_cubedsphere_mesh_generator_namelist( namelist_unit, vm, 0 )
  call close_file( namelist_unit )
  deallocate( filename )

  ! Create object to manipulate UGRID conforming NetCDF file
  allocate(ncdf_quad_type::ugrid_file)
  call ugrid_2d%set_file_handler(ugrid_file)

  ! Create object which can generate the cubedsphere mesh from
  ! specified inputs.
  csgen = gencube_ps_type( mesh_name, edge_cells, smooth_passes )
  call log_event( 'Generating cubed-sphere mesh, "' // trim(mesh_name)// &
                  '", with...', LOG_LEVEL_INFO )

  write(log_scratch_space, "(A,T17,I0)") "  Edge cells:",   edge_cells
  call log_event( log_scratch_space, LOG_LEVEL_INFO )
  write(log_scratch_space, "(A,T17,I0)") "  Smooth passes:", smooth_passes
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

  ! Mesh for the UGRID conforming NetCDF file is
  ! set by the generator passed to it
  call ugrid_2d%set_by_generator( csgen )

  ! Now the write out mesh to the NetCDF file
  call ugrid_2d%write_to_file(trim( mesh_filename) )

  call log_event( "...generation complete.", LOG_LEVEL_INFO )
  inquire(file=mesh_filename, size=fsize)
  write( log_scratch_space, '(A,I0,A)' )                       &
      'Writing ugrid mesh to '//trim(adjustl(mesh_filename))// &
      ' - ', fsize, ' bytes written.'
  call log_event( log_scratch_space, LOG_LEVEL_INFO )

end program cubedsphere_mesh_generator
