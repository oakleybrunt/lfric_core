!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!   Program to generate a biperiodic mesh and write this in ugrid format to
!   the specified file.
!   Passes command-line arguments to genbiperiodic_type and uses ncdf_quad_mod
!   to write resulting mesh.
!   Invocation without arguments, or omission of any one or more
!   argument, leads to the default output of:
!     -o ugrid_quads_2d.nc -nx 5 -ny 4 -dx 6000.0 -d 2000.0
!-------------------------------------------------------------------------------
program generate_biperiodic
!-------------------------------------------------------------------------------
use genbiperiodic_mod,       only : genbiperiodic_type
use generate_biperiodic_mod, only : parse_args
use ugrid_2d_mod,            only : ugrid_2d_type
use ugrid_file_mod,          only : ugrid_file_type
use ncdf_quad_mod,           only : ncdf_quad_type
use constants_mod,           only : i_def, r_def, str_def
use iso_fortran_env,         only : stdout => output_unit

implicit none
!-------------------------------------------------------------------------------
  type(genbiperiodic_type)               :: bpgen
  type(ugrid_2d_type)                    :: ugrid_2d
  class(ugrid_file_type), allocatable    :: ugrid_file
  character(len=str_def)                 :: filename, sztext
  integer(kind=i_def)                    :: nx, ny
  real(kind=r_def)                       :: dx, dy
  integer                                :: fsize

  call parse_args(filename, nx, ny, dx, dy)

  allocate(ncdf_quad_type::ugrid_file)
  call ugrid_2d%set_file_handler(ugrid_file)

  bpgen = genbiperiodic_type(nx, ny, dx ,dy)

  write(stdout, "(A)") "Generating biperiodic mesh with..."
  write(stdout, "(A,I5)") "  nx: ", nx
  write(stdout, "(A,I5)") "  ny: ", ny
  write(stdout, "(A,F6.1)") "  dx: ", dx
  write(stdout, "(A,F6.1)") "  dy: ", dy

  call ugrid_2d%set_by_generator(bpgen)
  write(stdout, "(A)") "...generation complete."

  write(stdout, "(A)", advance="NO") "Writing ugrid mesh to "//trim(adjustl(filename))//" ..."
  call ugrid_2d%write_to_file(trim(filename))
  inquire(file=filename, size=fsize)
  write(sztext, *) fsize
  write(stdout, "(A)") "... "//trim(adjustl(sztext))//" bytes written."

  stop

end program generate_biperiodic
