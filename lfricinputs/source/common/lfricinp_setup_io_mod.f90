! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_setup_io_mod

USE clock_mod,                     ONLY: clock_type
USE constants_mod,                 ONLY: i_def, str_max_filename
USE lfric_xios_file_mod,           ONLY: xios_file_type, &
                                         append_file_to_list
! Configuration modules
USE files_config_mod,              ONLY: ancil_directory,                      &
                                         checkpoint_stem_name,                 &
                                         land_area_ancil_path
USE io_config_mod,                 ONLY: checkpoint_frequency,                 &
                                         diagnostic_frequency,                 &
                                         checkpoint_write,                     &
                                         checkpoint_read,                      &
                                         write_diag

IMPLICIT NONE

CHARACTER(LEN=str_max_filename) :: checkpoint_write_fname = 'unset'
CHARACTER(LEN=str_max_filename) :: checkpoint_read_fname = 'unset'
CHARACTER(LEN=str_max_filename) :: ancil_fname = 'unset'

PRIVATE
PUBLIC :: checkpoint_write_fname, checkpoint_read_fname, ancil_fname,          &
          init_lfricinp_files

CONTAINS

SUBROUTINE init_lfricinp_files(files_list, clock)

! lfricinp modules
USE lfricinp_ancils_mod, ONLY: l_land_area_fraction

IMPLICIT NONE

type(xios_file_type), ALLOCATABLE, INTENT(OUT) :: files_list(:)
CLASS(clock_type),                 INTENT(IN)  :: clock

TYPE(xios_file_type) :: tmp_file

! Setup diagnostic output file
IF (write_diag) THEN
  CALL tmp_file%init_xios_file("lfric_diag", freq=diagnostic_frequency)
  call append_file_to_list(tmp_file, files_list)
END IF

IF (l_land_area_fraction) THEN
  ! Set land area ancil filename from namelist
  WRITE(ancil_fname,'(A)') TRIM(ancil_directory)//'/'//                        &
                           TRIM(land_area_ancil_path)
  CALL tmp_file%init_xios_file("land_area_ancil", path=ancil_fname)
  call append_file_to_list(tmp_file, files_list)
END IF

! Setup checkpoint writing context information
IF (checkpoint_write) THEN
  ! Create checkpoint filename from stem and first timestep
  WRITE(checkpoint_write_fname,'(A,A,I6.6)')                                   &
                       TRIM(checkpoint_stem_name),"_", clock%get_first_step()

  CALL tmp_file%init_xios_file("lfric_checkpoint_write",                       &
                               checkpoint_write_fname,                         &
                               freq=checkpoint_frequency)
  call append_file_to_list(tmp_file, files_list)
END IF

! Setup checkpoint reading context information
IF (checkpoint_read) THEN
  ! Create checkpoint filename from stem
  WRITE(checkpoint_read_fname,'(A)') TRIM(checkpoint_stem_name)

  CALL tmp_file%init_xios_file("lfric_checkpoint_read",                       &
                               checkpoint_read_fname)
  call append_file_to_list(tmp_file, files_list)
END IF

END SUBROUTINE init_lfricinp_files

END MODULE lfricinp_setup_io_mod