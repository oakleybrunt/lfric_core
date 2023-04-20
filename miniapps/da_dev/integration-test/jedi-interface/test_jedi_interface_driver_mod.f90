!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Drives the execution of the jedi interface tests
module test_jedi_interface_driver_mod

  use test_jedi_datetime_functions_mod, only : test_YYYYMMDD_to_JDN,         &
                                               test_JDN_to_YYYYMMDD_invalid, &
                                               test_hhmmss_to_seconds,       &
                                               test_seconds_to_hhmmss_large, &
                                               test_seconds_to_hhmmss_neg
  use test_jedi_datetime_mod,           only : test_init_lfric_calendar_start,     &
                                               test_init_lfric_calendar_start_err, &
                                               test_init_string_err,               &
                                               test_init_from_jedi_datetime_err

  implicit none

  private
  public test_jedi_interface_init,          &
         test_jedi_interface_final,         &
         run_init_lfric_calendar_start,     &
         run_init_lfric_calendar_start_err, &
         run_init_string_err,               &
         run_init_from_jedi_datetime_err,   &
         run_YYYYMMDD_to_JDN,               &
         run_JDN_to_YYYYMMDD_invalid,       &
         run_hhmmss_to_seconds,             &
         run_seconds_to_hhmmss_large,       &
         run_seconds_to_hhmmss_neg

contains

  !> @brief Initialise testing for the jedi-interface.
  subroutine test_jedi_interface_init()

    implicit none

  end subroutine test_jedi_interface_init

  !> @brief Finalises testing for the jedi-interface
  subroutine test_jedi_interface_final()

    implicit none

  end subroutine test_jedi_interface_final

  !> @brief Runs the init_lfric_calendar_start test
  subroutine run_init_lfric_calendar_start()

    implicit none

    call test_init_lfric_calendar_start()

  end subroutine run_init_lfric_calendar_start

  !> @brief Tests logging an error by init_lfric_calendar_start
  !!        if the calendar_start namelist variable isn't loaded
  subroutine run_init_lfric_calendar_start_err()

    implicit none

    call test_init_lfric_calendar_start_err()

  end subroutine run_init_lfric_calendar_start_err

  !> @brief Tests logging an error when initialising a datetime
  !!        with a bad string
  subroutine run_init_string_err()

    implicit none

    call test_init_string_err()

  end subroutine run_init_string_err

  !> @brief Tests logging an error when attempting to initialise
  !!        a datetime with another uninitialised datetime
  subroutine run_init_from_jedi_datetime_err()

    implicit none

    call test_init_from_jedi_datetime_err()

  end subroutine run_init_from_jedi_datetime_err

  !> @brief Runs the YYYYMMDD_to_JDN test
  subroutine run_YYYYMMDD_to_JDN()

    implicit none

    call test_YYYYMMDD_to_JDN()

  end subroutine run_YYYYMMDD_to_JDN

  !> @brief Runs the JDN_to_YYYYMMDD test with an invalid JDN
  subroutine run_JDN_to_YYYYMMDD_invalid()

    implicit none

    call test_JDN_to_YYYYMMDD_invalid()

  end subroutine run_JDN_to_YYYYMMDD_invalid

  !> @brief Runs the hhmmss_to_seconds test
  subroutine run_hhmmss_to_seconds()

    implicit none

    call test_hhmmss_to_seconds()

  end subroutine run_hhmmss_to_seconds

  !> @brief Runs the hhmmss_to_seconds test with a too large time
  subroutine run_seconds_to_hhmmss_large()

    implicit none

    call test_seconds_to_hhmmss_large()

  end subroutine run_seconds_to_hhmmss_large

  !> @brief Runs the hhmmss_to_seconds test with a negative time
  subroutine run_seconds_to_hhmmss_neg()

    implicit none

    call test_seconds_to_hhmmss_neg()

  end subroutine run_seconds_to_hhmmss_neg

end module test_jedi_interface_driver_mod
