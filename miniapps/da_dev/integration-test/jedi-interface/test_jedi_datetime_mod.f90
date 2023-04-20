!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Tests the jedi_datetime type init_lfric_calendar_start routine
module test_jedi_datetime_mod

  use constants_mod,                 only : i_def
  use log_mod,                       only : log_event, &
                                            LOG_LEVEL_INFO

  implicit none

  private
  public test_init_lfric_calendar_start,     &
         test_init_lfric_calendar_start_err, &
         test_init_string_err,               &
         test_init_from_jedi_datetime_err

  character(len=128) :: output

contains

  !> @brief Test initialising a datetime from
  !!        the lfric calendar_start namelist variable
  subroutine test_init_lfric_calendar_start()

    use jedi_datetime_mod, only : jedi_datetime_type

    implicit none

    type(jedi_datetime_type) :: jedi_datetime
    integer(i_def)           :: date, returned_date
    integer(i_def)           :: time, returned_time

    date = 2457389_i_def
    time = 54000_i_def

    call jedi_datetime%init_lfric_calendar_start()
    call jedi_datetime%get_date(returned_date)
    call jedi_datetime%get_time(returned_time)

    if ( (returned_date == date) .and. (returned_time == time) ) then
      call log_event( 'test PASS', LOG_LEVEL_INFO )
    else
      write ( output, '(4(A, I10))' ) 'test FAIL with values: date = ',  &
                                      returned_date, ' time = ',         &
                                      returned_time, 'expected date = ', &
                                      date, ' time = ', time
      call log_event( output, LOG_LEVEL_INFO )
    end if

  end subroutine test_init_lfric_calendar_start

  !> @brief Test logging an error initialising a datetime
  !!        from the lfric calendar_start namelist variable
  !!        if it hasn't been read yet
  subroutine test_init_lfric_calendar_start_err()

    use jedi_datetime_mod, only : jedi_datetime_type

    implicit none

    type(jedi_datetime_type) :: jedi_datetime

    call jedi_datetime%init_lfric_calendar_start()

  end subroutine test_init_lfric_calendar_start_err

  !> @brief Tests logging an error when initialising a datetime
  !!        with a bad string
  subroutine test_init_string_err()

    use jedi_datetime_mod, only : jedi_datetime_type

    implicit none

    type(jedi_datetime_type) :: jedi_datetime

    call jedi_datetime%init( '12-34-@5 11:34:23.6' )

  end subroutine test_init_string_err

  !> @brief Tests logging an error when attempting to initialise
  !!        a datetime with another uninitialised datetime
  subroutine test_init_from_jedi_datetime_err()

    use jedi_datetime_mod, only : jedi_datetime_type

    implicit none

    type(jedi_datetime_type) :: jedi_datetime
    type(jedi_datetime_type) :: jedi_datetime_2

    call jedi_datetime_2%init( jedi_datetime )

  end subroutine test_init_from_jedi_datetime_err

end module test_jedi_datetime_mod
