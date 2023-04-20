!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   A module providing JEDI/OOPS like datetime functions
!>
!> @details This module emulates the datetime functions used in the
!!          LFRic-JEDI interface and JEDI's OOPS
module jedi_datetime_functions_mod

  use constants_mod,           only : i_def, r_def, str_def, l_def
  use log_mod,                 only : log_event,         &
                                      log_scratch_space, &
                                      LOG_LEVEL_ERROR

  implicit none
  private

  public YYYYMMDD_to_JDN,   &
         JDN_to_YYYYMMDD,   &
         hhmmss_to_seconds, &
         seconds_to_hhmmss, &
         is_valid_datetime

! Make procedures public for unit / int testing
#ifdef UNIT_TEST
  public  is_valid_YYYYMMDD, &
          is_leap_year,      &
          is_valid_hhmmss
#elif defined INT_TEST
  public  is_valid_YYYYMMDD, &
          is_leap_year,      &
          is_valid_hhmmss
#else
  private is_valid_YYYYMMDD, &
          is_leap_year,      &
          is_valid_hhmmss
#endif

contains

  !> @brief   Converts year, month, day to a Julian Day Number (UTC)
  !> @details Uses the following formula from https://doi.org/10.1145/364096.364097
  !!          mirroring JEDI/OOPS https://github.com/JCSDA-internal/oops/blob/develop/src/oops/util/dateFunctions.cc
  !!          julian_day = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +
  !!                       ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -
  !!                       ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 +
  !!                       d - 32075
  !!
  !> @param [in] year  Year - YYYY
  !> @param [in] month Month - MM
  !> @param [in] day   Day - DD
  !!
  !> @result [inout] julian_day_number Julian Day Number
  subroutine YYYYMMDD_to_JDN( year, month, day, julian_day_number )

    implicit none

    integer(i_def), intent(in)    :: year    !< in YYYY format, eg 2020
    integer(i_def), intent(in)    :: month   !< in MM format, 01 (Jan) through 12 (Dec)
    integer(i_def), intent(in)    :: day     !< in DD format, 01 through 28, 29, 30, or 31
    integer(i_def), intent(inout) :: julian_day_number

    integer(i_def) :: l

    if ( .not. is_valid_YYYYMMDD( year, month, day ) ) then
      write ( log_scratch_space, '(A, I4, 2(A, I2))' ) &
        'Not a valid date: YYYY=', year, ' MM=', month, ' DD=', day
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    l = ( month - 14_i_def ) / 12_i_def

    julian_day_number = 1461_i_def * ( year + 4800_i_def + l ) / 4_i_def              &
                        + 367_i_def * ( month - 2_i_def - 12_i_def * l ) / 12_i_def   &
                        - 3_i_def * (( year + 4900_i_def + l ) / 100_i_def) / 4_i_def &
                        + day - 32075_i_def

  end subroutine YYYYMMDD_to_JDN

  !> @brief   Converts Julian Day Number (UTC) to year, month, day
  !> @details Uses the formula from https://doi.org/10.1145/364096.364097,
  !!          mirroring JEDI/OOPS https://github.com/JCSDA-internal/oops/blob/develop/src/oops/util/dateFunctions.cc
  !!
  !> @param [in] julian_day_number Julian Day Number
  !!
  !> @result [inout] year  Year - YYYY
  !> @result [inout] month Month - MM
  !> @result [inout] day   Day - DD
  subroutine JDN_to_YYYYMMDD( julian_day_number, year, month, day )

    implicit none

    integer(i_def), intent(in)    :: julian_day_number
    integer(i_def), intent(inout) :: year    !< in YYYY format, eg 2020
    integer(i_def), intent(inout) :: month   !< in MM format, 01 (Jan) through 12 (Dec)
    integer(i_def), intent(inout) :: day     !< in DD format, 01 through 28, 29, 30, or 31

    integer(i_def) :: l, n, i, j

    if ( julian_day_number < 0_i_def ) then
      write ( log_scratch_space, '(A)' ) &
        'Invalid Julian Date Passed (negative value)'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    l = julian_day_number + 68569_i_def
    n = ( 4_i_def * l ) / 146097_i_def
    l = l - (( (146097_i_def * n) + 3_i_def ) / 4_i_def)
    i = ( 4000_i_def * (l + 1_i_def) ) / 1461001_i_def
    l = l - (( 1461_i_def * i ) / 4_i_def) + 31_i_def
    j = ( 80_i_def * l ) / 2447_i_def
    day = l - (( 2447_i_def * j ) / 80_i_def)
    l = j / 11_i_def
    month = j + 2_i_def - ( 12_i_def * l )
    year = 100_i_def * ( n - 49_i_def ) + i + l

    if ( year >= huge(year) ) then
      write ( log_scratch_space, '(A)' ) &
        'Invalid Julian Date Passed'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end subroutine JDN_to_YYYYMMDD

  !> @brief Converts hours, minutes, and seconds to seconds (UTC)
  !!
  !> @param [in] hour   Hours - hh
  !> @param [in] minute Minutes - mm
  !> @param [in] second Seconds - ss
  !!
  !> @result [inout] time_seconds Time in seconds since the start of the day
  subroutine hhmmss_to_seconds( hour, minute, second, time_seconds )

    implicit none

    integer(i_def), intent(in)    :: hour    !< in hh format, 00 through 23
    integer(i_def), intent(in)    :: minute  !< in mm format, 00 through 59
    integer(i_def), intent(in)    :: second  !< in ss format, 00 through 59
    integer(i_def), intent(inout) :: time_seconds

    if ( .not. is_valid_hhmmss( hour, minute, second ) ) then
      write ( log_scratch_space, '(3(A, I2))' ) &
        'Not a valid time: hh=', hour, ' mm=', minute, ' ss=', second
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    time_seconds = hour * 3600_i_def + minute * 60_i_def + second

  end subroutine hhmmss_to_seconds

  !> @brief Converts seconds to hours, minutes, and seconds
  !!
  !> @param [in] time_seconds Time in seconds since the start of the day
  !!
  !> @param [inout] hour   Hours - hh
  !> @param [inout] minute Minutes - mm
  !> @param [inout] second Seconds - ss
  subroutine seconds_to_hhmmss( time_seconds, hour, minute, second )

    implicit none

    integer(i_def), intent(in)    :: time_seconds
    integer(i_def), intent(inout) :: hour    !< in hh format, 00 through 23
    integer(i_def), intent(inout) :: minute  !< in mm format, 00 through 59
    integer(i_def), intent(inout) :: second  !< in ss format, 00 through 59

    integer(i_def), parameter :: seconds_per_day = 86400_i_def
    integer(i_def)            :: local_sec

    if ( (time_seconds >= seconds_per_day) .or. (time_seconds < 0_i_def) ) then
      write ( log_scratch_space, '(A)' ) &
        'Not a valid time in seconds, 0 =< time =< 86400'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    local_sec = time_seconds
    hour = local_sec / 3600_i_def
    local_sec = mod( local_sec, 3600_i_def )
    minute = local_sec / 60_i_def
    local_sec = mod( local_sec, 60_i_def )
    second = local_sec

  end subroutine seconds_to_hhmmss

  !> @brief Tests if the year, month, day combination is valid
  !!
  !> @param [in] year  Year - YYYY
  !> @param [in] month Month - MM
  !> @param [in] day   Day - DD
  !!
  !> @result valid    .true. if the year, month, day combination is valid
  pure function is_valid_YYYYMMDD( year, month, day ) result(valid)

    implicit none

    integer(i_def), intent(in) :: year
    integer(i_def), intent(in) :: month
    integer(i_def), intent(in) :: day

    logical(l_def) :: valid
    valid = .false.

    if ( (year == 0_i_def) .and. (month == 0_i_def) .and. (day == 0_i_def) ) then
      valid = .true.
    else if ( (year >= 0_i_def)  .and. (year <= 9999_i_def) .and. &
              (month >= 1_i_def) .and. (month <= 12_i_def)  .and. &
              (day >= 1_i_def)   .and. (day <= 31_i_def) ) then

      if ( (month == 4_i_def) .or. (month == 6_i_def) .or. &
           (month == 9_i_def) .or. (month == 11_i_def) ) then
        if ( day <= 30_i_def ) valid = .true.
      else if ( month /= 2_i_def ) then
        if ( day <= 31_i_def ) valid = .true.
      else if ( is_leap_year(year) ) then
        if ( day <= 29_i_def ) valid = .true.
      else
        if ( day <= 28_i_def ) valid = .true.
      end if

    end if

  end function is_valid_YYYYMMDD

  !> @brief Tests if a year is a leap year
  !!
  !> @param [in] year          Year - YYYY
  !> @result     is_leap_year  .true. if the year is a leap year
  pure function is_leap_year( year )

    implicit none

    integer(i_def), intent(in) :: year

    logical(l_def) :: is_leap_year
    is_leap_year = .false.

    if ( ((mod(year, 4_i_def)==0_i_def) .and. ((mod(year, 100_i_def)/=0_i_def))) &
         .or. (mod(year, 400_i_def)==0_i_def)) then
      is_leap_year = .true.
    end if

  end function is_leap_year

  !> @brief Tests if the hour, minute, second combination is valid
  !!
  !> @param [in] hour   Hours - hh
  !> @param [in] minute Minutes - mm
  !> @param [in] second Seconds - ss
  !!
  !> @result valid .true. if the hour, minute, second combination is valid
  pure function is_valid_hhmmss( hour, minute, second ) result(valid)

    implicit none

    integer(i_def), intent(in) :: hour
    integer(i_def), intent(in) :: minute
    integer(i_def), intent(in) :: second

    logical(l_def) :: valid
    valid = .false.

    if ( (hour >= 0_i_def)   .and. (hour <= 23_i_def)   .and. &
         (minute >= 0_i_def) .and. (minute <= 59_i_def) .and. &
         (second >= 0_i_def) .and. (second <= 59_i_def) ) then
      valid = .true.
    end if

  end function is_valid_hhmmss

  !> @brief Tests if the passed date and time are valid
  !!
  !> @param [in] date Julian Day Number
  !> @param [in] time Seconds since start of day
  !!
  !> @result valid .true. if the date and time combination is valid
  function is_valid_datetime( date, time )

    implicit none

    integer(i_def), intent(in) :: date, time

    integer(i_def) :: year
    integer(i_def) :: month
    integer(i_def) :: day

    integer(i_def) :: hour
    integer(i_def) :: minute
    integer(i_def) :: second

    logical(l_def) :: valid_date, valid_time
    logical(l_def) :: is_valid_datetime

    is_valid_datetime = .false.

    call JDN_to_YYYYMMDD( date, year, month, day )
    call seconds_to_hhmmss( time, hour, minute, second )

    valid_date = is_valid_YYYYMMDD( year, month, day )
    valid_time = is_valid_hhmmss( hour, minute, second )

    if ( valid_date .and. valid_time ) is_valid_datetime = .true.

  end function is_valid_datetime

end module jedi_datetime_functions_mod
