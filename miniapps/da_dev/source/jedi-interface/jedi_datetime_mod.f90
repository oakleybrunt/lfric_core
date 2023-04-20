!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   A module providing JEDI/OOPS like date time functionality
!>
!> @details This module emulates the date time object used in the
!!          LFRic-JEDI interface
module jedi_datetime_mod

  use constants_mod,               only : i_def, r_def, str_def, l_def
  use jedi_datetime_functions_mod, only : YYYYMMDD_to_JDN,   &
                                          JDN_to_YYYYMMDD,   &
                                          hhmmss_to_seconds, &
                                          seconds_to_hhmmss, &
                                          is_valid_datetime
  use log_mod,                     only : log_event,         &
                                          log_scratch_space, &
                                          LOG_LEVEL_INFO,    &
                                          LOG_LEVEL_ERROR
  use time_config_mod,             only : calendar_start

  implicit none
  private

  !> @brief This type stores the JEDI date and time
  type, public :: jedi_datetime_type
    private

    integer(i_def) :: date    !< Julian day number
    integer(i_def) :: time    !< seconds since start of day

  contains

    procedure, public  :: init_lfric_calendar_start
    procedure, private :: init_iso_string
    procedure, private :: init_string
    procedure, private :: init_YYMMDD_hhmmss
    procedure, private :: init_YYYYMMDDhhmmss
    procedure, private :: init_from_jedi_datetime
    generic,   public  :: init => init_iso_string,     &
                                  init_YYMMDD_hhmmss,  &
                                  init_YYYYMMDDhhmmss, &
                                  init_from_jedi_datetime

    procedure, public  :: get_date
    procedure, public  :: get_time
    procedure, public  :: add_seconds
    procedure, public  :: seconds_between
    procedure, public  :: is_ahead
    procedure, public  :: to_string

  end type jedi_datetime_type

contains

  !> @brief Initialise a datetime using the lfric
  !!        calendar_start namelist variable
  subroutine init_lfric_calendar_start( self )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self

    call log_event( 'Creating JEDI datetime using LFRic calendar_start', LOG_LEVEL_INFO )

    if ( calendar_start == 'unset' ) then
      write ( log_scratch_space, '(A)' ) 'Creating JEDI datetime FAIL: &
                                        & calendar_start has not been loaded'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    call self%init_string( calendar_start )

  end subroutine init_lfric_calendar_start

  !> @brief   Initialise a datetime using an iso string (UTC)
  !> @details Initialise a jedi datetime using an iso string
  !!          of the form 2023-04-05T11:41:38Z. The T may be
  !!          a space and the timezone, Z is ignored
  !!
  !> @param [in] iso_datetime ISO datetime string
  subroutine init_iso_string( self, iso_datetime )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    character(*),                intent(in)    :: iso_datetime

    write ( log_scratch_space, '(2(A))' ) &
      'Creating JEDI datetime using iso datetime string: ', iso_datetime
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call self%init_string( iso_datetime )

  end subroutine init_iso_string

  !> @brief   Initialise a datetime using an iso string (UTC)
  !!
  !> @param [in] iso_datetime ISO datetime string
  subroutine init_string( self, iso_datetime )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    character(*),                intent(in)    :: iso_datetime

    integer(i_def) :: year
    integer(i_def) :: month
    integer(i_def) :: day

    integer(i_def) :: hour
    integer(i_def) :: minute
    integer(i_def) :: second

    integer(i_def) :: err

    read ( iso_datetime(1:4), '(I4)', iostat=err ) year
    read ( iso_datetime(6:7), '(I2)', iostat=err ) month
    read ( iso_datetime(9:10), '(I2)', iostat=err  ) day

    read ( iso_datetime(12:13), '(I2)', iostat=err  ) hour
    read ( iso_datetime(15:16), '(I2)', iostat=err  ) minute
    read ( iso_datetime(18:19), '(I2)', iostat=err  ) second

    if ( err /= 0_i_def) then
      write ( log_scratch_space, '(A)' ) &
              'Creating JEDI datetime FAIL: Failed to read iso_datetime string'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    call self%init_YYYYMMDDhhmmss( year, month, day, hour, minute, second )

  end subroutine init_string

  !> @brief   Initialise a datetime with two integers,
  !!          one for the date and one for time
  !> @details Initialise a datetime with a YYYYMMDD integer
  !!          for the year, month, and day eg. 20230405,
  !!          and a hhmmss integer for the hours, minutes,
  !!          and seconds eg. 130210.
  !!
  !> @param [in] YYYYMMDD Integer representing the year, month, day
  !> @param [in] hhmmss   Integer representing the hour, minute, second
  subroutine init_YYMMDD_hhmmss( self, YYYYMMDD, hhmmss )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_def), intent(in) :: YYYYMMDD !< year, month, and day eg 20230405
    integer(i_def), intent(in) :: hhmmss   !< hour, minute, second eg 130210

    integer(i_def) :: year
    integer(i_def) :: month
    integer(i_def) :: day

    integer(i_def) :: hour
    integer(i_def) :: minute
    integer(i_def) :: second

    integer(i_def) :: temp_int

    year     = YYYYMMDD / 10000_i_def       ! remove lower 4 digits, yields YYYY
    temp_int = mod( YYYYMMDD, 10000_i_def)  ! keep lower 4 digits, yields MMDD
    month    = temp_int / 100_i_def         ! remove lower 2 digits, yields MM
    day      = mod(temp_int, 100_i_def)     ! keep lower 2 digits, yields DD

    hour     = int(hhmmss, i_def) / 10000_i_def       ! remove lower 4 digits, yields hh
    temp_int = mod( int(hhmmss, i_def), 10000_i_def)  ! keep lower 4 digits, yields mmss
    minute   = temp_int / 100_i_def                   ! remove lower 2 digits, yields mm
    second   = mod( temp_int, 100_i_def )             ! keep lower 2 digits, yeilds ss

    call self%init_YYYYMMDDhhmmss( year, month, day, hour, minute, second )

  end subroutine init_YYMMDD_hhmmss

  !> @brief      Initialise a datetime with six integers for the year, month,
  !!             day, hour, minute, and second
  !!
  !> @param [in] year   Year   - YYYY
  !> @param [in] month  Month  - MM
  !> @param [in] day    Day    - DD
  !> @param [in] hour   Hour   - hh
  !> @param [in] minute Minute - mm
  !> @param [in] second Second - ss
  subroutine init_YYYYMMDDhhmmss( self, year, month, day, hour, minute, second )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_def), intent(in) :: year    !< in YYYY format, eg 2020
    integer(i_def), intent(in) :: month   !< in MM format, 01 (Jan) through 12 (Dec)
    integer(i_def), intent(in) :: day     !< in DD format, 01 through 28, 29, 30, or 31

    integer(i_def), intent(in) :: hour    !< in hh format, 00 through 23
    integer(i_def), intent(in) :: minute  !< in mm format, 00 through 59
    integer(i_def), intent(in) :: second  !< in ss format, 00 through 59

    call log_event( 'Initialising JEDI datetime', LOG_LEVEL_INFO )

    call YYYYMMDD_to_JDN( year, month, day, self%date )
    call hhmmss_to_seconds( hour, minute, second, self%time )

    write ( log_scratch_space, '(A, I8, A, I5)' )                   &
           'Initialising JEDI datetime SUCCESS: JDN = ', self%date, &
           ' Time (s) = ', self%time
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine init_YYYYMMDDhhmmss

  !> @brief Initialise a datetime from another datetime instance
  !!
  !> @param [inout] datetime A jedi_datetime instance to copy from
  subroutine init_from_jedi_datetime( self, datetime )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    type( jedi_datetime_type ),  intent(inout) :: datetime

    integer(i_def) :: date_to_copy, time_to_copy
    logical(l_def) :: valid

    call datetime%get_date( date_to_copy )
    call datetime%get_time( time_to_copy )

    valid = is_valid_datetime( date_to_copy, time_to_copy )

    if ( .not. valid ) then
      write ( log_scratch_space, '(A)' ) &
        'Cannot initialise this datetime from an invalid / unitiliased datetime'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    self%date = date_to_copy
    self%time = time_to_copy

  end subroutine init_from_jedi_datetime

  !> @brief Gets the date from the datetime instance
  !!
  !> @param [inout] date The Julian Day Number
  subroutine get_date( self, date )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_def),              intent(inout) :: date

    date = self%date

  end subroutine get_date

  !> @brief Gets the time from the datetime instance
  !!
  !> @param [inout] time Seconds since the start of the day
  subroutine get_time( self, time )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_def),              intent(inout) :: time

    time = self%time

  end subroutine get_time

  !> @brief Adds a time in seconds to the datetime instance
  !!
  !> @param [in] seconds Time in seconds to add, can be negative
  subroutine add_seconds( self, seconds )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_def),              intent(in)    :: seconds

    integer(i_def) :: days
    integer(i_def) :: new_time
    integer(i_def), parameter  :: seconds_per_day = 86400_i_def

    new_time = seconds
    days = new_time / seconds_per_day

    self%date = self%date + days
    new_time = mod( new_time, seconds_per_day )

    self%time = self%time + new_time

    if ( self%time < 0_i_def ) then
      self%date = self%date - 1_i_def
      self%time = self%time + seconds_per_day
    else if ( self%time >= seconds_per_day ) then
      self%date = self%date + 1_i_def
      self%time = self%time - seconds_per_day
    end if

  end subroutine

  !> @brief Calculates the number of seconds between two datetimes
  !!
  !> @param [in]    datetime     The jedi_datetime instance to compare with
  !> @param [inout] diff_seconds The difference between the two datetimes
  subroutine seconds_between( self, datetime, diff_seconds )

    implicit none

    class( jedi_datetime_type ), intent(in)    :: self
    class( jedi_datetime_type ), intent(in)    :: datetime
    integer(i_def),              intent(inout) :: diff_seconds

    integer(i_def), parameter :: seconds_per_day = 86400_i_def
    integer(i_def)            :: diff_date
    integer(i_def)            :: diff_time

    diff_date = datetime%date - self%date
    diff_time = datetime%time - self%time

    diff_seconds = diff_date * seconds_per_day + diff_time

  end subroutine seconds_between

  !> @brief Returns true if the datetime is ahead in time of the passed datetime
  !!
  !> @param [in] datetime The jedi_datetime instance to compare with
  function is_ahead( self, datetime )

    implicit none

    class(jedi_datetime_type), intent(in) :: self
    type(jedi_datetime_type),  intent(in) :: datetime

    logical(l_def) :: is_ahead

    integer(i_def) :: seconds_ahead

    is_ahead = .false.

    call self%seconds_between( datetime, seconds_ahead )
    if ( seconds_ahead < 0 ) is_ahead = .true.

  end function is_ahead

  !> @brief Returns the datetime as an iso string (UTC)
  !!
  !> @param [inout] iso_datetime The string to return
  subroutine to_string( self, iso_datetime )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    character(str_def),          intent(inout) :: iso_datetime

    integer(i_def) :: year
    integer(i_def) :: month
    integer(i_def) :: day

    integer(i_def) :: hour
    integer(i_def) :: minute
    integer(i_def) :: second

    character(len=4) :: temp_str_4
    character(len=2) :: temp_str_2
    character(len=1) :: dash, space, colon

    call JDN_to_YYYYMMDD( self%date, year, month, day )
    call seconds_to_hhmmss( self%time, hour, minute, second )

    dash  = '-'
    space = ' '
    colon = ':'

    write ( temp_str_4, '(I4)' ) year
    iso_datetime = temp_str_4 // dash
    write ( temp_str_2, '(I2.2)' ) month
    iso_datetime = trim(iso_datetime) // temp_str_2 // dash
    write ( temp_str_2, '(I2.2)' ) day
    iso_datetime = trim(iso_datetime) // temp_str_2

    write ( temp_str_2, '(I2.2)' ) hour
    iso_datetime = trim(iso_datetime) // space // temp_str_2 // colon
    write ( temp_str_2, '(I2.2)' ) minute
    iso_datetime = trim(iso_datetime) // temp_str_2 // colon
    write ( temp_str_2, '(I2.2)' ) second
    iso_datetime = trim(iso_datetime) // temp_str_2

  end subroutine to_string

end module jedi_datetime_mod
