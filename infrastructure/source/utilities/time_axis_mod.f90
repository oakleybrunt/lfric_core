!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief A module containing a time axis object
!>
!> @details Fields need to be updated at different times and frequencies. The
!>          time axis object can be linked to a field to provide information
!>          of how the field should be updated with time.
module time_axis_mod

  use constants_mod,        only: i_def, str_def, r_def, l_def
  use field_collection_mod, only: field_collection_type
  use log_mod,              only: log_event, log_scratch_space, LOG_LEVEL_ERROR

  implicit none

  private

  !> Time axis object type
  type, public :: time_axis_type

    private

    !> Name of the time_axis.
    character(str_def) :: name = 'unset'
    !> The data values of the time axis
    real(kind=r_def), allocatable :: time_data( : )
    !> The indices of the data points of the time axis
    integer(kind=i_def), allocatable :: index_data( : )
    !> The collection of fields associated with the time axis
    type(field_collection_type) :: fields
    !> Flag determining the cyclincal nature of time axis
    logical(l_def) :: cyclic = .true.
    !> String identifier for units of time axis
    character(str_def) :: time_units
  contains
    !> Initialiser for a field parent object
    procedure, public :: initialise
    !> Getter for time axis name
    procedure, public :: get_name
    !> Getter for time units
    procedure, public :: get_units
    !> Procedure for cycling through time data to find the correct entry
    procedure, public :: shift_forward
    !> Returns start/end times for active time window
    procedure, public :: get_time_window
    !> Returns start/end indices for active time window
    procedure, public :: get_index_window
    !> Aligns the active time window with the current model time
    procedure, public :: align
    !> Returns in time axis is cyclic or not
    procedure, public :: is_cyclic
  end type time_axis_type

contains

  !> Initialise a <code>time_axis_type</code> object.
  !>
  !> @param [in] input_data The input time data array
  !> @param [in] input_index_data The corresponding indices for the data
  !> @param [in] name The time axis name
  !> @param [in] input_fields The field collection to contain associated fields
  !> @param [in] cyclic Flag determining if axis is cyclic
  subroutine initialise(self, input_data, input_index_data, name, input_fields, input_units, cyclic)

    implicit none

    class(time_axis_type),        intent(inout) :: self
    real(r_def),    dimension(:), intent(in)    :: input_data
    integer(i_def), dimension(:), intent(in)    :: input_index_data
    character(*),                 intent(in)    :: name
    type(field_collection_type),  intent(in)    :: input_fields
    character(*),                 intent(in)    :: input_units
    logical(l_def), optional,     intent(in)    :: cyclic

    self%time_data = input_data

    self%index_data = input_index_data

    self%name = name

    self%fields = input_fields

    self%time_units = input_units

    if ( present(cyclic) ) self%cyclic = cyclic

  end subroutine initialise

  !> Returns start and end time data for active time window
  !> @result output_name The time axis name
  function get_name(self) result(output_name)

    implicit none

    class(time_axis_type), intent(inout) :: self

    character(str_def) :: output_name

    output_name = self%name

  end function get_name

  !> Returns start and end time data for active time window
  !> @result output_unitsThe time axis units
  function get_units(self) result(output_units)

    implicit none

    class(time_axis_type), intent(inout) :: self

    character(str_def) :: output_units

    output_units = self%time_units

  end function get_units

  !> Performs a cshift on the data and index data arrays.
  subroutine shift_forward(self)

    implicit none

    class(time_axis_type), intent(inout) :: self

    self%time_data = cshift(self%time_data, 1)
    self%index_data = cshift(self%index_data, 1)

  end subroutine shift_forward

  !> Returns start and end time data for active time window
  !> @result time_window The active time window of the time axis
  function get_time_window(self) result(time_window)

    implicit none

    class(time_axis_type), intent(inout) :: self

    real(r_def), dimension(2) :: time_window

    time_window = self%time_data(1:2)

  end function get_time_window

  !> Returns start and end time indices for active time window
  !> @result i_window The indices of the active time window
  function get_index_window(self) result(i_window)

    implicit none

    class(time_axis_type), intent(inout) :: self

    integer(i_def), dimension(2) :: i_window

    i_window = self%index_data(1:2)

  end function get_index_window

  !> Takes model time and shifts forward through time axis so active time
  !> window is aligned with model time
  !> @param [in] input_time The current time of the model
  subroutine align(self, input_time)

    implicit none

    class(time_axis_type), intent(inout) :: self
    real(r_def),           intent(in)    :: input_time

    real(r_def),    dimension(2) :: t_window
    integer(i_def), dimension(2) :: i_window
    real(r_def) :: year_length
    integer(i_def) :: n

    ! Until final clock/calendar implementation, the time is currently in the
    ! form of a real number of days (future work will introduce time units)
    t_window = self%get_time_window()
    do n=1,size(self%time_data)
      if ( (t_window(1) <= input_time) .and. (input_time < t_window(2)) ) then
        return
      else
        call self%shift_forward()
        t_window = self%get_time_window()
      end if
    end do

    ! If we're still going then the correct time is either on the boundary
    ! between years or out of the axis bounds
    if ( self%is_cyclic() ) then
      i_window = self%get_index_window()
      do n=1,size(self%index_data)
        if ( (i_window(1) > i_window(2)) ) then
          return
        else
          call self%shift_forward()
          i_window = self%get_index_window()
        end if
      end do
    else
      write( log_scratch_space, '(A,E16.8,A)' ) "Model time ", input_time, &
                                      "out of bounds for non-cyclic time axis"
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end subroutine align

  !> Identifies if the time axis is cyclic in nature
  !> @result cyclic_flag Flag determining if time axis is cyclic
  function is_cyclic(self) result(cyclic_flag)

    implicit none

    class(time_axis_type), intent(inout) :: self

    logical(l_def) :: cyclic_flag

    cyclic_flag = self%cyclic

  end function is_cyclic

end module time_axis_mod