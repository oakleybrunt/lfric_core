!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the stepping of the diagnostics miniapp
!>
module diagnostics_step_mod

    use clock_mod,                only : clock_type
    use constants_mod,            only : i_def, str_def, r_def
    use diagnostics_alg_mod,      only : diagnostics_alg
    use field_mod,                only : field_type
    use field_parent_mod,         only : field_parent_type
    use field_collection_mod,     only : field_collection_type, &
                                         field_collection_iterator_type
    use io_config_mod,            only : write_diag
    use gungho_model_data_mod,    only : model_data_type
    use colours__prognostics__meta_mod, only : colours__prognostics__meta_type
    use colours__diagnostics__meta_mod, only : colours__diagnostics__meta_type
    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    private
    public diagnostics_step


contains

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !> Performs time steps.
    !>
    subroutine diagnostics_step( mesh_id,      &
                                 twod_mesh_id, &
                                 model_data,   &
                                 clock )

        implicit none

        ! Model run working data set
        integer(i_def),                  intent(in)    :: mesh_id ! included for consistency with gungho
        integer(i_def),                  intent(in)    :: twod_mesh_id ! included for consistency with gungho
        type( model_data_type ), target, intent(inout) :: model_data
        class(clock_type),               intent(in)    :: clock ! included for consistency with gungho

        type(field_type), pointer :: red => null()
        type(field_type), pointer :: green => null()
        type(field_type), pointer :: blue => null()
        type(field_type), pointer :: hex => null()

        type(colours__prognostics__meta_type) :: prognostics_meta
        type(colours__diagnostics__meta_type) :: diagnostics_meta
        type(field_collection_type), pointer :: prognostic_fields => null()
        character(str_def) :: hex_id

        ! Demonstrate here that fields can be obtained from their field
        ! collection in model_data, or directly from the depository

        ! Use the metadata to get field collection from model_data
        ! Get individual fields using their unique IDs
        prognostics_meta = colours__prognostics__meta_type()
        prognostic_fields => model_data%get_field_collection( prognostics_meta%name )
        red => prognostic_fields%get_field( prognostics_meta%red%get_unique_id() )
        green => prognostic_fields%get_field( prognostics_meta%green%get_unique_id() )
        blue => prognostic_fields%get_field( prognostics_meta%blue%get_unique_id() )

        ! Get an individual field from depository
        diagnostics_meta = colours__diagnostics__meta_type()
        hex_id = diagnostics_meta%hex%get_unique_id()
        hex => model_data%depository%get_field( hex_id )

        ! Call an algorithm
        call diagnostics_alg(&
                red, &
                green, &
                blue, &
                hex &
                )

        ! Write the fields
        if (write_diag) then
            call log_event("Writing " // red%get_name(), LOG_LEVEL_INFO)
            call red%write_field(red%get_name())
            call log_event("Writing " // green%get_name(), LOG_LEVEL_INFO)
            call green%write_field(green%get_name())
            call log_event("Writing " // blue%get_name(), LOG_LEVEL_INFO)
            call blue%write_field(blue%get_name())
            call log_event("Writing " // hex%get_name(), LOG_LEVEL_INFO)
            call hex%write_field(hex%get_name())
        end if

    end subroutine diagnostics_step

end module diagnostics_step_mod
