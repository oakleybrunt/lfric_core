!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> A simple algorithm for testing the psy layer.
!>
module dynamo_algorithm_mod

  use lfric
  use log_mod,    only: log_event, log_scratch_space, LOG_LEVEL_INFO
  use psy,        only: invoke_rhs_v3, invoke_v3_solver_kernel

  implicit none

  private
  public :: dynamo_algorithm

contains

  !> A simple algorithm which calls two kernels.
  !>
  subroutine dynamo_algorithm( pressure_density, rhs )

    implicit none

    type( field_type ), intent( inout ) :: pressure_density
    type( field_type ), intent( inout ) :: rhs

    !Construct PSy layer given a list of kernels. This is the line the code
    !generator may parse and do its stuff.

    call log_event( "Dynamo: calling 1st kernel", LOG_LEVEL_INFO )
    !PSY call invoke ( v3_rhs_kernel_type(rhs) )
    call invoke_rhs_v3( rhs )

    call log_event( "Dynamo:calling 2nd kernel", LOG_LEVEL_INFO )
    !PSY call invoke ( v3_solver_kernel_type(pressure_density,rhs) )
    call invoke_v3_solver_kernel( pressure_density, rhs )

    call print_field( 'RHS field...', rhs )
    call print_field( 'LHS field...', pressure_density )

  end subroutine dynamo_algorithm

  !> Send a field to the log.
  !>
  subroutine print_field( title, field )

    use lfric
    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    character( * ),     intent( in ) :: title
    type( field_type ), intent( in ) :: field

    integer                   :: cell
    integer                   :: layer
    integer                   :: df
    integer,          pointer :: map( : )

    call log_event( title, LOG_LEVEL_INFO )

    do cell=1,field%vspace%get_ncell()
      call field%vspace%get_cell_dofmap(cell,map)
      do df=1,field%vspace%get_ndf()
        do layer=0,field%get_nlayers()-1
          write( log_scratch_space, '( I4, I4, I4, F8.2 )' ) &
              cell, df, layer+1, field%data( map( df ) + layer )
          call log_event( log_scratch_space, LOG_LEVEL_INFO )
        end do
      end do
    end do

  end subroutine print_field

end module dynamo_algorithm_mod
