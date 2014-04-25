!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @mainpage Dynamo
!> PsyKAl is the architecture for Gung Ho. Whlist the computational and optimisation
!> infrastructure is being developed, the science code is being developed using 
!> a hand-rolled Psy layer, Psy-lite. A PsyKAl-lite needs a dynamo!
!> Eventually, PsyKAlite is replaced with the real Psy and Dynamo becomes Gung Ho.

!> @brief Main program used to illustrate dynamo functionality.

!> @details Creates the function spaces and alls the set_up to populate them (
!> either read or compute) then individual calls to the psy-layer with kernels
!> as if the code has been pre-processed by Psyclone.
!> Comments starting with !PSY are what the code would lool like before Psyclone
!> generated anything.

program dynamo

  use lfric
  use log_mod,              only : log_event, log_scratch_space, LOG_LEVEL_INFO
  use set_up_mod,           only : set_up
  use dynamo_algorithm_mod, only : dynamo_algorithm

  implicit none

  type( function_space_type )      :: v3_function_space, v2_function_space, &
                                      v1_function_space, v0_function_space
  type( field_type )               :: pressure_density, rhs
  type( gaussian_quadrature_type ) :: gq
  integer                          :: num_layers

  call log_event( 'Dynamo running...', LOG_LEVEL_INFO )

  call set_up( v0_function_space, v1_function_space, v2_function_space, &
               v3_function_space, num_layers )

  gq = gaussian_quadrature_type( )

  pressure_density = field_type( vector_space = v3_function_space, &
                                 gq = gq,                          &
                                 num_layers = num_layers)

  rhs = field_type( vector_space = v3_function_space, &
                    gq = gq,                          &
                    num_layers = num_layers )

  call dynamo_algorithm( pressure_density, rhs )

  call log_event( 'Dynamo completed', LOG_LEVEL_INFO )

end program dynamo
