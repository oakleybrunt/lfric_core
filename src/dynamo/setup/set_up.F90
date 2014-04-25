!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Sets up the function spaces

!> @details Thie code generates a mesh and determines the dofmaps. This will
!> be replaced with code that reads this in a mesh generation and paritioning
!> will be a pre-processor stage. The constructurs for the functions spaces
!> are called here once the information needed to populate them has been
!> created. 
!> There are no tests for this code as this will be replaced.

module set_up_mod

  use constants_mod, only: dp
  use function_space_mod, only : function_space_type
  use reference_element_mod, only : reference_cube

  use mesh_generator_mod,    only: mesh_generator_init,        &
                                   mesh_generator_cubedsphere, &
                                   mesh_generator_biperiodic,  &
                                   mesh_connectivity
  use num_dof_mod,           only: num_dof_init
  use compute_basis_function_mod, only : compute_basis
  use dofmap_mod,           only: get_dofmap
  implicit none
contains 

  subroutine set_up(v0, v1, v2, v3, num_layers)

    use log_mod, only : log_event, LOG_LEVEL_INFO

    implicit none

    type(function_space_type), intent(inout) :: v0, v1, v2, v3
    integer, intent(out)                     :: num_layers

    integer                                  :: num_cells, element_order
    logical                                  :: l_spherical
    real(kind=dp), parameter                 :: delta=1.0_dp

    integer                                  :: v_unique_dofs(4,2)
    integer                                  :: v_dof_entity(4,0:3)    

    character(len = 100)                     :: filename

    ! hard-coded these numbers are
    num_cells = 9
    num_layers = 3
    element_order = 0
    l_spherical = .false.
    filename = '../../data/Cubegrid.dat' 
    call log_event( "set_up: generating/reading the mesh", LOG_LEVEL_INFO )

!  ----------------------------------------------------------
!  Mesh generation, really a preprocessor step for reading
! -----------------------------------------------------------

    ! Setup reference cube  
    call reference_cube()
    ! Initialise mesh
    call mesh_generator_init(num_cells,num_layers)
    ! Genereate mesh  
    if ( l_spherical ) then
       call mesh_generator_cubedsphere(filename,num_cells,num_layers,delta)
    else
       call mesh_generator_biperiodic(num_cells,3,3,num_layers,delta,delta,delta)
    end if
    ! Extend connectivity ( cells->faces, cells->edges )  
    call mesh_connectivity(num_cells)    


! -----------------------------------------------------------
! Initialise FE elements on the mesh constructed above
! really another pre-processor step
! ----------------------------------------------------------
    call log_event( "set_up: building function spaces", LOG_LEVEL_INFO )
    ! initialise numbers of dofs    
    call num_dof_init(num_cells,num_layers,element_order,v_unique_dofs,v_dof_entity)
    ! call the constructors for the vspaces
    v0 = function_space_type( &
         num_cells = num_cells ,num_dofs = v_unique_dofs(1,2), &
         num_unique_dofs = v_unique_dofs(1,1) ,  &
         dim_space = 1, dim_space_p1 = 3,  &
         ngp = 3 )

    v1 = function_space_type( &
         num_cells = num_cells ,num_dofs = v_unique_dofs(2,2), &
         num_unique_dofs = v_unique_dofs(2,1) ,  &
         dim_space = 3, dim_space_p1 = 3,  &
         ngp = 3 )

    v2 = function_space_type( &
         num_cells = num_cells ,num_dofs = v_unique_dofs(3,2), &
         num_unique_dofs = v_unique_dofs(3,1) ,  &
         dim_space = 3, dim_space_p1 = 1,  &
         ngp = 3 )

    v3 = function_space_type( &
         num_cells = num_cells ,num_dofs = v_unique_dofs(4,2), &
         num_unique_dofs = v_unique_dofs(4,1) ,  &
         dim_space = 1, dim_space_p1 = 1,  &
         ngp = 3 )

    call log_event( "set_up: computing basis functions", LOG_LEVEL_INFO )

    ! compute the value of the basis functions and populate the
    ! basis_function_type
    call compute_basis(element_order,v0,v1,v2,v3,v_unique_dofs,v_dof_entity)  
    call log_event( "set_up: computing the dof_map", LOG_LEVEL_INFO )
    ! compute the dof maps for each function space
    call get_dofmap(num_layers,v0,v_dof_entity(1,:))
    call get_dofmap(num_layers,v1,v_dof_entity(2,:))
    call get_dofmap(num_layers,v2,v_dof_entity(3,:))
    call get_dofmap(num_layers,v3,v_dof_entity(4,:))

    return

  end subroutine set_up

end module set_up_mod
