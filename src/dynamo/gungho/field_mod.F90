!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Module for the field

!> @detail This module contains the type definition for the field_type


module field_mod
use function_space_mod, only: function_space_type
use gaussian_quadrature_mod, only : gaussian_quadrature_type, ngp
use constants_mod, only: dp
implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public :: field_type
  private
  integer :: nlayers, undf
  !> Each field has a pointer to the function space on which it lives
  type(function_space_type), pointer, public :: vspace
  !> Each field has a pointer to the gaussian quadrature rule which will be
  !! used to integrate over its values
  type(gaussian_quadrature_type), pointer, public :: gaussian_quadrature
  !> allocatable array of type real which hold the values of the field
  real(kind=dp), allocatable, public :: data(:)
contains
  !> function wrapper to call to get_ncell on the underlying function space
  !! @param[in] self the calling field
  !! @return an integer, the number of cells in a single layer (columns)
  procedure :: get_ncell

  !> function accessor function
  !! @param[in] self the calling field
  !! @return an integer, the number of layers
  procedure :: get_nlayers
end type field_type

!overload the default structure constructor for field
!> The interface used to overload the constructor
!! Sets the pointers and the number of layers
!! allocates the size of the data array by calling the procedure
!! get_undf(), the number of unique degrees of freedom on the underlying 
!! function space
!! @param vector_space the function space that the field lives on
!! @param gq the gaussian quadrature rule
!! @param num_layers integer number of layers for the field
!! @return self the field
  interface field_type
     module procedure field_constructor
  end interface

  public get_ncell
  public get_nlayers
contains
!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

  type(field_type) function field_constructor( &
       vector_space,                                   &
       gq,                                             &
       num_layers)                                     &
       result(self)
    ! constructor
    type(function_space_type), target, intent(in) :: vector_space
    type(gaussian_quadrature_type), target, intent(in) :: gq
    integer, intent(in) :: num_layers
    
    self%vspace => vector_space
    self%gaussian_quadrature => gq
    self%nlayers = num_layers
    self%undf = self%vspace%get_undf()

    ! allocate the array in memory
    allocate(self%data(self%undf))
    
  end function field_constructor

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
  !> function wrapper to call to get_ncell on the underlying function space
  !! @param[in] self the calling field
  !! @return an integer, the number of cells in a single layer (columns)
  integer function get_ncell(self)
    class(field_type) :: self
    get_ncell=self%vspace%get_ncell()
    return
  end function get_ncell

  !> function accessor function
  !! @param[in] self the calling field
  !! @return an integer, the number of layers
  integer function get_nlayers(self)
    class(field_type) :: self
    get_nlayers=self%nlayers
    return
  end function get_nlayers

end module field_mod

