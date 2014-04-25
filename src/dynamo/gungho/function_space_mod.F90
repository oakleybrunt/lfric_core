!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief A type which holds information about the function space.

!> @details A container which holds data about the function space and has a
!> a basis function and accessor functions.
!> Private functions popoulate the values of the basis and differential basis 
!> functions and get an individual value for the index specified.
!> @param dofmap is a two-dim array, which is of size number_of_dofs_per_cell 
!> and number_of_cells
!> @param ndf private how many degrees of freedom per cell
!> @param ncell private  how many cells (in a layer)
!> @param undf private how many unique degrees of freedom
!> @param npg private number of gaussian quadrature points


 

module function_space_mod

use constants_mod, only:dp
use basis_function_mod, only: basis_function_type

implicit none
private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
  
type, public :: function_space_type
  private
  integer              :: ndf, ncell, undf, ngp
  integer              :: dim_space, dim_space_p1
  !> A two dimensional, allocatable array which holds the indirection map 
  !! or dofmap for the whole function space over the bottom level of the domain.
  integer, allocatable :: dofmap(:,:)
  !> The basis functions for this function space
  type(basis_function_type) :: basis_function
  ! accessor functions go here
contains
  !final :: destructor

!> Function to get the total unique degrees of freedom for this space
!! returns an integer
!! @param[in] self the calling function space
  procedure :: get_undf

!> Function Returns the number of cells in the function space
!> @param[in] self the calling function space.
!> @return Integer the number of cells
  procedure :: get_ncell

!> Subroutine Returns a pointer to the dofmap for the cell 
!! @param[in] self The calling functions_space
!! @param[in] cell Which cell
!! @param[out] map The pointer which points to a slice of the dofmap
  procedure :: get_cell_dofmap

!> Subroutine which populates the dofmap with data
!! @param[in] self The calling function space
!! @param[in] cell which cell
!! @param[in] map a 1-d array of integers which is the dofmap for this cell
  procedure :: populate_cell_dofmap

!> Function which obtains the number of dofs per cell
!! @param[in] self The calling functions space
!! return an integer, the number of dofs per cell
  procedure :: get_ndf

!>  Subroutine Wrapper to the accessor procedure for the basis function
!! @param[in] self the calling function space
!! @param[out] basis a pointer to the array to hold the values of the basis function 
  procedure :: get_basis

!> Subroutine Wrapper to the accessor procedure for the basis function
!! to set an invididual value
!! @param[inout] self the calling basis function
!! @param[in] p1 integer specifing element of first array dimension
!! @param[in] p2 integer specifing element of second array dimension
!! @param[in] p3 integer specifing element of third array dimension
!! @param[in] p4 integer specifing element of fourth array dimension
!! @param[in] p5 integer specifing element of fifth array dimension
!! @param[in] basis a real value to set the basis function
  procedure :: set_basis

!>  Subroutine Wrapper to the accessor procedure for the basis function for the
!! differential or "next" function space
!! @param[in] self the calling function space
!! @param[out] basis a pointer to the real array to hold the values of the
!!  basis function 
  procedure :: get_diff_basis

!> Subroutine Wrapper to the accessor procedure for the "next" basis function
!! to set an invididual value 
!! @param[inout] self the calling function space
!! @param[in] p1 integer specifing element of first array dimension
!! @param[in] p2 integer specifing element of second array dimension
!! @param[in] p3 integer specifing element of third array dimension
!! @param[in] p4 integer specifing element of fourth array dimension
!! @param[in] p5 integer specifing element of fifth array dimension
!! @param[in] basis a real value to set the basis function
  procedure :: set_diff_basis
  

end type function_space_type

!-------------------------------------------------------------------------------
! Constructors
!-------------------------------------------------------------------------------

! overload the default structure constructor for function space
!> Constructor for the function space. Allocates the memory for the dofmap
!! and calls the constructor on for the basis function.
!! @param[in] num_cells
!! @param[in] num_dofs
!! @param[in] num_unique_dofs
!! @param[in] dim_space The dimension of this function space
!! @param[in] dim_space_p1 The dimension of the next function space
!! @param[in] ngp The number of guassian quadrature points
interface function_space_type
   module procedure constructor
end interface

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public get_ncell, get_cell_dofmap
contains

type(function_space_type) function constructor(num_cells,num_dofs, &
                                               num_unique_dofs,  &
                                               dim_space, dim_space_p1,  &
                                               ngp ) &
     result(self)
  !-----------------------------------------------------------------------------
  ! Constructor
  !-----------------------------------------------------------------------------

  !Arguments
  integer, intent(in) :: num_cells, num_dofs, num_unique_dofs
  integer, intent(in) :: dim_space, dim_space_p1
  integer, intent(in) :: ngp

  self%ncell        = num_cells
  self%ndf          = num_dofs
  self%undf         = num_unique_dofs
  self%dim_space    = dim_space
  self%dim_space_p1 = dim_space_p1
  self%ngp          = ngp
  
  ! allocate some space
  allocate(self%dofmap(0:num_cells,num_dofs))
  ! this would need populating 
  ! call the contructor for the basis function
  self%basis_function = basis_function_type(ndofs_per_cell = num_dofs, &
                              vdim = dim_space, vdim_diff = dim_space_p1, ngp = ngp)
  return
end function constructor

!subroutine destructor()
!  !-----------------------------------------------------------------------------
!  ! Destructor. Allocatables are handled by any F2003-compliant compiler
!  ! anyway.
!  !-----------------------------------------------------------------------------
!  implicit none
!
!  type(function_space_type) :: self
!
!  !deallocate( self%v3dofmap)
!  !deallocate( self%Rv3)
!
!  return
!end subroutine final_lfric

!-----------------------------------------------------------------------------
! Get total unique dofs for this space
!-----------------------------------------------------------------------------

!> Function to get the total unique degrees of freedom for this space
!! returns an integer
!! @param[in] self the calling function space
integer function get_undf(self)
  class(function_space_type) :: self
  get_undf=self%undf
  return
end function get_undf

!-----------------------------------------------------------------------------
! Get the number of cells for this function space
!-----------------------------------------------------------------------------
!> Function Returns the number of cells in the function space
!> @param[in] self the calling function space.
!> @return Integer the number of cells
integer function get_ncell(self)
  class(function_space_type) :: self
  get_ncell=self%ncell
  return
end function get_ncell


!-----------------------------------------------------------------------------
! Get the number of dofs for a single cell 
!-----------------------------------------------------------------------------
integer function get_ndf(self)
  class(function_space_type) :: self
  get_ndf=self%ndf
  return
end function get_ndf

!-----------------------------------------------------------------------------
! Get the dofmap for a single cell
!-----------------------------------------------------------------------------
!> Subroutine Returns a pointer to the dofmap for the cell 
!! @param[in] self The calling functions_space
!! @param[in] cell Which cell
!! @param[out] map The pointer which points to a slice of the dofmap
subroutine get_cell_dofmap(self,cell,map)
  implicit none
  class(function_space_type), target, intent(in) :: self
  integer,                            intent(in) :: cell
  integer, pointer,                   intent(out) :: map(:)

  map => self%dofmap(cell,:)
  return
end subroutine get_cell_dofmap

!-----------------------------------------------------------------------------
! Copy data in the dofmap
!-----------------------------------------------------------------------------
subroutine populate_cell_dofmap(self,cell,map)

  implicit none

  class(function_space_type), intent(inout) :: self
  integer, intent(in) :: cell
  integer, intent(in) :: map(self%ndf)

  integer        :: dof

  do dof = 1,self%ndf
     self%dofmap(cell,dof) = map(dof)
  end do
  return 
end subroutine populate_cell_dofmap

!-----------------------------------------------------------------------------
! Get the basis function
!-----------------------------------------------------------------------------
subroutine get_basis(self,basis)
  implicit none
  class(function_space_type), intent(in)    :: self  
  real(kind=dp),   pointer , intent(out) :: basis(:,:,:,:,:)
  call self%basis_function%get_basis(basis)
  return
end subroutine get_basis

!-----------------------------------------------------------------------------
! Set the basis function
!-----------------------------------------------------------------------------
subroutine set_basis(self,basis,p1,p2,p3,p4,p5)
  implicit none
  class(function_space_type), intent(inout) :: self  
  real(kind=dp), intent(in) :: basis
  integer,intent(in) :: p1,p2,p3,p4,p5
  call self%basis_function%set_basis(basis,p1,p2,p3,p4,p5)
  return
end subroutine set_basis

!-----------------------------------------------------------------------------
! Get the differential of the basis function
!-----------------------------------------------------------------------------
subroutine get_diff_basis(self,diff_basis)
  implicit none
  class(function_space_type), intent(in)  :: self  
  real(kind=dp), pointer,     intent(out) :: diff_basis(:,:,:,:,:)
  call self%basis_function%get_diff_basis(diff_basis)
  return
end subroutine get_diff_basis

!-----------------------------------------------------------------------------
! Set the differential of the basis function
!-----------------------------------------------------------------------------
subroutine set_diff_basis(self,basis,p1,p2,p3,p4,p5)
  implicit none
  class(function_space_type), intent(inout) :: self  
  real(kind=dp), intent(in) :: basis
  integer,intent(in) :: p1,p2,p3,p4,p5

  call self%basis_function%set_diff_basis(basis,p1,p2,p3,p4,p5)
  return
end subroutine set_diff_basis


end module function_space_mod
