!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
! Computes the basis functions on quadrature points for 4 element spaces 
!-------------------------------------------------------------------------------
!> @brief The basis function type
!> @detail Constainer type for arrays containing the values of the basis
!! contains integers for array sizes. Two arrays of reals for the basis function
!! for this space and the "next" space. Accessor functions to get and set the 
!! array values.

module basis_function_mod

use constants_mod, only: dp

implicit none
private


!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public :: basis_function_type
   private

   integer :: dim, diff_dim,ndf,ngp   
   !> 5-dim allocatable array of reals which hold the values of the basis function
   real(kind=dp), allocatable :: basis(:,:,:,:,:)
   !> 5-dim allocatable array of reals which hold the values of the basis function
   !! for the differential or "next" functions space
   real(kind=dp), allocatable :: diff_basis(:,:,:,:,:)

 contains
!> subroutine Set values of the basis function 
!! @param[inout] self the calling basis function
!! @param[in] p1 the first array dimension element
!! @param[in] p2 the second array dimension element
!! @param[in] p3 the third array dimension element
!! @param[in] p4 the fourth array dimension element
!! @param[in] p5 the fifth array dimension element
!! @param[in] basis real value of basis element
   procedure :: set_basis

!> subroutine Set values of the differential or next basis function 
!! @param[inout] self the calling basis function
!! @param[in] p1 the first array dimension element
!! @param[in] p2 the second array dimension element
!! @param[in] p3 the third array dimension element
!! @param[in] p4 the fourth array dimension element
!! @param[in] p5 the fifth array dimension element
!! @param[in] basis real value of basis element
   procedure :: set_diff_basis

!> subroutine returns a pointer to the array of reals which are the values 
!! of the basis functions
!! @param[in] self the calling basis function
!! @param[out] basis the pointer which points to the basis function array
   procedure :: get_basis

!> subroutine returns a pointer to the array of reals which are the values
!!  of the differented basis functions
!! @param[in] self the calling basis function
!! @param[out] basis the pointer which points to the diff basis array
   procedure :: get_diff_basis
   
end type basis_function_type

! ------------------------------------
! Constructors   
! -------------------------------------

! overload the default structure construcor for the function spacs
!> Constructor for the basis_function. Sets the integers and allocates
!! the two arrays for the basis and diff basis.
!! @param[in] ndofs_per_cell number of degrees of freedom per cell
!! @param[in] vdim the dimension of the function_space to which this belongs
!! @param[in] vdim_diff the dimension of the differentiaed or "next" function space
!! @param[in] ngp The number of gaussian quadrature points
!! @return self The basis_function
interface basis_function_type
   module procedure basis_init
end interface


!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

contains 


!-----------------------------------------------------------------------------
! The constructor
!-----------------------------------------------------------------------------
  type(basis_function_type) function basis_init(ndofs_per_cell, vdim, vdim_diff,ngp) & 
       result(self)
    integer, intent(in) :: ndofs_per_cell, vdim, vdim_diff,ngp

    self%dim = vdim
    self%diff_dim = vdim_diff
    self%ndf = ndofs_per_cell
    self%ngp = ngp
    allocate( self%basis(ndofs_per_cell,self%ngp,self%ngp,self%ngp,vdim) )
    allocate( self%diff_basis(ndofs_per_cell,self%ngp,self%ngp,self%ngp,vdim_diff) )

  end function basis_init

  subroutine set_basis(self,basis,p1,p2,p3,p4,p5)

    implicit none
    class(basis_function_type), intent(inout) :: self
    real(kind=dp), intent(in) :: basis
    integer, intent(in) :: p1,p2,p3,p4,p5

    self%basis(p1,p2,p3,p4,p5) = basis
    
  end subroutine set_basis

  subroutine set_diff_basis(self,basis,p1,p2,p3,p4,p5)

    implicit none
    class(basis_function_type), intent(inout) :: self
    real(kind=dp), intent(in) :: basis
    integer, intent(in) :: p1,p2,p3,p4,p5
    self%diff_basis(p1,p2,p3,p4,p5) = basis
    
  end subroutine set_diff_basis

  subroutine get_basis(self,basis)
    class(basis_function_type), target, intent(in)   :: self
    real(kind=dp), pointer,             intent(out)  :: basis(:,:,:,:,:)
    basis => self%basis
    return
  end subroutine get_basis

  subroutine get_diff_basis(self,basis)
    class(basis_function_type), target, intent(in)  :: self
    real(kind=dp), pointer,             intent(out) :: basis(:,:,:,:,:)

    basis => self%diff_basis
  end subroutine get_diff_basis


end module basis_function_mod
