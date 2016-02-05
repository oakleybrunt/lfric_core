!------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------
!> @brief A type which holds information about the master dofmap.
!> @details Type for storing the master_dofmap ( dofmap for a single cell)
module master_dofmap_mod

use constants_mod, only : i_def

implicit none

private
type, public :: master_dofmap_type
  private 
  integer(i_def), allocatable :: dofmap(:,:) 
contains
  procedure :: get_master_dofmap
end type master_dofmap_type

interface master_dofmap_type
  module procedure master_dofmap_constructor
end interface

contains 

!-----------------------------------------------------------------------------
! Construct the master dofmap
!-----------------------------------------------------------------------------
!> Function to construct the master (cell) dofmap
!> @param[in] master_dofmap
!> @return The master dofmap object
function master_dofmap_constructor( master_dofmap ) result(self)

  implicit none

  integer(i_def), intent(in) :: master_dofmap(:,:)
  type(master_dofmap_type) :: self

  integer(i_def) :: dim1, dim2

  dim1 = size(master_dofmap,1)
  dim2 = size(master_dofmap,2)-1

  allocate( self%dofmap(dim1,0:dim2) )
  self%dofmap(:,:) = master_dofmap(:,:)

  return
end function master_dofmap_constructor

!-----------------------------------------------------------------------------
! Get the master dofmap for a single cell
!-----------------------------------------------------------------------------
!> Subroutine Returns a pointer to the dofmap for the cell 
!! @param[in] self The calling function_space
!! @param[in] cell Which cell
!! @return The pointer which points to a slice of the dofmap
function get_master_dofmap(self,cell) result(map)
  implicit none
  class(master_dofmap_type), target, intent(in) :: self
  integer(i_def),                    intent(in) :: cell
  integer(i_def), pointer                       :: map(:) 

  map => self%dofmap(:,cell)
  return
end function get_master_dofmap
 
end module master_dofmap_mod

