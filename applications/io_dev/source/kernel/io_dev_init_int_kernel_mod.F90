!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Kernel to initialise IO_Dev integer field data to dofmap value
module io_dev_init_int_kernel_mod

  use kernel_mod,              only : kernel_type
  use argument_mod,            only : arg_type, GH_FIELD, GH_INTEGER, GH_INC, &
                                      ANY_SPACE_1, CELL_COLUMN
  use constants_mod,           only : i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: io_dev_init_int_kernel_type
private
type(arg_type) :: meta_args(1) = (/                         &
     arg_type(GH_FIELD,   GH_INTEGER, GH_INC,  ANY_SPACE_1) &
     /)
integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: io_dev_init_int_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public io_dev_init_int_code
contains

!> @brief   Initialise IO_Dev field data to its dof_id
!> @param[in]  nlayers   Number of layers
!> @param[out] int_data  Integer field data
!> @param[in]  ndf_x     Number of degrees of freedom per cell for the output field
!> @param[in]  undf_x    Number of unique degrees of freedom for the output field
!> @param[in]  map_x     Dofmap for the cell at the base of the column for the output field
!>                       nodal points of the x function space
subroutine io_dev_init_int_code( nlayers,                    &
                                 int_data,                   &
                                 ndf_x, undf_x, map_x        )

  implicit none

  ! Arguments
  integer(kind=i_def),                    intent(in)    :: nlayers
  integer(kind=i_def),                    intent(in)    :: ndf_x, undf_x
  integer(kind=i_def), dimension(ndf_x),  intent(in)    :: map_x
  integer(kind=i_def), dimension(undf_x), intent(inout) :: int_data

  ! Internal variables
  integer :: df_x, k

  ! Set field data equal to DoF value
  do k = 0, nlayers-1
    do df_x = 1,ndf_x
      int_data(map_x(df_x)+k) = int( map_x(df_x)+k, i_def )
    end do
  end do

end subroutine io_dev_init_int_code

end module io_dev_init_int_kernel_mod