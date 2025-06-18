!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes weights for the scalar inter element-order mapping.
!> @details Compute the weights for the mapping higher to lower-order as the
!!          basis functions of the higher-order space evaluated at the dofs
!!          of the lower-order space inside the coarse mesh cells. Invert this
!!          matrix to get the weights for the mapping lower to higher-order.
!!          The segment centre quadrature rule gives the fine dof locations.

module sci_weights_scalar_inter_element_order_kernel_mod

use argument_mod,      only: arg_type, func_type,       &
                             GH_FIELD, GH_REAL,         &
                             GH_READ, GH_WRITE,         &
                             GH_COARSE, GH_FINE,        &
                             GH_BASIS,                  &
                             ANY_DISCONTINUOUS_SPACE_1, &
                             ANY_DISCONTINUOUS_SPACE_2, &
                             ANY_DISCONTINUOUS_SPACE_3, &
                             CELL_COLUMN,               &
                             GH_QUADRATURE_XYoZ
use constants_mod,     only: i_def, r_def
use kernel_mod,        only: kernel_type
use matrix_invert_mod, only: matrix_invert_plu

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>

type, public, extends(kernel_type) :: weights_scalar_inter_element_order_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/ &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1, mesh_arg=GH_FINE  ), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2, mesh_arg=GH_COARSE), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3, mesh_arg=GH_COARSE), &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3, mesh_arg=GH_COARSE)  &
       /)
  type(func_type) :: meta_funcs(1) = (/func_type(ANY_DISCONTINUOUS_SPACE_2, GH_BASIS)/)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: weights_scalar_inter_element_order_code
end type weights_scalar_inter_element_order_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------

public :: weights_scalar_inter_element_order_code

contains

  !> @brief Computes weights for the scalar inter element-order mapping
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in]     dummy_fine               An unused fine field supplied
  !!                                         to ensure that appropriate fine
  !!                                         field properties are provided.
  !> @param[in]     dummy_coarse             An unused coarse field supplied
  !!                                         to ensure that appropriate coarse
  !!                                         field properties are provided.
  !> @param[in,out] weights_fe_to_fv         The coarse grid 2D field that
  !!                                         contains the weights to be
  !!                                         computed here for mapping high
  !!                                         order to low order
  !> @param[in,out] weights_fv_to_fe         The coarse grid 2D field that
  !!                                         contains the weights to be
  !!                                         computed here for mapping low
  !!                                         order to high order
  !> @param[in]     ndf_fine                 Num of DoFs per cell on the fine
  !!                                         grid
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     ndf_coarse               Num of DoFs per cell on the coarse
  !!                                         grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     basis_coarse             Basis functions on the coarse grid
  !> @param[in]     undf_weights             Total num of DoFs on the weights
  !!                                         field for this mesh partition
  !> @param[in]     map_weights              DoFmap of cells on the weights
  !!                                         fields
  !> @param[in]     nqp_xy                   Number of quadrature points in the
  !!                                         horizontal
  !> @param[in]     nqp_z                    Number of quadrature points in the
  !!                                         vertical
  !> @param[in]     weights_xy               Weights for the horizontal
  !!                                         quadrature points
  !> @param[in]     weights_z                Weights for the vertical
  !!                                         quadrature points
  subroutine weights_scalar_inter_element_order_code(nlayers,                 &
                                                     cell_map,                &
                                                     ncell_fine_per_coarse_x, &
                                                     ncell_fine_per_coarse_y, &
                                                     ncell_fine,              &
                                                     dummy_fine,              &
                                                     dummy_coarse,            &
                                                     weights_fe_to_fv,        &
                                                     weights_fv_to_fe,        &
                                                     ndf_fine,                &
                                                     undf_fine,               &
                                                     map_fine,                &
                                                     ndf_coarse,              &
                                                     undf_coarse,             &
                                                     map_coarse,              &
                                                     basis_coarse,            &
                                                     undf_weights,            &
                                                     map_weights,             &
                                                     nqp_xy,                  &
                                                     nqp_z,                   &
                                                     weights_xy,              &
                                                     weights_z)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in) :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in) :: ncell_fine
    integer(kind=i_def), intent(in) :: ndf_fine, ndf_coarse
    integer(kind=i_def), intent(in) :: undf_fine, undf_coarse
    integer(kind=i_def), intent(in) :: undf_weights
    integer(kind=i_def), intent(in) :: nqp_xy, nqp_z

    real(kind=r_def),    intent(in) :: basis_coarse(1, ndf_coarse, nqp_xy, nqp_z)
    real(kind=r_def),    intent(in) :: weights_xy(nqp_xy)
    real(kind=r_def),    intent(in) :: weights_z(nqp_z)

    ! Fields
    real(kind=r_def),    intent(inout) :: weights_fe_to_fv(undf_weights)
    real(kind=r_def),    intent(inout) :: weights_fv_to_fe(undf_weights)
    real(kind=r_def),    intent(in)    :: dummy_fine(undf_fine)
    real(kind=r_def),    intent(in)    :: dummy_coarse(undf_coarse)

    ! Maps
    integer(kind=i_def), intent(in) :: map_fine(ncell_fine, ncell_fine)
    integer(kind=i_def), intent(in) :: map_coarse(ndf_coarse)
    integer(kind=i_def), intent(in) :: map_weights(ndf_coarse)
    integer(kind=i_def), intent(in) :: cell_map(ncell_fine_per_coarse_x, &
                                                ncell_fine_per_coarse_y)

    ! Internal arguments
    integer(kind=i_def) :: df, x_idx, y_idx, j

    real(kind=r_def) :: weights_matrix(ndf_coarse, ndf_coarse)
    real(kind=r_def) :: inv_weights_matrix(ndf_coarse, ndf_coarse)

    !---------------------------------------------------------------------------
    ! Calculate the weights for the high to low mapping, store these to invert
    ! mapping
    !---------------------------------------------------------------------------
    ! If f(x) is a scalar field, d_i is the i-th dof of the higher-order
    ! space and b_i is the corresponding basis function, then
    !
    !   f(x) = sum_i f(d_i)*b_i(x)
    !
    ! Taking x to be the location of a lowest-order dof, the weights for
    ! mapping from f on the higher-order space to f on the lowest-order space
    ! are the higher-order basis functions evaluated at the lowest-order dof
    ! locations
    !
    ! We call the higher-order function space the "coarse space", since it is
    ! stored on a coarse mesh, and likewise the lowest-order function space is
    ! called the "fine space".

    do df = 1, ndf_coarse
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          j = x_idx + (y_idx - 1)*ncell_fine_per_coarse_x
          weights_matrix(df, j) = basis_coarse(1,df,j,1)
          weights_fe_to_fv(map_weights(df) + j - 1) = weights_matrix(df, j)
        end do
      end do
    end do

    !---------------------------------------------------------------------------
    ! Invert the weights matrix on each coarse cell
    !---------------------------------------------------------------------------
    ! Matrix produced does not typically have an LU decomposition, so use
    ! matrix inversion with pivoting
    call matrix_invert_plu(weights_matrix, inv_weights_matrix, ndf_coarse)

    !---------------------------------------------------------------------------
    ! Flatten the produced inverse weights matrix
    !---------------------------------------------------------------------------
    do df = 1, ndf_coarse
      do y_idx = 1, ncell_fine_per_coarse_y
        do x_idx = 1, ncell_fine_per_coarse_x
          j = x_idx + (y_idx - 1)*ncell_fine_per_coarse_x
          weights_fv_to_fe(map_weights(df) + j - 1) = inv_weights_matrix(df, j)
        end do
      end do
    end do

  end subroutine weights_scalar_inter_element_order_code

end module sci_weights_scalar_inter_element_order_kernel_mod
