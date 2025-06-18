!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief   The weighted mapping from a lower to a higher-order scalar space.
!> @details Take a scalar space which is lowest-order in the horizontal but
!!          shares the same number of levels, and map it to a higher-order
!!          scalar space on a more coarse mesh.

module sci_map_scalar_fv_to_fe_kernel_mod

use argument_mod,              only: arg_type,                  &
                                     GH_FIELD, GH_REAL,         &
                                     GH_READ, GH_WRITE,         &
                                     ANY_DISCONTINUOUS_SPACE_1, &
                                     ANY_DISCONTINUOUS_SPACE_2, &
                                     ANY_DISCONTINUOUS_SPACE_3, &
                                     GH_COARSE, GH_FINE, CELL_COLUMN
use constants_mod,             only: i_def, r_def
use kernel_mod,                only: kernel_type

implicit none

private

type, public, extends(kernel_type) :: map_scalar_fv_to_fe_kernel_type
   private
   type(arg_type) :: meta_args(3) = (/                                       &
        arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1,     &
                 mesh_arg=GH_COARSE),                                        &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2,     &
                 mesh_arg=GH_FINE),                                          &
        arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3,     &
                 mesh_arg=GH_COARSE)                                         &
        /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: map_scalar_fv_to_fe_code
end type map_scalar_fv_to_fe_kernel_type

public :: map_scalar_fv_to_fe_code

contains

  !> @brief Performs the mapping between scalar fields of different orders
  !> @param[in]     nlayers                  Number of layers in a model column
  !> @param[in]     cell_map                 A 2D index map of which fine grid
  !!                                         cells lie in the coarse grid cell
  !> @param[in]     ncell_fine_per_coarse_x  Number of fine cells per coarse
  !!                                         cell in the horizontal x-direction
  !> @param[in]     ncell_fine_per_coarse_y  Number of fine cells per coarse
  !!                                         cell in the horizontal y-direction
  !> @param[in]     ncell_fine               Number of cells in the partition
  !!                                         for the fine grid
  !> @param[in,out] coarse_field             Higher-order coarse grid field to
  !!                                         map to
  !> @param[in]     fine_field               Lowest-order fine grid field to
  !!                                         map from
  !> @param[in]     weights                  Weights for the mapping
  !> @param[in]     ndf_coarse               Num of DoFs per cell on the coarse
  !!                                         grid
  !> @param[in]     undf_coarse              Total num of DoFs on the coarse
  !!                                         grid for this mesh partition
  !> @param[in]     map_coarse               DoFmap of cells on the coarse grid
  !> @param[in]     ndf_fine                 Num of DoFs per cell on the fine
  !!                                         grid
  !> @param[in]     undf_fine                Total num of DoFs on the fine grid
  !!                                         for this mesh partition
  !> @param[in]     map_fine                 DoFmap of cells on the fine grid
  !> @param[in]     undf_weights             Num of DoFs on the weights field
  !> @param[in]     map_weights              DoFmap of cells on the weights
  !!                                         field

  subroutine map_scalar_fv_to_fe_code(nlayers,                 &
                                      cell_map,                &
                                      ncell_fine_per_coarse_x, &
                                      ncell_fine_per_coarse_y, &
                                      ncell_fine,              &
                                      coarse_field,            &
                                      fine_field,              &
                                      weights,                 &
                                      ndf_coarse,              &
                                      undf_coarse,             &
                                      map_coarse,              &
                                      ndf_fine,                &
                                      undf_fine,               &
                                      map_fine,                &
                                      undf_weights,            &
                                      map_weights)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_x
    integer(kind=i_def), intent(in)    :: ncell_fine_per_coarse_y
    integer(kind=i_def), intent(in)    :: ncell_fine
    integer(kind=i_def), intent(in)    :: ndf_fine, ndf_coarse
    integer(kind=i_def), intent(in)    :: undf_fine
    integer(kind=i_def), intent(in)    :: undf_coarse
    integer(kind=i_def), intent(in)    :: undf_weights

    ! Fields
    real(kind=r_def), intent(inout) :: coarse_field(undf_coarse)
    real(kind=r_def), intent(in)    :: fine_field(undf_fine)
    real(kind=r_def), intent(in)    :: weights(undf_weights)

    ! Maps
    integer(kind=i_def), intent(in)    :: cell_map(ncell_fine_per_coarse_x, &
                                                   ncell_fine_per_coarse_y)
    integer(kind=i_def), intent(in)    :: map_fine(ndf_fine, ncell_fine)
    integer(kind=i_def), intent(in)    :: map_coarse(ndf_coarse)
    integer(kind=i_def), intent(in)    :: map_weights(ndf_coarse)

    ! Internal variables
    integer(kind=i_def) :: k, x_idx, y_idx, df, j, df_f, top_layer, top_df

    ! Assume lowest-order W3 or Wtheta space for fine mesh
    df_f = 1
    ! Loop is 0 -> nlayers-1 for W3 fields, but 0 -> nlayers for Wtheta fields
    top_layer = nlayers - 2 + ndf_fine
    ! Need to loop over volume dofs of W3 or one face of Wtheta on coarse space
    top_df = int(ndf_coarse/ndf_fine, i_def)

    do y_idx = 1, ncell_fine_per_coarse_y
      do x_idx = 1, ncell_fine_per_coarse_x
        j = x_idx + (y_idx - 1)*ncell_fine_per_coarse_x - 1
        do k = 0, top_layer
          do df = 1, top_df
            coarse_field(map_coarse(df)+k) = coarse_field(map_coarse(df)+k) &
              + fine_field(map_fine(df_f,cell_map(x_idx,y_idx))+k) &
                * weights(map_weights(df)+ j)
          end do
        end do
      end do
    end do

  end subroutine map_scalar_fv_to_fe_code


end module sci_map_scalar_fv_to_fe_kernel_mod
