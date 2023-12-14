!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the vertical departure distances, to be consistent
!!        with the transport of dry density.
!> @details This code calculates the distance which is swept through a cell
!!          in the z direction during one time step of consistent transport. The
!!          arrival point is the cell face and the departure point is calculated
!!          and stored as a field.
!!          The departure points are a dimensionless displacement, corresponding
!!          to the number of cells moved by a fluid parcel. The part of the dry
!!          flux corresponding to the fractional part of the departure points is
!!          also computed.

module consistent_vertical_deppt_kernel_mod

use argument_mod,                only : arg_type,              &
                                        GH_FIELD, GH_REAL,     &
                                        GH_WRITE, GH_READ,     &
                                        CELL_COLUMN
use fs_continuity_mod,           only : W3, W2v
use constants_mod,               only : r_tran, i_def
use kernel_mod,                  only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: consistent_vertical_deppt_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                  &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2v),    & ! dep_pts
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2v),    & ! frac_dry_flux
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2v),    & ! dry_flux
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),     & ! rho_d
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3)      & ! detj_at_w3
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: consistent_vertical_deppt_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: consistent_vertical_deppt_code

contains

!> @brief Kernel to computes the vertical departure distances, to be consistent
!!        with the transport of dry density.
!> @param[in]     nlayers        The number of layers in the mesh
!> @param[in,out] dep_pts_z      Field with vertical departure distances
!> @param[in,out] frac_dry_flux  Fractional part of the vertical dry flux, to
!!                               be computed here
!> @param[in]     dry_flux       The vertical mass flux used in a vertical
!!                               transport step for the dry density
!> @param[in]     rho_d          The dry density field before the vertical
!!                               transport step
!> @param[in]     detj_at_w3     det(J) at W3 points
!> @param[in]     ndf_w2v        Number of DoFs per cell for W2V
!> @param[in]     undf_w2v       Number of W2V DoFs in memory for this partition
!> @param[in]     map_w2v        Map of lowest-cell W2V DoFs
!> @param[in]     ndf_w3         Number of DoFs per cell for W3
!> @param[in]     undf_w3        Number of W3 DoFs in memory for this partition
!> @param[in]     map_w3         Map of lowest-cell W3 DoFs
subroutine consistent_vertical_deppt_code( nlayers,             &
                                           dep_pts_z,           &
                                           frac_dry_flux,       &
                                           dry_flux,            &
                                           rho_d,               &
                                           detj_at_w3,          &
                                           ndf_w2v,             &
                                           undf_w2v,            &
                                           map_w2v,             &
                                           ndf_w3,              &
                                           undf_w3,             &
                                           map_w3 )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2v, ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w2v, undf_w3
  integer(kind=i_def), intent(in)    :: map_w2v(ndf_w2v)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  real(kind=r_tran),   intent(in)    :: rho_d(undf_w3)
  real(kind=r_tran),   intent(in)    :: dry_flux(undf_w2v)
  real(kind=r_tran),   intent(in)    :: detj_at_w3(undf_w3)
  real(kind=r_tran),   intent(inout) :: dep_pts_z(undf_w2v)
  real(kind=r_tran),   intent(inout) :: frac_dry_flux(undf_w2v)

  integer(kind=i_def) :: j, k, cell_idx, max_num_cells, num_int_cells
  integer(kind=i_def) :: sign_flux, offset
  real(kind=r_tran)   :: flux_face_k, running_flux_k, flux_cell_j
  real(kind=r_tran)   :: frac_dep_dist, frac_flux_face_k

  ! Set the bottom values
  frac_dry_flux(map_w2v(1)) = 0.0_r_tran
  dep_pts_z(map_w2v(1)) = 0.0_r_tran

  do k = 1, nlayers - 1

    flux_face_k = dry_flux(map_w2v(1)+k)

    ! Get an upper limit on the number of cells to step through
    if (flux_face_k >= 0.0_r_tran) then
      max_num_cells = k
      offset = 0
      sign_flux = 1
    else
      max_num_cells = nlayers - k
      offset = -1
      sign_flux = -1
    end if

    num_int_cells = 0
    running_flux_k = 0.0_r_tran

    ! Step backwards through flux to find departure point
    do j = 1, max_num_cells
      ! Index of cell to look at
      cell_idx = map_w3(1) + k - sign_flux*j + offset
      flux_cell_j = rho_d(cell_idx) * detj_at_w3(cell_idx)

      ! We have found integer number of cells, so exit do-loop
      if (running_flux_k + flux_cell_j > abs(flux_face_k)) EXIT

      ! Increment running values, if we have not exited loop
      running_flux_k = running_flux_k + flux_cell_j
      num_int_cells = num_int_cells + 1
    end do

    ! Set fractional flux. running_flux_k is now the integer flux
    frac_flux_face_k = sign_flux*(abs(flux_face_k) - running_flux_k)

    ! Determine fractional distance
    cell_idx = map_w3(1) + k - sign_flux*(num_int_cells + 1) + offset
    frac_dep_dist = frac_flux_face_k / (rho_d(cell_idx) * detj_at_w3(cell_idx))

    ! Set the values of the output fields
    frac_dry_flux(map_w2v(1)+k) = frac_flux_face_k
    dep_pts_z(map_w2v(1)+k) = real(sign_flux*num_int_cells, r_tran) + frac_dep_dist

  end do

  ! Set the top values
  frac_dry_flux(map_w2v(1)+nlayers) = 0.0_r_tran
  dep_pts_z(map_w2v(1)+nlayers) = 0.0_r_tran

end subroutine consistent_vertical_deppt_code

end module consistent_vertical_deppt_kernel_mod
