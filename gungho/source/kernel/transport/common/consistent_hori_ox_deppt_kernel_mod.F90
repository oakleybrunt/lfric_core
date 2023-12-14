!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel to compute the horizontal departure distances in the local x
!!        direction, to be consistent with the transport of dry density.
!> @details This code calculates the distance which is swept through a cell
!!          in the x direction during one time step of consistent transport. The
!!          arrival point is the cell face and the departure point is calculated
!!          and stored as a field.
!!          The departure points are a dimensionless displacement, corresponding
!!          to the number of cells moved by a fluid parcel. The part of the dry
!!          flux corresponding to the fractional part of the departure points is
!!          also computed.

module consistent_hori_ox_deppt_kernel_mod

use argument_mod,                only : arg_type,                  &
                                        GH_FIELD, GH_REAL,         &
                                        GH_WRITE, GH_READ,         &
                                        GH_SCALAR, GH_INTEGER,     &
                                        STENCIL, X1D, CELL_COLUMN, &
                                        ANY_DISCONTINUOUS_SPACE_1, &
                                        GH_LOGICAL
use fs_continuity_mod,           only : W3, W2h
use constants_mod,               only : r_tran, i_def, l_def
use kernel_mod,                  only : kernel_type
use log_mod,                     only : log_event,         &
                                        log_scratch_space, &
                                        LOG_LEVEL_ERROR
use reference_element_mod,       only : E, W

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: consistent_hori_ox_deppt_kernel_type
  private
  type(arg_type) :: meta_args(10) = (/                                       &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2h),                          & ! dep_pts
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W2h),                          & ! frac_dry_flux
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2h),                          & ! dry_flux
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3, STENCIL(X1D)),             & ! rho_d_y
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3, STENCIL(X1D)),             & ! rho_d_x
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3, STENCIL(X1D)),             & ! detj_at_w3
       arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! i_start
       arg_type(GH_FIELD,  GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! i_end
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             & ! stencil_extent
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              & ! cap_dep_points
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: consistent_hori_ox_deppt_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: consistent_hori_ox_deppt_code

contains

!> @brief Kernel to computes the horizontal departure distances in the local x
!!        direction, to be consistent with the transport of dry density.
!> @param[in]     nlayers           The number of layers in the mesh
!> @param[in,out] dep_pts_x         Field with dep distances in x-direction
!> @param[in,out] frac_dry_flux     Fractional part of the dry flux in the local
!!                                  x direction, to be computed here
!> @param[in]     dry_flux          The mass flux used in a x-direction
!!                                  transport step for the dry density
!> @param[in]     rho_d_y           The dry density field before the y-direction
!!                                  transport step (following a y-step)
!> @param[in]     stencil_size_y    Size of the stencil for the density field
!> @param[in]     stencil_map_y     Map of DoFs in the stencil for density field
!> @param[in]     rho_d_x           The dry density field before the y-direction
!!                                  transport step (following an x-step)
!> @param[in]     stencil_size_x    Size of the stencil for the density field
!> @param[in]     stencil_map_x     Map of DoFs in the stencil for density field
!> @param[in]     detj_at_w3        det(J) at W3 points
!> @param[in]     stencil_size_detj Size of the stencil for the detj field
!> @param[in]     stencil_map_detj  Map of DoFs in the stencil for detj field
!> @param[in]     i_start           Start index for change in panel ID orientation
!> @param[in]     i_end             End index for change in panel ID orientation
!> @param[in]     stencil_extent    Max number of stencil cells in one direction
!> @param[in]     cap_dep_points    Flag for whether departure points should be
!!                                  capped if they exceed the stencil depth
!> @param[in]     ndf_w2h           Number of DoFs per cell for W2H
!> @param[in]     undf_w2h          Num of W2H DoFs in memory for this partition
!> @param[in]     map_w2h           Map of lowest-cell W2H DoFs
!> @param[in]     ndf_w3            Number of DoFs per cell for W3
!> @param[in]     undf_w3           Num of W3 DoFs in memory for this partition
!> @param[in]     map_w3            Map of lowest-cell W3 DoFs
!> @param[in]     ndf_wp            Number of degrees of freedom for panel ID
!!                                  index function space per cell
!> @param[in]     undf_wp           Number of unique degrees of freedom for
!!                                  panel ID index function space
!> @param[in]     map_wp            Map for panel ID index function space
subroutine consistent_hori_ox_deppt_code( nlayers,             &
                                          dep_pts_x,           &
                                          frac_dry_flux,       &
                                          dry_flux,            &
                                          rho_d_y,             &
                                          stencil_size_y,      &
                                          stencil_map_y,       &
                                          rho_d_x,             &
                                          stencil_size_x,      &
                                          stencil_map_x,       &
                                          detj_at_w3,          &
                                          stencil_size_detj,   &
                                          stencil_map_detj,    &
                                          i_start,             &
                                          i_end,               &
                                          stencil_extent,      &
                                          cap_dep_points,      &
                                          ndf_w2h,             &
                                          undf_w2h,            &
                                          map_w2h,             &
                                          ndf_w3,              &
                                          undf_w3,             &
                                          map_w3,              &
                                          ndf_wp,              &
                                          undf_wp,             &
                                          map_wp )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  logical(kind=l_def), intent(in)    :: cap_dep_points
  integer(kind=i_def), intent(in)    :: ndf_w2h, ndf_w3, ndf_wp
  integer(kind=i_def), intent(in)    :: undf_w2h, undf_w3, undf_wp
  integer(kind=i_def), intent(in)    :: stencil_size_x, stencil_size_y
  integer(kind=i_def), intent(in)    :: stencil_size_detj
  integer(kind=i_def), intent(in)    :: stencil_extent
  integer(kind=i_def), intent(in)    :: map_w2h(ndf_w2h)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)
  integer(kind=i_def), intent(in)    :: map_wp(ndf_wp)
  integer(kind=i_def), intent(in)    :: stencil_map_x(ndf_w3, stencil_size_x)
  integer(kind=i_def), intent(in)    :: stencil_map_y(ndf_w3, stencil_size_y)
  integer(kind=i_def), intent(in)    :: stencil_map_detj(ndf_w3, stencil_size_detj)
  real(kind=r_tran),   intent(in)    :: rho_d_x(undf_w3)
  real(kind=r_tran),   intent(in)    :: rho_d_y(undf_w3)
  integer(kind=i_def), intent(in)    :: i_start(undf_wp)
  integer(kind=i_def), intent(in)    :: i_end(undf_wp)
  real(kind=r_tran),   intent(in)    :: dry_flux(undf_w2h)
  real(kind=r_tran),   intent(in)    :: detj_at_w3(undf_w3)
  real(kind=r_tran),   intent(inout) :: dep_pts_x(undf_w2h)
  real(kind=r_tran),   intent(inout) :: frac_dry_flux(undf_w2h)

  integer(kind=i_def) :: df, j, k, cell_idx, num_int_cells
  integer(kind=i_def) :: sign_flux, offset, df_offset
  real(kind=r_tran)   :: flux_face_i, running_flux_i, flux_cell_j
  real(kind=r_tran)   :: frac_dep_dist, frac_flux_face_i
  integer(kind=i_def) :: half_level, df_ctr, df_idx
  integer(kind=i_def) :: stencil_size, stencil_half, lam_edge_size
  integer(kind=i_def) :: possible_dfs(2), dfs_for_this_column(2)
  real(kind=r_tran)   :: rho_d_local(stencil_size_y)
  real(kind=r_tran)   :: rho_d_x_local(stencil_size_x)
  real(kind=r_tran)   :: rho_d_y_local(stencil_size_y)
  real(kind=r_tran)   :: detj_local(stencil_size_detj)

  ! ========================================================================== !
  ! Determine which horizontal DoFs to loop over
  ! ========================================================================== !
  possible_dfs = (/ W, E /)
  dfs_for_this_column = (/ 0, 0 /)
  half_level = nlayers / 2

  ! If calculation has already happened for this set of faces, don't do it again
  ! This is determined by whether the flux is  still zero
  df_ctr = 0
  do df_idx = 1, SIZE(possible_dfs)
    ! Check if frac_flux values are non-zero:
    ! As fluxes are on shared dofs and have been initialized to zero,
    ! if any flux in the column on the given dof is non-zero then the
    ! fluxes have already been computed and don't need to be computed again. To save
    ! time we only check 2 fluxes - the lowest level and the half domain level.
    if ( frac_dry_flux(map_w2h(possible_dfs(df_idx)) ) == 0.0_r_tran .AND. &
         frac_dry_flux(map_w2h(possible_dfs(df_idx)) + half_level) == 0.0_r_tran ) then
      df_ctr = df_ctr + 1
      dfs_for_this_column(df_ctr) = possible_dfs(df_idx)
    end if
  end do

  ! Set stencil info -----------------------------------------------------------
  ! Use stencil_size_y as each stencil size should be equal
  stencil_size = stencil_size_y
  stencil_half = (stencil_size + 1_i_def) / 2_i_def
  lam_edge_size = 2_i_def*stencil_extent+1_i_def

  ! Calculation near LAM boundaries --------------------------------------------
  if (lam_edge_size > stencil_size) then
    ! Set output to zero
    do df_idx = 1, df_ctr
      df = dfs_for_this_column(df_idx)
      do k = 0, nlayers - 1
        frac_dry_flux(map_w2h(df) + k) = 0.0_r_tran
        dep_pts_x(map_w2h(df) + k) = 0.0_r_tran
      end do
    end do

  ! ========================================================================== !
  ! Calculation in domain interior
  ! ========================================================================== !
  else

    ! Not at edge of LAM so compute fluxes
    ! Loop through horizontal DoFs
    do df_idx = 1, df_ctr
      df = dfs_for_this_column(df_idx)
      df_offset = (df+1)/2 - 1 ! 0 for W, 1 for E

      ! Loop through layers
      do k = 0, nlayers - 1

        flux_face_i = dry_flux(map_w2h(df)+k)

        ! Determine how to step through the stencil
        if (flux_face_i >= 0.0_r_tran) then
          offset = df_offset
          sign_flux = 1
        else
          offset = df_offset-1
          sign_flux = -1
        end if

        ! -------------------------------------------------------------------- !
        ! Fill local arrays
        ! -------------------------------------------------------------------- !
        ! Stencil has order e.g.        | 5 | 4 | 3 | 2 | 1 | 6 | 7 | 8 | 9 | for extent 4
        ! Local fields have order e.g.  | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | for extent 4
        do j = 1, stencil_half
          rho_d_x_local(j) = rho_d_x(stencil_map_x(1,stencil_half+1-j) + k)
          rho_d_y_local(j) = rho_d_y(stencil_map_y(1,stencil_half+1-j) + k)
          detj_local(j) = detj_at_w3(stencil_map_detj(1,stencil_half+1-j) + k)
        end do
        do j = stencil_half+1, stencil_size_y
          rho_d_x_local(j) = rho_d_x(stencil_map_x(1,j) + k)
          rho_d_y_local(j) = rho_d_y(stencil_map_y(1,j) + k)
          detj_local(j) = detj_at_w3(stencil_map_detj(1,j) + k)
        end do

        ! Correction for when a panel boundary has been crossed
        rho_d_local(:) = rho_d_y_local(:)
        do j = i_start(map_wp(1)), i_end(map_wp(1))
          rho_d_local(j) = rho_d_x_local(j)
        end do

        ! -------------------------------------------------------------------- !
        ! Find departure point
        ! -------------------------------------------------------------------- !
        num_int_cells = 0
        running_flux_i = 0.0_r_tran

        ! Step backwards through flux to find the dimensionless departure point:
        ! First, the integer part of departure point (num_int_cells)
        do j = 1, stencil_extent
          ! Index of cell to look at
          cell_idx = stencil_half + offset - sign_flux*j
          flux_cell_j = rho_d_local(cell_idx) * detj_local(cell_idx)

          ! Check if the mass swept up to this cell now exceeds the flux:
          ! if so, the we have found integer number of cells and the departure
          ! point will lie in cell (num_int_cells + 1), so exit do-loop
          if (running_flux_i + flux_cell_j > abs(flux_face_i)) EXIT

          ! Increment running values, if we have not exited loop
          running_flux_i = running_flux_i + flux_cell_j
          num_int_cells = num_int_cells + 1
        end do

        ! Second, the fractional part of departure point
        if (num_int_cells < stencil_extent) then
          ! Departure point is within stencil: determine fractional departure
          ! distance and fractional flux

          ! Set fractional flux. running_flux_i is now the integer flux
          frac_flux_face_i = sign_flux*(abs(flux_face_i) - running_flux_i)

          ! Determine fractional distance
          cell_idx = stencil_half + offset - sign_flux*(num_int_cells + 1)
          frac_dep_dist = frac_flux_face_i / (rho_d_local(cell_idx) * detj_local(cell_idx))

        else if (cap_dep_points) then
          ! Departure point has exceeded stencil, but the user has specified
          ! to cap the departure point so set fractional parts to be zero
          frac_flux_face_i = 0.0_r_tran
          frac_dep_dist = 0.0_r_tran

        else
          ! Departure point has exceeded stencil depth so throw an error
          write(log_scratch_space, '(A,I8)') 'consistent_hori_ox: Consistent ' // &
            'tracer departure points have exceeded stencil depth of ', stencil_extent
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        ! Set the values of the output fields
        frac_dry_flux(map_w2h(df)+k) = frac_flux_face_i
        dep_pts_x(map_w2h(df)+k) = real(sign_flux*num_int_cells, r_tran) + frac_dep_dist
      end do
    end do
  end if

end subroutine consistent_hori_ox_deppt_code

end module consistent_hori_ox_deppt_kernel_mod
