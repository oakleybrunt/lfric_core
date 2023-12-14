!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Computes a correction to horizontal fluxes for consistent transport
!> @details The shifted horizontal fluxes used in consistent conservative
!!          transport equation need an adjustment to correctly capture the
!!          dispersion relation. This kernel computes them.
!!          Only implemented for the lowest-order elements.

module consistent_dispersion_kernel_mod
use argument_mod,            only : arg_type, GH_FIELD, GH_READ,     &
                                    CELL_COLUMN, GH_REAL, GH_WRITE,  &
                                    ANY_SPACE_2
use fs_continuity_mod,       only : W2
use constants_mod,           only : r_tran, r_def, i_def
use kernel_mod,              only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------

type, public, extends(kernel_type) :: consistent_dispersion_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                      &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_SPACE_2), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W2),          &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_SPACE_2)  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: consistent_dispersion_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public consistent_dispersion_code

contains

!> @brief Computes a correction to horizontal fluxes for consistent transport
!> @param[in]     nlayers_shifted    Number of layers in the shifted mesh
!> @param[in,out] flux_X_correction  The flux correction to be computed in
!!                                   W2 on the shifted mesh
!> @param[in]     dry_flux_prime     The dry flux in W2 on the prime mesh
!> @param[in]     theta_factor       0.25*dtheta_dz*dz in shifted W2
!> @param[in]     ndf_w2_shifted     Num of DoFs per cell for shifted W2
!> @param[in]     undf_w2_shifted    Num of DoFs per partition for shifted W2
!> @param[in]     map_w2_shifted     Base cell DoF-map for shifted W2
!> @param[in]     ndf_w2_prime       Num of DoFs per cell for prime W2
!> @param[in]     undf_w2_prime      Num of DoFs per partition for prime W2
!> @param[in]     map_w2_prime       Base cell DoF-map for prime W2
subroutine consistent_dispersion_code( nlayers_shifted,    &
                                       flux_X_correction,  &
                                       dry_flux_prime,     &
                                       theta_factor,       &
                                       ndf_w2_shifted,     &
                                       undf_w2_shifted,    &
                                       map_w2_shifted,     &
                                       ndf_w2_prime,       &
                                       undf_w2_prime,      &
                                       map_w2_prime )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers_shifted
  integer(kind=i_def), intent(in)    :: undf_w2_prime, ndf_w2_prime
  integer(kind=i_def), intent(in)    :: undf_w2_shifted, ndf_w2_shifted
  integer(kind=i_def), intent(in)    :: map_w2_prime(ndf_w2_prime)
  integer(kind=i_def), intent(in)    :: map_w2_shifted(ndf_w2_shifted)
  real(kind=r_tran),   intent(inout) :: flux_X_correction(undf_w2_shifted)
  real(kind=r_tran),   intent(in)    :: dry_flux_prime(undf_w2_prime)
  real(kind=r_tran),   intent(in)    :: theta_factor(undf_w2_shifted)

  ! Internal variables
  integer(kind=i_def) :: k, face

  ! -------------------------------------------------------------------------- !
  ! Average dz/4 * dtheta/dz * delta_F to each face
  ! -------------------------------------------------------------------------- !

  do face = 1, 4
    ! Faces in bottom layer: use dry fluxes from lowest two layers
    k = 0
    flux_X_correction(map_w2_shifted(face)+k) =                                &
      theta_factor(map_w2_shifted(face)+k)                                     &
      * (dry_flux_prime(map_w2_prime(face)+k+1) - dry_flux_prime(map_w2_prime(face)+k))

    ! No contributions to bottom or top layers, so loop over internal layers
    do k = 1, nlayers_shifted - 2
      flux_X_correction(map_w2_shifted(face)+k) =                              &
        theta_factor(map_w2_shifted(face)+k)                                   &
        * (dry_flux_prime(map_w2_prime(face)+k) - dry_flux_prime(map_w2_prime(face)+k-1))
    end do

    ! Faces in top layer: use dry fluxes from highest two layers
    k = nlayers_shifted - 1
    flux_X_correction(map_w2_shifted(face)+k) =                                &
      theta_factor(map_w2_shifted(face)+k)                                     &
      * (dry_flux_prime(map_w2_prime(face)+k-1) - dry_flux_prime(map_w2_prime(face)+k-2))
  end do


end subroutine consistent_dispersion_code

end module consistent_dispersion_kernel_mod