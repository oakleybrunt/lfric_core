!-------------------------------------------------------------------------------
! (c) The copyright relating to this work is owned jointly by the Crown, 
! Met Office and NERC 2014. 
! However, it has been created with the help of the GungHo Consortium, 
! whose members are identified at https://puma.nerc.ac.uk/trac/GungHo/wiki
!-------------------------------------------------------------------------------

!> @brief Collects LFRic methods into a single module.

!> @details Various parts of the LFRic infrastructure are collected together into
!>            a high-level module.

module lfric
  use constants_mod,           only: dp
  use function_space_mod,      only: function_space_type
  use field_mod,               only: field_type
  use kernel_mod,              only: kernel_type
  use basis_function_mod,      only: basis_function_type
  use gaussian_quadrature_mod, only: gaussian_quadrature_type, ngp

end module lfric
