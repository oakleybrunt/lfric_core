!----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Provides an implementation of the Psy layer for physics

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_phys_mod

  use constants_mod,         only : i_def, r_def
  use field_mod,             only : field_type, field_proxy_type
  use mesh_mod,              only : mesh_type

  implicit none
  public

contains
  !---------------------------------------------------------------------
  !> LFRic and PSyclone currently do not have a mechanism to loop over a subset
  !> of cells in a horizontal domain. The relevant PSyclone ticket relating to
  !> this is #487.
  !> The orographic drag kernel only needs to be applied to a subset of land
  !> points where the standard deviation of subgrid orography is more than zero.
  !>
  !> invoke_orographic_drag_kernel: Invokes the kernel which calls the UM
  !> orographic drag scheme only on points where the standard deviation
  !> of subgrid orography is more than zero.
  subroutine invoke_orographic_drag_kernel(                            &
                      du_blk, dv_blk, du_orog_gwd, dv_orog_gwd,        &
                      dtemp_blk, dtemp_orog_gwd, u1_in_w3, u2_in_w3,   &
                      wetrho_in_w3, theta, exner_in_w3, sd_orog,       &
                      grad_xx_orog, grad_xy_orog, grad_yy_orog,        &
                      mr_v, mr_cl, mr_ci,                              &
                      height_w3, height_wth)

    use orographic_drag_kernel_mod, only: orographic_drag_kernel_code
    use mesh_mod, only: mesh_type
    implicit none

    ! Increments from orographic drag
    type(field_type), intent(inout) :: du_blk, dv_blk,           & ! Winds
                                       du_orog_gwd, dv_orog_gwd, &
                                       dtemp_blk, dtemp_orog_gwd   ! Temperature

    ! Inputs to orographic drag scheme
    type(field_type), intent(in) :: u1_in_w3, u2_in_w3,         & ! Winds
                                    wetrho_in_w3, theta,        & ! Density, Temperature
                                    exner_in_w3,                & ! Exner pressure
                                    sd_orog, grad_xx_orog,      & ! Orography ancils
                                    grad_xy_orog, grad_yy_orog, & !
                                    mr_v, mr_cl, mr_ci,         & ! mixing ratios
                                    height_w3, height_wth         ! Heights

    integer :: cell

    ! Number of degrees of freedom
    integer :: ndf_w3, undf_w3, ndf_wtheta, undf_wtheta

    ! These are currently in ANY_SPACE_1
    ! but need to change to DISCONTINOUS_SPACE_1
    ! here and throughout code once LFRic ticket #1968 is on trunk
    integer :: ndf_any_space_1_sd_orog, undf_any_space_1_sd_orog

    integer :: nlayers

    type(field_proxy_type) :: du_blk_proxy, dv_blk_proxy,             &
                              du_orog_gwd_proxy, dv_orog_gwd_proxy,   &
                              dtemp_blk_proxy, dtemp_orog_gwd_proxy,  &
                              u1_in_w3_proxy, u2_in_w3_proxy,         &
                              wetrho_in_w3_proxy, theta_proxy,        &
                              exner_in_w3_proxy,                      &
                              sd_orog_proxy, grad_xx_orog_proxy,      &
                              grad_xy_orog_proxy, grad_yy_orog_proxy, &
                              mr_v_proxy, mr_cl_proxy, mr_ci_proxy,   &
                              height_w3_proxy, height_wth_proxy

    integer, pointer :: map_any_space_1_sd_orog(:,:) => null(), &
                        map_w3(:,:) => null(),                  &
                        map_wtheta(:,:) => null()

    TYPE(mesh_type), pointer :: mesh => null()

    ! Initialise field and/or operator proxies
    du_blk_proxy = du_blk%get_proxy()
    dv_blk_proxy = dv_blk%get_proxy()
    du_orog_gwd_proxy = du_orog_gwd%get_proxy()
    dv_orog_gwd_proxy = dv_orog_gwd%get_proxy()
    dtemp_blk_proxy = dtemp_blk%get_proxy()
    dtemp_orog_gwd_proxy = dtemp_orog_gwd%get_proxy()
    u1_in_w3_proxy = u1_in_w3%get_proxy()
    u2_in_w3_proxy = u2_in_w3%get_proxy()
    wetrho_in_w3_proxy = wetrho_in_w3%get_proxy()
    theta_proxy = theta%get_proxy()
    exner_in_w3_proxy = exner_in_w3%get_proxy()
    sd_orog_proxy = sd_orog%get_proxy()
    grad_xx_orog_proxy = grad_xx_orog%get_proxy()
    grad_xy_orog_proxy = grad_xy_orog%get_proxy()
    grad_yy_orog_proxy = grad_yy_orog%get_proxy()
    mr_v_proxy = mr_v%get_proxy()
    mr_cl_proxy = mr_cl%get_proxy()
    mr_ci_proxy = mr_ci%get_proxy()
    height_w3_proxy = height_w3%get_proxy()
    height_wth_proxy = height_wth%get_proxy()

    ! Initialise number of layers
    nlayers = du_blk_proxy%vspace%get_nlayers()

    ! Create a mesh object
    mesh => du_blk%get_mesh()

    ! Look-up dofmaps for each function space
    map_w3 => du_blk_proxy%vspace%get_whole_dofmap()
    map_wtheta => dtemp_blk_proxy%vspace%get_whole_dofmap()
    map_any_space_1_sd_orog => sd_orog_proxy%vspace%get_whole_dofmap()

    ! Initialise number of DoFs for w3
    ndf_w3 = du_blk_proxy%vspace%get_ndf()
    undf_w3 = du_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for wtheta
    ndf_wtheta = dtemp_blk_proxy%vspace%get_ndf()
    undf_wtheta = dtemp_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for any_space_1_sd_orog
    ndf_any_space_1_sd_orog = sd_orog_proxy%vspace%get_ndf()
    undf_any_space_1_sd_orog = sd_orog_proxy%vspace%get_undf()

    ! Call kernels and communication routines
    if (sd_orog_proxy%is_dirty(depth=1)) then
      call sd_orog_proxy%halo_exchange(depth=1)
    end if

    if (grad_xx_orog_proxy%is_dirty(depth=1)) then
      call grad_xx_orog_proxy%halo_exchange(depth=1)
    end if

    if (grad_xy_orog_proxy%is_dirty(depth=1)) then
      call grad_xy_orog_proxy%halo_exchange(depth=1)
    end if

    if (grad_yy_orog_proxy%is_dirty(depth=1)) then
      call grad_yy_orog_proxy%halo_exchange(depth=1)
    end if

    ! Loop over cells
    do cell=1,mesh%get_last_edge_cell()
      ! Only call orographic_drag_kernel_code at points where the
      ! standard deviation of the subgrid orography is more than zero.
      if ( sd_orog_proxy%data(map_any_space_1_sd_orog(1, cell)) > 0.0_r_def ) then

        call orographic_drag_kernel_code(                                  &
                      nlayers, du_blk_proxy%data, dv_blk_proxy%data,       &
                      du_orog_gwd_proxy%data, dv_orog_gwd_proxy%data,      &
                      dtemp_blk_proxy%data, dtemp_orog_gwd_proxy%data,     &
                      u1_in_w3_proxy%data, u2_in_w3_proxy%data,            &
                      wetrho_in_w3_proxy%data, theta_proxy%data,           &
                      exner_in_w3_proxy%data,                              &
                      sd_orog_proxy%data, grad_xx_orog_proxy%data,         &
                      grad_xy_orog_proxy%data, grad_yy_orog_proxy%data,    &
                      mr_v_proxy%data, mr_cl_proxy%data, mr_ci_proxy%data, &
                      height_w3_proxy%data, height_wth_proxy%data,         &
                      ndf_w3, undf_w3, map_w3(:,cell),                     &
                      ndf_wtheta, undf_wtheta, map_wtheta(:,cell),         &
                      ndf_any_space_1_sd_orog, undf_any_space_1_sd_orog,   &
                      map_any_space_1_sd_orog(:,cell) )

      end if ! sd_orog_proxy%data(map_any_space_1_sd_orog(1, cell)) > 0.0_r_def

    end do


    ! Set halos dirty/clean for fields modified in the above loop
    call du_blk_proxy%set_dirty()
    call dv_blk_proxy%set_dirty()
    call du_orog_gwd_proxy%set_dirty()
    call dv_orog_gwd_proxy%set_dirty()
    call dtemp_blk_proxy%set_dirty()
    call dtemp_orog_gwd_proxy%set_dirty()

  end subroutine invoke_orographic_drag_kernel

end module psykal_lite_phys_mod
