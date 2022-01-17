!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief UKCA initialisation subroutine for UM science configuration

module um_ukca_init_mod

  use log_mod, only : log_event,                                               &
                      log_scratch_space,                                       &
                      LOG_LEVEL_ERROR, LOG_LEVEL_INFO

  use constants_mod, only : r_um, i_um

  ! LFRic namelists which have been read
  use aerosol_config_mod,        only: glomap_mode,                            &
                                       glomap_mode_ukca
  use section_choice_config_mod, only: aerosol, aerosol_um

  ! UM modules used
  use nlsizes_namelist_mod, only: bl_levels, row_length, rows, model_levels
  use timestep_mod,         only: timestep
  use cv_run_mod,           only: l_param_conv

  ! JULES modules used 
  use jules_surface_types_mod, only: ntype, npft,                              &
                                     brd_leaf, ndl_leaf,                       &
                                     c3_grass, c4_grass,                       &
                                     shrub, urban, lake, soil, ice
  use jules_soil_mod,          only: dzsoil_io

  ! UKCA API module
  use ukca_api_mod, only: ukca_setup,                                          &
                          ukca_get_ntp_varlist,                                &
                          ukca_get_tracer_varlist,                             &
                          ukca_get_envgroup_varlists,                          &
                          ukca_get_emission_varlist,                           &
                          ukca_register_emission,                              &
                          ukca_chem_offline,                                   &
                          ukca_age_reset_by_level,                             &
                          ukca_activation_arg,                                 &
                          ukca_maxlen_fieldname,                               &
                          ukca_maxlen_message,                                 &
                          ukca_maxlen_procname,                                &
                          ukca_maxlen_emiss_long_name,                         &
                          ukca_maxlen_emiss_tracer_name,                       &
                          ukca_maxlen_emiss_var_name

  implicit none

  private
  public :: um_ukca_init

  ! Number of entrainment levels considered by UM routine tr_mix
  integer(i_um), parameter, public :: nlev_ent_tr_mix = 3

  ! Number of emission field entries required for the UKCA configuration
  integer(i_um) :: n_emiss_slots

  ! UKCA field names for environmental drivers (subset used in this module):
  ! GLOMAP-specific drivers in full-height group
  character(len=*), parameter, public :: fldname_autoconv = 'autoconv'
  character(len=*), parameter, public :: fldname_accretion = 'accretion'
  character(len=*), parameter, public :: fldname_rim_cry = 'rim_cry'
  character(len=*), parameter, public :: fldname_rim_agg = 'rim_agg'
  ! Activate-specific drivers in boundary levels group
  character(len=*), parameter, public :: fldname_bl_tke = 'bl_tke'

  ! List of tracers required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public :: tracer_names(:) =>  &
                                                           null()

  ! List of non-transported prognostics required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public :: ntp_names(:) =>     &
                                                           null()

  ! Lists of environmental driver fields required for the UKCA configuration
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_scalar_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_integer(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flat_logical(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_flatpft_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullht_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullht0_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_fullhtp1_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_bllev_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_entlev_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_land_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landtile_real(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landtile_logical(:) => null()
  character(len=ukca_maxlen_fieldname), pointer, public ::                     &
    env_names_landpft_real(:) => null()

  ! Lists of emissions for the UKCA configuration
  character(len=ukca_maxlen_emiss_tracer_name), pointer ::                     &
    emiss_names(:) => null() 
  character(len=ukca_maxlen_emiss_tracer_name), allocatable, public ::         &
    emiss_names_flat(:)
  character(len=ukca_maxlen_emiss_tracer_name), allocatable, public ::         &
    emiss_names_fullht(:)


contains

  subroutine um_ukca_init()
  !> @brief Set up the UKCA model

  implicit none

    if ( aerosol == aerosol_um .and. glomap_mode == glomap_mode_ukca ) then

        call aerosol_ukca_init( row_length, rows, model_levels, bl_levels,     &
                                ntype, npft, brd_leaf, ndl_leaf,               &
                                c3_grass, c4_grass, shrub,                     &
                                urban, lake, soil, ice,                        &
                                dzsoil_io(1), timestep, l_param_conv )

    end if

  end subroutine um_ukca_init


  !> @brief Set up UKCA for GLOMAP-mode with Offline Oxidants chemistry
  !> @details Configure UKCA using options generally consistent with a GA7 run.
  !>          Where possible, all values are taken from GA7 science settings.
  !>          Also register each emission field to be supplied to UKCA and
  !>          obtain lists of UKCA tracers, non-transported prognostics and
  !>          environmental drivers required for the selected configuration.
  !> @param[in] row_length      X dimension of model grid 
  !> @param[in] rows            Y dimension of model grid
  !> @param[in] model_levels    Z dimension of model grid
  !> @param[in] bl_levels       No. of Z levels in boundary layer
  !> @param[in] ntype           No. of surface types considered in interactive dry deposition
  !> @param[in] npft            No. of plant functional types
  !> @param[in] i_brd_leaf      Index of surface type 'broad-leaf tree'
  !> @param[in] i_ndl_leaf      Index of surface type 'needle-leaf tree'
  !> @param[in] i_c3_grass      Index of surface type 'c3 grass'
  !> @param[in] i_c4_grass      Index of surface type 'c4 grass'
  !> @param[in] i_shrub         Index of surface type 'shrub'
  !> @param[in] i_urban         Index of surface type 'urban'
  !> @param[in] i_lake          Index of surface type 'lake'
  !> @param[in] i_soil          Index of surface type 'soil'
  !> @param[in] i_ice           Index of surface type 'ice'
  !> @param[in] dzsoil_layer1   Thickness of surface soil layer (m)
  !> @param[in] timestep        Model time step (s)
  !> @param[in] l_param_conv    True if convection is parameterized
  
  subroutine aerosol_ukca_init( row_length, rows, model_levels, bl_levels,     &
                                ntype, npft, i_brd_leaf, i_ndl_leaf,           &
                                i_c3_grass, i_c4_grass, i_shrub,               &
                                i_urban, i_lake, i_soil, i_ice,                &
                                dzsoil_layer1, timestep, l_param_conv )

    ! UM module used
    use dms_flux_mod_4a, only: i_liss_merlivat

    ! UM modules containing thing that needs setting
    use mphys_diags_mod, ONLY: l_praut_diag, l_pracw_diag, l_piacw_diag,       &
                               l_psacw_diag
    use bl_diags_mod, ONLY: bl_diag

    implicit none

    integer, intent(in) :: row_length      ! X dimension of model grid 
    integer, intent(in) :: rows            ! Y dimension of model grid 
    integer, intent(in) :: model_levels    ! Z dimension of model grid
    integer, intent(in) :: bl_levels       ! No. of Z levels in boundary layer 
    integer, intent(in) :: ntype           ! No. of surface types considered in
                                           ! interactive dry deposition 
    integer, intent(in) :: npft            ! No. of plant functional types
    integer, intent(in) :: i_brd_leaf      ! Index of type 'broad-leaf tree'
    integer, intent(in) :: i_ndl_leaf      ! Index of type 'needle-leaf tree'  
    integer, intent(in) :: i_c3_grass      ! Index of type 'c3 grass'
    integer, intent(in) :: i_c4_grass      ! Index of type 'c4 grass'
    integer, intent(in) :: i_shrub         ! Index of type 'shrub'
    integer, intent(in) :: i_urban         ! Index of type 'urban'
    integer, intent(in) :: i_lake          ! Index of type 'lake'
    integer, intent(in) :: i_soil          ! Index of type 'soil'
    integer, intent(in) :: i_ice           ! Index of type 'ice'
    real(r_um), intent(in) :: dzsoil_layer1! Thickness of surface soil layer (m)
    real(r_um), intent(in) :: timestep     ! Model time step (s)
    logical, intent(in) :: l_param_conv    ! True if convection is parameterized

    ! Local variables

    integer(i_um) :: emiss_id
    integer :: n_emissions
    integer :: n
    integer :: n_previous
    integer :: i

    logical :: l_three_dim

    character(len=ukca_maxlen_emiss_long_name) :: long_name 
    character(len=ukca_maxlen_emiss_tracer_name), allocatable :: tmp_names(:)

    character(len=ukca_maxlen_emiss_var_name) :: field_varname
    character(len=*), parameter  :: emiss_units = 'kg m-2 s-1' 

    ! Variables for UKCA error handling
    integer :: ukca_errcode
    character(len=ukca_maxlen_message) :: ukca_errmsg
    character(len=ukca_maxlen_procname) :: ukca_errproc

    ! Set up proto-GA configuration based on GA7.
    ! The ASAD Newton-Raphson Offline Oxidants scheme (ukca_chem_offline)
    ! is substituted for the Explicit backward-Euler scheme used in GA7
    ! which cannot be called by columns. Hence, configuration values for
    ! nrsteps and l_ukca_asad_columns are included.
    ! In lieu of plume scavenging, which is not yet available in LFRic,
    ! convective rainout is switched on (l_cv_rainout = .true.).
    ! Other configuration values, with the exception of temporary logicals
    ! specifiying fixes, are set to match GA7 science settings or taken
    ! from the LFRic context. Unlike GA7, all fixes that are controlled by
    ! temporary logicals will be on. (i.e. the defaults for these
    ! logicals, .true. by convention in UKCA, are not overridden.)

    call ukca_setup( ukca_errcode,                                             &
           ! Context information
           row_length=row_length,                                              &
           rows=rows,                                                          &
           model_levels=model_levels,                                          &
           bl_levels=bl_levels,                                                &
           nlev_ent_tr_mix=nlev_ent_tr_mix,                                    &
           ntype=ntype,                                                        &
           npft=npft,                                                          &
           i_brd_leaf=i_brd_leaf,                                              &
           i_ndl_leaf=i_ndl_leaf,                                              &
           i_c3_grass=i_c3_grass,                                              &
           i_c4_grass=i_c4_grass,                                              &
           i_shrub=i_shrub,                                                    &
           i_urban=i_urban,                                                    &
           i_lake=i_lake,                                                      &
           i_soil=i_soil,                                                      &
           i_ice=i_ice,                                                        &
           dzsoil_layer1=dzsoil_layer1,                                        &
           timestep=timestep,                                                  &
           ! General UKCA configuration options
           i_ukca_chem=ukca_chem_offline,                                      &
           l_ukca_mode=.true.,                                                 &
           l_fix_tropopause_level=.true.,                                      &
           l_ukca_persist_off=.true.,                                          &
           ! Chemistry configuration options
           i_ukca_chem_version=117,                                            &
           chem_timestep=3600,                                                 &
           nrsteps=45,                                                         &
           l_ukca_asad_columns=.true.,                                         &
           l_ukca_intdd=.true.,                                                &
           l_ukca_ddep_lev1=.false.,                                           &
           l_ukca_ddepo3_ocean=.false.,                                        &
           l_ukca_dry_dep_so2wet=.false.,                                      &
           ! UKCA emissions configuration options
           mode_parfrac=2.5_r_um,                                              &
           l_ukca_enable_seadms_ems=.true.,                                    &
           i_ukca_dms_flux=i_liss_merlivat,                                    &
           l_ukca_scale_seadms_ems=.false.,                                    &
           l_ukca_scale_soa_yield=.true.,                                      &
           soa_yield_scaling=2.0_r_um,                                         &
           l_support_ems_vertprof=.true.,                                      &
           ! UKCA environmental driver configuration options
           l_param_conv=l_param_conv,                                          &
           l_ctile=.true.,                                                     &
           ! General GLOMAP configuration options
           i_mode_nzts=15,                                                     &
           i_mode_setup=8,                                                     &
           l_mode_bhn_on=.true.,                                               &
           l_mode_bln_on=.false.,                                              &
           i_mode_nucscav=3,                                                   &
           mode_activation_dryr=37.5_r_um,                                     &
           mode_incld_so2_rfrac=0.25_r_um,                                     &
           l_cv_rainout=.true.,                                                &
           ! GLOMAP emissions configuration options
           l_ukca_primsu=.true.,                                               &
           l_ukca_primss=.true.,                                               &
           l_ukca_primdu=.true.,                                               &
           l_ukca_primbcoc=.true.,                                             &
           l_bcoc_bf=.true.,                                                   &
           l_bcoc_bm=.true.,                                                   &
           l_bcoc_ff=.true.,                                                   &
           l_ukca_scale_biom_aer_ems=.true.,                                   &
           biom_aer_ems_scaling=2.0_r_um,                                      &
           ! GLOMAP feedback configuration options
           l_ukca_radaer=.true.,                                               &
           i_ukca_activation_scheme=ukca_activation_arg,                       &
           i_ukca_nwbins=20,                                                   &
           ! Return status information
           error_message=ukca_errmsg,                                          &
           error_routine=ukca_errproc)

    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Retrieve the lists of required fields for the configuration
    CALL get_ukca_field_lists()

    ! Switch on optional UM microphysics diagnostics required by UKCA
    if (any(env_names_fullht_real(:) == fldname_autoconv))                     &
      l_praut_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_accretion))                    &
      l_pracw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_cry))                      &
      l_piacw_diag = .true.
    if (any(env_names_fullht_real(:) == fldname_rim_agg))                      &
      l_psacw_diag = .true.

    ! Switch on optional UM boundary layer diagnostics required by UKCA
    if (any(env_names_bllev_real(:) == fldname_bl_tke))                        &
      bl_diag%l_request_tke = .true.

    ! Emissions registration:
    ! One emission field is provided for each active emission species (as for
    ! GA7).
    ! Register 2D emissions first, followed by 3D emissions; emissions
    ! must be registered in order of dimensionality for compatibility with
    ! the UKCA time step call which expects an array of 2D emission fields
    ! and an array of 3D emission fields. These arrays are expected to
    ! contain fields corresponding to consecutive emission id numbers, 
    ! starting at 1, with the 3D fields having the highest numbers.
    ! The names of the emission species corresponding to each field array
    ! are given in separate reference arrays created below.

    n_emissions = size(emiss_names)
    allocate(tmp_names(n_emissions)) 

    ! Register 2D emissions
    n = 0    
    do i = 1, n_emissions
      if (.NOT.( emiss_names(i) == 'SO2_nat' .or.                              &
                 emiss_names(i) == 'BC_biomass' .or.                           &
                 emiss_names(i) == 'OM_biomass' )) then
        field_varname = 'emissions_' // emiss_names(i)
        l_three_dim = .false.
        ! Assign long name as in GA7 ancil file.
        ! *** WARNING ***
        ! This is currently hard-wired in the absence of functionality
        ! to handle the NetCDF variable attributes that provide the long names
        ! in the UM. There is therefore no guarantee of compatibility with
        ! arbitrary ancil files and care must be taken to ensure that the
        ! data in the ancil files provided are consistent with the names
        ! defined here.
        select case(emiss_names(i))
        case('BC_biofuel')
          long_name = 'BC biofuel surf emissions'
        case('BC_fossil')
          long_name = 'BC fossil fuel surf emissions'
        case('DMS')
          long_name = 'DMS emissions expressed as sulfur'
        case('Monoterp')
          long_name = 'Monoterpene surf emissions expressed as carbon'
        case('OM_biofuel')
          long_name = 'OC biofuel surf emissions expressed as carbon'
        case('OM_fossil')
          long_name = 'OC fossil fuel surf emissions expressed as carbon'
        case('SO2_low')
          long_name = 'SO2 low level emissions expressed as sulfur'
        case('SO2_high')
          long_name = 'SO2 high level emissions expressed as sulfur'
        end select
        if (emiss_names(i) == 'SO2_high') then
          call ukca_register_emission( n_emiss_slots, field_varname,           &
                                       emiss_names(i), emiss_units,            &
                                       l_three_dim, emiss_id,                  &
                                       long_name=long_name,                    &
                                       vert_fact='high_level',                 &
                                       lowest_lev=8, highest_lev=8 )
        else
          call ukca_register_emission( n_emiss_slots, field_varname,           &
                                       emiss_names(i), emiss_units,            &
                                       l_three_dim, emiss_id,                  &
                                       long_name=long_name,                    &
                                       vert_fact='surface' )
        end if
        n = n + 1
        if (emiss_id /= int( n, i_um )) then
          write( log_scratch_space, '(A,I0,A,I0)' )                            &
            'Unexpected id (', emiss_id,                                       &
            ') assigned on registering a 2D UKCA emission. Expected ', n
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        tmp_names(n) = emiss_names(i)  
      end if
    end do

    ! Create reference array of 2D emission names
    allocate(emiss_names_flat(n))
    do i = 1, n
      emiss_names_flat(i) = tmp_names(i)
    end do

    ! Register 3D emissions
    n_previous = n
    n = 0    
    do i = 1, n_emissions 
      if ( emiss_names(i) == 'SO2_nat' .or.                                    &
           emiss_names(i) == 'BC_biomass' .or.                                 &
           emiss_names(i) == 'OM_biomass' ) then
        field_varname = 'emissions_' // emiss_names(i)
        l_three_dim = .true.
        ! Assign long name as in GA7 ancil file.
        ! *** WARNING ***
        ! This is currently hard-wired in the absence of functionality
        ! to handle the NetCDF variable attributes. There is therefore
        ! no guarantee of compatibility with arbitrary ancil files.
        select case(emiss_names(i))
        case('BC_biomass')
          long_name = 'BC biomass 3D emissions'
        case('OM_biomass')
          long_name = 'OC biomass 3D emissions expressed as carbon'
        case('SO2_nat')
          long_name = 'SO2 natural emissions expressed as sulfur'
        end select
        call ukca_register_emission( n_emiss_slots, field_varname,             &
                                     emiss_names(i), emiss_units, l_three_dim, &
                                     emiss_id, long_name=long_name,            &
                                     vert_fact='all_levels' )
        n = n + 1
        if (emiss_id /= int( n + n_previous, i_um )) then
          write( log_scratch_space, '(A,I0,A,I0)' )                            &
            'Unexpected id (', emiss_id,                                       &
            ') assigned on registering 3D UKCA emission. Expected ',n 
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        tmp_names(n) = emiss_names(i)  
      end if
    end do

    ! Create reference array of 3D emission names
    allocate(emiss_names_fullht(n))
    do i = 1, n
      emiss_names_fullht(i) = tmp_names(i)
    end do

  end subroutine aerosol_ukca_init


  !>@brief Get the lists of fields required for the UKCA configuration
  subroutine get_ukca_field_lists()

    implicit none

    ! Local variables

    integer :: i
    integer :: n 
    integer :: n_tot 

    ! Number of emission entries for each emitted species
    integer, pointer :: n_slots(:) => null()   

    ! Variables for UKCA error handling
    integer :: ukca_errcode
    character(len=ukca_maxlen_message) :: ukca_errmsg
    character(len=ukca_maxlen_procname) :: ukca_errproc

    ! Get list of tracers required by the current UKCA configuration
    CALL ukca_get_tracer_varlist( tracer_names, ukca_errcode,                  &
                                  error_message=ukca_errmsg,                   &
                                  error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'Tracers required (', size(tracer_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(tracer_names)
      write( log_scratch_space, '(A)' ) tracer_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Get list of NTPs required by the current UKCA configuration
    CALL ukca_get_ntp_varlist( ntp_names, ukca_errcode,                        &
                               error_message=ukca_errmsg,                      &
                               error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'NTPs required (', size(ntp_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(ntp_names)
      write( log_scratch_space, '(A)' ) ntp_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Get lists of environmental drivers required by the current UKCA
    ! configuration. 
    ! Note that these group lists derived from UKCA's master list must
    ! be requested here to force UKCA to set up the group lists prior to
    ! their use in kernel calls. (Setup within kernels is not thread-safe.) 

    CALL ukca_get_envgroup_varlists(                                           &
           ukca_errcode,                                                       &
           varnames_scalar_real_ptr=env_names_scalar_real,                     &
           varnames_flat_integer_ptr=env_names_flat_integer,                   &
           varnames_flat_real_ptr=env_names_flat_real,                         &
           varnames_flat_logical_ptr=env_names_flat_logical,                   &
           varnames_flatpft_real_ptr=env_names_flatpft_real,                   &
           varnames_fullht_real_ptr=env_names_fullht_real,                     &
           varnames_fullht0_real_ptr=env_names_fullht0_real,                   &
           varnames_fullhtp1_real_ptr=env_names_fullhtp1_real,                 &
           varnames_bllev_real_ptr=env_names_bllev_real,                       &
           varnames_entlev_real_ptr=env_names_entlev_real,                     &
           varnames_land_real_ptr=env_names_land_real,                         &
           varnames_landtile_real_ptr=env_names_landtile_real,                 &
           varnames_landtile_logical_ptr=env_names_landtile_logical,           &
           varnames_landpft_real_ptr=env_names_landpft_real,                   &
           error_message=ukca_errmsg, error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    n_tot = 0

    n = size(env_names_scalar_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in scalar real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_scalar_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_integer)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat integer group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_flat_integer(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_flat_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flat_logical)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat logical group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_flat_logical(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_flatpft_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in flat PFT real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_flatpft_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullht_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_fullht_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullht0_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height + level 0 real group required (',  &
      n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_fullht0_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_fullhtp1_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in full-height + 1 real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_fullhtp1_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_bllev_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in boundary layer real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_bllev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_entlev_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in entrainment levels real group required (',     &
      n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_entlev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_land_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_land_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landtile_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land tile real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_landtile_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landtile_logical)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land tile logical group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_landtile_logical(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    n = size(env_names_landpft_real)
    write( log_scratch_space, '(A,I0,A)' )                                     &
      'Environmental drivers in land PFT real group required (', n, '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, n 
      write( log_scratch_space, '(A)' ) env_names_landpft_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do
    n_tot = n_tot + n

    write( log_scratch_space, '(A,I0)' )                                      &
      'Total number of environmental drivers required: ', n_tot
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Get list of active offline emissions species in the current UKCA
    ! configuration
    CALL ukca_get_emission_varlist( emiss_names, n_slots, ukca_errcode,        &
                                    error_message=ukca_errmsg,                 &
                                    error_routine=ukca_errproc )
    if (ukca_errcode /= 0) then
      write( log_scratch_space, '(A,I0,A,A,A,A)' )                             &
           'UKCA error ', ukca_errcode, ' in ', ukca_errproc, ': ',            &
           ukca_errmsg
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if
    write( log_scratch_space, '(A,I0,A)' )                                     &
         'Offline emissions active (', size(emiss_names), '):'
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    do i = 1, size(emiss_names)
      write( log_scratch_space, '(A)' ) emiss_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    end do

    ! Determine number of offline emission entries required.
    ! Allow for one emission field for each active emission species
    ! (consistent with GA7 requirements).
    n_emiss_slots = 0
    do i = 1, size(emiss_names)
      n_emiss_slots = n_emiss_slots + n_slots(i) 
    end do

  end subroutine get_ukca_field_lists

end module um_ukca_init_mod
