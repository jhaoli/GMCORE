module types_mod

  use mesh_mod
  use parallel_mod

  implicit none

  private

  public allocate_data
  public deallocate_data
  public coef_type
  public state_type
  public static_type
  public tend_type

  type coef_type
    ! Coriolis coefficient at full/half meridional grids
    real, allocatable :: full_f(:)
    real, allocatable :: half_f(:)
    ! Curvature coefficient at full/half meridional grids 
    real, allocatable :: full_c(:)
    real, allocatable :: half_c(:)
    ! Zonal difference coefficient at full/half meridional grids
    real, allocatable :: full_dlon(:)
    real, allocatable :: half_dlon(:)
    ! Meridional difference coefficient at full/half meridional grids
    real, allocatable :: full_dlat(:)
    real, allocatable :: half_dlat(:)
  end type coef_type

  type state_type
    real, allocatable :: u(:,:)
    real, allocatable :: v(:,:)
    real, allocatable :: gd(:,:) ! Geopotential depth
  end type state_type

  type static_type
    real, allocatable :: ghs(:,:) ! Surface geopotential
  end type static_type

  type dia_type
    real, allocatable :: gd_corner(:,:)
    real, allocatable :: pot_vor(:,:)
    real, allocatable :: kinetic_energy(:,:)
    real, allocatable :: energy(:,:)
    real, allocatable :: gd_lon(:,:)
    real, allocatable :: gd_lat(:,:)
    real, allocatable :: normal_lon_flux(:,:)
    real, allocatable :: normal_lat_flux(:,:)
    real, allocatable :: mass_flux_lon_t(:,:)
    real, allocatable :: mass_flux_lat_t(:,:) 
    real, allocatable :: pv_lon(:,:), pv_lat(:,:)
    !
    real, allocatable :: tangent_wind_lon(:,:)
    real, allocatable :: tangent_wind_lat(:,:)
    real, allocatable :: gradient_pv_lon(:,:)
    real, allocatable :: gradient_pv_lat(:,:)
    real, allocatable :: tangent_gradient_pv_lon(:,:)
    real, allocatable :: tangent_gradient_pv_lat(:,:)
  end type dia_type

  type tend_type
    real, allocatable :: u_nonlinear(:,:)
    real, allocatable :: v_nonlinear(:,:)
    real, allocatable :: u_pgf(:,:)
    real, allocatable :: v_pgf(:,:)
    real, allocatable :: mass_div(:,:)
    real, allocatable :: mass_div_lon(:,:)
    real, allocatable :: mass_div_lat(:,:)
    real, allocatable :: du(:,:)
    real, allocatable :: dv(:,:)
    real, allocatable :: dgd(:,:)
    type(dia_type) diag
  end type tend_type

  interface allocate_data
    module procedure allocate_coef_data
    module procedure allocate_static_data
    module procedure allocate_state_data
    module procedure allocate_tend_data
  end interface allocate_data

  interface deallocate_data
    module procedure deallocate_coef_data
    module procedure deallocate_static_data
    module procedure deallocate_state_data
    module procedure deallocate_tend_data
  end interface deallocate_data

contains

  subroutine allocate_coef_data(coef)

    type(coef_type), intent(out) :: coef

    allocate(coef%full_f(mesh%num_full_lat))
    allocate(coef%half_f(mesh%num_half_lat))
    allocate(coef%full_c(mesh%num_full_lat))
    allocate(coef%half_c(mesh%num_half_lat))
    allocate(coef%full_dlon(mesh%num_full_lat))
    allocate(coef%half_dlon(mesh%num_half_lat))
    allocate(coef%full_dlat(mesh%num_full_lat))
    allocate(coef%half_dlat(mesh%num_half_lat))

  end subroutine allocate_coef_data

  subroutine allocate_static_data(static)

    type(static_type), intent(out) :: static

    if (.not. allocated(static%ghs))        call parallel_allocate(static%ghs)

  end subroutine allocate_static_data

  subroutine allocate_state_data(state)

    type(state_type), intent(out) :: state

    if (.not. allocated(state%u))           call parallel_allocate(state%u,             half_lon=.true.)
    if (.not. allocated(state%v))           call parallel_allocate(state%v,             half_lat=.true.)
    if (.not. allocated(state%gd))          call parallel_allocate(state%gd)

  end subroutine allocate_state_data

  subroutine allocate_tend_data(tend)

    type(tend_type), intent(out) :: tend

    if (.not. allocated(tend%u_nonlinear))    call parallel_allocate(tend%u_nonlinear,      half_lon=.true.)
    if (.not. allocated(tend%u_pgf))        call parallel_allocate(tend%u_pgf,          half_lon=.true.)
    if (.not. allocated(tend%du))           call parallel_allocate(tend%du,             half_lon=.true.)
    if (.not. allocated(tend%v_nonlinear))    call parallel_allocate(tend%v_nonlinear,      half_lat=.true.)
    if (.not. allocated(tend%v_pgf))        call parallel_allocate(tend%v_pgf,          half_lat=.true.)
    if (.not. allocated(tend%dv))           call parallel_allocate(tend%dv,             half_lat=.true.)
    if (.not. allocated(tend%mass_div))     call parallel_allocate(tend%mass_div)
    if (.not. allocated(tend%mass_div_lon))     call parallel_allocate(tend%mass_div_lon)
    if (.not. allocated(tend%mass_div_lat))     call parallel_allocate(tend%mass_div_lat)
    if (.not. allocated(tend%dgd))          call parallel_allocate(tend%dgd)
    if (.not. allocated(tend%diag%gd_corner)) call parallel_allocate(tend%diag%gd_corner, half_lon=.true., half_lat=.true.)
    if (.not. allocated(tend%diag%gd_lon)) call parallel_allocate(tend%diag%gd_lon, half_lon=.true.)
    if (.not. allocated(tend%diag%gd_lat)) call parallel_allocate(tend%diag%gd_lat, half_lat=.true.)
    if (.not. allocated(tend%diag%pot_vor)) call parallel_allocate(tend%diag%pot_vor, half_lon=.true., half_lat=.true.)
    if (.not. allocated(tend%diag%kinetic_energy))  call parallel_allocate(tend%diag%kinetic_energy)
    if (.not. allocated(tend%diag%energy))  call parallel_allocate(tend%diag%energy)
    if (.not. allocated(tend%diag%normal_lon_flux))  call parallel_allocate(tend%diag%normal_lon_flux, half_lon=.true.)
    if (.not. allocated(tend%diag%normal_lat_flux))  call parallel_allocate(tend%diag%normal_lat_flux, half_lat=.true.)
    if (.not. allocated(tend%diag%mass_flux_lon_t))  call parallel_allocate(tend%diag%mass_flux_lon_t, half_lon=.true.)
    if (.not. allocated(tend%diag%mass_flux_lat_t))  call parallel_allocate(tend%diag%mass_flux_lat_t, half_lat=.true.)
    if (.not. allocated(tend%diag%pv_lon))    call parallel_allocate(tend%diag%pv_lon, half_lon=.true.)
    if (.not. allocated(tend%diag%pv_lat))    call parallel_allocate(tend%diag%pv_lat, half_lat=.true.) 
  !!
    if (.not. allocated(tend%diag%tangent_wind_lon)) call parallel_allocate(tend%diag%tangent_wind_lon, half_lon=.true.)
    if (.not. allocated(tend%diag%tangent_wind_lat)) call parallel_allocate(tend%diag%tangent_wind_lat, half_lat=.true.)
    if (.not. allocated(tend%diag%gradient_pv_lon)) call parallel_allocate(tend%diag%gradient_pv_lon, half_lon=.true.)
    if (.not. allocated(tend%diag%gradient_pv_lat)) call parallel_allocate(tend%diag%gradient_pv_lat, half_lat=.true.)
    if (.not. allocated(tend%diag%tangent_gradient_pv_lon)) call parallel_allocate(tend%diag%tangent_gradient_pv_lon, half_lon=.true.)
    if (.not. allocated(tend%diag%tangent_gradient_pv_lat)) call parallel_allocate(tend%diag%tangent_gradient_pv_lat, half_lat=.true.)
  end subroutine allocate_tend_data

  subroutine deallocate_coef_data(coef)

    type(coef_type), intent(inout) :: coef

    if (allocated(coef%full_f))    deallocate(coef%full_f)
    if (allocated(coef%half_f))    deallocate(coef%half_f)
    if (allocated(coef%full_c))    deallocate(coef%full_c)
    if (allocated(coef%half_c))    deallocate(coef%half_c)
    if (allocated(coef%full_dlon)) deallocate(coef%full_dlon)
    if (allocated(coef%half_dlon)) deallocate(coef%half_dlon)
    if (allocated(coef%full_dlat)) deallocate(coef%full_dlat)
    if (allocated(coef%half_dlat)) deallocate(coef%half_dlat)

  end subroutine deallocate_coef_data

  subroutine deallocate_static_data(static)

    type(static_type), intent(inout) :: static

    if (allocated(static%ghs)) deallocate(static%ghs)

  end subroutine deallocate_static_data

  subroutine deallocate_state_data(state)

    type(state_type), intent(inout) :: state

    if (allocated(state%u))  deallocate(state%u)
    if (allocated(state%v))  deallocate(state%v)
    if (allocated(state%gd)) deallocate(state%gd)

  end subroutine deallocate_state_data

  subroutine deallocate_tend_data(tend)

    type(tend_type), intent(inout) :: tend

    if (allocated(tend%u_nonlinear))    deallocate(tend%u_nonlinear)
    if (allocated(tend%v_nonlinear))    deallocate(tend%v_nonlinear)
    if (allocated(tend%u_pgf))        deallocate(tend%u_pgf)
    if (allocated(tend%v_pgf))        deallocate(tend%v_pgf)
    if (allocated(tend%mass_div))     deallocate(tend%mass_div)
    if (allocated(tend%mass_div_lon))     deallocate(tend%mass_div_lon)
    if (allocated(tend%mass_div_lat))     deallocate(tend%mass_div_lat)
    if (allocated(tend%du))           deallocate(tend%du)
    if (allocated(tend%dv))           deallocate(tend%dv)
    if (allocated(tend%dgd))          deallocate(tend%dgd)
    if (allocated(tend%diag%gd_corner)) deallocate(tend%diag%gd_corner)
    if (allocated(tend%diag%gd_lon)) deallocate(tend%diag%gd_lon)
    if (allocated(tend%diag%gd_lat)) deallocate(tend%diag%gd_lat)
    if (allocated(tend%diag%pot_vor)) deallocate(tend%diag%pot_vor)
    if (allocated(tend%diag%energy))  deallocate(tend%diag%energy)
    if (allocated(tend%diag%kinetic_energy))  deallocate(tend%diag%kinetic_energy)
    if (allocated(tend%diag%normal_lon_flux))  deallocate(tend%diag%normal_lon_flux)
    if (allocated(tend%diag%normal_lat_flux))  deallocate(tend%diag%normal_lat_flux)
    if (allocated(tend%diag%mass_flux_lon_t))  deallocate(tend%diag%mass_flux_lon_t)
    if (allocated(tend%diag%mass_flux_lat_t))  deallocate(tend%diag%mass_flux_lat_t)
    if (allocated(tend%diag%pv_lon))  deallocate(tend%diag%pv_lon)
    if (allocated(tend%diag%pv_lat))  deallocate(tend%diag%pv_lat)
    !!
    if (allocated(tend%diag%tangent_wind_lon)) deallocate(tend%diag%tangent_wind_lon)
    if (allocated(tend%diag%tangent_wind_lat)) deallocate(tend%diag%tangent_wind_lat)
    if (allocated(tend%diag%gradient_pv_lon)) deallocate(tend%diag%gradient_pv_lon)
    if (allocated(tend%diag%gradient_pv_lat)) deallocate(tend%diag%gradient_pv_lat)
    if (allocated(tend%diag%tangent_gradient_pv_lon)) deallocate(tend%diag%tangent_gradient_pv_lon)
    if (allocated(tend%diag%tangent_gradient_pv_lat)) deallocate(tend%diag%tangent_gradient_pv_lat)
    end subroutine deallocate_tend_data


end module types_mod
