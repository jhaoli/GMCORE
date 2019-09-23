module history_mod

  use io_mod
  use log_mod
  use mesh_mod
  use params_mod
  use parallel_mod
  use types_mod
  use diag_mod
  use string_mod

  implicit none

  private

  public history_init
  public history_final
  public history_write

  ! A-grid velocity
  real, allocatable :: u(:,:)
  real, allocatable :: v(:,:)
  real, allocatable :: gh(:,:)

  interface history_write
    module procedure history_write_state
    module procedure history_write_tendency
  end interface history_write

contains

  subroutine history_init()

    call io_create_dataset(desc=case_desc, file_prefix=case_name // '.h0')
    call io_add_meta('time_step_size', time_step_size)
    call io_add_meta('time_scheme', time_scheme)
    call io_add_meta('split_scheme', split_scheme)
    call io_add_meta('subcycles', subcycles)
    call io_add_dim('time', add_var=.true.)
    call io_add_dim('lon',  size=mesh%num_full_lon, add_var=.true.)
    call io_add_dim('lat',  size=mesh%num_full_lat, add_var=.true.)
    call io_add_dim('ilon', size=mesh%num_half_lon, add_var=.true.)
    call io_add_dim('ilat', size=mesh%num_half_lat, add_var=.true.)
    call io_add_var('u',    long_name='u wind component',     units='m s-1',  dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('v',    long_name='v wind component',     units='m s-1',  dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('gh',   long_name='geopotential height',  units='m2 s-2', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('ghs',  long_name='surface geopotential', units='m2 s-2', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('vor',  long_name='relative vorticity',   units='s-1',    dim_names=['ilon', 'ilat', 'time'])
    call io_add_var('div',  long_name='divergence',           units='s-1',    dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('pv',   long_name='potential vorticity',  units='s-1 m-1',dim_names=['ilon', 'ilat', 'time'])
    call io_add_var('te',   long_name='total energy',         units='m4 s-4', dim_names=['time'])
    call io_add_var('tm',   long_name='total mass',           units='m2 s-2', dim_names=['time'])
    call io_add_var('tes',  long_name='total potential enstrophy',units='m-1 s-2',dim_names=['time'])
    call io_add_var('tav',  long_name='total absolute vorticity',units='m s-2',dim_names=['time'])

    call io_create_dataset(name='debug', desc=case_desc, file_prefix=case_name // '.debug')
    call io_add_dim('time', 'debug', add_var=.true.)
    call io_add_dim('lon',  'debug', size=mesh%num_full_lon, add_var=.true.)
    call io_add_dim('lat',  'debug', size=mesh%num_full_lat, add_var=.true.)
    call io_add_dim('ilon', 'debug', size=mesh%num_half_lon, add_var=.true.)
    call io_add_dim('ilat', 'debug', size=mesh%num_half_lat, add_var=.true.)
!     call io_add_var('u_adv_lon',    'debug', long_name='u_adv_lon',           units='', dim_names=['ilon', 'lat ', 'time'])
!     call io_add_var('u_adv_lat',    'debug', long_name='u_adv_lat',           units='', dim_names=['ilon', 'lat ', 'time'])
!     call io_add_var('v_adv_lon',    'debug', long_name='v_adv_lon',           units='', dim_names=['lon ', 'ilat', 'time'])
!     call io_add_var('v_adv_lat',    'debug', long_name='v_adv_lat',           units='', dim_names=['lon ', 'ilat', 'time'])
!     call io_add_var('fv',           'debug', long_name='fv',                  units='', dim_names=['ilon', 'lat ', 'time'])
!     call io_add_var('fu',           'debug', long_name='fu',                  units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('dEdlon',  'debug', long_name='zonal energy gradient',    units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('dEdlat', 'debug', long_name='meridional energy gradient',units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('mass_div',     'debug', long_name='mass_div',            units='', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('du',           'debug', long_name='du',                  units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('dv',           'debug', long_name='dv',                  units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('dgd',          'debug', long_name='dgd',                 units='', dim_names=['lon ', 'lat ', 'time'])
!     call io_add_var('iap_u',        'debug', long_name='IAP u',               units='', dim_names=['ilon', 'lat ', 'time'])
!     call io_add_var('iap_v',        'debug', long_name='IAP v',               units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('gd_corner',    'debug', long_name='potential thickness corner', units='', dim_names=['ilon', 'ilat', 'time'])
    call io_add_var('u_nonlinear',  'debug', long_name='u_none_linear_force', units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('v_nonlinear',  'debug', long_name='v_none_linear_force', units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('normal_lon_flux', 'debug', long_name='normal u flux',    units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('normal_lat_flux', 'debug', long_name='normal v flux',    units='', dim_names=['lon ', 'ilat', 'time'])
    ! call io_add_var('center_energy', 'debug', long_name='cell center energy', units='', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('PV',           'debug', long_name='potential vorticity', units='', dim_names=['ilon', 'ilat', 'time'])
    call io_add_var('PV_on_lon',    'debug', long_name='PV on longitude',     units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('PV_on_lat',    'debug', long_name='PV on latitude',      units='', dim_names=['lon ', 'ilat', 'time'])
    call io_add_var('kinetic_energy','debug', long_name='kinetic energy',     units='', dim_names=['lon ', 'lat ', 'time'])
    call io_add_var('mass_flux_lon_t','debug', long_name='tangential_flux on u',units='', dim_names=['ilon', 'lat ', 'time'])
    call io_add_var('mass_flux_lat_t','debug', long_name='tangential_flux on v',units='', dim_names=['lon ', 'ilat', 'time'])
    !!
!     call io_add_var('tangent_wind_lon', 'debug', long_name='tangent wind lon', units='ms-1', dim_names=['ilon', 'lat ', 'time'])
!     call io_add_var('tangent_wind_lat', 'debug', long_name='tangent wind lat', units='ms-1', dim_names=['lon ', 'ilat', 'time'])
!     call io_add_var('gradient_pv_lon', 'debug', long_name='gradient pv lon', units='m-1s-1', dim_names=['ilon', 'lat ', 'time'])
!     call io_add_var('gradient_pv_lat', 'debug', long_name='gradient pv lat', units='m-1s-1', dim_names=['lon ', 'ilat', 'time'])
!     call io_add_var('tangent_gradient_pv_lon', 'debug', long_name='tangent_gradient_pv_lon',units='', dim_names=['ilon', 'lat ', 'time'])
!     call io_add_var('tangent_gradient_pv_lat', 'debug', long_name='tangent_gradient_pv_lat',units='', dim_names=['lon ', 'ilat', 'time'])
    if (.not. allocated(u)) call parallel_allocate(u)
    if (.not. allocated(v)) call parallel_allocate(v)
    if (.not. allocated(gh)) call parallel_allocate(gh)

    call log_notice('History module is initialized.')

  end subroutine history_init

  subroutine history_final()

    if (allocated(u)) deallocate(u)
    if (allocated(v)) deallocate(v)
    if (allocated(gh)) deallocate(gh)

    call log_notice('History module is finalized.')

  end subroutine history_final

  subroutine history_write_state(state, static, diag)

    type(state_type), intent(in) :: state
    type(static_type), intent(in) :: static
    type(diag_type), intent(in) :: diag

    integer i, j

    ! Convert wind from C grid to A grid.
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        u(i,j) = 0.5 * (state%u(i,j) + state%u(i-1,j))
        v(i,j) = 0.5 * (state%v(i,j) + state%v(i,j-1))
        gh(i,j) = state%gd(i,j) + static%ghs(i,j)
      end do
    end do

    call io_start_output()
    call io_output('lon',   mesh%full_lon_deg(:))
    call io_output('lat',   mesh%full_lat_deg(:))
    call io_output('ilon',  mesh%half_lon_deg(:))
    call io_output('ilat',  mesh%half_lat_deg(:))
    call io_output('u',     u(1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('v',     v(1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('gh',    gh(1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('ghs',   static%ghs(1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('vor',   diag%vor(1:mesh%num_half_lon,1:mesh%num_half_lat))
    call io_output('div',   diag%div(1:mesh%num_full_lon,1:mesh%num_full_lat))
    call io_output('pv',    diag%pv(1:mesh%num_half_lon,1:mesh%num_half_lat))
    call io_output('te',    diag%total_energy)
    call io_output('tm',    diag%total_mass)
    call io_output('tes',   diag%total_potential_enstrophy)
    call io_output('tav',   diag%total_absolute_vorticity)
    call io_end_output()

  end subroutine history_write_state

  subroutine history_write_tendency(state, tend, tag)
    type(state_type), intent(in):: state 
    type(tend_type), intent(in) :: tend
    integer, intent(in) :: tag
 
    call io_start_output('debug', trim(to_string(tag)))
    call io_output('lon',           mesh%full_lon_deg(:),   'debug')
    call io_output('lat',           mesh%full_lat_deg(:),   'debug')
    call io_output('ilon',          mesh%half_lon_deg(:),   'debug')
    call io_output('ilat',          mesh%half_lat_deg(:),   'debug')

    call io_output('u_nonlinear',     tend%u_nonlinear(1:mesh%num_half_lon,1:mesh%num_full_lat),    'debug')
    call io_output('v_nonlinear',     tend%v_nonlinear(1:mesh%num_full_lon,1:mesh%num_half_lat),    'debug')
    call io_output('dEdlon',         tend%dEdlon(1:mesh%num_half_lon,1:mesh%num_full_lat),        'debug')
    call io_output('dEdlat',         tend%dEdlat(1:mesh%num_full_lon,1:mesh%num_half_lat),        'debug')
    call io_output('mass_div',      tend%mass_div(1:mesh%num_full_lon,1:mesh%num_full_lat), 'debug')
    call io_output('du',            tend%du(1:mesh%num_half_lon,1:mesh%num_full_lat),           'debug')
    call io_output('dv',            tend%dv(1:mesh%num_full_lon,1:mesh%num_half_lat),           'debug')
    call io_output('dgd',           tend%dgd(1:mesh%num_full_lon,1:mesh%num_full_lat),          'debug')
    call io_output('normal_lon_flux', tend%diag%normal_lon_flux(1:mesh%num_half_lon,1:mesh%num_full_lat), 'debug')
    call io_output('normal_lat_flux', tend%diag%normal_lat_flux(1:mesh%num_full_lon,1:mesh%num_half_lat), 'debug')
    ! call io_output('center_energy', tend%diag%energy(1:mesh%num_full_lon,1:mesh%num_full_lat),  'debug')
    call io_output('kinetic_energy', tend%diag%kinetic(1:mesh%num_full_lon,1:mesh%num_full_lat),  'debug')
    call io_output('PV',            tend%diag%pot_vor(1:mesh%num_half_lon,1:mesh%num_half_lat),  'debug')
    call io_output('PV_on_lon',     tend%diag%pv_lon(1:mesh%num_half_lon,1:mesh%num_full_lat), 'debug')
    call io_output('PV_on_lat',     tend%diag%pv_lat(1:mesh%num_full_lon,1:mesh%num_half_lat), 'debug')
    call io_output('gd_corner',     tend%diag%gd_corner(1:mesh%num_half_lon,1:mesh%num_half_lat),  'debug')
    call io_output('mass_flux_lon_t', tend%diag%mass_flux_lon_t(1:mesh%num_full_lon,1:mesh%num_half_lat), 'debug')
    call io_output('mass_flux_lat_t', tend%diag%mass_flux_lat_t(1:mesh%num_half_lon,1:mesh%num_full_lat), 'debug')
!!
!     call io_output('tangent_wind_lon', tend%diag%tangent_wind_lon(1:mesh%num_half_lon,1:mesh%num_full_lat), 'debug')
!     call io_output('tangent_wind_lat', tend%diag%tangent_wind_lat(1:mesh%num_full_lon,1:mesh%num_half_lat), 'debug')
!     call io_output('gradient_pv_lon', tend%diag%gradient_pv_lon(1:mesh%num_half_lon,1:mesh%num_full_lat), 'debug')
!     call io_output('gradient_pv_lat', tend%diag%gradient_pv_lat(1:mesh%num_full_lon,1:mesh%num_half_lat), 'debug')
!     call io_output('tangent_gradient_pv_lon', tend%diag%tangent_gradient_pv_lon(1:mesh%num_half_lon,1:mesh%num_full_lat), 'debug')
!     call io_output('tangent_gradient_pv_lat', tend%diag%tangent_gradient_pv_lat(1:mesh%num_full_lon,1:mesh%num_half_lat), 'debug')
    call io_end_output('debug')

  end subroutine history_write_tendency

end module history_mod
