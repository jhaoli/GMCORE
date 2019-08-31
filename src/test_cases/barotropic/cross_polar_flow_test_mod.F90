module cross_polar_flow_test_mod

  use mesh_mod
  use parallel_mod
  use params_mod
  use data_mod

  implicit none

  private

  public cross_polar_flow_test_set_initial_condition

  real, parameter :: v0 = 20 !m s-1
  real, parameter :: gd0 = 5.7684e4 ! m2 s-2

contains

  subroutine cross_polar_flow_test_set_initial_condition()

    real cos_lat, sin_lat, cos_lon, sin_lon
    integer i, j

    write(6, *) '[Notice]: Use cross polar flow initial condition.'

    static%ghs(:,:) = 0.0

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        sin_lon = mesh%half_sin_lon(i)
        state(1)%u(i,j) = -v0 * sin_lon * sin_lat * (4.0 * cos_lat**2 - 1.0)
      end do
    end do

    call parallel_fill_halo(state(1)%u, all_halo=.true.)

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      sin_lat = mesh%half_sin_lat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        cos_lon = mesh%full_cos_lon(i)
        state(1)%v(i,j) = v0 * sin_lat**2 * cos_lon
      end do
    end do

    call parallel_fill_halo(state(1)%v, all_halo=.true.)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      cos_lat = mesh%full_cos_lat(j)
      sin_lat = mesh%full_sin_lat(j)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sin_lon = mesh%full_sin_lon(i)
        state(1)%gd(i,j) = gd0 + 2 * radius * omega * v0 * sin_lat**3 * cos_lat * sin_lon
      end do
    end do

    call parallel_fill_halo(state(1)%gd, all_halo=.true.)

  end subroutine cross_polar_flow_test_set_initial_condition

end module cross_polar_flow_test_mod
