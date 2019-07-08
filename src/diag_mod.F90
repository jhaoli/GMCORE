module diag_mod

  use ieee_arithmetic
  use parallel_mod
  use mesh_mod
  use data_mod
  use types_mod
  use log_mod
  use params_mod

  implicit none

  private

  public diag_init
  public diag_run
  public diag_final
  public diag_total_energy
  public diag_type
  public diag

  type diag_type
    real total_mass
    real total_energy
    real total_enstrophy
    real, allocatable :: vor(:,:)
    real, allocatable :: div(:,:)
    real, allocatable :: pv(:,:)
  end type diag_type

  type(diag_type) diag

contains

  subroutine diag_init()

    if (.not. allocated(diag%vor)) call parallel_allocate(diag%vor, half_lon=.true., half_lat=.true.)
    if (.not. allocated(diag%pv))  call parallel_allocate(diag%pv, half_lon=.true., half_lat=.true.)
    if (.not. allocated(diag%div)) call parallel_allocate(diag%div)

    call log_notice('Diag module is initialized.')

  end subroutine diag_init

  subroutine diag_run(state)

    type(state_type), intent(in) :: state

    real vm1, vp1, um1, up1
    integer i, j
    real sp, np, area_pole

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        diag%div(i,j) = ((state%u(i,j) - state%u(i-1,j)) * mesh%vertex_lat_distance(j) + &
                         (state%v(i,j+1) * mesh%vertex_lon_distance(j+1) -&
                          state%v(i,j) * mesh%vertex_lon_distance(j))) / mesh%cell_area(j)  
      end do
    end do
    call parallel_fill_halo(diag%div, all_halo=.true.)

    diag%total_mass = 0.0
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        diag%total_mass = diag%total_mass + state%gd(i,j) * mesh%cell_area(j) 
      end do
    end do

    if (ieee_is_nan(diag%total_mass)) then
      call log_error('Total mass is NaN!')
    end if

    diag%total_energy = diag_total_energy(state)

    if (ieee_is_nan(diag%total_energy)) then
      call log_error('Total energy is NaN!')
    end if

    if (ieee_is_nan(diag%total_enstrophy)) then
      call log_error('Total potential enstrophy is NaN!')
    end if

  end subroutine diag_run

  subroutine diag_final()

    if (allocated(diag%vor)) deallocate(diag%vor)
    if (allocated(diag%div)) deallocate(diag%div)
    if (allocated(diag%pv)) deallocate(diag%pv)
  end subroutine diag_final
    

  real function diag_total_energy(state) result(res)

    type(state_type), intent(in) :: state

    integer i, j

    res = 0.0
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        res = res + (mesh%lon_edge_left_area(j) * state%gd(i,j) + mesh%lon_edge_left_area(j) * state%gd(i+1,j)) /&
              (mesh%lon_edge_area(j) * g) * state%u(i,j)**2 * mesh%lon_edge_area(j) 
      end do
    end do
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + (mesh%lat_edge_up_area(j) * state%gd(i,j) + mesh%lat_edge_down_area(j) * state%gd(i,j-1)) /&
             (mesh%lat_edge_area(j) * g) * state%v(i,j)**2 * mesh%lat_edge_area(j) 
      end do
    end do
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        res = res + (state%gd(i,j)**2 / g * 0.5 + state%gd(i,j) * static%ghs(i,j) / g) * mesh%cell_area(j)
      end do
    end do

  end function diag_total_energy

end module diag_mod
