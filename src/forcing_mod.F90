module forcing_mod

  use ieee_arithmetic
  use parallel_mod
  use mesh_mod
  use data_mod
  use types_mod
  use log_mod
  use params_mod

implicit none

private

public force_run

contains

  subroutine force_run(state, tend)
    type(state_type), intent(in)    :: state
    type(tend_type),  intent(inout) :: tend
    real, parameter    :: u0 = 120.0
    integer, parameter :: m = 12
    real, parameter    :: dh = 8000.0 / g
    real, parameter    :: tau_u = 100*86400.0, tau_h = 100*86400.0   !s 

    real :: ueqm, heqm, h0
    integer :: i,j

    select case (test_case)
    case ('held_suarez')
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        ueqm = u0 * sin(2*mesh%full_lat(j))**2 * (sin(m*mesh%full_lat(j))**2 - 0.5)
      	do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%force%u(i,j) = (ueqm - state%u(i,j)) / tau_u
      	end do
      end do
      call parallel_fill_halo(tend%force%u, all_halo=.true.)
  
      tend%force%v(:,:) = 0.0
      call parallel_fill_halo(tend%force%v, all_halo=.true.)
  
      h0 = 10**3 - 2 * dh / 3.
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        heqm = h0 + dh * (1 - mesh%full_sin_lat(j)**2)
      	do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%force%gd(i,j) = (heqm * g - state%gd(i,j)) / tau_h
      	end do
      end do
      call parallel_fill_halo(tend%force%gd, all_halo=.true.)
    case default
      tend%force%u = 0.0
      tend%force%v = 0.0
      tend%force%gd = 0.0
    end select

  end subroutine force_run

end module forcing_mod
