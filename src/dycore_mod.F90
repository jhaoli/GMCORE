module dycore_mod

  use params_mod, split_scheme_in => split_scheme
  use data_mod
  use log_mod
  use types_mod
  use mesh_mod
  use time_mod, old => old_time_idx, new => new_time_idx
  use parallel_mod
  use io_mod
  use diag_mod
  use history_mod
  use restart_mod

  implicit none

  private

  public dycore_init
  public dycore_restart
  public dycore_run
  public dycore_final
  
  ! 1: csp1: first order conservative split
  ! 2: csp2: second order conservative split
  ! 3: isp: inproved second order conservative split
  integer split_scheme
  integer, parameter :: csp1 = 1
  integer, parameter :: csp2 = 2
  integer, parameter :: isp  = 3
  integer, parameter :: all_pass = 0
  integer, parameter :: fast_pass = 1
  integer, parameter :: slow_pass = 2

  integer, parameter :: center_diff = 0
  integer, parameter :: upwind = 1
  integer, parameter :: weno = 2

  interface
    subroutine integrator_interface(time_step_size, old, new, pass, qcon_modified_)
      real, intent(in) :: time_step_size
      integer, intent(in) :: old
      integer, intent(in) :: new
      integer, intent(in) :: pass 
      logical, intent(in), optional :: qcon_modified_
    end subroutine
  end interface

  procedure(integrator_interface), pointer :: integrator

  ! For debug
  integer :: tag = 0

contains

  subroutine dycore_init()

    if (case_name == '') then
      call log_error('case_name is not set!')
    end if

    call log_init()
    call mesh_init()
    call time_init()
    call parallel_init()
    call io_init()
    call diag_init()
    call history_init()
    call restart_init()
    call data_init()

    select case (time_scheme)
    case ('predict_correct')
      integrator => predict_correct
    case default
      call log_error('Unknown time_scheme ' // trim(time_scheme) // '!')
    end select

    select case (split_scheme_in)
    case ('csp1')
      split_scheme = csp1
    case ('csp2')
      split_scheme = csp2
    case ('isp')
      split_scheme = isp
    case default
      split_scheme = 0
      call log_notice('No fast-slow split.')
    end select

    call log_notice('Dycore module is initialized.')

  end subroutine dycore_init

  subroutine dycore_restart()

    call restart_read(state(old), static)

  end subroutine dycore_restart

  subroutine dycore_run()

    call reset_cos_lat_at_poles()
    call diag_run(state(old))
    call output(state(old))
    call log_add_diag('total_mass', diag%total_mass)
    call log_add_diag('total_energy', diag%total_energy)
    call log_step()


    do while (.not. time_is_finished())
      tag = 1
      call time_integrate()
      call time_advance()
      call diag_run(state(old))
      call output(state(old))
      call log_add_diag('total_mass', diag%total_mass)
      call log_add_diag('total_energy', diag%total_energy)
      call log_step()
    end do

  end subroutine dycore_run

  subroutine dycore_final()

    call mesh_final()
    call parallel_final()
    call diag_final()
    call history_final()
    call data_final()

    call log_notice('Dycore module is finalized.')

  end subroutine dycore_final

  subroutine reset_cos_lat_at_poles()

    integer j

    j = parallel%full_lat_south_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_south_pole_idx) * 0.25
    coef%full_dlon(j) = radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = radius * mesh%dlat * mesh%full_cos_lat(j)

    j = parallel%full_lat_north_pole_idx
    mesh%full_cos_lat(j) = mesh%half_cos_lat(parallel%half_lat_north_pole_idx) * 0.25
    coef%full_dlon(j) = radius * mesh%dlon * mesh%full_cos_lat(j)
    coef%full_dlat(j) = radius * mesh%dlat * mesh%full_cos_lat(j)

  end subroutine reset_cos_lat_at_poles

  subroutine output(state)

    type(state_type), intent(in) :: state 

    if (time_is_alerted('hist0.output')) call history_write(state, static, diag)
!     call history_write(state, static, diag)
 
!     if (time_is_alerted('restart.output')) call restart_write(state, static)
!     if (time_is_alerted('debug.output')) call history_write(state, tend(old), tag)
  end subroutine output

  subroutine space_operators(state, tend, pass)
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pass 
    integer i, j

    call tend_diag_operator(state, tend)  
    select case (pass)
    case (all_pass)
      call nonlinear_coriolis_operator(state, tend)
      call energy_gradient_operator(state, tend)
      call mass_divergence_operator(state, tend)

      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%u_pgf(i,j) + tend%u_nonlinear(i,j) 
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = tend%v_pgf(i,j) + tend%v_nonlinear(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx  
          tend%dgd(i,j) = tend%mass_div(i,j) 
        end do
      end do
    case (slow_pass)
#ifndef NDEBUG
      tend%u_pgf = 0.0
      tend%v_pgf = 0.0
      tend%mass_div = 0.0
#endif 
      call nonlinear_coriolis_operator(state, tend)

      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%u_nonlinear(i,j) 
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = tend%v_nonlinear(i,j)
        end do
      end do
      tend%dgd = 0.0
    case (fast_pass)
#ifndef NDEBUG
      tend%u_nonlinear = 0.0
      tend%v_nonlinear = 0.0
#endif
      call energy_gradient_operator(state, tend)
      call mass_divergence_operator(state, tend)

      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%u_pgf(i,j) 
        end do
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dv(i,j) = tend%v_pgf(i,j)
        end do
      end do

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = tend%mass_div(i,j)
        end do
      end do
    end select

!     call check_spaceoperator(tend, state)
!     stop 'check!'
  end subroutine space_operators

  subroutine tend_diag_operator(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real um1, up1, vm1, vp1, r1, r2, r3, r4, hd_corner, area_pole, sp, np
    integer i, j
  
    !! normal thickness flux
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        if (j == parallel%full_lat_start_idx .or. j == parallel%full_lat_end_idx) then
          tend%diag%normal_u_flux(i,j) = 0.0
        else 
          tend%diag%normal_u_flux(i,j) = 0.5 * (state%gd(i,j) + state%gd(i+1,j)) / g * state%u(i,j) 
        end if 
      end do
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%normal_v_flux(i,j) = 0.5 * (state%gd(i,j) + state%gd(i,j+1)) / g * state%v(i,j)
      end do 
    end do 

    !! vertex potential vorticity 
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      r1 = radius**2 * mesh%dlon * 0.5 * (mesh%full_sin_lat(j+1) - mesh%half_sin_lat(j)) / mesh%full_area(j+1)
      r2 = radius**2 * mesh%dlon * 0.5 * (mesh%half_sin_lat(j) - mesh%full_sin_lat(j)) / mesh%full_area(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        um1 = state%u(i,j) * mesh%full_cos_lat(j)
        up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
        vm1 = state%v(i,j)
        vp1 = state%v(i+1,j)
        tend%diag%hd_corner(i,j) = 1 / mesh%half_area(j) *(r1 * mesh%full_area(j+1) * (state%gd(i,j+1) + state%gd(i+1,j+1))+&
                                                           r2 * mesh%full_area(j  ) * (state%gd(i,j  ) + state%gd(i+1,j  ))) / g 
!         tend%diag%pot_vor(i,j) = ((vp1 - vm1) / coef%half_dlon(j) - (up1 - um1) / coef%half_dlat(j) + coef%half_f(j)) / tend%diag%hd_corner(i,j)
       
        tend%diag%pot_vor(i,j) = (1 / mesh%half_area(j) * (state%u(i,j) * radius * mesh%full_cos_lat(j) * mesh%dlon +&
                                                          state%v(i+1,j) * radius * mesh%dlat -&
                                                          state%u(i,j+1) * radius * mesh%full_cos_lat(j+1) * mesh%dlon -&
                                                          state%v(i,j) * radius * mesh%dlat) + coef%half_f(j)) / tend%diag%hd_corner(i,j)
        
      end do
    end do
  
    if (parallel%has_south_pole) then
      j = parallel%half_lat_south_pole_idx
      r1 = 0.5 
      r2 = radius**2 * mesh%dlon * 0.5 * (mesh%full_sin_lat(j+1) - mesh%half_sin_lat(j)) / mesh%full_area(j+1)

      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        um1 = state%u(i,j) * mesh%full_cos_lat(j)
        up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
        vm1 = state%v(i,j)
        vp1 = state%v(i+1,j)
        tend%diag%hd_corner(i,j) = 1 / mesh%half_area(j) * (r1 * mesh%full_area(j  ) * (state%gd(i,j  ) + state%gd(i+1,j  )) +&
                                                            r2 * mesh%full_area(j+1) * (state%gd(i,j+1) + state%gd(i+1,j+1))) / g 
!         tend%diag%pot_vor(i,j) = ((vp1 - vm1) / coef%half_dlon(j) - (up1 - um1) / coef%half_dlat(j) + coef%half_f(j)) / tend%diag%hd_corner(i,j)
        
        tend%diag%pot_vor(i,j) = (1 / mesh%half_area(j) * (state%v(i+1,j) * radius * mesh%dlat -&
                                                          state%u(i,j+1) * radius * mesh%full_cos_lat(j+1) * mesh%dlon -&
                                                          state%v(i,j) * radius * mesh%dlat) + coef%half_f(j))/ tend%diag%hd_corner(i,j) 
      end do                                       
    end if
     
    if (parallel%has_north_pole) then
      j = parallel%half_lat_north_pole_idx 
      r1 = 0.5 
      r2 = radius**2 * mesh%dlon * 0.5 * (mesh%half_sin_lat(j) - mesh%full_sin_lat(j)) / mesh%full_area(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        um1 = state%u(i,j) * mesh%full_cos_lat(j)
        up1 = state%u(i,j+1) * mesh%full_cos_lat(j+1)
        vm1 = state%v(i,j)
        vp1 = state%v(i+1,j)
        tend%diag%hd_corner(i,j) = 1 / mesh%half_area(j) * (r1 * mesh%full_area(j+1) * (state%gd(i,j+1) + state%gd(i+1,j+1)) +&
                                                            r2 * mesh%full_area(j) * (state%gd(i,j) + state%gd(i+1,j))) / g 
!         tend%diag%pot_vor(i,j) = ((vp1 - vm1) / coef%half_dlon(j) - (up1 - um1) / coef%half_dlat(j) + coef%half_f(j)) / tend%diag%hd_corner(i,j)
        tend%diag%pot_vor(i,j) = (1 / mesh%half_area(j) * (state%u(i,j) * radius * mesh%full_cos_lat(j) * mesh%dlon +&
                                                          state%v(i+1,j) * radius * mesh%dlat -&
                                                          state%v(i,j) * radius * mesh%dlat) + coef%half_f(j)) / tend%diag%hd_corner(i,j)
      end do                                       
    end if
     
    !! kinetic energy on primal cell center 
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
!       r1 = 0.25 * radius**2 * (mesh%half_sin_lat(j) - mesh%half_sin_lat(j-1)) * mesh%dlon 
!       r2 = r1
!       r3 = 0.5 * radius**2 * (mesh%half_sin_lat(j) - mesh%full_sin_lat(j  )) * mesh%dlon 
!       r4 = 0.5 * radius**2 * (mesh%full_sin_lat(j) - mesh%half_sin_lat(j-1)) * mesh%
      r1 = radius**2 * 0.5 * mesh%dlon * (mesh%half_sin_lat(j) - mesh%full_sin_lat(j)) * mesh%full_cos_lat(j) / (mesh%full_cos_lat(j) + mesh%half_cos_lat(j)) +&
           radius**2 * 0.5 * mesh%dlon * (mesh%full_sin_lat(j) - mesh%half_sin_lat(j-1)) * mesh%full_cos_lat(j) / (mesh%full_cos_lat(j) + mesh%half_cos_lat(j-1))
      r2 = r1
      r3 = radius**2 * mesh%dlon * (mesh%half_sin_lat(j) - mesh%full_sin_lat(j)) * mesh%half_cos_lat(j) / (mesh%half_cos_lat(j) + mesh%full_cos_lat(j))
      r4 = radius**2 * mesh%dlon * (mesh%full_sin_lat(j) - mesh%half_sin_lat(j-1)) * mesh%half_cos_lat(j-1) / (mesh%full_cos_lat(j) + mesh%half_cos_lat(j-1)) 
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        um1 = state%u(i-1,j) 
        up1 = state%u(i,j)
        vm1 = state%v(i,j-1)
        vp1 = state%v(i,j)
        tend%diag%kinetic_energy(i,j) = 1 / mesh%full_area(j) * (r1 * state%u(i,j)**2 + r2 * state%u(i-1,j)**2 +&
                                                                 r3 * state%v(i,j)**2 + r4 * state%v(i,j-1)**2) 
      end do 
    end do

    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + state%v(i,j)**2 
      end do
      call parallel_zonal_sum(sp)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%kinetic_energy(i,j) = sp / mesh%num_full_lon 
      end do
    end if

    if (parallel%has_north_pole) then
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np + state%v(i,j-1)**2 
      end do
      call parallel_zonal_sum(np)
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%kinetic_energy(i,j) = np / mesh%num_full_lon
      end do
    end if

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%energy(i,j) = tend%diag%kinetic_energy(i,j) + state%gd(i,j) + static%ghs(i,j)
      end do 
    end do 

    !! tangential velocity on primal cell eges
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        r1 = mesh%half_cos_lat(j) / mesh%full_cos_lat(j)
        r2 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%tangential_u_flux(i,j) = 0.125 *(r1 * ((state%gd(i,j) + state%gd(i,j+1)) * state%v(i,j  ) + (state%gd(i+1,j) + state%gd(i+1,j+1)) * state%v(i+1,j  )) + &
                                                   r2 * ((state%gd(i,j) + state%gd(i,j-1)) * state%v(i,j-1) + (state%gd(i+1,j) + state%gd(i+1,j-1)) * state%v(i+1,j-1))) / g  
      end do 
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%tangential_v_flux(i,j) = 0.125 * ((state%gd(i-1,j) + state%gd(i,j  )) * state%u(i-1,j) + (state%gd(i-1,j+1) + state%gd(i,j+1  )) * state%u(i-1,j+1) +&
                                                   (state%gd(i,j  ) + state%gd(i+1,j)) * state%u(i,j  ) + (state%gd(i,j+1  ) + state%gd(i+1,j+1)) * state%u(i,j+1 )) / g  
      end do 
    end do 

    call parallel_fill_halo(tend%diag%normal_u_flux, all_halo=.true.)
    call parallel_fill_halo(tend%diag%normal_v_flux, all_halo=.true.)
    call parallel_fill_halo(tend%diag%pot_vor, all_halo=.true.)
    call parallel_fill_halo(tend%diag%energy, all_halo=.true.)
    call parallel_fill_halo(tend%diag%tangential_u_flux, all_halo=.true.)
    call parallel_fill_halo(tend%diag%tangential_v_flux, all_halo=.true.)

  end subroutine tend_diag_operator

  subroutine nonlinear_coriolis_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j
    real, allocatable :: q_u(:,:), q_v(:,:)
    real :: r1, r2
    real :: dt, dx, dy, hh0, alpha
    real, parameter :: lat0= 90
    real :: full_beta(parallel%full_lat_start_idx_no_pole: parallel%full_lat_end_idx_no_pole)
    real :: half_beta(parallel%half_lat_start_idx: parallel%half_lat_end_idx)

    if (.not. allocated(q_u))    call parallel_allocate(q_u, half_lon=.true.)
    if (.not. allocated(q_v))    call parallel_allocate(q_v, half_lat=.true.)

 
!！（1） average nethod 
!     do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
!       do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
!         q_v(i,j) = 0.5 * (tend%diag%pot_vor(i-1,j) + tend%diag%pot_vor(i,j))
!       end do 
!     end do 
    
!     do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
!       do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
!         q_u(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i,j-1))
!       end do 
!     end do 
!!(2) replaced by upwind  
!     do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
!       do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
!         if (mesh%half_lat_deg(j) < -1.0 * lat0 .or. mesh%half_lat_deg(j) > lat0) then
!           if (tend%diag%tangential_v_flux(i,j) > 0) then 
!             q_v(i,j) = tend%diag%pot_vor(i-1,j)
!           elseif (tend%diag%tangential_v_flux(i,j) < 0) then
!             q_v(i,j) = tend%diag%pot_vor(i,j)
!           else
!             q_v(i,j) = 0.5 * (tend%diag%pot_vor(i-1,j) + tend%diag%pot_vor(i,j))
!           endif 
!         else 
!           q_v(i,j) = 0.5 * (tend%diag%pot_vor(i-1,j) + tend%diag%pot_vor(i,j))
!         end if 
!       end do 
!     end do 
!     do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
!       do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
!         if (mesh%full_lat_deg(j) < -1.0 * lat0 .or. mesh%full_lat_deg(j) > lat0) then
!           if (tend%diag%tangential_u_flux(i,j) > 0) then
!             q_u(i,j) = tend%diag%pot_vor(i,j-1)
!           elseif (tend%diag%tangential_u_flux(i,j) < 0) then
!             q_u(i,j) = tend%diag%pot_vor(i,j)
!           else
!             q_u(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i,j-1))
!           end if
!         else
!           q_u(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i,j-1))
!         end if 
!       end do 
!     end do 
!! (3) upwind intepolation with parameter-beta
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        full_beta(j) = 4 / pi**2 * mesh%full_lat(j)**2
      end do
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        half_beta(j) = 4 / pi**2 * mesh%half_lat(j)**2
      end do

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          q_v(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i-1,j)) - &
                    half_beta(j) * 0.5 * sign(1.0, tend%diag%tangential_v_flux(i,j)) * (tend%diag%pot_vor(i,j) - tend%diag%pot_vor(i-1,j)) 
        end do 
      end do 

      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole 
       do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          if (j == parallel%full_lat_start_idx .or. j == parallel%full_lat_end_idx) then
            q_u(i,j) = 0.0
          else
            q_u(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i,j-1)) - &
                   full_beta(j) * 0.5 * sign(1.0, tend%diag%tangential_u_flux(i,j)) * (tend%diag%pot_vor(i,j) - tend%diag%pot_vor(i,j-1))
          end if 
       end do 
      end do
    
    call parallel_fill_halo(q_u, all_halo = .true.)
    call parallel_fill_halo(q_v, all_halo = .true.)
    call parallel_fill_halo(tend%diag%normal_u_flux, all_halo=.true.)

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      r1 = mesh%half_cos_lat(j) / mesh%full_cos_lat(j)
      r2 = mesh%half_cos_lat(j-1) / mesh%full_cos_lat(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_nonlinear(i,j) = 0.125 * (r1 * (tend%diag%normal_v_flux(i  ,j  ) * (q_u(i,j) + q_v(i,  j))  +&
                                               tend%diag%normal_v_flux(i+1,j  ) * (q_u(i,j) + q_v(i+1,j))) +&
                                         r2 * (tend%diag%normal_v_flux(i,j-1  ) * (q_u(i,j) + q_v(i,j-1))  +&
                                               tend%diag%normal_v_flux(i+1,j-1) * (q_u(i,j) + q_v(i+1,j-1))) )
        
      end do
    end do
    
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_nonlinear(i,j) = -0.125 * (tend%diag%normal_u_flux(i-1,j+1) * (q_v(i,j) + q_u(i-1,j+1)) + &
                                         tend%diag%normal_u_flux(i,j+1  ) * (q_v(i,j) + q_u(i,  j+1)) + &
                                         tend%diag%normal_u_flux(i-1,j  ) * (q_v(i,j) + q_u(i-1,j  )) + &
                                         tend%diag%normal_u_flux(i,j    ) * (q_v(i,j) + q_u(i,j    )) )
      end do 
    end do  

    if (allocated(q_u))    deallocate(q_u)
    if (allocated(q_v))    deallocate(q_v)
  end subroutine nonlinear_coriolis_operator

  subroutine energy_gradient_operator(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_pgf(i,j) = -(tend%diag%energy(i+1,j) - tend%diag%energy(i,j)) / coef%full_dlon(j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf(i,j) = -(tend%diag%energy(i,j+1) - tend%diag%energy(i,j)) / coef%full_dlat(j) * mesh%full_cos_lat(j)
      end do
    end do
  end subroutine energy_gradient_operator

  subroutine mass_divergence_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i, j
    real sp, np

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div(i,j) = -1.0 / mesh%full_area(j) * (tend%diag%normal_v_flux(i,j  ) * radius * mesh%half_cos_lat(j  ) * mesh%dlon -&
                                                        tend%diag%normal_v_flux(i,j-1) * radius * mesh%half_cos_lat(j-1) * mesh%dlon +&
                                                        tend%diag%normal_u_flux(i,j  ) * radius * mesh%dlat - &
                                                        tend%diag%normal_u_flux(i-1,j) * radius * mesh%dlat) 
      end do 
    end do 
    
    if (parallel%has_south_pole) then
      j = parallel%full_lat_south_pole_idx
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + tend%diag%normal_v_flux(i,j) * radius * mesh%half_cos_lat(j) * mesh%dlon
      end do
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div(i,j) = -sp / (mesh%num_full_lon * mesh%full_area(j))
      end do 
    end if 
    if (parallel%has_north_pole) then 
      j = parallel%full_lat_north_pole_idx
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np - tend%diag%normal_v_flux(i,j-1) * radius * mesh%half_cos_lat(j-1) * mesh%dlon
      end do
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div(i,j) = -np / (mesh%num_full_lon * mesh%full_area(j))
      end do 
    end if 
    

  end subroutine mass_divergence_operator

  subroutine update_state(dt, tend, old_state, new_state)

    real, intent(in) :: dt
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: old_state
    type(state_type), intent(inout) :: new_state

    integer i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%gd(i,j) = (old_state%gd(i,j) / g  + dt * tend%dgd(i,j)) * g
      end do
    end do

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u(i,j) = old_state%u(i,j) + dt * tend%du(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        new_state%v(i,j) = old_state%v(i,j) + dt * tend%dv(i,j)
      end do
    end do

    call parallel_fill_halo(new_state%gd(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%u(:,:), all_halo=.true.)
    call parallel_fill_halo(new_state%v(:,:), all_halo=.true.)

  end subroutine update_state

  subroutine time_integrate()

    select case (split_scheme)
    case(csp2)
      call csp2_splitting()
    case default
      call integrator(time_step_size, old, new, all_pass)
    end select

  end subroutine time_integrate

  subroutine csp2_splitting()
    real fast_dt
    integer subcycle, t1, t2

    fast_dt = time_step_size / subcycles
    t1 = 0 
    t2 = old

    call integrator(0.5 * time_step_size, old, t1, slow_pass)
    do subcycle = 1, subcycles
      call integrator(fast_dt, t1, t2, fast_pass)
      call time_swap_indices(t1,t2)
    end do 
    call integrator(0.5 * time_step_size, t1, new, slow_pass)

  end subroutine csp2_splitting

  subroutine predict_correct(time_step_size, old, new, pass, qcon_modified_)

    real, intent(in) :: time_step_size
    integer, intent(in) :: old
    integer, intent(in) :: new
    integer, intent(in) :: pass 
    logical, intent(in), optional :: qcon_modified_

    real dt, ip1, ip2

    dt = time_step_size * 0.5

    ! Do first predict step.
    call space_operators(state(old), tend(old), pass)
    call update_state(dt, tend(old), state(old), state(new))
     
    ! Do second predict step.
    call space_operators(state(new), tend(old), pass)
    call update_state(dt, tend(old), state(old), state(new))

    ! Do correct stepe
    call space_operators(state(new), tend(new), pass)
    dt = time_step_size 
    call update_state(dt, tend(new), state(old), state(new))

  end subroutine predict_correct

  subroutine check_spaceoperator(tend, state)
    type(tend_type), intent(in) :: tend
    type(state_type), intent(in) :: state
    integer :: i, j
    real hh0, hh1

    real ip_egf
    real ip_fv 
    real ip_fu
    real ip_energy_div
    real ip_mass_div

    ip_egf = 0.0
    ip_fv = 0.0
    ip_fu = 0.0
    ip_energy_div = 0.0
    ip_mass_div = 0.0

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
!         hh0 = (state%gd(i,j) + state%gd(i+1,j)) * 0.5 / g
!         ip_egf = ip_egf + tend%u_pgf(i,j) * hh0 * state%u(i,j) * mesh%full_cos_lat(j)
        ip_egf = ip_egf + tend%u_pgf(i,j) * tend%diag%normal_u_flux(i,j) * mesh%full_cos_lat(j)
        ip_fv = ip_fv + tend%u_nonlinear(i,j) * tend%diag%normal_u_flux(i,j) * mesh%full_cos_lat(j)
      end do 
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
!         hh1 = (state%gd(i,j) + state%gd(i,j+1)) * 0.5 / g
!         ip_egf = ip_egf + tend%v_pgf(i,j) * hh1 * state%v(i,j) * mesh%half_cos_lat(j)
        ip_egf = ip_egf + tend%v_pgf(i,j) * tend%diag%normal_v_flux(i,j) * mesh%half_cos_lat(j)
        ip_fu = ip_fu + tend%v_nonlinear(i,j) * tend%diag%normal_v_flux(i,j) * mesh%half_cos_lat(j)
      end do 
    end do 

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        ip_energy_div = ip_energy_div + tend%diag%energy(i,j) * tend%mass_div(i,j) * mesh%full_cos_lat(j)
        ip_mass_div = ip_mass_div + tend%mass_div(i,j) * mesh%full_cos_lat(j) 
      end do
    end do
    print*, 'nonlinear Coriolis:    ' ,' total energy transfer tendency:', '      total mass tendency:'
    print*, ip_fv + ip_fu, ip_egf + ip_energy_div, ip_mass_div

  end subroutine check_spaceoperator  
end module dycore_mod
