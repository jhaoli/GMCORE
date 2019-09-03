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
  use forcing_mod

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

    call diag_run(state(old))
    call calc_total_vorticity_enstrophy(state(old),tend(old))
    call history_write(state(old), static, diag)
    call log_add_diag('total_mass', diag%total_mass)
    call log_add_diag('total_energy', diag%total_energy)
    call log_add_diag('total_potential_enstrophy', diag%total_enstrophy)
    call log_add_diag('total_absolute_vorticity', diag%total_absolute_vorticity)
    call log_step()

    do while (.not. time_is_finished())
      tag = 1
      call time_integrate()
      call time_advance()
      call diag_run(state(old))
      call calc_total_vorticity_enstrophy(state(old),tend(old))
      call output(state(old))
      call log_add_diag('total_mass', diag%total_mass)
      call log_add_diag('total_energy', diag%total_energy)
      call log_add_diag('total_potential_enstrophy', diag%total_enstrophy)
      call log_add_diag('total_absolute_vorticity', diag%total_absolute_vorticity)
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

  subroutine output(state)

    type(state_type), intent(in) :: state 

    if (time_is_alerted('hist0.output')) call history_write(state, static, diag)
 
!     if (time_is_alerted('restart.output')) call restart_write(state, static)
!     if (time_is_alerted('debug.output')) call history_write(state, tend(old), tag)
  end subroutine output

  subroutine space_operators(state, tend, pass)
    type(state_type), intent(inout) :: state
    type(tend_type), intent(inout) :: tend
    integer, intent(in) :: pass 
    integer i, j

    call calc_normal_mass_flux(state, tend)  
    call calc_pv_on_vertex(state, tend)
    call calc_energy_on_center(state, tend)
    call calc_tangent_mass_flux(state, tend)

    select case (pass)
    case (all_pass)
      
      call nonlinear_coriolis_operator(state, tend)
      call energy_gradient_operator(state, tend)
      call mass_divergence_operator(state, tend)
      call force_run(state, tend)

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%u_pgf(i,j) + tend%u_nonlinear(i,j) + tend%force%u(i,j)
        end do
      end do
      if (test_case == "held_suarez") then
        tend%dv(:,:) = 0.0  
      else
        do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
          do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
            tend%dv(i,j) = tend%v_pgf(i,j) + tend%v_nonlinear(i,j) + tend%force%v(i,j)
          end do
        end do  
      end if 

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%dgd(i,j) = tend%mass_div(i,j) + tend%force%gd(i,j)
        end do
      end do
    case (slow_pass)
#ifndef NDEBUG
      tend%u_pgf = 0.0
      tend%v_pgf = 0.0
      tend%mass_div = 0.0
#endif
      
      call nonlinear_coriolis_operator(state, tend)
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%u_nonlinear(i,j) 
        end do
      end do
      do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
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
      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%du(i,j) = tend%u_pgf(i,j) 
        end do
      end do

      do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
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

  subroutine calc_normal_mass_flux(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j
    call calc_gd_on_edge(state, tend)

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%normal_lon_flux(i,j) = tend%diag%gd_lon(i,j) / g * state%u(i,j) 
      end do
    end do 
    
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        if (j == parallel%half_lat_start_idx .or. j == parallel%half_lat_end_idx) then
          tend%diag%normal_lat_flux(i,j) = 0.0
        else
          tend%diag%normal_lat_flux(i,j) = tend%diag%gd_lat(i,j) / g * state%v(i,j)
        end if 
      end do 
    end do 

    call parallel_fill_halo(tend%diag%normal_lon_flux, all_halo=.true.)
    call parallel_fill_halo(tend%diag%normal_lat_flux, all_halo=.true.)

  end subroutine calc_normal_mass_flux

  subroutine calc_gd_on_edge(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%gd_lon(i,j) = (mesh%lon_edge_left_area(j) * state%gd(i,j  ) +&
                                 mesh%lon_edge_right_area(j) * state%gd(i+1,j)) / mesh%lon_edge_area(j) 
      end do 
    end do 
    
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%gd_lat(i,j) = (mesh%lat_edge_up_area(j) * state%gd(i,j) +&
                                  mesh%lat_edge_down_area(j) * state%gd(i,j-1)) / mesh%lat_edge_area(j)  
      end do 
    end do 
  end subroutine

  subroutine calc_pv_on_vertex(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real :: r1, r2, area_pole, sp, np 
    integer :: i, j
    
    diag%total_enstrophy = 0.0
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%gd_corner(i,j) = 1.0 / mesh%vertex_area(j) * (mesh%subcell_area(1,j  ) * (state%gd(i,j  ) + state%gd(i+1,j  )) +&
                                                                mesh%subcell_area(2,j-1) * (state%gd(i,j-1) + state%gd(i+1,j-1))) 
        diag%vor(i,j) = 1.0 / mesh%vertex_area(j) * (state%u(i,j-1) * mesh%cell_lon_distance(j-1) -& 
                                                     state%u(i,j  ) * mesh%cell_lon_distance(j) +& 
                                                     state%v(i+1,j) * mesh%cell_lat_distance(j) -& 
                                                     state%v(i,j  ) * mesh%cell_lat_distance(j))  
        tend%diag%pot_vor(i,j) = (diag%vor(i,j) + coef%half_f(j)) / (tend%diag%gd_corner(i,j) / g )
        diag%pv(i,j) = tend%diag%pot_vor(i,j)
      end do
    end do

    if (parallel%has_south_pole) then
      j = parallel%half_lat_south_pole_idx
      r1 = 1.0 / mesh%num_full_lon
      sp = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        sp = sp + r1 * state%gd(i,j)
      end do 
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%gd_corner(i,j) = sp 
      end do 
      sp = 0.0
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        sp = sp - state%u(i,j) * mesh%cell_lon_distance(j) 
      end do 
      area_pole = mesh%num_half_lon * mesh%vertex_area(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        diag%vor(i,j) = sp / area_pole
        tend%diag%pot_vor(i,j) = (diag%vor(i,j) + coef%half_f(j)) / (tend%diag%gd_corner(i,j) / g)
        diag%pv(i,j) = tend%diag%pot_vor(i,j)
      end do                              
    end if
    
    if (parallel%has_north_pole) then
      j = parallel%half_lat_north_pole_idx 
      r1 = 1.0 / mesh%num_full_lon
      np = 0.0
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        np = np + r1 * state%gd(i,j-1)
      end do 
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%gd_corner(i,j) = np 
      end do
      np = 0.0
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        np = np + state%u(i,j-1) * mesh%cell_lon_distance(j-1) 
      end do  
      area_pole = mesh%num_half_lon * mesh%vertex_area(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        diag%vor(i,j) = np / area_pole
        tend%diag%pot_vor(i,j) = (diag%vor(i,j) + coef%half_f(j)) / (tend%diag%gd_corner(i,j) / g)
        diag%pv(i,j) = tend%diag%pot_vor(i,j)
      end do
    end if
    call parallel_fill_halo(diag%vor, all_halo=.true.)
    call parallel_fill_halo(diag%pv, all_halo=.true.)
    call parallel_fill_halo(tend%diag%pot_vor, all_halo=.true.)
  end subroutine calc_pv_on_vertex

  subroutine calc_energy_on_center(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer :: i, j
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx 
        tend%diag%kinetic_energy(i,j) = 1.0 / mesh%cell_area(j) * (mesh%lon_edge_left_area(j) * state%u(i,j)**2 +&
                                                                   mesh%lon_edge_right_area(j) * state%u(i-1,j)**2 +&
                                                                   mesh%lat_edge_down_area(j+1) * state%v(i,j+1)**2 +&
                                                                   mesh%lat_edge_up_area(j) * state%v(i,j)**2)
      end do 
    end do 
     
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%energy(i,j) = tend%diag%kinetic_energy(i,j) + state%gd(i,j) + static%ghs(i,j)
      end do 
    end do
    
    call parallel_fill_halo(tend%diag%energy, all_halo=.true.)

  end subroutine calc_energy_on_center

  subroutine calc_tangent_mass_flux(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    real :: r1, r2
    integer :: i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%diag%mass_flux_lat_t(i,j) = 1.0 / mesh%cell_lon_distance(j) *&
                                         (mesh%vertex_lon_distance(j+1) * mesh%half_tangent_wgt(1,j+1) * tend%diag%normal_lat_flux(i,j+1  ) +&
                                          mesh%vertex_lon_distance(j+1) * mesh%half_tangent_wgt(1,j+1) * tend%diag%normal_lat_flux(i+1,j+1) +&
                                          mesh%vertex_lon_distance(j  ) * mesh%half_tangent_wgt(2,j  ) * tend%diag%normal_lat_flux(i,j    ) +&
                                          mesh%vertex_lon_distance(j  ) * mesh%half_tangent_wgt(2,j  ) * tend%diag%normal_lat_flux(i+1,j  ))  
      end do 
    end do 

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        if (j == parallel%half_lat_start_idx .or. j == parallel%half_lon_end_idx) then
          tend%diag%mass_flux_lon_t(i,j) = 0.0
        else 
          tend%diag%mass_flux_lon_t(i,j) = 1.0 / mesh%cell_lat_distance(j) *&
                                          (mesh%vertex_lat_distance(j  ) * mesh%full_tangent_wgt(1,j  ) * tend%diag%normal_lon_flux(i-1,j  ) +&
                                           mesh%vertex_lat_distance(j  ) * mesh%full_tangent_wgt(1,j  ) * tend%diag%normal_lon_flux(i,j    ) +&
                                           mesh%vertex_lat_distance(j-1) * mesh%full_tangent_wgt(2,j-1) * tend%diag%normal_lon_flux(i-1,j-1) +&
                                           mesh%vertex_lat_distance(j-1) * mesh%full_tangent_wgt(2,j-1) * tend%diag%normal_lon_flux(i,j-1  )) 
        end if
      end do 
    end do 

    call parallel_fill_halo(tend%diag%mass_flux_lat_t, all_halo=.true.)
    call parallel_fill_halo(tend%diag%mass_flux_lon_t, all_halo=.true.)
  end subroutine calc_tangent_mass_flux

  subroutine nonlinear_coriolis_operator(state, tend)

    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j
    real :: r1, r2

    select case(pv_scheme)
    case (1)  
      call calc_pv_on_edge_midpoint(state, tend)
    case (2) 
      call calc_pv_on_edge_upwind1D(state, tend)
    case (3)
      call calc_pv_on_edge_apvm(state, tend)
    case (4)
      call calc_pv_on_edge_upwind2D(state, tend)
    case default
      call log_error('Unknown PV scheme.')
    end select

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_nonlinear(i,j) = 0.5 / mesh%cell_lon_distance(j) *&
                                (mesh%vertex_lon_distance(j+1) * mesh%half_tangent_wgt(1,j+1) * tend%diag%normal_lat_flux(i,  j+1) * (tend%diag%pv_lon(i,j) + tend%diag%pv_lat(i, j+1 )) +&
                                 mesh%vertex_lon_distance(j+1) * mesh%half_tangent_wgt(1,j+1) * tend%diag%normal_lat_flux(i+1,j+1) * (tend%diag%pv_lon(i,j) + tend%diag%pv_lat(i+1,j+1)) +&
                                 mesh%vertex_lon_distance(j  ) * mesh%half_tangent_wgt(2,j  ) * tend%diag%normal_lat_flux(i,j    ) * (tend%diag%pv_lon(i,j) + tend%diag%pv_lat(i,j    )) +&
                                 mesh%vertex_lon_distance(j  ) * mesh%half_tangent_wgt(2,j  ) * tend%diag%normal_lat_flux(i+1,j  ) * (tend%diag%pv_lon(i,j) + tend%diag%pv_lat(i+1,j  )))
      
      end do
    end do      

    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_nonlinear(i,j) = -0.5 / mesh%cell_lat_distance(j) *&
                               (mesh%vertex_lat_distance(j  ) * mesh%full_tangent_wgt(1,j  ) * tend%diag%normal_lon_flux(i-1,j  ) * (tend%diag%pv_lat(i,j) + tend%diag%pv_lon(i-1,j  )) + &
                                mesh%vertex_lat_distance(j  ) * mesh%full_tangent_wgt(1,j  ) * tend%diag%normal_lon_flux(i,j    ) * (tend%diag%pv_lat(i,j) + tend%diag%pv_lon(i,  j  )) + &
                                mesh%vertex_lat_distance(j-1) * mesh%full_tangent_wgt(2,j-1) * tend%diag%normal_lon_flux(i-1,j-1) * (tend%diag%pv_lat(i,j) + tend%diag%pv_lon(i-1,j-1)) + &
                                mesh%vertex_lat_distance(j-1) * mesh%full_tangent_wgt(2,j-1) * tend%diag%normal_lon_flux(i,j-1  ) * (tend%diag%pv_lat(i,j) + tend%diag%pv_lon(i,j-1  )) )
      
      end do 
    end do  
 
  end subroutine nonlinear_coriolis_operator
  
  subroutine calc_pv_on_edge_midpoint(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j  

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%diag%pv_lat(i,j) = 0.5 * (tend%diag%pot_vor(i-1,j) + tend%diag%pot_vor(i,j))
        end do 
      end do 

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
        do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%diag%pv_lon(i,j) = 0.5 * (tend%diag%pot_vor(i,j+1) + tend%diag%pot_vor(i,j))
        end do 
      end do 
    call parallel_fill_halo(tend%diag%pv_lat, all_halo = .true.)
    call parallel_fill_halo(tend%diag%pv_lon, all_halo = .true.)
  end subroutine calc_pv_on_edge_midpoint

  subroutine calc_pv_on_edge_upwind1D(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j  
    real :: full_beta(parallel%full_lat_start_idx: parallel%full_lat_end_idx)
    real :: half_beta(parallel%half_lat_start_idx: parallel%half_lat_end_idx)

      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          tend%diag%pv_lat(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i-1,j)) - &
                    mesh%half_upwind_beta(j) * 0.5 * sign(1.0, tend%diag%mass_flux_lon_t(i,j)) * (tend%diag%pot_vor(i,j) - tend%diag%pot_vor(i-1,j)) 
        end do 
      end do 

      do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
       do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
          tend%diag%pv_lon(i,j) = 0.5 * (tend%diag%pot_vor(i,j+1) + tend%diag%pot_vor(i,j)) - &
                  mesh%full_upwind_beta(j) * 0.5 * sign(1.0, tend%diag%mass_flux_lat_t(i,j)) * (tend%diag%pot_vor(i,j+1) - tend%diag%pot_vor(i,j))
       end do 
      end do

    call parallel_fill_halo(tend%diag%pv_lat, all_halo = .true.)
    call parallel_fill_halo(tend%diag%pv_lon, all_halo = .true.)
  end subroutine calc_pv_on_edge_upwind1D

  subroutine calc_pv_on_edge_upwind2D(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type) , intent(inout) :: tend
    integer :: i, j

    call calc_dpv_on_edge(state, tend)

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%pv_lat(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i-1,j)) - &
                                mesh%half_upwind_beta(j) * 0.5 * (sign(1.0, tend%diag%mass_flux_lon_t(i,j)) * tend%diag%dpv_lon_t(i,j) +&
                                                                  sign(1.0, state%v(i,j)) * tend%diag%dpv_lat_n(i,j))
      end do 
    end do 
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%pv_lon(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i,j-1)) - &
                                mesh%full_upwind_beta(j) * 0.5 * (sign(1.0, state%u(i,j)) * tend%diag%dpv_lon_n(i,j) + &
                                                                  sign(1.0, tend%diag%mass_flux_lat_t(i,j)) * tend%diag%dpv_lat_t(i,j))
      end do 
    end do 
    call parallel_fill_halo(tend%diag%pv_lon, all_halo = .true.)
    call parallel_fill_halo(tend%diag%pv_lat, all_halo = .true.)
  end subroutine calc_pv_on_edge_upwind2D

  subroutine calc_dpv_on_edge(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j
    ! tangent pv difference
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%dpv_lon_t(i,j) = tend%diag%pot_vor(i,j) - tend%diag%pot_vor(i-1,j)
      end do 
    end do 

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%dpv_lat_t(i,j) = tend%diag%pot_vor(i,j+1) - tend%diag%pot_vor(i,j)
      end do 
    end do 
    call parallel_fill_halo(tend%diag%dpv_lon_t, all_halo = .true.)
    call parallel_fill_halo(tend%diag%dpv_lat_t, all_halo = .true.)

    ! normal pv difference
    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%diag%dpv_lat_n(i,j) = 0.25 * (tend%diag%dpv_lat_t(i-1,j  ) + tend%diag%dpv_lat_t(i,j  ) +&
                                           tend%diag%dpv_lat_t(i-1,j-1) + tend%diag%dpv_lat_t(i,j-1))

      end do
    end do 
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%diag%dpv_lon_n(i,j) = 0.25 * (tend%diag%dpv_lon_t(i,j+1) + tend%diag%dpv_lon_t(i+1,j+1) +&
                                           tend%diag%dpv_lon_t(i,j  ) + tend%diag%dpv_lon_t(i+1,j  ))
      end do  
    end do
    call parallel_fill_halo(tend%diag%dpv_lat_n, all_halo = .true.)
    call parallel_fill_halo(tend%diag%dpv_lon_n, all_halo = .true.)
  end subroutine calc_dpv_on_edge

  subroutine calc_pv_on_edge_apvm(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    real le, de, un, ut, vn, vt 
    integer :: i, j

    call calc_dpv_on_edge(state, tend)
    call calc_gd_on_edge(state, tend)

    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
        le = mesh%vertex_lon_distance(j)
        de = mesh%cell_lat_distance(j)
        do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
          if (j == parallel%half_lat_start_idx .or. j == parallel%half_lat_end_idx) then
            tend%diag%pv_lat(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i-1,j)) 
          else
            ut = tend%diag%mass_flux_lon_t(i,j) / tend%diag%gd_lat(i,j) / g 
            vn = state%v(i,j)
            tend%diag%pv_lat(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i-1,j)) - &
                                    mesh%half_upwind_beta(j) * 0.5 * (ut * tend%diag%dpv_lon_t(i,j) / le + vn * tend%diag%dpv_lat_n(i,j) / de) * time_step_size
          end if
        end do  
    end do 

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      le = mesh%vertex_lat_distance(j)
      de = mesh%cell_lon_distance(j)
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        un = state%u(i,j)
        vt = tend%diag%mass_flux_lat_t(i,j) / tend%diag%gd_lon(i,j) / g
        tend%diag%pv_lon(i,j) = 0.5 * (tend%diag%pot_vor(i,j) + tend%diag%pot_vor(i,j+1)) - &
                                mesh%full_upwind_beta(j) * 0.5 * (un * tend%diag%dpv_lon_n(i,j) / de + vt * tend%diag%dpv_lat_t(i,j) / le) * time_step_size
      end do 
    end do 
    call parallel_fill_halo(tend%diag%pv_lat, all_halo = .true.)
    call parallel_fill_halo(tend%diag%pv_lon, all_halo = .true.)
  end subroutine calc_pv_on_edge_apvm
  

  subroutine calc_total_vorticity_enstrophy(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer i, j 

    call calc_gd_on_edge(state, tend)
    call calc_pv_on_vertex(state, tend)
    
    diag%total_enstrophy = 0.0
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        diag%total_enstrophy = diag%total_enstrophy + 0.5 * (tend%diag%gd_corner(i,j) / g) * tend%diag%pot_vor(i,j)**2 * mesh%vertex_area(j)
      end do 
    end do 

    diag%total_absolute_vorticity = 0.0
    do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        diag%total_absolute_vorticity = diag%total_absolute_vorticity + (tend%diag%gd_corner(i,j) / g) * tend%diag%pot_vor(i,j) * mesh%vertex_area(j)
      end do 
    end do 
  end subroutine calc_total_vorticity_enstrophy

  subroutine energy_gradient_operator(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend

    integer i,j
    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tend%u_pgf(i,j) = -(tend%diag%energy(i+1,j) - tend%diag%energy(i,j)) / mesh%cell_lon_distance(j) 
      end do
    end do

    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%v_pgf(i,j) = -(tend%diag%energy(i,j) - tend%diag%energy(i,j-1)) / mesh%cell_lat_distance(j) 
      end do
    end do
  end subroutine energy_gradient_operator

  subroutine mass_divergence_operator(state, tend)
    type(state_type), intent(in) :: state
    type(tend_type), intent(inout) :: tend
    integer :: i, j

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tend%mass_div(i,j) = -1.0 / mesh%cell_area(j) * (tend%diag%normal_lon_flux(i,j  ) * mesh%vertex_lat_distance(j  ) -&
                                                         tend%diag%normal_lon_flux(i-1,j) * mesh%vertex_lat_distance(j  ) +&
                                                         tend%diag%normal_lat_flux(i,j+1) * mesh%vertex_lon_distance(j+1)-&
                                                         tend%diag%normal_lat_flux(i,j  ) * mesh%vertex_lon_distance(j  ))
      end do 
    end do 
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

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        new_state%u(i,j) = old_state%u(i,j) + dt * tend%du(i,j)
      end do
    end do

    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
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

    real dt, ip1, ip2, beta

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
    real hh0, hh1, r1, r3, r4

    real ip_egf
    real ip_fv 
    real ip_fu
    real ip_energy_div
    real ip_mass_div
    real ip_mass_div_dual
    real tmp 

    ip_egf = 0.0
    ip_fv = 0.0
    ip_fu = 0.0
    ip_energy_div = 0.0
    ip_mass_div = 0.0
  
    ip_mass_div_dual = 0.0 

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      do i = parallel%half_lon_start_idx, parallel%half_lon_end_idx
        tmp = ((state%gd(i+1,j) + static%ghs(i+1,j)) - (state%gd(i,j) + static%ghs(i,j))) / mesh%cell_lon_distance(j) 
        ip_egf = ip_egf + tmp * tend%diag%normal_lon_flux(i,j) * mesh%lon_edge_area(j) * 2 /radius**2 
        ip_fv = ip_fv + tend%u_nonlinear(i,j) * tend%diag%normal_lon_flux(i,j) * mesh%lon_edge_area(j) / radius**2 
      end do 
    end do 

    do j = parallel%half_lat_start_idx_no_pole, parallel%half_lat_end_idx_no_pole
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tmp = ((state%gd(i,j) + static%ghs(i,j)) - (state%gd(i,j-1) + static%ghs(i,j-1))) / mesh%cell_lat_distance(j)
        ip_egf = ip_egf + tmp * tend%diag%normal_lat_flux(i,j) * mesh%lat_edge_area(j) * 2 / radius**2
        ip_fu = ip_fu + tend%v_nonlinear(i,j) * tend%diag%normal_lat_flux(i,j) * mesh%lat_edge_area(j) / radius**2 
      end do 
    end do 

    do j = parallel%full_lat_start_idx, parallel%full_lat_end_idx
      r1 = radius * mesh%full_cos_lat(j) * mesh%dlon * radius * mesh%dlat
      do i = parallel%full_lon_start_idx, parallel%full_lon_end_idx
        tmp = state%gd(i,j) + static%ghs(i,j)
        ip_energy_div = ip_energy_div + tmp * tend%mass_div(i,j) * mesh%cell_area(j) / radius**2 
        ip_mass_div = ip_mass_div + tend%mass_div(i,j)  * mesh%cell_area(j) / radius**2 
      end do 
    end do
    print*, 'nonlinear Coriolis:    ' ,' total energy conversion:', '      total mass tendency:'
    print*, (ip_fv + ip_fu),& 
            ip_egf - ip_energy_div,& 
            ip_mass_div

  end subroutine check_spaceoperator  
end module dycore_mod
