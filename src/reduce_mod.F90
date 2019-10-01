module reduce_mod
 
  use mesh_mod
  use log_mod
  use params_mod
  use parallel_mod
  use types_mod
  ! use time_mod
  use sphere_geometry_mod

  implicit none
 
  private
 
  public reduce_init
  public reduce_run
  public reduce_final

  public reduced_start_idx_at_full_lat
  public reduced_end_idx_at_full_lat
  public reduced_start_idx_at_half_lat
  public reduced_end_idx_at_half_lat

  public append_reduced_tend_to_raw_tend_at_full_lat
  public append_reduced_tend_to_raw_tend_at_half_lat

  public full_reduce_factor
  public half_reduce_factor

  public reduced_mesh_type
  public reduced_static_type
  public reduced_state_type
  public reduced_diag_type

  public full_reduced_state
  ! public half_reduced_state
  public full_reduced_static
  ! public half_reduced_static
  public full_reduced_diag
  ! public half_reduced_diag

  public full_reduced_mesh
  ! public calc_reduced_kinetic

  integer, allocatable :: full_reduce_factor(:)
  integer, allocatable :: half_reduce_factor(:)
  real, allocatable :: full_reduce_weight(:,:)
  real, allocatable :: half_reduce_weight(:,:)


  type reduced_static_type
    real, allocatable :: ghs(:,:,:)
  end type reduced_static_type
  
  type reduced_state_type
    real, allocatable :: u(:,:,:)
    real, allocatable :: v(:,:,:)
    real, allocatable :: gd(:,:,:)
  end type reduced_state_type
  
  type reduced_diag_type
    real, allocatable :: kinetic(:,:,:)
    real, allocatable :: normal_lat_flux(:,:,:)
    real, allocatable :: normal_lon_flux(:,:,:)
  end type reduced_diag_type

  type reduced_mesh_type
    integer num_full_lon
    integer num_half_lon
    real dlon
    real, allocatable :: full_lon(:)   ! only for computaion of area, the first level is enogh
    real, allocatable :: full_lon_deg(:)
    real, allocatable :: half_lon(:)
    real, allocatable :: half_lon_deg(:)
    real full_lat
    real full_sin_lat 
    real full_cos_lat
    real, allocatable :: half_lat(:)
    real, allocatable :: half_sin_lat(:)
    real, allocatable :: half_cos_lat(:)
    real cell_area
    real, allocatable :: vertex_area(:) 
    real, allocatable :: subcell_area(:)

    real lon_edge_area 
    real lon_edge_right_area
    real lon_edge_left_area
    real, allocatable :: lat_edge_area(:) 
    real lat_edge_up_area
    real lat_edge_down_area

    real vertex_lat_distance
    real cell_lon_distance
    real, allocatable :: vertex_lon_distance(:)
    real, allocatable :: cell_lat_distance(:)
  end type reduced_mesh_type

  type(reduced_mesh_type), allocatable :: full_reduced_mesh(:,:)
  ! type(reduced_mesh_type), allocatable :: half_reduced_mesh(:)

  type(reduced_state_type), allocatable :: full_reduced_state(:)
  ! type(reduced_state_type), allocatable :: half_reduced_state(:)
  type(reduced_static_type), allocatable :: full_reduced_static(:)
  ! type(reduced_static_type), allocatable :: half_reduced_static(:)
  type(reduced_diag_type), allocatable :: full_reduced_diag(:)
  ! type(reduced_diag_type), allocatable :: half_reduced_diag(:)
  
contains
  
  subroutine reduce_init()
    call reduce_mesh_init()
    call reduce_array_init()
    call log_notice('Reduce module is initialized.')
  end subroutine reduce_init

  subroutine reduce_mesh_init()
    integer :: i, j, k
    real x(3), y(3), z(3)
    integer :: start_idx, end_idx
    real :: total_area, tmp
    
    if(.not. allocated(full_reduce_factor)) allocate(full_reduce_factor(mesh%num_full_lat))
    if(.not. allocated(half_reduce_factor)) allocate(half_reduce_factor(mesh%num_half_lat))

    full_reduce_factor(:) = 1
    half_reduce_factor(:) = 1
    if (use_zonal_reduce) then
      if (parallel%has_south_pole) then
        do j = 1, size(zonal_reduce_factors)
          if (zonal_reduce_factors(j) == 0) exit
          if (mod(mesh%num_full_lon, zonal_reduce_factors(j)) /=0)then
            call log_error('Zonal reduce factor ' // trim(to_string(zonal_reduce_factors(j))) // 'cannot divide zonal grid number ' // trim(to_string(mesh%num_full_lon)) //'!')
          end if 
          full_reduce_factor(parallel%full_lat_start_idx+j)   = zonal_reduce_factors(j)
          half_reduce_factor(parallel%half_lat_start_idx+j-1) = zonal_reduce_factors(j)
        end do 
      end if 
      if (parallel%has_north_pole) then
        do j = 1, size(zonal_reduce_factors)
          if (zonal_reduce_factors(j) == 0) exit
          if (mod(mesh%num_full_lon, zonal_reduce_factors(j)) /=0)then
            call log_error('Zonal reduce factor ' // trim(to_string(zonal_reduce_factors(j))) // 'cannot divide zonal grid number ' // trim(to_string(mesh%num_full_lon)) //'!')
          end if 
          full_reduce_factor(parallel%full_lat_end_idx-j)   = zonal_reduce_factors(j)
          half_reduce_factor(parallel%half_lat_end_idx-j+1) = zonal_reduce_factors(j)
        end do 
      end if 
    end if 

    if (.not. allocated(full_reduced_mesh))  allocate(full_reduced_mesh(parallel%full_lat_start_idx_no_pole:parallel%full_lat_end_idx_no_pole,3))

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) /= 1) then
        if(.not. allocated(full_reduced_mesh(j,1)%full_lon))     allocate(full_reduced_mesh(j,1)%full_lon(    reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        if(.not. allocated(full_reduced_mesh(j,2)%full_lon))     allocate(full_reduced_mesh(j,2)%full_lon(    reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        if(.not. allocated(full_reduced_mesh(j,3)%full_lon))     allocate(full_reduced_mesh(j,3)%full_lon(    reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        if(.not. allocated(full_reduced_mesh(j,2)%full_lon_deg)) allocate(full_reduced_mesh(j,2)%full_lon_deg(reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        
        if(.not. allocated(full_reduced_mesh(j,1)%half_lon))     allocate(full_reduced_mesh(j,1)%half_lon(    reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        if(.not. allocated(full_reduced_mesh(j,2)%half_lon))     allocate(full_reduced_mesh(j,2)%half_lon(    reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        if(.not. allocated(full_reduced_mesh(j,3)%half_lon))     allocate(full_reduced_mesh(j,3)%half_lon(    reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        if(.not. allocated(full_reduced_mesh(j,2)%half_lon_deg)) allocate(full_reduced_mesh(j,2)%half_lon_deg(reduced_start_idx_at_full_lat(j):reduced_end_idx_at_full_lat(j)))
        
        if(.not. allocated(full_reduced_mesh(j,1)%half_lat))     allocate(full_reduced_mesh(j,1)%half_lat(2))
        if(.not. allocated(full_reduced_mesh(j,2)%half_lat))     allocate(full_reduced_mesh(j,2)%half_lat(2))
        if(.not. allocated(full_reduced_mesh(j,3)%half_lat))     allocate(full_reduced_mesh(j,3)%half_lat(2))

        if(.not. allocated(full_reduced_mesh(j,1)%half_sin_lat)) allocate(full_reduced_mesh(j,1)%half_sin_lat(2))
        if(.not. allocated(full_reduced_mesh(j,2)%half_sin_lat)) allocate(full_reduced_mesh(j,2)%half_sin_lat(2))
        if(.not. allocated(full_reduced_mesh(j,3)%half_sin_lat)) allocate(full_reduced_mesh(j,3)%half_sin_lat(2))

        if(.not. allocated(full_reduced_mesh(j,1)%half_cos_lat)) allocate(full_reduced_mesh(j,1)%half_cos_lat(2))
        if(.not. allocated(full_reduced_mesh(j,2)%half_cos_lat)) allocate(full_reduced_mesh(j,2)%half_cos_lat(2))
        if(.not. allocated(full_reduced_mesh(j,3)%half_cos_lat)) allocate(full_reduced_mesh(j,3)%half_cos_lat(2))
        
        if(.not. allocated(full_reduced_mesh(j,2)%subcell_area)) allocate(full_reduced_mesh(j,2)%subcell_area(2))
        
        if(.not. allocated(full_reduced_mesh(j,2)%lat_edge_area))allocate(full_reduced_mesh(j,2)%lat_edge_area(2))
        if(.not. allocated(full_reduced_mesh(j,2)%vertex_lon_distance)) allocate(full_reduced_mesh(j,2)%vertex_lon_distance(2))
        if(.not. allocated(full_reduced_mesh(j,2)%cell_lat_distance))   allocate(full_reduced_mesh(j,2)%cell_lat_distance(2))
      end if 
    end do

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) /= 1) then
        full_reduced_mesh(j,2)%num_full_lon = mesh%num_full_lon / full_reduce_factor(j)
        full_reduced_mesh(j,2)%num_half_lon = full_reduced_mesh(j,2)%num_full_lon
        
        full_reduced_mesh(j,2)%dlon = mesh%dlon * full_reduce_factor(j)

        do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
          start_idx = 1 + (i-1) * full_reduce_factor(j)
          end_idx = start_idx + full_reduce_factor(j) -1
          ! print*, start_idx, end_idx 
          tmp = sum(mesh%full_lon(start_idx:end_idx)) / full_reduce_factor(j)
          full_reduced_mesh(j,2)%full_lon(i) = tmp
          ! print*, tmp*rad_to_deg
          
          full_reduced_mesh(j,2)%half_lon(i) = full_reduced_mesh(j,2)%full_lon(i) + 0.5 * full_reduced_mesh(j,2)%dlon
           
          full_reduced_mesh(j,2)%full_lon_deg(i) = full_reduced_mesh(j,2)%full_lon(i) * rad_to_deg
          full_reduced_mesh(j,2)%half_lon_deg(i) = full_reduced_mesh(j,2)%half_lon(i) * rad_to_deg
          ! print*, full_reduced_mesh(j)%full_lon_deg(i, k), full_reduced_mesh(j)%half_lon_deg(i, k)
          full_reduced_mesh(j,3)%full_lon(i) = full_reduced_mesh(j,2)%full_lon(i)
          full_reduced_mesh(j,1)%full_lon(i) = full_reduced_mesh(j,2)%full_lon(i)
          full_reduced_mesh(j,3)%half_lon(i) = full_reduced_mesh(j,2)%half_lon(i)
          full_reduced_mesh(j,1)%half_lon(i) = full_reduced_mesh(j,2)%half_lon(i)
        end do
        

        full_reduced_mesh(j,1)%full_lat = mesh%full_lat(j-1) 
        full_reduced_mesh(j,2)%full_lat = mesh%full_lat(j)
        full_reduced_mesh(j,3)%full_lat = mesh%full_lat(j+1)
        
        full_reduced_mesh(j,1)%half_lat(1) = mesh%half_lat(j-2)
        full_reduced_mesh(j,1)%half_lat(2) = mesh%half_lat(j-1)

        full_reduced_mesh(j,2)%half_lat(1) = mesh%half_lat(j-1)
        full_reduced_mesh(j,2)%half_lat(2) = mesh%half_lat(j)

        full_reduced_mesh(j,3)%half_lat(1) = mesh%half_lat(j)
        full_reduced_mesh(j,3)%half_lat(2) = mesh%half_lat(j+1)

       
        ! print*, full_reduced_mesh(j)%full_lat * rad_to_deg, full_reduced_mesh(j)%half_lat(1)*rad_to_deg, full_reduced_mesh(j)%half_lat(2) * rad_to_deg
        
        full_reduced_mesh(j,2)%full_sin_lat    = sin(full_reduced_mesh(j,2)%full_lat)
        full_reduced_mesh(j,2)%half_sin_lat(1) = sin(full_reduced_mesh(j,2)%half_lat(1))
        full_reduced_mesh(j,2)%half_sin_lat(2) = sin(full_reduced_mesh(j,2)%half_lat(2))

        full_reduced_mesh(j,2)%full_cos_lat    = cos(full_reduced_mesh(j,2)%full_lat)
        full_reduced_mesh(j,2)%half_cos_lat(1) = cos(full_reduced_mesh(j,2)%half_lat(1))
        full_reduced_mesh(j,2)%half_cos_lat(2) = cos(full_reduced_mesh(j,2)%half_lat(2))
      end if 
    end do
    
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) /= 1) then
        
        full_reduced_mesh(j,2)%cell_area = radius**2 * full_reduced_mesh(j,2)%dlon * (full_reduced_mesh(j,2)%half_sin_lat(2) - full_reduced_mesh(j,2)%half_sin_lat(1)) 
        full_reduced_mesh(j,2)%subcell_area(1) = radius**2 * full_reduced_mesh(j,2)%dlon * 0.5 * (full_reduced_mesh(j,2)%full_sin_lat - full_reduced_mesh(j,2)%half_sin_lat(1))
        full_reduced_mesh(j,2)%subcell_area(2) = radius**2 * full_reduced_mesh(j,2)%dlon * 0.5 * (full_reduced_mesh(j,2)%half_sin_lat(2) - full_reduced_mesh(j,2)%full_sin_lat) 

        call cartesian_transform(full_reduced_mesh(j,2)%full_lon(1), full_reduced_mesh(j,2)%full_lat, x(1), y(1), z(1))
        call cartesian_transform(full_reduced_mesh(j,2)%half_lon(1), full_reduced_mesh(j,2)%half_lat(1), x(2), y(2), z(2))
        call cartesian_transform(full_reduced_mesh(j,2)%half_lon(1), full_reduced_mesh(j,2)%half_lat(2), x(3), y(3), z(3))
        
        full_reduced_mesh(j,2)%lon_edge_left_area  = calc_area(x, y, z)
        full_reduced_mesh(j,2)%lon_edge_right_area = full_reduced_mesh(j,2)%lon_edge_left_area
        full_reduced_mesh(j,2)%lon_edge_area       = full_reduced_mesh(j,2)%lon_edge_left_area + full_reduced_mesh(j,2)%lon_edge_right_area
        
        call cartesian_transform(full_reduced_mesh(j,2)%full_lon(2), full_reduced_mesh(j,2)%full_lat, x(1), y(1), z(1))
        call cartesian_transform(full_reduced_mesh(j,2)%half_lon(2), full_reduced_mesh(j,2)%half_lat(2), x(2), y(2), z(2))
        call cartesian_transform(full_reduced_mesh(j,2)%half_lon(1), full_reduced_mesh(j,2)%half_lat(2), x(3), y(3), z(3))
        full_reduced_mesh(j,2)%lat_edge_down_area = calc_area_with_last_small_arc(x, y, z)

        call cartesian_transform(full_reduced_mesh(j,2)%full_lon(2), full_reduced_mesh(j,2)%full_lat, x(1), y(1), z(1))
        call cartesian_transform(full_reduced_mesh(j,2)%half_lon(1), full_reduced_mesh(j,2)%half_lat(1), x(2), y(2), z(2))
        call cartesian_transform(full_reduced_mesh(j,2)%half_lon(2), full_reduced_mesh(j,2)%half_lat(1), x(3), y(3), z(3))
        full_reduced_mesh(j,2)%lat_edge_up_area = calc_area_with_last_small_arc(x, y, z)

        call cartesian_transform(full_reduced_mesh(j,3)%full_lon(2), full_reduced_mesh(j,3)%full_lat, x(1), y(1), z(1))
        call cartesian_transform(full_reduced_mesh(j,3)%half_lon(1), full_reduced_mesh(j,3)%half_lat(1), x(2), y(2), z(2))
        call cartesian_transform(full_reduced_mesh(j,3)%half_lon(2), full_reduced_mesh(j,3)%half_lat(1), x(3), y(3), z(3))
        full_reduced_mesh(j,3)%lat_edge_up_area = calc_area_with_last_small_arc(x, y, z)

        full_reduced_mesh(j,2)%lat_edge_area(2) = full_reduced_mesh(j,2)%lat_edge_down_area + full_reduced_mesh(j,3)%lat_edge_up_area
        
        call cartesian_transform(full_reduced_mesh(j,1)%full_lon(2), full_reduced_mesh(j,1)%full_lat, x(1), y(1), z(1))
        call cartesian_transform(full_reduced_mesh(j,1)%half_lon(2), full_reduced_mesh(j,1)%half_lat(2), x(2), y(2), z(2))
        call cartesian_transform(full_reduced_mesh(j,1)%half_lon(1), full_reduced_mesh(j,1)%half_lat(2), x(3), y(3), z(3))
        full_reduced_mesh(j,1)%lat_edge_down_area = calc_area_with_last_small_arc(x, y, z)

        full_reduced_mesh(j,2)%lat_edge_area(1) = full_reduced_mesh(j,2)%lat_edge_up_area + full_reduced_mesh(j,1)%lat_edge_down_area
      end if 
      ! print*, 4 * mesh%cell_area(j), full_reduced_mesh(j,2)%cell_area, full_reduced_mesh(j,2)%lon_edge_area+ full_reduced_mesh(j,2)%lat_edge_down_area +full_reduced_mesh(j,2)%lat_edge_up_area
    end do
    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) /= 1) then
        full_reduced_mesh(j,2)%vertex_lat_distance = radius * mesh%dlat
        full_reduced_mesh(j,2)%cell_lon_distance = 2.0 * full_reduced_mesh(j,2)%lon_edge_area / full_reduced_mesh(j,2)%vertex_lat_distance
        ! print*, full_reduced_mesh(j,2)%cell_lon_distance, mesh%cell_lon_distance(j) * full_reduce_factor(j) 
      end if 
    end do 
    

    do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      if (full_reduce_factor(j) /= 1) then
        full_reduced_mesh(j,2)%vertex_lon_distance(1) = radius * full_reduced_mesh(j,2)%half_cos_lat(1) * full_reduced_mesh(j,2)%dlon
        full_reduced_mesh(j,2)%cell_lat_distance(1) = 2.0 * full_reduced_mesh(j,2)%lat_edge_area(1) / full_reduced_mesh(j,2)%vertex_lon_distance(1)

        full_reduced_mesh(j,2)%vertex_lon_distance(2) = radius * full_reduced_mesh(j,2)%half_cos_lat(2) * full_reduced_mesh(j,2)%dlon
        full_reduced_mesh(j,2)%cell_lat_distance(2) = 2.0 * full_reduced_mesh(j,2)%lat_edge_area(2) / full_reduced_mesh(j,2)%vertex_lon_distance(2)
      end if 
    end do 

    !! Allocate half_reduced_mesh
    ! if (.not. allocated(half_reduced_mesh))  allocate(half_reduced_mesh(parallel%half_lat_start_idx:parallel%half_lat_end_idx))
    
    ! do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
    !   if (half_reduce_factor(j) /= 1) then
    !     half_reduced_mesh(j)%num_full_lon = mesh%num_full_lon / half_reduce_factor(j)
    !     half_reduced_mesh(j)%num_half_lon = half_reduced_mesh(j)%num_full_lon
        
    !     half_reduced_mesh(j)%dlon = mesh%dlon * half_reduce_factor(j)

    !     allocate(half_reduced_mesh(j)%full_lon(    reduced_start_idx_at_half_lat(j):reduced_end_idx_at_half_lat(j)))
    !     allocate(half_reduced_mesh(j)%full_lon_deg(reduced_start_idx_at_half_lat(j):reduced_end_idx_at_half_lat(j)))
    !     allocate(half_reduced_mesh(j)%half_lon(    reduced_start_idx_at_half_lat(j):reduced_end_idx_at_half_lat(j)))
    !     allocate(half_reduced_mesh(j)%half_lon_deg(reduced_start_idx_at_half_lat(j):reduced_end_idx_at_half_lat(j)))
    !     allocate(half_reduced_mesh(j)%half_lat(2))
    !     allocate(half_reduced_mesh(j)%half_sin_lat(2))
    !     do i = reduced_start_idx_at_half_lat(j), reduced_end_idx_at_half_lat(j)
    !       half_reduced_mesh(j)%full_lon(i) = - 0.5 * mesh%dlon + (i-1) * half_reduced_mesh(j)%dlon + 0.5 * half_reduced_mesh(j)%dlon
    !       ! print*,i, full_reduced_mesh(j)%full_lon(i)* rad_to_deg
    !       if (half_reduced_mesh(j)%full_lon(i) > 2 * pi) then
    !         half_reduced_mesh(j)%full_lon(i) = half_reduced_mesh(j)%full_lon(i) - 2 * pi
    !       end if

    !       half_reduced_mesh(j)%half_lon(i) = half_reduced_mesh(j)%full_lon(i) + 0.5 * half_reduced_mesh(j)%dlon
    !       if (half_reduced_mesh(j)%half_lon(i) > 2 * pi) then
    !         half_reduced_mesh(j)%half_lon(i) = half_reduced_mesh(j)%half_lon(i) - 2 * pi
    !       end if 
    !       half_reduced_mesh(j)%full_lon_deg(i) = half_reduced_mesh(j)%full_lon(i) * rad_to_deg
    !       half_reduced_mesh(j)%half_lon_deg(i) = half_reduced_mesh(j)%half_lon(i) * rad_to_deg
    !       ! print*, full_reduced_mesh(j)%full_lon_deg(i), full_reduced_mesh(j)%half_lon_deg(i)
    !     end do
    !     half_reduced_mesh(j)%full_lat = mesh%half_lat(j)
    !     half_reduced_mesh(j)%half_lat(1) = mesh%full_lat(j)
    !     half_reduced_mesh(j)%half_lat(2) = mesh%full_lat(j+1)
    !     ! print*, full_reduced_mesh(j)%full_lat * rad_to_deg, full_reduced_mesh(j)%half_lat(1)*rad_to_deg, full_reduced_mesh(j)%half_lat(2) * rad_to_deg
        
    !     half_reduced_mesh(j)%full_sin_lat    = sin(half_reduced_mesh(j)%full_lat)
    !     half_reduced_mesh(j)%half_sin_lat(1) = sin(half_reduced_mesh(j)%half_lat(1))
    !     half_reduced_mesh(j)%half_sin_lat(2) = sin(half_reduced_mesh(j)%half_lat(2))

    !   end if 
    ! end do 

    ! do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
    !   if (half_reduce_factor(j) /= 1) then
    !     ! full_reduced_mesh%vertex_area(j) = radius**2 * full_reduced_mesh%dlon * (half_reduced_mesh%full_sin_lat(j+1) - half_reduced_mesh%full_sin_lat(j))
    !     call cartesian_transform(half_reduced_mesh(j)%full_lon(2,1), half_reduced_mesh(j)%half_lat(1), x(1), y(1), z(1))
    !     call cartesian_transform(half_reduced_mesh(j)%half_lon(2,1), half_reduced_mesh(j)%full_lat, x(2), y(2), z(2))
    !     call cartesian_transform(half_reduced_mesh(j)%half_lon(1,1), half_reduced_mesh(j)%full_lat, x(3), y(3), z(3))
    !     half_reduced_mesh(j)%lat_edge_down_area = calc_area_with_last_small_arc(x, y, z)
    !     call cartesian_transform(half_reduced_mesh(j)%full_lon(2,1), half_reduced_mesh(j)%half_lat(2), x(1), y(1), z(1))
    !     call cartesian_transform(half_reduced_mesh(j)%half_lon(1,1), half_reduced_mesh(j)%full_lat, x(2), y(2), z(2))
    !     call cartesian_transform(half_reduced_mesh(j)%half_lon(2,1), half_reduced_mesh(j)%full_lat, x(3), y(3), z(3))
    !     half_reduced_mesh(j)%lat_edge_up_area = calc_area_with_last_small_arc(x, y, z)
    !     half_reduced_mesh(j)%lat_edge_area = half_reduced_mesh(j)%lat_edge_down_area + half_reduced_mesh(j)%lat_edge_up_area
    !   end if  
    ! end do
    
  end subroutine reduce_mesh_init

  subroutine reduce_array_init()

    integer i,j,k

    if(.not. allocated(full_reduce_weight)) allocate(full_reduce_weight(maxval(zonal_reduce_factors), mesh%num_full_lat))
    if(.not. allocated(half_reduce_weight)) allocate(half_reduce_weight(maxval(zonal_reduce_factors), mesh%num_half_lat))
    if(.not. allocated(full_reduced_state)) allocate(full_reduced_state(parallel%full_lat_start_idx_no_pole:parallel%full_lat_end_idx_no_pole))
    ! if(.not. allocated(half_reduced_state)) allocate(half_reduced_state(parallel%half_lat_start_idx:parallel%half_lat_end_idx))
    if(.not. allocated(full_reduced_static)) allocate(full_reduced_static(parallel%full_lat_start_idx_no_pole:parallel%full_lat_end_idx_no_pole))
    ! if(.not. allocated(half_reduced_static)) allocate(half_reduced_static(parallel%half_lat_start_idx:parallel%half_lat_end_idx))
    
    if(.not. allocated(full_reduced_diag))  allocate(full_reduced_diag(parallel%full_lat_start_idx_no_pole:parallel%full_lat_end_idx_no_pole))
    ! if(.not. allocated(half_reduced_diag))  allocate(half_reduced_diag(parallel%half_lat_start_idx:parallel%half_lat_end_idx))

    if (use_zonal_reduce) then
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      	if (full_reduce_factor(j) /= 1) then
          full_reduce_weight(:,j) = 1.0 / full_reduce_factor(j)
      	end if 
      end do  
      do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      	if (half_reduce_factor(j) /= 1) then
          half_reduce_weight(:,j) = 1.0 / half_reduce_factor(j)
      	end if 
      end do 

      !Allocate reduced data arrays.
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      	if (full_reduce_factor(j) /= 1) then
          call parallel_allocate(full_reduced_state(j)%u,      dim=[2,3], size=[full_reduce_factor(j),3], half_lon=.true., extended_halo=.true.)
          call parallel_allocate(full_reduced_state(j)%v,      dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true., extended_halo=.true.) 
          call parallel_allocate(full_reduced_state(j)%gd,     dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true., extended_halo=.true.)
          call parallel_allocate(full_reduced_static(j)%ghs,   dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true., extended_halo=.true.)

          call parallel_allocate(full_reduced_diag(j)%kinetic,         dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true., extended_halo=.true.)
          call parallel_allocate(full_reduced_diag(j)%normal_lat_flux, dim=[2,3], size=[full_reduce_factor(j),3], full_lon=.true., extended_halo=.true.)
          call parallel_allocate(full_reduced_diag(j)%normal_lon_flux, dim=[2,3], size=[full_reduce_factor(j),3], half_lon=.true., extended_halo=.true.)
          
        end if 
      end do 
      ! print*, size(full_reduced_state(j)%u, dim=1)
      ! stop "reduce_mod"
      ! do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx
      ! 	if (half_reduce_factor(j) /= 1) then
      !     call parallel_allocate(half_reduced_state(j)%u,      dim=[2,3], size=[half_reduce_factor(j),3], half_lon=.true.)
      !     call parallel_allocate(half_reduced_state(j)%v,      dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
      !     call parallel_allocate(half_reduced_state(j)%gd,     dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
      !     ! call parallel_allocate(half_reduced_static(j)%ghs,   dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
      !     ! call parallel_allocate(half_reduced_diag(j)%kinetic, dim=[2,3], size=[half_reduce_factor(j),3], full_lon=.true.)
      ! 	end if 
      ! end do 
    end if 
 	end subroutine reduce_array_init


 	subroutine reduce_final()
		if (allocated(full_reduce_factor)) deallocate(full_reduce_factor)
    if (allocated(full_reduce_weight)) deallocate(full_reduce_weight)
    ! if (allocated(half_reduce_factor)) deallocate(half_reduce_factor)
    ! if (allocated(half_reduce_weight)) deallocate(half_reduce_weight)
    if (allocated(full_reduced_state)) deallocate(full_reduced_state)
    ! if (allocated(half_reduced_state)) deallocate(half_reduced_state)
    if (allocated(full_reduced_static)) deallocate(full_reduced_static)
    if (allocated(full_reduced_diag))  deallocate(full_reduced_diag)
    ! if (allocated(half_reduced_diag))  deallocate(half_reduced_diag)
 	end subroutine reduce_final
  
  integer function reduced_start_idx_at_full_lat(j) result(start_idx)
    
    integer, intent(in) :: j 
    start_idx = 1

  end function reduced_start_idx_at_full_lat

  integer function reduced_end_idx_at_full_lat(j) result(end_idx)
  
    integer, intent(in) :: j

    end_idx = mesh%num_full_lon / full_reduce_factor(j)
      
  end function reduced_end_idx_at_full_lat

  integer function reduced_start_idx_at_half_lat(j) result(start_idx)
    
    integer, intent(in) :: j
    
    start_idx = 1
      
  end function reduced_start_idx_at_half_lat

  integer function reduced_end_idx_at_half_lat(j) result(end_idx)

    integer, intent(in) :: j
    
    end_idx = mesh%num_full_lon / half_reduce_factor(J)
      
  end function reduced_end_idx_at_half_lat

 	subroutine reduce_run(state, static)
    type(state_type),  intent(in) :: state
    type(static_type), intent(in) :: static

    integer :: j, k, i, m, n
    n = parallel%lon_halo_width_for_reduce

    if (use_zonal_reduce) then
    	do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
    		if (full_reduce_factor(j) > 1) then
    			do k = 1, full_reduce_factor(j)
           
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%u(:,j-1),  full_reduced_state(j)%u(:,k,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%u(:,j  ),  full_reduced_state(j)%u(:,k,2))

            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%u(:,j+1),  full_reduced_state(j)%u(:,k,3))

            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%v(:,j-1),  full_reduced_state(j)%v(:,k,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%v(:,j  ),  full_reduced_state(j)%v(:,k,2))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%v(:,j+1),  full_reduced_state(j)%v(:,k,3))

            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%gd(:,j-1), full_reduced_state(j)%gd(:,k,1))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%gd(:,j  ), full_reduced_state(j)%gd(:,k,2))
            call average_raw_array_to_reduced_array_at_full_lat(j, k, state%gd(:,j+1), full_reduced_state(j)%gd(:,k,3))

            call average_raw_array_to_reduced_array_at_full_lat(j, k, static%ghs(:,j), full_reduced_static(j)%ghs(:,k,2))      
    			end do
    		end if 
      end do
   
      !(1) diagnose kinetic energy at reduced cell centers
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (full_reduce_factor(j) > 1) then
          do k = 1, full_reduce_factor(j)
            do i = reduced_start_idx_at_full_lat(j)-1, reduced_end_idx_at_full_lat(j)
              full_reduced_diag(j)%kinetic(i,k,2) = 1.0 / full_reduced_mesh(j,2)%cell_area * (full_reduced_mesh(j,2)%lon_edge_left_area  * full_reduced_state(j)%u(i,  k,2)**2 +&
                                                                                              full_reduced_mesh(j,2)%lon_edge_right_area * full_reduced_state(j)%u(i-1,k,2)**2 +&
                                                                                              full_reduced_mesh(j,2)%lat_edge_down_area  * full_reduced_state(j)%v(i,  k,2)**2 +&
                                                                                              full_reduced_mesh(j,2)%lat_edge_up_area    * full_reduced_state(j)%v(i,  k,1)**2 )
              ! print*, full_reduced_state(j)%u(i,k,2)**2, full_reduced_diag(j)%kinetic(i,k)
            end do
            ! print*, full_reduced_diag(j)%kinetic(1:5,k,2)
            ! Fill halo for reduced array.
            ! m = (size(full_reduced_diag(j)%kinetic(:,k,2)) - 2 * parallel%lon_halo_width_for_reduce) / full_reduce_factor(j) 
            ! full_reduced_diag(j)%kinetic(1-n:0,k,2) = full_reduced_diag(j)%kinetic(m-n+1:m, k, 2)
            ! full_reduced_diag(j)%kinetic(m+1:m+n, k,2) = full_reduced_diag(j)%kinetic(1:n, k,2)
            ! print*, full_reduced_diag(j)%kinetic(:,k,2)
            
          end do 
        end if 
      end do  
      !(2) diagnose normal mass flux at reduced cell edges
      ! do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
      !   if (full_reduce_factor(j) > 1) then
      !     do k = 1, full_reduce_factor(j)
      !       do i = reduced_start_idx_at_full_lat(j), reduced_end_idx_at_full_lat(j)
      !       full_reduced_diag(j)%normal_lat_flux(i,k,1) = (full_reduced_state(j)%gd(i,k,1) * full_reduced_mesh(j,1)%lat_edge_down_area +&
      !                                                      full_reduced_state(j)%gd(i,k,2) * full_reduced_mesh(j,2)%lat_edge_up_area ) /&
      !                                                      full_reduced_mesh(j,2)%lat_edge_area(1) * full_reduced_state(j)%v(i,k,1)
      !       full_reduced_diag(j)%normal_lat_flux(i,k,2) = (full_reduced_state(j)%gd(i,k,3) * full_reduced_mesh(j,3)%lat_edge_up_area +&
      !                                                      full_reduced_state(j)%gd(i,j,2) * full_reduced_mesh(j,2)%lat_edge_down_area ) /&
      !                                                      full_reduced_mesh(j,2)%lat_edge_area(2) * full_reduced_state(j)%v(i,k,2)                                      
      !       end do 
      !     end do 
      !   end if 
      ! end do 
      do j = parallel%full_lat_start_idx_no_pole, parallel%full_lat_end_idx_no_pole
        if (full_reduce_factor(j) > 1) then
          do k = 1, full_reduce_factor(j)
            do i = reduced_start_idx_at_full_lat(j)-1, reduced_end_idx_at_full_lat(j)
            full_reduced_diag(j)%normal_lon_flux(i,k,2) = (full_reduced_state(j)%gd(i,  k,2) * full_reduced_mesh(j,2)%lon_edge_left_area +&
                                                           full_reduced_state(j)%gd(i+1,k,2) * full_reduced_mesh(j,2)%lon_edge_right_area ) /&
                                                           full_reduced_mesh(j,2)%lon_edge_area * full_reduced_state(j)%u(i,k,2)                                          
            end do
          end do 
        end if 
      end do 
  
     ! do j = parallel%half_lat_start_idx, parallel%half_lat_end_idx 
     !    if (half_reduce_factor(j) > 1) then
     !    	do k =1, half_reduce_factor(j)
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%u(:,j-1),  half_reduced_state(j)%u(:,k,1))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%u(:,j  ),  half_reduced_state(j)%u(:,k,2))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%u(:,j+1),  half_reduced_state(j)%u(:,k,3))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%v(:,j-1),  half_reduced_state(j)%v(:,k,1))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%v(:,j  ),  half_reduced_state(j)%v(:,k,2))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%v(:,j+1),  half_reduced_state(j)%v(:,k,3))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%gd(:,j-1), half_reduced_state(j)%gd(:,k,1))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%gd(:,j  ), half_reduced_state(j)%gd(:,k,2))
     !        ! call average_raw_array_to_reduced_array_at_half_lat(j, k, state%gd(:,j+1), half_reduced_state(j)%gd(:,k,3))
     !      end do 
     !    end if 
     ! end do 
    end if

 	end subroutine reduce_run


 	subroutine average_raw_array_to_reduced_array_at_full_lat(j, k, raw_array, reduced_array)
    integer,intent(in) :: j
    integer,intent(in) :: k
    real, intent(in) :: raw_array(:)
    real, intent(out) :: reduced_array(:)

    integer :: i, l, count, m, n 
    
    n = parallel%lon_halo_width_for_reduce

    reduced_array(:) = 0.0
    l = n + 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_array) - parallel%lon_halo_width_for_reduce
      count = count + 1
      reduced_array(l) = reduced_array(l) + raw_array(i) * full_reduce_weight(count,j)
      if (count == full_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if 
    end do 

    ! Fill halo for reduced_array.
    m = (size(raw_array) - 2 * parallel%lon_halo_width_for_reduce) / full_reduce_factor(j) + 2 * n 

    reduced_array(1:n) = reduced_array(m-2*n+1:m-n)
    reduced_array(m-n+1:m) = reduced_array(1+n:2*n)
   
 	end subroutine average_raw_array_to_reduced_array_at_full_lat

  subroutine average_raw_array_to_reduced_array_at_half_lat(j, k, raw_array, reduced_array)
    integer, intent(in) :: j
    integer, intent(in) :: k
    real, intent(in) :: raw_array(:)
    real, intent(out) :: reduced_array(:)

    integer :: i, l, count, m, n 

    n = parallel%lon_halo_width_for_reduce
    reduced_array(:) = 0.0
    l = n + 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_array) - parallel%lon_halo_width_for_reduce
      count = count + 1
      reduced_array(l) = reduced_array(l) + raw_array(i) * half_reduce_weight(count,j)
      if (count == half_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if 
    end do 

    ! Fill halo for reduced_array
    m = (size(raw_array) - 2 * parallel%lon_halo_width_for_reduce) / half_reduce_factor(j) + 2 * n 
    reduced_array(1:n) = reduced_array(m-2*n+1:m-n)
    reduced_array(m-n+1:m) = reduced_array(1+n:2*n)

  end subroutine average_raw_array_to_reduced_array_at_half_lat

  subroutine append_reduced_tend_to_raw_tend_at_full_lat(j, k, reduced_tend, raw_tend)
    integer, intent(in) :: j 
    integer, intent(in) :: k 
    real, intent(in) :: reduced_tend(:)
    real, intent(inout) :: raw_tend(:)

    integer :: i, l, count

    l = 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_tend) - parallel%lon_halo_width_for_reduce
      count = count + 1
      raw_tend(i) = raw_tend(i) + reduced_tend(l) * full_reduce_weight(count,j)
      if (count == full_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if 
    end do 
  end subroutine append_reduced_tend_to_raw_tend_at_full_lat

  subroutine append_reduced_tend_to_raw_tend_at_half_lat(j, k, reduced_tend, raw_tend)
    integer, intent(in) :: j
    integer, intent(in) :: k 
    real, intent(in) :: reduced_tend(:)
    real, intent(inout) :: raw_tend(:)

    integer :: i, l, count

    l = 1
    count = 0
    do i = k + parallel%lon_halo_width_for_reduce, k - 1 + size(raw_tend) - parallel%lon_halo_width_for_reduce
      count = count + 1
      raw_tend(i) = raw_tend(i) + reduced_tend(l) * half_reduce_weight(count,j)
      if (count == half_reduce_factor(j)) then
        l = l + 1
        count = 0
      end if  
    end do 
  end subroutine append_reduced_tend_to_raw_tend_at_half_lat

end module reduce_mod
