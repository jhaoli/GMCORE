module mesh_mod

  use log_mod
  use params_mod
  use sphere_geometry_mod

  implicit none

  private

  public mesh
  public mesh_init
  public mesh_final

  type mesh_type
    integer num_full_lon
    integer num_half_lon
    integer num_full_lat
    integer num_half_lat
    real dlon
    real dlat
    real, allocatable :: full_lon(:)
    real, allocatable :: half_lon(:)
    real, allocatable :: full_lat(:)
    real, allocatable :: half_lat(:)
    real, allocatable :: full_cos_lon(:)
    real, allocatable :: half_cos_lon(:)
    real, allocatable :: full_sin_lon(:)
    real, allocatable :: half_sin_lon(:)
    real, allocatable :: full_cos_lat(:)
    real, allocatable :: half_cos_lat(:)
    real, allocatable :: full_sin_lat(:)
    real, allocatable :: half_sin_lat(:)
    ! For output
    real, allocatable :: full_lon_deg(:)
    real, allocatable :: half_lon_deg(:)
    real, allocatable :: full_lat_deg(:)
    real, allocatable :: half_lat_deg(:)
    ! Area for weighting
    real, allocatable :: cell_area(:) 
    real, allocatable :: vertex_area(:)  
    real, allocatable :: lon_edge_area(:)
    real, allocatable :: lon_edge_left_area(:)
    real, allocatable :: lon_edge_right_area(:)
    real, allocatable :: lat_edge_area(:)
    real, allocatable :: lat_edge_up_area(:)
    real, allocatable :: lat_edge_down_area(:)
    real, allocatable :: subcell_area(:,:)
    ! cell distance and  vertex distance
    real, allocatable :: cell_lon_distance(:)
    real, allocatable :: cell_lat_distance(:)
    real, allocatable :: vertex_lon_distance(:)
    real, allocatable :: vertex_lat_distance(:)
    ! weights for reconstructing tangent wind
    real, allocatable :: full_tangent_wgt(:,:)
    real, allocatable :: half_tangent_wgt(:,:)
  end type mesh_type

  type(mesh_type) mesh

contains

  subroutine mesh_init()

    integer i, j
    real x(3), y(3), z(3), total_area, total_distance

    mesh%num_full_lon = num_lon
    mesh%num_half_lon = num_lon
    mesh%num_full_lat = num_lat
    mesh%num_half_lat = num_lat - 1

    allocate(mesh%full_lon(mesh%num_full_lon))
    allocate(mesh%half_lon(mesh%num_half_lon))
    allocate(mesh%full_lat(mesh%num_full_lat))
    allocate(mesh%half_lat(mesh%num_half_lat))
    allocate(mesh%full_cos_lon(mesh%num_full_lon))
    allocate(mesh%half_cos_lon(mesh%num_half_lon))
    allocate(mesh%full_sin_lon(mesh%num_full_lon))
    allocate(mesh%half_sin_lon(mesh%num_half_lon))
    allocate(mesh%full_cos_lat(mesh%num_full_lat))
    allocate(mesh%half_cos_lat(mesh%num_half_lat))
    allocate(mesh%full_sin_lat(mesh%num_full_lat))
    allocate(mesh%half_sin_lat(mesh%num_half_lat))
    allocate(mesh%full_lon_deg(mesh%num_full_lon))
    allocate(mesh%half_lon_deg(mesh%num_half_lon))
    allocate(mesh%full_lat_deg(mesh%num_full_lat))
    allocate(mesh%half_lat_deg(mesh%num_half_lat))
    allocate(mesh%cell_area(mesh%num_full_lat))
    allocate(mesh%vertex_area(mesh%num_half_lat))
    allocate(mesh%lon_edge_area(mesh%num_full_lat))
    allocate(mesh%lon_edge_left_area(mesh%num_full_lat))
    allocate(mesh%lon_edge_right_area(mesh%num_full_lat))
    allocate(mesh%lat_edge_area(mesh%num_half_lat))
    allocate(mesh%lat_edge_up_area(mesh%num_half_lat))
    allocate(mesh%lat_edge_down_area(mesh%num_half_lat))
    allocate(mesh%subcell_area(2, mesh%num_full_lat))
    allocate(mesh%cell_lon_distance(mesh%num_full_lat))
    allocate(mesh%cell_lat_distance(mesh%num_half_lat))
    allocate(mesh%vertex_lon_distance(mesh%num_half_lat))
    allocate(mesh%vertex_lat_distance(mesh%num_full_lat))
    allocate(mesh%full_tangent_wgt(2,mesh%num_full_lat))
    allocate(mesh%half_tangent_wgt(2,mesh%num_half_lat))

    mesh%dlon = 2 * pi / mesh%num_full_lon
    do i = 1, mesh%num_full_lon
      mesh%full_lon(i) = (i - 1) * mesh%dlon
      mesh%half_lon(i) = mesh%full_lon(i) + 0.5 * mesh%dlon
      mesh%full_lon_deg(i) = mesh%full_lon(i) * rad_to_deg
      mesh%half_lon_deg(i) = mesh%half_lon(i) * rad_to_deg
    end do

    mesh%dlat = pi / mesh%num_half_lat
    do j = 1, mesh%num_half_lat
      mesh%full_lat(j) = - 0.5 * pi + (j - 1) * mesh%dlat
      mesh%half_lat(j) = mesh%full_lat(j) + 0.5 * mesh%dlat
      mesh%full_lat_deg(j) = mesh%full_lat(j) * rad_to_deg
      mesh%half_lat_deg(j) = mesh%half_lat(j) * rad_to_deg
    end do
    mesh%full_lat(num_lat) = 0.5 * pi
    mesh%full_lat_deg(num_lat) = 90.0

    do i = 1, mesh%num_full_lon
      mesh%full_cos_lon(i) = cos(mesh%full_lon(i))
      mesh%full_sin_lon(i) = sin(mesh%full_lon(i))
    end do

    do i = 1, mesh%num_half_lon
      mesh%half_cos_lon(i) = cos(mesh%half_lon(i))
      mesh%half_sin_lon(i) = sin(mesh%half_lon(i))
    end do

    do j = 1, mesh%num_half_lat
      mesh%half_cos_lat(j) = cos(mesh%half_lat(j))
      mesh%half_sin_lat(j) = sin(mesh%half_lat(j))
    end do

    do j = 1, mesh%num_full_lat
      mesh%full_cos_lat(j) = cos(mesh%full_lat(j))
      mesh%full_sin_lat(j) = sin(mesh%full_lat(j))
    end do
    mesh%full_cos_lat(1) = 0.0
    mesh%full_cos_lat(mesh%num_full_lat) = 0.0
    mesh%full_sin_lat(1) = -1.0
    mesh%full_sin_lat(mesh%num_full_lat) = 1.0


    do j = 1, mesh%num_full_lat
      if (j == 1) then
        mesh%cell_area(j) = radius**2 * mesh%dlon * ( mesh%half_sin_lat(1) + 1.0)
        mesh%subcell_area(1,j) = 0.0
        mesh%subcell_area(2,j) = radius**2 * mesh%dlon * 0.5 * (mesh%half_sin_lat(j) +1.0)
      else if (j == mesh%num_full_lat) then 
        mesh%cell_area(j) = radius**2 * mesh%dlon * (1 - mesh%half_sin_lat(j-1))
        mesh%subcell_area(1,j) = radius**2 * mesh%dlon * 0.5 * (1.0 - mesh%half_sin_lat(j-1))
        mesh%subcell_area(2,j) = 0.0
      else
        mesh%cell_area(j) = radius**2 * mesh%dlon * (mesh%half_sin_lat(j) - mesh%half_sin_lat(j-1)) 
        mesh%subcell_area(1,j) = radius**2 * mesh%dlon * 0.5 * (mesh%full_sin_lat(j) - mesh%half_sin_lat(j-1))
        mesh%subcell_area(2,j) = radius**2 * mesh%dlon * 0.5 * (mesh%half_sin_lat(j) - mesh%full_sin_lat(j)) 
        call cartesian_transform(mesh%full_lon(1), mesh%full_lat(j  ), x(1), y(1), z(1))
        call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j-1), x(2), y(2), z(2))
        call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j  ), x(3), y(3), z(3))
        mesh%lon_edge_left_area(j) = calc_area(x, y, z)
        mesh%lon_edge_right_area(j) = mesh%lon_edge_left_area(j)
        mesh%lon_edge_area(j) = mesh%lon_edge_left_area(j) + mesh%lon_edge_right_area(j)
      end if 
    end do 

    do j = 1, mesh%num_half_lat
      mesh%vertex_area(j) = radius**2 * mesh%dlon * (mesh%full_sin_lat(j+1) - mesh%full_sin_lat(j))
      call cartesian_transform(mesh%full_lon(2), mesh%full_lat(j+1), x(1), y(1), z(1))
      call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j  ), x(2), y(2), z(2))
      call cartesian_transform(mesh%half_lon(2), mesh%half_lat(j  ), x(3), y(3), z(3))
      mesh%lat_edge_up_area(j) = calc_area_with_last_small_arc(x, y, z)
      call cartesian_transform(mesh%full_lon(2), mesh%full_lat(j), x(1), y(1), z(1))
      call cartesian_transform(mesh%half_lon(2), mesh%half_lat(j), x(2), y(2), z(2))
      call cartesian_transform(mesh%half_lon(1), mesh%half_lat(j), x(3), y(3), z(3))
      mesh%lat_edge_down_area(j) = calc_area_with_last_small_arc(x, y, z)
      mesh%lat_edge_area(j) = mesh%lat_edge_up_area(j) + mesh%lat_edge_down_area(j)
    end do 

    total_area = 0.0
    do j = 1, mesh%num_full_lat
      total_area = total_area + mesh%cell_area(j) * mesh%num_full_lon
    end do 
    if ( abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0E-11) then
      call log_notice('Failed to calculate cell area!')
    endif

    total_area = 0.0
    do j = 1, mesh%num_half_lat
      total_area = total_area + mesh%vertex_area(j) * mesh%num_half_lon
    end do 
    if ( abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0E-11) then
      call log_notice('Failed to calculate vertex area!')
    endif

    total_area = 0.0
    do j = 1, mesh%num_full_lat
      total_area = total_area + mesh%lon_edge_area(j) * mesh%num_full_lon
    end do 
    do j = 1, mesh%num_half_lat
      total_area = total_area + mesh%lat_edge_area(j) * mesh%num_full_lon
    end do 
    if (abs((4 * pi * radius**2 - total_area) / (4 * pi * radius**2)) > 1.0e-10) then
      call log_error('Failed to calculate edge area!')
    end if

    do j = 1, mesh%num_full_lat
      if (j == 1) then 
        total_area = mesh%subcell_area(2,j) * 2
      elseif(j == mesh%num_full_lat) then
        total_area = mesh%subcell_area(1,j) * 2
      else
        total_area = (mesh%subcell_area(1,j) + mesh%subcell_area(2,j)) * 2
      endif 
      if (abs(total_area - mesh%cell_area(j)) / mesh%cell_area(j) > 1.0E-12) then
        call log_error('Failed to calculate subcell area!')
      end if 
    end do 


    do j = 1, mesh%num_full_lat
      if (j == 1 .or. j == mesh%num_full_lat) then
        mesh%vertex_lat_distance(j) = radius * mesh%dlat * 0.5
        mesh%cell_lon_distance(j) = 0.0
      else
        mesh%vertex_lat_distance(j) = radius * mesh%dlat
        mesh%cell_lon_distance(j) = 2.0 * mesh%lon_edge_area(j) / mesh%vertex_lat_distance(j)
      end if 
    end do 

    do j = 1, mesh%num_half_lat
      mesh%vertex_lon_distance(j) = radius * mesh%half_cos_lat(j) * mesh%dlon
      mesh%cell_lat_distance(j) = 2.0 * mesh%lat_edge_area(j) / mesh%vertex_lon_distance(j)
    end do 

    
    total_distance = 0.0
    do j = 1, mesh%num_half_lat
      total_distance = total_distance + mesh%cell_lat_distance(j)
    end do
    if (abs((pi * radius - total_distance) / ( pi * radius)) > 1.0e-4) then
      call log_error('Failed to calculate cell_lat_distance!', __FILE__, __LINE__)
    end if
    
    do j = 1, mesh%num_full_lat
      mesh%full_tangent_wgt(1,j) = mesh%subcell_area(2,j) / mesh%cell_area(j)
      mesh%full_tangent_wgt(2,j) = mesh%subcell_area(1,j) / mesh%cell_area(j)
!       print*, mesh%full_lat_deg(j), mesh%subcell_area(2,j), mesh%full_tangent_wgt(1,j), mesh%full_tangent_wgt(2,j)
    end do 
    do j = 1, mesh%num_half_lat
      mesh%half_tangent_wgt(1,j) = mesh%subcell_area(1,j) / mesh%cell_area(j)
      mesh%half_tangent_wgt(2,j) = mesh%subcell_area(2,j+1) / mesh%cell_area(j+1)
!       print*, mesh%half_lat_deg(j), mesh%half_tangent_wgt(1,j), mesh%half_tangent_wgt(2,j)
    end do 
!     open(11,file='full_tangent_wgt.txt')
!       do j = 1, mesh%num_full_lat
!          write(11,*) mesh%full_lat_deg(j), mesh%full_tangent_wgt(1,j), mesh%full_tangent_wgt(2,j)
!       end do 
!     close(11) 
!     open(11,file='half_tangent_wgt.txt')
!       do j = 1, mesh%num_half_lat
!          write(11,*) mesh%full_lat_deg(j), mesh%half_tangent_wgt(1,j), mesh%half_tangent_wgt(2,j)
!       end do 
!     close(11)    
    call log_notice('Mesh module is initialized.')
! stop  
  end subroutine mesh_init

  subroutine mesh_final()

    if (allocated(mesh%full_lon)) deallocate(mesh%full_lon)
    if (allocated(mesh%full_lat)) deallocate(mesh%full_lat)
    if (allocated(mesh%half_lon)) deallocate(mesh%half_lon)
    if (allocated(mesh%half_lat)) deallocate(mesh%half_lat)
    if (allocated(mesh%full_cos_lat)) deallocate(mesh%full_cos_lat)
    if (allocated(mesh%half_cos_lat)) deallocate(mesh%half_cos_lat)
    if (allocated(mesh%full_sin_lat)) deallocate(mesh%full_sin_lat)
    if (allocated(mesh%half_sin_lat)) deallocate(mesh%half_sin_lat)
    if (allocated(mesh%full_lon_deg)) deallocate(mesh%full_lon_deg)
    if (allocated(mesh%half_lon_deg)) deallocate(mesh%half_lon_deg)
    if (allocated(mesh%full_lat_deg)) deallocate(mesh%full_lat_deg)
    if (allocated(mesh%half_lat_deg)) deallocate(mesh%half_lat_deg)
    if (allocated(mesh%cell_area)) deallocate(mesh%cell_area)
    if (allocated(mesh%vertex_area)) deallocate(mesh%vertex_area)
    if (allocated(mesh%lon_edge_area)) deallocate(mesh%lon_edge_area)
    if (allocated(mesh%lon_edge_left_area)) deallocate(mesh%lon_edge_left_area)
    if (allocated(mesh%lon_edge_right_area)) deallocate(mesh%lon_edge_right_area)
    if (allocated(mesh%lat_edge_area)) deallocate(mesh%lat_edge_area)
    if (allocated(mesh%lat_edge_up_area)) deallocate(mesh%lat_edge_up_area)
    if (allocated(mesh%lat_edge_down_area)) deallocate(mesh%lat_edge_down_area)
    if (allocated(mesh%subcell_area)) deallocate(mesh%subcell_area)
    if (allocated(mesh%cell_lon_distance)) deallocate(mesh%cell_lon_distance)
    if (allocated(mesh%cell_lat_distance)) deallocate(mesh%cell_lat_distance)
    if (allocated(mesh%vertex_lon_distance)) deallocate(mesh%vertex_lon_distance)
    if (allocated(mesh%vertex_lat_distance)) deallocate(mesh%vertex_lat_distance)
    if (allocated(mesh%full_tangent_wgt)) deallocate(mesh%full_tangent_wgt)
    if (allocated(mesh%half_tangent_wgt)) deallocate(mesh%half_tangent_wgt)

    call log_notice('Mesh module is finalized.')

  end subroutine mesh_final

end module mesh_mod
