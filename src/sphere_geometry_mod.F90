module sphere_geometry_mod
  
  use params_mod
  use types_mod
  use log_mod
  use mesh_mod

  implicit none

  private
  public calc_edge_area

  type point_type
    real lon
    real lat
    real x
    real y
    real z
  end type point_type

  interface calc_edge_area
!     module procedure calc_lon_edge_area
    module procedure calc_lat_edge_area
  end interface calc_edge_area
contains

  subroutine calc_lon_edge_area
  	integer :: i, j
  	type(point_type) pointA, pointB, pointC

    do j = int(mesh%num_full_lat / 2)+1, mesh%num_full_lat-1

    	pointA%lon = mesh%full_lon(2)
    	pointA%lat = mesh%full_lat(j)

    	pointB%lon = mesh%half_lon(2)
			pointB%lat = mesh%half_lat(j)

			pointC%lon = mesh%half_lon(2)
			pointC%lat = mesh%half_lat(j-1)

			call cartesian_transform(pointA)
			call cartesian_transform(pointB)
			call cartesian_transform(pointC)
			
			
      mesh%lon_edge_area(j) = radius * calc_area([pointA%x, pointA%y, pointA%z],&
      																					 [pointB%x, pointB%y, pointB%z],&
      																					 [pointC%x, pointC%y, pointC%z])
      print*, j, mesh%lon_edge_area(j)
    end do
    print*, sum(mesh%lon_edge_area)*2 * mesh%num_full_lon / (4 * pi * radius**2)
    stop
  end subroutine calc_lon_edge_area

  subroutine calc_lat_edge_area
   type(point_type) pointA, pointB, pointC
   	
  end subroutine calc_lat_edge_area

  subroutine cartesian_transform(point)

    class(point_type), intent(inout) :: point

    real cos_lat

    cos_lat = cos(point%lat)
    point%x = radius * cos_lat * cos(point%lon)
    point%y = radius * cos_lat * sin(point%lon)
    point%z = radius * sin(point%lat)

  end subroutine cartesian_transform

  subroutine cartesian_transform_2(lon, lat, x, y, z)

    real, intent(in)  :: lon, lat
    real, intent(out) :: x, y, z

    real cos_lat

    cos_lat = cos(lat)
    x = radius * cos_lat * cos(lon)
    y = radius * cos_lat * sin(lon)
    z = radius * sin(lat)

  end subroutine cartesian_transform_2


  real function calc_area(x, y, z) result(res)
  	
    real, intent(in) :: x(:)
    real, intent(in) :: y(:)
    real, intent(in) :: z(:)

    integer n, im1, i, ip1
    real angle
     
    res = 0.0
    n = size(x)
#ifndef NDEBUG
    if (n < 3) then
      call log_error('Spherical polygon number is less than 3!')
    end if
#endif
    res = 0.0
    do i = 1, n
      im1 = merge(i - 1, n, i /= 1)
      ip1 = merge(i + 1, 1, i /= n)
      angle = calc_sphere_angle([x(im1),y(im1),z(im1)], [x(i),y(i),z(i)], [x(ip1),y(ip1),z(ip1)])
      res = res + angle
    end do
    res = radius**2 * (res - (n - 2) * pi)
! #ifndef NDEBUG
!     if (abs(res) > 1.0 .or. res <= 0.0) then
!       call log_error('Encounter bad spherical polygon to calculate area!')
!     end if
! #endif
  end function calc_area

	real function calc_sphere_angle(a, b, c) result(res)

    real, intent(in) :: a(3)
    real, intent(in) :: b(3)
    real, intent(in) :: c(3)

    real nab(3) ! Normal vector of plane AB
    real nbc(3) ! Normal vector of plane BC

    nab = norm_vector(cross_product(a, b))
    nbc = norm_vector(cross_product(b, c))
    res = acos(- max(min(dot_product(nab, nbc), 1.0d0), -1.0d0))

    ! Judge the cyclic direction with respect to point A to handle obtuse angle.
    if (dot_product(cross_product(nab, nbc), a) < 0.0) res = 2 * pi - res

  end function calc_sphere_angle

  function norm_vector(x) result(res)

    real(real_kind), intent(in) :: x(:)
    real(real_kind) res(size(x))

    real(real_kind) n

    n = sqrt(sum(x * x))
    if (n /= 0) then
      res = x / n
    else
      res = x
    end if

  end function norm_vector

  function cross_product(x, y) result(res)
    real, intent(in) :: x(3)
    real, intent(in) :: y(3)
  	real res(3)
    res(1) = x(2) * y(3) - x(3) * y(2)
    res(2) = x(3) * y(1) - x(1) * y(3)
    res(3) = x(1) * y(2) - x(2) * y(1)

  end function cross_product

end module sphere_geometry_mod
