module kohonen_layer_utilities
!
use kohonen_prototype_utilities;
!
implicit none;
!
type,extends(kohonen_layer_base) ::  kohonen_layer 
  private
    character(len=50) :: layer_type
    !class(kohonen_container),allocatable :: grid(:,:,:)
    type(kohonen_prototype),allocatable :: grid(:,:)
    integer,allocatable :: number_patterns(:,:),cells_index(:,:)
    real(kind=8),allocatable :: u_matrix(:,:),current_distances(:,:)
    type(kohonen_map_parameters) :: parameters
    type(factory_distance) :: factory
    class(distance_base),allocatable :: distance_function
  contains
    procedure,public :: create => create_layer
    procedure,public :: destroy => destroy_layer
!    procedure,public :: train (function for maps)
    procedure,public :: calculate_distance => calculate_distance_layer
    procedure,public :: get => get_container_layer
    procedure,public :: set => set_container_layer
    procedure,public :: get_distances => get_distances_layer
end type kohonen_layer
!
 subroutine calculate_distance_layer(layer,input,dist)
 class(kohonen_layer) :: layer
 class(kohonen_container_base),allocatable :: input
 !scaled distances
 real(kind=8),dimension(:,:),intent(inout) :: dist
!

!
 end subroutine calculate_distance_layer
! 
end module kohonen_layer_utilities