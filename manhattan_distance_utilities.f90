module manhattan_distance_utilities
!
use distance_base_utilities;
!
implicit none
!
type,extends(distance_base) :: manhattan_distance
  contains
  procedure,public :: calculate => calculate_manhattan_distance
end type manhattan_distance
!
 contains
! 
  function calculate_manhattan_distance(distance,vector1,vector2) result(d)
!  
    class(manhattan_distance) :: distance
    real(kind=8),dimension(:,:),intent(inout) :: vector1,vector2
    real(kind=8) :: d
 !
    d=sum(abs(vector1-vector2));
 !
  end function calculate_manhattan_distance
!  
end module manhattan_distance_utilities