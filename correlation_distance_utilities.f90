module correlation_distance_utilities

use distance_base_utilities;

implicit none;

type,extends(distance_base) :: correlation_distance 
 contains
   procedure,public :: calculate => calculate_correlation_distance
end type correlation_distance

 contains
! 
 function calculate_correlation_distance(distance,vector1,vector2) result(d)
! 
   class(correlation_distance) :: distance
   real(kind=8),dimension(:,:),intent(inout) :: vector1,vector2
   real(kind=8) :: d
!
   real(kind=8),dimension(size(vector1,1)*size(vector1,2)) :: v1,v2
   real(kind=8) :: correlation,m1,m2,s1,s2
   integer :: n11,n21
!
   n11=size(vector1,1);
   n21=size(vector1,2);
   v1=reshape(vector1,(/n11*n21/))
   v2=reshape(vector2,(/n11*n21/))
   m1=sum(v1)/float(n11*n21);
   m1=sum(v2)/float(n11*n21);
   s1=sum((v1-m1)**2);
   s2=sum((v2-m2)**2);
   correlation=sum((v1-m1)*(v2-m2))/(s1*s2);
   d=sqrt((2.0*(1.0d0-correlation)));
!
 end function calculate_correlation_distance
 

end module correlation_distance_utilities