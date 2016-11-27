module dtw_distance_utilities
!
use distance_base_utilities;
!
type,extends(distance_base) :: dtw_distance 
 contains
   procedure,public :: calculate => calculate_dtw_distance
   procedure,private,nopass :: calculate_warping_function
end type dtw_distance
!
 contains
!=============================================================
 function calculate_dtw_distance(distance,vector1,vector2) result(d)
!=============================================================
   class(dtw_distance) :: distance
   real(kind=8),dimension(:,:),intent(inout) :: vector1,vector2
   real(kind=8) :: d
!
   real(kind=8),dimension(0:size(vector1,1),0:size(vector2,1)) ::g
   real(kind=8),dimension(size(vector1,1),3) :: warping
!   integer :: i
!
   call calculate_warping_function(vector1,vector2,g,warping);
!    do i=1,size(warping,1)
!       write(*,*) warping(i,:);
!    enddo
   d=warping(size(vector1,1),3) !/dble(size(vector1,1))
!
 end function calculate_dtw_distance
!=============================================================
 subroutine calculate_warping_function(vector1,vector2,g,warping)
!=============================================================
!   class(dtw_distance) :: distance
   real(kind=8),dimension(:,:),intent(inout) :: vector1,vector2
   real(kind=8),dimension(0:size(vector1,1),0:size(vector2,1)),intent(out) :: g
   real(kind=8),dimension(:,:),intent(out) :: warping
!
   integer :: i,j,p
   integer,dimension(1) :: pos
   real(kind=8) :: distance1,test1,test2,test3
   real(kind=8),dimension(3) :: test_val
!
   g=1.0d99;
   g(0,0)=0.0d0;
   do j=1,size(vector2,1);
      do i=1,size(vector1,1);
         distance1=dabs(vector1(i,1)-vector2(j,1));
         test_val(1)=g(i,j-1)+distance1;
         test_val(2)=g(i-1,j-1)+distance1;
         test_val(3)=g(i-1,j)+distance1;
         g(i,j)=minval(test_val);         
      enddo
   enddo
!    open(1,file='warping1.out',status='unknown');
!    do i=1,size(vector1,1)
!       write(1,*) (g(i,j),j=1,size(vector2,1));
!    enddo
!    close(1)
!  write(*,*) g(size(vector1,1),size(vector2,1))
!
   warping(size(vector1,1),1)=dble(size(vector1,1));
   warping(size(vector1,1),2)=dble(size(vector2,1));
   warping(size(vector1,1),3)=g(size(vector1,1),size(vector2,1));
   i=size(vector1,1);j=size(vector2,1);
   do p=size(vector1,1)-1,1,-1
      test_val(1)=g(i,j-1);
      test_val(2)=g(i-1,j-1);
      test_val(3)=g(i-1,j);
      warping(p,3)=minval(test_val);
      pos=minloc(test_val);
      select case(pos(1));
        case(1)
         i=i;j=j-1;
        case(2)
         i=i-1;j=j-1;
        case(3)
        i=i-1;j=j;
      end select
      warping(p,1)=dble(i);warping(p,2)=dble(j);
   enddo   
!
 end subroutine calculate_warping_function 
 
end module dtw_distance_utilities