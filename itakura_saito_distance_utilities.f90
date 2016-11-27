module itakura_saito_distance_utilities
!
use distance_base_utilities;
!
implicit none
!
real(kind=8),parameter :: pi=4.D0*DATAN(1.D0)
!
type,extends(distance_base) :: itakura_saito_distance 
 contains
   procedure,public :: calculate => calculate_itakura_distance
   procedure,private,nopass :: calculate_itakura_divergence
end type itakura_saito_distance

 contains 
!=============================================================
 function calculate_itakura_distance(distance,vector1,vector2) result(d)
!=============================================================
   class(itakura_saito_distance) :: distance
   real(kind=8),dimension(:,:),intent(inout) :: vector1,vector2
   real(kind=8) :: d
!
   d=0.5d0*(calculate_itakura_divergence(vector1,vector2)+&
            calculate_itakura_divergence(vector2,vector1));
!
 end function calculate_itakura_distance
!=============================================================
 function calculate_itakura_divergence(vector1,vector2) result(d)
!=============================================================
!   class(itakura_saito_distance) :: distance
   real(kind=8),dimension(:,:),intent(inout) :: vector1,vector2
   real(kind=8) :: d
!
   integer :: nrow,ncol,i
   real(kind=8),dimension(size(vector1,1)) :: term1,term2,term3,freq,v1,v2
!
   nrow=size(vector1,1);
   ncol=size(vector1,2);
   v1(1:nrow)=vector1(1:nrow,1);
   v2(1:nrow)=vector2(1:nrow,1);
!
   do i=1,nrow
      freq(i)=dble(i)*(pi/dble(nrow));
!      write(*,*) freq(i)
   enddo
!
   where(abs(v1) .lt. 1.0e-5)
     v1=1.0e-5
   end where
!
   where(abs(v2) .lt. 1.0e-5)
     v2=1.0e-5
   end where

!
   if(nrow .eq. 1 .or. ncol .eq. 1) then
     term1=v1/v2;
     term2=log(term1);
     term3=term1-term2-1.0d0;
     d=trapezoidal_real(freq,term3);
     d=d/pi;
   else
     write(6,*) 'ERROR: 1D vector (spectra) is required in this case'
      stop;
   endif
!
 end function calculate_itakura_divergence
!=============================================================
  function trapezoidal_real(t,f) result(i)
!=============================================================
  real(kind=8) :: i
  real(kind=8),dimension(:),intent(inout) :: t
  real(kind=8),dimension(:),intent(inout) :: f
!
  integer :: ndat,ncol,nrow
  real(kind=8) :: h
!
  ndat=size(t,1);nrow=size(f,1);
  if(ndat .ne. nrow) stop 'error in data matrix'
  h=(t(ndat)-t(1))/dble(ndat);
  i=0.0d0;
  i=0.5d0*(f(1)+2.0d0*sum(f(2:ndat-1))+f(ndat));
  i=h*i;
!
  end function trapezoidal_real  

end module itakura_saito_distance_utilities