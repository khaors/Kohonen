module wavelet_distance_utilities
!
use distance_base_utilities;
!
implicit none
!
type,extends(distance_base) :: wavelet_distance 
 contains
   procedure,public :: calculate => calculate_wavelet_distance
end type wavelet_distance
!
 contains
! 
 function calculate_wavelet_distance(distance,vector1,vector2) result(d)
! 
   class(wavelet_distance) :: distance
   real(kind=8),dimension(:,:),intent(inout) :: vector1,vector2
   real(kind=8) :: d
!
   real(kind=8),dimension(size(vector1,2),size(vector1,2)) :: cov,u,s,v
   real(kind=8),dimension(size(vector1,2)) :: w
   integer :: nrow,ncol,i,j,num_k
   real(kind=8) :: total_variance,dd
   real(kind=8),dimension(size(vector1,2)) :: eigenvalues
   real(kind=8),dimension(size(vector1,1),6) :: L1,L2
   real(kind=8),dimension(size(vector1,1)-1,6) :: SL
   real(kind=8),dimension(size(vector1,2)-1,6) :: SU
   real(kind=8),dimension(6) :: weights,dl,du
   real(kind=8),dimension(1) :: d1
!
   nrow=size(vector1,1);
   ncol=size(vector1,2);
   cov=matmul(transpose(vector1),vector2);
   u=cov;
   call svdcmp(u,ncol,ncol,w,v);
!   call r8mat_svd_lapack(nrow, ncol, cov, u, w, v )
   eigenvalues=w;
   total_variance=sum(w);
! !
   num_k=6;
   L1=matmul(vector1,u(:,1:num_k));
   L2=matmul(vector2,v(:,1:num_k));
! !
! !   W1N=matmul(L1,transpose(v(:,1:num_k)));
! !   W1N=matmul(L2,transpose(u(:,1:num_k)));
! !   
  SL=0.0d0;
  SU=0.0d0;
  do j=1,nrow-1;
    do i=1,num_k;
        SL(j,i)=SL(j,i)+atan(abs((L1(j,i)-L2(j,i))-(L1(j+1,i)-L2(j+1,i))));        
    enddo;
  enddo;
! 
  do j=1,ncol-1;
     do i=1,num_k;
        SU(j,i)=SU(j,i)+atan(abs((u(j,i)-v(j,i))-(u(j+1,i)-v(j+1,i))));
    enddo;
  enddo;
! !
   weights=eigenvalues(1:num_k)/sum(eigenvalues(1:num_k))
! !
  dl=sum(SL,1);
  du=sum(SU,1);
  d=dot_product(weights,(dl+du));
!
!    d1=sum(sum((vector1-vector2)**2,2),1)
!    d=d1(1)
 end function calculate_wavelet_distance
!*****************************************************************************80
subroutine r8mat_svd_lapack ( m, n, a, u, s, v )
!*****************************************************************************80
!
!! R8MAT_SVD_LAPACK gets the SVD of a matrix using a call to LAPACK.
!
!  Discussion:
!
!    The singular value decomposition of a real MxN matrix A has the form:
!
!      A = U * S * V'
!
!    where
!
!      U is MxM orthogonal,
!      S is MxN, and entirely zero except for the diagonal;
!      V is NxN orthogonal.
!
!    Moreover, the nonzero entries of S are positive, and appear
!    in order, from largest magnitude to smallest.
!
!    This routine calls the LAPACK routine DGESVD to compute the
!    factorization.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    08 February 2007
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, integer ( kind = 4 ) M, N, the number of rows and columns
!    in the matrix A.
!
!    Input, real ( kind = 8 ) A(M,N), the matrix whose singular value
!    decomposition we are investigating.
!
!    Output, real ( kind = 8 ) U(M,M), S(M,N), V(N,N), the factors
!    that form the singular value decomposition of A.
!
  implicit none

  integer ( kind = 4 ) m
  integer ( kind = 4 ) n

  real ( kind = 8 ) a(m,n)
  real ( kind = 8 ) a_copy(m,n)
  integer ( kind = 4 ) i
  integer ( kind = 4 ) info
  integer ( kind = 4 ) lda
  integer ( kind = 4 ) ldu
  integer ( kind = 4 ) ldv
  character jobu
  character jobv
  integer ( kind = 4 ) lwork
  real ( kind = 8 ) sdiag(min(m,n))
  real ( kind = 8 ) s(m,n)
  real ( kind = 8 ) u(m,m)
  real ( kind = 8 ) v(n,n)
  real ( kind = 8 ), allocatable, dimension ( : ) :: work
!
  external dgesvd
!
  lwork = max ( 3 * min ( m, n ) + max ( m, n ), 5 * min ( m, n ) )

  allocate ( work(1:lwork) )
!
!  Compute the eigenvalues and eigenvectors.
!
  jobu = 'A'
  jobv = 'A'
  lda = m
  ldu = m
  ldv = n
!
!  The input matrix is destroyed by the routine.  Since we need to keep
!  it around, we only pass a copy to the routine.
!
  a_copy(1:m,1:n) = a(1:m,1:n)

  call dgesvd ( jobu, jobv, m, n, a_copy, lda, sdiag, u, ldu, v, ldv, work, &
    lwork, info )

  if ( info /= 0 ) then
    write ( *, '(a)' ) ' '
    write ( *, '(a)' ) 'R8MAT_SVD_LAPACK - Failure!'
    write ( *, '(a)' ) '  The SVD could not be calculated.'
    write ( *, '(a)' ) '  LAPACK routine DGESVD returned a nonzero'
    write ( *, '(a,i8)' ) '  value of the error flag, INFO = ', info
    return
  end if
!
!  Make the MxN matrix S from the diagonal values in SDIAG.
!
  s(1:m,1:n) = 0.0D+00
  do i = 1, min ( m, n )
    s(i,i) = sdiag(i)
  end do
!
!  Transpose V.
!
  v = transpose ( v )

  deallocate ( work )

  return
end subroutine r8mat_svd_lapack

SUBROUTINE svdcmp(a,m,n,w,v) 
integer :: m,n,NMAX 
real(kind=8)  :: a(m,n),v(n,n),w(n) 
PARAMETER (NMAX=500)  !Maximum anticipated value of n. 
!-------------------------------------------------------------------------------------- 
! Given a matrix A(1:m,1:n), this routine computes its singular value decomposition, 
! A = U · W · Vt. The matrix U replaces A on output. The diagonal matrix of singular 
! values W is output as a vector W(1:n). The matrix V (not the transpose Vt) is output 
! as V(1:n,1:n). 
!--------------------------------------------------------------------------------------
integer :: i,its,j,jj,k,l,nm 
real(kind=8)  :: anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX)!,pythag 
  g=0.d0  !Householder reduction to bidiagonal form. 
  scale=0.d0 
  anorm=0.d0 
do i=1,n 
  l=i+1 
  rv1(i)=scale*g 
  g=0.d0 
  s=0.d0 
  scale=0.d0 
  if(i.le.m)then 
  do k=i,m 
    scale=scale+abs(a(k,i)) 
  end do 
  if(scale.ne.0.d0)then 
  do k=i,m 
    a(k,i)=a(k,i)/scale 
    s=s+a(k,i)*a(k,i) 
  end do 
  f=a(i,i) 
  g=-dsign(dsqrt(s),f) 
  h=f*g-s 
  a(i,i)=f-g 
  do j=l,n 
    s=0.d0 
    do k=i,m 
      s=s+a(k,i)*a(k,j) 
    end do 
    f=s/h 
    do k=i,m 
      a(k,j)=a(k,j)+f*a(k,i) 
    end do
  end do 
  do k=i,m 
    a(k,i)=scale*a(k,i) 
  end do 
  endif 
  endif 
  w(i)=scale *g 
  g=0.d0 
  s=0.d0 
  scale=0.d0 
  if((i.le.m).and.(i.ne.n))then 
  do k=l,n 
    scale=scale+abs(a(i,k)) 
  end do 
  if(scale.ne.0.d0)then 
  do k=l,n 
    a(i,k)=a(i,k)/scale 
    s=s+a(i,k)*a(i,k) 
  end do 
  f=a(i,l) 
  g=-sign(sqrt(s),f) 
  h=f*g-s 
  a(i,l)=f-g 
  do k=l,n 
    rv1(k)=a(i,k)/h 
  end do 
  do j=l,m 
    s=0.d0 
    do k=l,n 
      s=s+a(j,k)*a(i,k) 
    end do 
    do k=l,n 
      a(j,k)=a(j,k)+s*rv1(k) 
    end do 
  end do 
  do k=l,n 
    a(i,k)=scale*a(i,k) 
  end do 
  endif 
  endif 
  anorm=max(anorm,(abs(w(i))+abs(rv1(i)))) 
end do !do i=1,n
 
do i=n,1,-1 !Accumulation of right-hand transformations. 
  if(i.lt.n)then 
  if(g.ne.0.d0)then 
  do j=l,n       !Double division to avoid possible underflow. 
    v(j,i)=(a(i,j)/a(i,l))/g 
  end do 
  do j=l,n 
    s=0.d0 
    do k=l,n 
      s=s+a(i,k)*v(k,j) 
    end do 
    do k=l,n 
      v(k,j)=v(k,j)+s*v(k,i) 
    end do 
  end do 
  endif 
  do j=l,n 
    v(i,j)=0.d0 
    v(j,i)=0.d0 
  end do 
  endif 
  v(i,i)=1.d0
  g=rv1(i) 
  l=i 
end do 

do i=min(m,n),1,-1 !Accumulation of left-hand transformations. 
  l=i+1 
  g=w(i) 
  do j=l,n 
    a(i,j)=0.d0 
  end do 
  if(g.ne.0.d0)then 
  g=1.d0/g 
  do j=l,n 
    s=0.d0 
    do k=l,m 
      s=s+a(k,i)*a(k,j) 
    end do 
    f=(s/a(i,i))*g 
    do k=i,m 
      a(k,j)=a(k,j)+f*a(k,i) 
    end do 
  end do 
  do j=i,m 
    a(j,i)=a(j,i)*g 
  end do 
  else
  do j= i,m 
    a(j,i)=0.d0 
  end do 
  endif 
  a(i,i)=a(i,i)+1.d0 
end do 

do k=n,1,-1 !Diagonalization of the bidiagonal form: Loop over 
            !singular values, and over allowed iterations. 
do its=1,30 
do l=k,1,-1 !Test for splitting. 
  nm=l-1 !Note that rv1(1) is always zero. 
  if((abs(rv1(l))+anorm).eq.anorm) goto 2 
  if((abs(w(nm))+anorm).eq.anorm) goto 1 
end do 
1 c=0.d0 !Cancellation of rv1(l), if l > 1. 
s=1.d0 
do i=l,k 
  f=s*rv1(i) 
  rv1(i)=c*rv1(i) 
  if((abs(f)+anorm).eq.anorm) goto 2 
  g=w(i) 
  h=pythag(f,g) 
  w(i)=h 
  h=1.d0/h 
  c= (g*h) 
  s=-(f*h) 
  do j=1,m 
    y=a(j,nm) 
    z=a(j,i) 
    a(j,nm)=(y*c)+(z*s) 
    a(j,i)=-(y*s)+(z*c) 
  end do 
end do 
2 z=w(k) 
if(l.eq.k)then   !Convergence. 
if(z.lt.0.d0)then !Singular value is made nonnegative. 
w(k)=-z 
do j=1,n 
  v(j,k)=-v(j,k) 
end do
endif 
goto 3 
endif 
if(its.eq.30) stop 'no convergence in svdcmp' 
x=w(l) !Shift from bottom 2-by-2 minor. 
nm=k-1 
y=w(nm) 
g=rv1(nm) 
h=rv1(k) 
f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y) 
g=pythag(f,1.d0) 
f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x 
c=1.d0 !Next QR transformation: 
s=1.d0 
do j=l,nm 
  i=j+1 
  g=rv1(i) 
  y=w(i) 
  h=s*g 
  g=c*g 
  z=pythag(f,h) 
  rv1(j)=z 
  c=f/z 
  s=h/z 
  f= (x*c)+(g*s) 
  g=-(x*s)+(g*c) 
  h=y*s 
  y=y*c 
  do jj=1,n 
    x=v(jj,j) 
    z=v(jj,i) 
    v(jj,j)= (x*c)+(z*s) 
    v(jj,i)=-(x*s)+(z*c) 
  end do 
  z=pythag(f,h) 
  w(j)=z !Rotation can be arbitrary if z = 0. 
  if(z.ne.0.d0)then 
  z=1.d0/z 
  c=f*z 
  s=h*z 
  endif 
  f= (c*g)+(s*y) 
  x=-(s*g)+(c*y) 
  do jj=1,m 
    y=a(jj,j) 
    z=a(jj,i) 
    a(jj,j)= (y*c)+(z*s) 
    a(jj,i)=-(y*s)+(z*c) 
  end do 
end do !j=l;nm 
rv1(l)=0.d0 
rv1(k)=f 
w(k)=x 
end do !its=1,30
3 continue 
end do !k=n,1,-1 
return 
END SUBROUTINE svdcmp
 
FUNCTION pythag(a,b) result(p)
real(kind=8)  :: a,b 
real(kind=8) :: p
!Computes sqrt(a**2 + b**2) without destructive underflow or overflow.
real(kind=8)  :: absa,absb 
  absa=abs(a) 
  absb=abs(b) 
  if(absa.gt.absb)then 
    p=absa*sqrt(1.+(absb/absa)**2) 
  else 
    if(absb.eq.0.)then 
      p=0. 
    else
      p=absb*sqrt(1.+(absa/absb)**2) 
    endif 
  endif 
  return 
END function pythag


end module wavelet_distance_utilities