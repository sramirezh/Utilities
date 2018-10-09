!========================================
!> Class vector and matrix computations
!> \author Hiroaki Yoshida 
!> \details
!> last updated 2009 09 15 \n
!> last updated 2010 09 16 \n
!> last updated 2011 01 18 encapsulation \n
!========================================
module class_vecmat
!
!-----------------------------------------------------------------------------
! Declarations of structure variables
!-----------------------------------------------------------------------------
!
 integer,parameter :: Vecmat_nmax=16 !< max dimension
!
!> structure for subroutines named lu_***
 type luprop 
    private
    real(8) :: lu(Vecmat_nmax,Vecmat_nmax)   !< n*n array: lu decomposition of a
    integer :: id(Vecmat_nmax)               !< n array: index for pivot selection
    integer :: nn                            !< dimension of matrix
    integer :: err                           !< error flag (return 1 if error occurs, 0 otherwise)
 end type luprop

contains
!
!========================================
!> take inner product of x and y
 subroutine vm_inner(c,x,y,nn)

  integer,intent(in) :: nn    !< dimension
  real(8),intent(out) :: c    !< scalor output
  real(8),intent(in) :: x(nn) !< input vector
  real(8),intent(in) :: y(nn) !< input vector


  integer :: i

  c=0d0
  do i=1,nn
   c=c+x(i)*y(i)
  end do

  return
 end subroutine
!
!========================================
!> take product b=ax
 subroutine vm_product(y,a,x,nn)

  integer,intent(in) :: nn !< dimension
  real(8),intent(out) :: y(nn) !< output vector
  real(8),intent(in) :: a(nn,nn) !< input matrix
  real(8),intent(in) :: x(nn) !< output vector


  integer :: i,j

  do i=1,nn
   y(i)=0d0
   do j=1,nn
    y(i)=y(i)+a(j,i)*x(j)
   end do
  end do
!
  return
 end subroutine
!
!========================================
!> take product z=x times y
 subroutine vm_cross(z,x,y)
!
  real(8),intent(out) :: z(3)!< 3d array, output vector
  real(8),intent(in) :: x(3) !< 3d array, input vector
  real(8),intent(in) :: y(3) !< 3d array, input vector
!
  z(1)=x(2)*y(3)-x(3)*y(2)
  z(2)=x(3)*y(1)-x(1)*y(3)
  z(3)=x(1)*y(2)-x(2)*y(1)
!
  return
 end subroutine
!
!========================================
!> compute lu decomposition of matrix a and put into lu \n
!> pivot is taken into account
 subroutine lu_decomposition(a,nn,lup)
!
  implicit none
!
  integer,intent(in) :: nn !< dimension
  real(8),intent(in) :: a(nn,nn) !< nn*nn array, input matrix: a(j,i) j=column, i=row
!
  type(luprop),intent(out) :: lup !< structure of type luprop, output
!
  real(8) :: lu(Vecmat_nmax,Vecmat_nmax)
  integer :: id(Vecmat_nmax)
!
  integer :: i,j !coution i:row, j:column
  integer :: k
  real(8) :: weight(Vecmat_nmax)
  real(8) :: u,tmp,dummy,sum
  integer :: imax,idummy
!
!
  lup%err=0
!
!
  do i=1,nn
   u=0d0
   do j=1,nn
    tmp=dabs(a(j,i))
    if(tmp>u)then
     u=tmp
    end if
   end do
!
   if(u==0d0)then
    lup%lu(:,:)=0d0
    lup%err=1
    return
   else
    weight(i)=1d0/u
   end if
!
  end do
!
!
  do i=1,nn
   do j=1,nn
    lu(i,j)=a(i,j)
   end do
  end do
!
!
  do i=1,nn
   id(i)=i
  end do
!
  do j=1,nn
   do i=1,j-1
    sum=lu(j,i)
    do k=1,i-1
     sum=sum-lu(k,i)*lu(j,k)
    end do
    lu(j,i)=sum
   end do
   u=0d0
   do i=j,nn
    sum=lu(j,i)
    do k=1,j-1
     sum=sum-lu(k,i)*lu(j,k)
    end do
    lu(j,i)=sum
!
    dummy=weight(i)*dabs(sum)
    if(dummy>=u)then
     u=dummy
     imax=i
    end if
   end do
!
   if(j /= imax)then
    do k=1,nn
     dummy=lu(k,imax)
     lu(k,imax)=lu(k,j)
     lu(k,j)=dummy
     idummy=id(j)
     id(j)=id(imax)
     id(imax)=idummy
    end do
    weight(imax)=weight(j)
   end if
!
   if(lu(j,j)==0d0)lu(j,j)=1d-16
   if(j<nn)then
    dummy=1d0/lu(j,j)
    do i=j+1,nn
     lu(j,i)=lu(j,i)*dummy
    end do
   end if
!
  end do !index j
!
  lup%id(:)=id(:)
  lup%lu(:,:)=lu(:,:)
  lup%nn=nn
!
  return
!
 end subroutine
!
!========================================
!> compute determinant of a()
 subroutine lu_det(det,lup)
!
  implicit none
!
  real(8),intent(out) :: det !< determinant, output
  type(luprop),intent(in) :: lup!< structure of type luprop, input
!
  real(8) :: lu(Vecmat_nmax,Vecmat_nmax)
  integer :: nn
!
  integer :: idd(Vecmat_nmax)
  integer :: i,j,sign
!
!
!-----count number of exchanges
  nn=lup%nn
  idd(:)=lup%id(:)
  lu(:,:)=lup%lu(:,:)
!
  sign=1
  do i=1,nn
   if(i /= idd(i))then
    do j=i+1,nn
     if(idd(j)==i)then
      idd(j)=idd(i)
      sign=-sign
      exit
     end if
    end do
   end if
  end do
!
!-----compute determinant
  det=dble(sign)
  do i=1,nn
   det=det*lu(i,i)
  end do
!
  return
!
  end subroutine
!
!========================================
!> compute inverse matrix of a()
 subroutine lu_inverse(inv,lup)
!
  implicit none
!
  type(luprop),intent(in) :: lup !< structure of type luprop, input
  real(8),intent(out) :: inv(lup%nn,lup%nn) !< inverse matrix, output
!
  real(8) :: lu(Vecmat_nmax,Vecmat_nmax)
  integer :: id(Vecmat_nmax)
  integer :: nn
!
  integer :: i,j !caution i:row, j:column
  integer :: k
  real(8) :: sum
  integer :: pivot
!
  nn=lup%nn
  lu(:,:)=lup%lu(:,:)
  id(:)=lup%id(:)
!
  inv(:,:)=0d0
!
  do k=1,nn
   do i=1,nn
    pivot=id(i)
    if(pivot==k)then
     sum=1d0
    else
     sum=0d0
    end if
!
    do j=1,i-1
     sum=sum-lu(j,i)*inv(k,j)
    end do
    inv(k,i)=sum
   end do
!
   do i=nn,1,-1
    pivot=id(i)
    sum=inv(k,i)
    do j=i+1,nn
     sum=sum-lu(j,i)*inv(k,j)
    end do
    inv(k,i)=sum/lu(i,i)
   end do
  end do
!
  return
!
  end subroutine
!
!========================================
!> solve ax=b with respect to x
 subroutine lu_solver(x,b,lup)
!
  implicit none
!
  type(luprop),intent(in) :: lup !< structure of type luprop, input
  real(8),intent(out) :: x(lup%nn) !< solution x, output
  real(8),intent(in) :: b(lup%nn)  !< vector b
  real(8) :: lu(Vecmat_nmax,Vecmat_nmax)
  integer :: id(Vecmat_nmax)
  integer :: nn
!
  integer :: i,j !caution i:row, j:column
  integer :: k
  real(8) :: sum
  integer :: pivot
!
!
  nn=lup%nn
  lu(:,:)=lup%lu(:,:)
  id(:)=lup%id(:)
!
  x(:)=0d0
  do i=1,nn
   pivot=id(i)
   sum=b(pivot)
!
   do j=1,i-1
    sum=sum-lu(j,i)*x(j)
   end do
   x(i)=sum
  end do
!
  do i=nn,1,-1
   sum=x(i)
   do j=i+1,nn
    sum=sum-lu(j,i)*x(j)
   end do
   x(i)=sum/lu(i,i)
  end do
!
  return
!
  end subroutine
!
!
!
end module
