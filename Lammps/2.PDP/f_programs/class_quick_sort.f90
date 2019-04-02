!--------------------------------------------
!>\author Hiroaki Yoshida
!>\details
!> last updated 2008 09 18 \n
!> last updated 2009 10 16 \n
!> last updated 2011 01 21 \n
!--------------------------------------------
 module class_quick_sort
!
! Declarations of structure variables
!-----------------------------------------------------------------------------
!>structure for quick-sort subroutines
!>\details
!>all members are PRIVATE
 type quickSortProp
  private
  real(8),allocatable :: x(:) !<original data
  real(8),allocatable :: y(:) !<sorted data
  integer,allocatable :: id(:) !<id of original data
  integer :: np !<number of data
 end type quickSortProp
!
 private
  real(8),allocatable :: a(:) !<array for sort
  integer,allocatable :: b(:) !<array for index renumbering
 public::quickSortProp,quick_sort,quick_sort_id,quick_sort_des
!
 contains
!------------------------------------------------------------
!> (constructor) sort array x(n) in increasing order and put into y(n)
!------------------------------------------------------------
 subroutine quick_sort(x,y,n,qProp)
  implicit none
!
  real(8),intent(out) :: y(n) !<sorted data
  type(quickSortProp),intent(out) :: qProp !<structure of type 'quickSortProp'
  integer,intent(in) :: n !<number of data
  real(8),intent(in) :: x(n) !<original data
  integer :: i
!
  allocate(a(n),b(n))
!
  do i=1,n
   b(i)=i
  end do
!
  a(:)=x(:)
!
  qProp%np=n
  call qs_recursion(1,n)
!
  allocate(qProp%x(1:n),qProp%y(1:n),qProp%id(1:n))
  do i=1,n
   qProp%x(i)=x(i)
   qProp%y(i)=a(i)
   y(i)=a(i)
   qProp%id(i)=b(i)
  end do
!
  deallocate(a,b)
  return
 end subroutine
!
!
!--------------------------------------------------------------
!> return the original id number of is-th value in sorted order
 function quick_sort_id(is,qProp)
  integer quick_sort_id !<output id number
  integer,intent(in) :: is !<id in the sorted order
  type(quickSortProp),intent(in) :: qProp !<structure of type 'quickSortProp'
!
   quick_sort_id=qProp%id(is)
  return
 end function
!
!--------------------------------------------------------------
!> destructor for quickSortProp structure
 subroutine quick_sort_des(qProp)
  type(quickSortProp),intent(in) :: qProp !<structure of type 'quickSortProp'
!
!   deallocate(qProp%x,qProp%y,qProp%id)
  return
 end subroutine

!--------------------------------------------------------------
!> <<<PRIVATE>>>recursion quick sort algorithm
  recursive subroutine qs_recursion(i,j)
   integer,intent(in) :: i !first index
   integer,intent(in) :: j !last index
   integer :: k,th
   integer :: ii
   real(8) :: dth
!
   if(j<=i)return
   th=thres(i,j)
   if(th/=-1)then
    dth=a(th)
    k=partition(i,j,dth)
    call qs_recursion(i,k-1)
    call qs_recursion(k,j)
   end if
   return
  end subroutine
!
!---------------------------------------
!> <<<PRIVATE>>> get index of larger value from first two different values
  function thres(i,j)
   integer :: thres !<threshold
   integer,intent(in) :: i !<first index
   integer,intent(in) :: j !<last index
   integer :: k
   k=i+1
    do
     if(a(i)/=a(k))exit
     k=k+1
     if(k>j)exit
    end do 
    if(k>j)then
     thres=-1
    else if(a(i)>=a(k))then
     thres=i
    else 
     thres=k
    end if
   return
  end function
!
!---------------------------------------
!>  <<<PRIVATE>>> 
!> \details
!> Find the value larger than x in increasing order from i
!> and also find the value smaller than x in decreasing order from j.
!> Search is terminated when two searching points meet.
!> Return the index where two searches meet.
  function partition(i,j,x)
   integer :: partition
   integer,intent(in) :: i!<first index
   integer,intent(in) :: j!<last index
   real(8),intent(in) :: x!<threshold
   real(8) :: t,u
   integer :: l,r
!
   l=i
   r=j
   do
!   !find value larger than leading one
    do
     if(a(l)>=x)exit
     l=l+1
     if(l>j)exit
    end do
!   !find value smaller than leading one
    do
     if(a(r)<x)exit
     r=r-1
     if(r<i)exit
    end do
!
    if(l>r)exit
    t=a(l)
    a(l)=a(r)
    a(r)=t
!
    u=b(l)
    b(l)=b(r)
    b(r)=u
!
    l=l+1
    r=r-1
   end do
 !
   partition=l
  end function
end module
