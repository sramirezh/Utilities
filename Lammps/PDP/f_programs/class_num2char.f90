!-------------------------------------------------------------
!>converting a number into a string character
!>\author Hiroaki Yoshida
!>\details
!>last up date 2011 01 18
!-------------------------------------------------------------
 module class_num2char
!
 implicit none
!
 private::construct_character
!
 contains
!
!-------------------------------------------------------------
!>convert integer into character
!
!>\details
!>if the digit number of n is smaller than i, 0 is put before converted n
!>argument "i" can be omitted (default=2) \n
!> \n
!>-caution- \n
!> series of spaces are put behind the converted number, \n
!> so that total length of int2char is 11 \n
!> please use "trim" function to remove those spaces \n
!> \n
!>-example-
!>trim(int2char(4))="04" \n
!>trim(int2char(4,1))="4" \n
!>trim(int2char(134,5))="00134" \n
!>trim(int2char(-67,3))="-067"
!
!-------------------------------------------------------------
 function int2char(n,i)
  character :: int2char*11
  integer,intent(in) :: n !<integer number to be converted into character
  integer,intent(in),optional :: i !<number of minimum digit
  integer :: nu,j,a,m
  character,parameter :: c*11="0123456789-"
  character :: c1*11
  logical :: fl_minus
!
  fl_minus=.false.
  if(n<0)fl_minus=.true.
  nu=abs(n)
!
  c1=""
  do 
   a=mod(nu,10)+1
   c1=c(a:a)//trim(c1)
   nu=nu/10
   if(nu==0)exit
  end do
!
   m=2
   if (present(i)) then
    if (i>0 .and. i<11) m=i
   end if
!
   do j=1,m-1
    if (len(trim(c1))<m) c1=c(1:1)//trim(c1)
   end do
!
   if (fl_minus) c1=c(11:11)//trim(c1)
!
   int2char=adjustl(c1)
   return
 end function
!
!
!
!-------------------------------------------------------------
!>convert a real number into character
!
!>\details
!>-caution- \n
!>unlike int2char, total length of real2char is 16 \n
!>if the result exceeds this limit, 16 characters from the left are returned \n
!>\n
!>-example- \n
!>trim(dble2char(7.2))="07.20" \n
!>trim(dble2char(7.2,1))="7.20" \n
!>trim(dble2char(728.0505,4,3))="0728.051" \n
!>trim(dble2char(-728.0505,4,3))="-0728.051"
!-------------------------------------------------------------
 function real2char(r,i,d)
  character :: real2char*16 
  real,intent(in) :: r !<real number to be converted into character
  integer,intent(in),optional :: i !<number of minimum digit of integer part (default=2)
  integer,intent(in),optional :: d !<number of digit of decimal part (default=2)
  integer :: mi,md
  integer :: nu_int,nu_dec
  real :: nu
  character :: c1*16
  logical :: fl_minus
!
  mi=2
  if (present(i)) then
   if (i>0 .and. i<11) mi=i
  end if
!
  md=2
  if (present(d)) then
   if (d>=0 .and. d<11) md=d
  end if
!
  fl_minus=.false.
  if(r<0)fl_minus=.true.
  nu=abs(r)
!
  nu=nu+0.1e0**(md+1)
!
  nu_int=int(nu)
!
  nu_dec=nint((nu-real(nu_int))*10e0**md)
!
  if(nu_dec==10**md)then
     nu_int=nu_int+1
     nu_dec=0
  end if
!
  real2char=adjustl(construct_character(nu_int,nu_dec,mi,md,fl_minus))
!
  return
 end function
!
!
!
!-------------------------------------------------------------
!>convert double precision real number into character
!
!>\details
!>-caution- \n
!>unlike int2char, total length of dble2char is 16 \n
!>if the result exceeds this limit, 16 characters from the left are returned \n
!>\n
!>-example- \n
!>trim(dble2char(7.2))="07.20" \n
!>trim(dble2char(7.2,1))="7.20" \n
!>trim(dble2char(728.0505,4,3))="0728.051" \n
!>trim(dble2char(-728.0505,4,3))="-0728.051"
!-------------------------------------------------------------
 function dble2char(r,i,d)
  character :: dble2char*16
  real(8),intent(in) :: r !<double number to be converted into character
  integer,intent(in),optional :: i !<number of minimum digit of integer part (default=2)
  integer,intent(in),optional :: d !<number of digit of decimal part (default=2)
  integer :: mi,md
  integer :: nu_int,nu_dec
  real(8) :: nu
  character :: c1*16
  logical :: fl_minus
!
  mi=2
  if (present(i)) then
   if (i>0 .and. i<11) mi=i
  end if
!
  md=2
  if (present(d)) then
   if (d>=0 .and. d<11) md=d
  end if
!
  fl_minus=.false.
  if(r<0)fl_minus=.true.
  nu=dabs(r)
!
  nu=nu+0.1d0**(md+1)
!
  nu_int=int(nu)
!
  nu_dec=nint((nu-dble(nu_int))*10d0**md)
!
  if(nu_dec==10**md)then
     nu_int=nu_int+1
     nu_dec=0
  end if
!
  dble2char=adjustl(construct_character(nu_int,nu_dec,mi,md,fl_minus))
!
  return
 end function
!

!-------------------------------------------------------------
!> <<<PRIVATE>>> construct character
!-------------------------------------------------------------
 function construct_character(in,de,i,d,fl)
  character :: construct_character*16
  integer,intent(in) :: in !<integer for integer part
  integer,intent(in) :: de !<integer for decimal part
  integer,intent(in) :: i !<number of minimum digit of integer part (default=2)
  integer,intent(in) :: d !<number of digit of decimal part (default=2)
  logical,intent(in) :: fl !<flag for sign
  character :: c*16

  if(d>0)then
   c=trim(int2char(in,i))//"."//trim(int2char(de,d))
  else
   c=trim(int2char(in,i))
  end if   

  if (fl) c="-"//trim(c)

  construct_character=c

  return
 end function

 end module
