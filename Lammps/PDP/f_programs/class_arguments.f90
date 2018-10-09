!>class for arguments
!>created by Hiroaki Yoshida
!>2012/06/20
!>last updated 2012/03/26
 module class_arguments
!
   implicit none
!
  integer,parameter :: ArgLen=128 !Max length of an argument
  integer,parameter ::  nargMax=128 !Max number of arguments
  type ArgProp
     integer :: narg !<number of arguments
     character(len=ArgLen) :: varg(nargMax) !value of arguments
  end type ArgProp
!
 contains
!-----------------------------------------------------------------------------
!>Initializing
 subroutine argInit(Prop)
!
    implicit none
!
    type(ArgProp),intent(inout) :: Prop !structure of type ArgProp
    integer :: i
!
    Prop%narg=command_argument_count()
    do i=1,Prop%narg
     call get_command_argument(i,Prop%varg(i))
    end do
!
   return
 end subroutine
!
!-----------------------------------------------------------------------------
!>Display
 subroutine argDisplay(Prop)
!
    implicit none
!
    type(ArgProp),intent(inout) :: Prop !structure of type ArgProp
    integer :: i
!
    write(*,*)"---Arguments are---------"
    do i=1,Prop%narg
     write(*,*) trim(Prop%varg(i))
    end do
    write(*,*)"---End of Arguments------"
!
   return
 end subroutine
!
!-----------------------------------------------------------------------------
!>exist?
 subroutine argExist(opt,exi,Prop,eini)
!
    implicit none
! 
    logical,intent(inout) :: exi
    type(ArgProp),intent(in) :: Prop !structure of type ArgProp
    character(len=*),intent(in) :: opt
    logical,intent(in),optional :: eini
!
    integer :: i
!
    if(present(eini))then
     exi=eini
    end if
!
    if(Prop%narg>0)then
     do i=1,Prop%narg
      if(Prop%varg(i)==opt)then
       exi=.true.
      end if
     end do
    end if
!
   return
 end subroutine
!
!-----------------------------------------------------------------------------
!>read integer value
 subroutine argGetInt(opt,n,Prop,nini)
!
    implicit none
! 
    type(ArgProp),intent(in) :: Prop !structure of type ArgProp
    character(len=*),intent(in) :: opt
    integer,intent(inout) :: n
    integer,intent(in),optional :: nini
!
    integer :: i
!
    if(present(nini))then
     n=nini
    end if
!
    if(Prop%narg>0)then
     do i=1,Prop%narg
      if(Prop%varg(i)==opt)then
       if(i<Prop%narg)then
        read(Prop%varg(i+1),*)n
       end if
      end if
     end do
    end if
!
   return
 end subroutine
!
!-----------------------------------------------------------------------------
!>read real value
 subroutine argGetReal(opt,r,Prop,rini)
!
    implicit none
! 
    type(ArgProp),intent(in) :: Prop !structure of type ArgProp
    character(len=*),intent(in) :: opt
    real(8),intent(inout) :: r
    real(8),intent(in),optional :: rini
!
    integer :: i
!
    if(present(rini))then
     r=rini
    end if
!
    if(Prop%narg>0)then
     do i=1,Prop%narg
      if(Prop%varg(i)==opt)then
       if(i<Prop%narg)then
        read(Prop%varg(i+1),*)r
       end if
      end if
     end do
    end if
!
   return
 end subroutine
!
!-----------------------------------------------------------------------------
!>read character
 subroutine argGetChar(opt,c,Prop,cini)
! 
    implicit none
! 
    type(ArgProp),intent(in) :: Prop !structure of type ArgProp
    character(len=*),intent(in) :: opt
    character(len=*),intent(inout) :: c
    character(len=*),intent(in),optional :: cini
!
    integer :: i
!
    if(present(cini))then
     c=cini
    end if
!
    if(Prop%narg>0)then
     do i=1,Prop%narg
      if(Prop%varg(i)==opt)then
       if(i<Prop%narg)then
        read(Prop%varg(i+1),'(a)')c
       end if
      end if
     end do
    end if
!
   return
 end subroutine
!
!-----------------------------------------------------------------------------
 end module
