!>class for reading files
!>created by Hiroaki Yoshida
!>2013/07/06 in Grenoble
 module class_readfiles
!
   implicit none
!
   integer,parameter :: ReadFilesMaxRow=1000000 !maximum number of rows
   integer,parameter :: ReadFilesMaxAtom=256000 !maximum number of atoms
   integer,parameter :: ReadFilesMaxMol=128000 !maximum number of molecules
   integer,parameter :: ReadFilesMaxAtomType=128 !maximum number of atom types
   type ReadFilesProp
      integer ::ThStep(ReadFilesMaxRow)
      real(8) ::ThVal(20,ReadFilesMaxRow)
      integer :: Th1Step
      real(8) ::Th1Val(20)
      integer ::Nrow
!
      integer ::atomID(ReadFilesMaxAtom)
      integer ::atomType(ReadFilesMaxAtom)
      integer ::molID(ReadFilesMaxAtom)
      real(8) ::mass(ReadFilesMaxAtom),charge(ReadFilesMaxAtom)
      real(8) ::rx(ReadFilesMaxAtom),ry(ReadFilesMaxAtom),rz(ReadFilesMaxAtom)
      real(8) ::molx(ReadFilesMaxAtom),moly(ReadFilesMaxAtom),molz(ReadFilesMaxAtom)
      real(8) ::vx(ReadFilesMaxAtom),vy(ReadFilesMaxAtom),vz(ReadFilesMaxAtom)
      real(8) ::rx0(ReadFilesMaxAtom),ry0(ReadFilesMaxAtom),rz0(ReadFilesMaxAtom)
      real(8) ::dmx(ReadFilesMaxAtom),dmy(ReadFilesMaxAtom),dmz(ReadFilesMaxAtom)
!
      integer ::molType(ReadFilesMaxMol),atomPerMol(ReadFilesMaxMol)
      real(8) ::cmass(ReadFilesMaxMol)
      real(8) ::crx(ReadFilesMaxMol),cry(ReadFilesMaxMol),crz(ReadFilesMaxMol)
      real(8) ::cvx(ReadFilesMaxMol),cvy(ReadFilesMaxMol),cvz(ReadFilesMaxMol)
!
      real(8) ::massflux(3)
      real(8) ::chargeflux(3)
      real(8) ::dipole_moment(3)
!
      integer :: Np
      integer :: Nptype(ReadFilesMaxAtomType)
      integer :: Mp
      real(8) :: xhigh,xlow,yhigh,ylow,zhigh,zlow
      real(8) :: lx,ly,lz
   end type ReadFilesProp
!
 contains
!
! ****************************************************************** 
! *** This subroutine open data file ***
 subroutine read_line_ini(cnfile,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128 !<file name to be read
 integer, intent(in) :: iunit !<unit number
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
! 
 integer :: i
 real*8 :: dum
 character :: dummy
!
! ******************************************************************
!
 open(unit=iunit,file=cnfile,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
!
 return
 end subroutine
!
!
! ****************************************************************** 
! *** This subroutine read one line data ***
 subroutine read_line(ieof,cnfile,iunit,Prop)
!
 implicit none
!
 integer,intent(out) :: ieof
 character, intent(in) :: cnfile*128 !<file name to be read
 integer, intent(in) :: iunit !<unit number
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
 integer :: idum
 character :: line*2048
! 
!
! ******************************************************************
!
  read(iunit,"(a)",iostat=ieof) line
  if(ieof<0)line=""
  read(line,*,iostat=idum) Prop%Th1Step, Prop%Th1Val(1), Prop%Th1Val(2), &
        Prop%Th1Val(3), Prop%Th1Val(4), Prop%Th1Val(5), Prop%Th1Val(6), &
        Prop%Th1Val(7), Prop%Th1Val(8), Prop%Th1Val(9), Prop%Th1Val(10), &
        Prop%Th1Val(11), Prop%Th1Val(12), Prop%Th1Val(13), Prop%Th1Val(14), &
        Prop%Th1Val(15), Prop%Th1Val(16), Prop%Th1Val(17), Prop%Th1Val(18), &
        Prop%Th1Val(19), Prop%Th1Val(20)
!
 return
 end subroutine
!
! ****************************************************************** 
! *** This subroutine read one line data ***
 subroutine read_line_close(iunit,Prop)
!
 implicit none
!
 integer, intent(in) :: iunit !<unit number
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
!
  close(iunit)
!
 return
 end subroutine
!
! ****************************************************************** 
! *** This subroutine reads thermo data ***
 subroutine read_thermo(cnfile,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128 !<file name to be read
 integer, intent(in) :: iunit !<unit number
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
! 
 integer :: i
 integer :: ieof
 real*8 :: dum
 character :: dummy,line*2048
!
! ******************************************************************
!
 open(unit=iunit,file=cnfile,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  Prop%nrow=0
  do i = 1, ReadFilesMaxRow
   read(iunit,"(a)",iostat=ieof) line
   if(ieof<0)exit
   read(line,*,iostat=ieof) Prop%ThStep(i), Prop%ThVal(1,i), Prop%ThVal(2,i), &
        Prop%ThVal(3,i), Prop%ThVal(4,i), Prop%ThVal(5,i), Prop%ThVal(6,i), &
        Prop%ThVal(7,i), Prop%ThVal(8,i), Prop%ThVal(9,i), Prop%ThVal(10,i)
   Prop%nrow=Prop%nrow+1
  end do
 close(iunit)
!
 return
 end subroutine
!
!-----------------------------------------------------------------
! calculate center of mass of molecules
 subroutine mass_center(Prop)
!
 implicit none
!
  type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
  integer i,mid
  real(8) :: x,y,z,cxtem,cytem,cztem
!
! ******************************************************************
!
 Prop%crx(:)=0d0
 Prop%cry(:)=0d0
 Prop%crz(:)=0d0
 Prop%cvx(:)=0d0
 Prop%cvy(:)=0d0
 Prop%cvz(:)=0d0
 Prop%cmass(:)=0d0
 Prop%mp=1
 Prop%molType(:)=0
 Prop%AtomPerMol(:)=0
!
 do i=1,Prop%Np
  mid=Prop%molID(i)
  x=Prop%rx(i)
  y=Prop%ry(i)
  z=Prop%rz(i)
  Prop%AtomPerMol(mid)=Prop%AtomPerMol(mid)+1
!  if(Prop%atomType(i)==1)then
  Prop%crx(mid)=Prop%crx(mid)+Prop%mass(i)*x
  Prop%cry(mid)=Prop%cry(mid)+Prop%mass(i)*y
  Prop%crz(mid)=Prop%crz(mid)+Prop%mass(i)*z
!  Prop%cvx(mid)=Prop%cvx(mid)+Prop%mass(i)*Prop%vx(i)
!  Prop%cvy(mid)=Prop%cvy(mid)+Prop%mass(i)*Prop%vy(i)
!  Prop%cvz(mid)=Prop%cvz(mid)+Prop%mass(i)*Prop%vz(i)
  Prop%cmass(mid)=Prop%cmass(mid)+Prop%mass(i)
!  end if
  if(Prop%AtomType(i)>Prop%molType(mid))Prop%molType(mid)=Prop%AtomType(i)
!
  if(Prop%mp<mid)Prop%mp=mid
 end do
!
 do i=1,Prop%mp
  Prop%crx(i)=Prop%crx(i)/Prop%cmass(i)
  Prop%cry(i)=Prop%cry(i)/Prop%cmass(i)
  Prop%crz(i)=Prop%crz(i)/Prop%cmass(i)
! print*,i,Prop%cmass(i)
! print*,i,Prop%AtomPerMol(i)
!  Prop%cvx(i)=Prop%cvx(i)/Prop%cmass(i)
!  Prop%cvy(i)=Prop%cvy(i)/Prop%cmass(i)
!  Prop%cvz(i)=Prop%cvz(i)/Prop%cmass(i)
!print*,Prop%moltype(i),Prop%crz(i)
 end do
!
! write(*,*)"number of molecules is",Prop%mp
!stop
!
 return
 end subroutine
!
!-----------------------------------------------------------------
! calculate  mass flux
 subroutine mass_flux(Prop)
!
 implicit none
!
  type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
  integer i
!
! ******************************************************************
!
 Prop%massflux=0d0
 do i=1,Prop%np
  Prop%massflux(1) = Prop%massflux(1) + Prop%mass(i)*Prop%vx(i)
  Prop%massflux(2) = Prop%massflux(2) + Prop%mass(i)*Prop%vy(i)
  Prop%massflux(3) = Prop%massflux(3) + Prop%mass(i)*Prop%vz(i)
 end do
!
 return
 end subroutine
!
!-----------------------------------------------------------------
! calculate  charge flux
 subroutine charge_flux(Prop)
!
 implicit none
!
  type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
  integer i
!
! ******************************************************************
!
!
 Prop%chargeflux=0d0
 do i=1,Prop%np
!
  if(Prop%atomid(i)==1)Prop%charge(i)=-0.8476d0
  if(Prop%atomid(i)==2)Prop%charge(i)=-0.4238d0
  if(Prop%atomid(i)==3)Prop%charge(i)=1d0
  if(Prop%atomid(i)==4)Prop%charge(i)=-1d0
!
  Prop%chargeflux(1) = Prop%chargeflux(1) + Prop%charge(i)*Prop%vx(i)
  Prop%chargeflux(2) = Prop%chargeflux(2) + Prop%charge(i)*Prop%vy(i)
  Prop%chargeflux(3) = Prop%chargeflux(3) + Prop%charge(i)*Prop%vz(i)
 end do
!
 return
 end subroutine
!
!-----------------------------------------------------------------
! calculate  charge flux
 subroutine dipole_moment(Prop)
!
 implicit none
!
  type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
  integer :: i,j,m,mid,midn
!
! ******************************************************************
!
!
 Prop%dipole_moment(:)=0d0
 do i=1,Prop%np
  Prop%dipole_moment(1) = Prop%dipole_moment(1) + Prop%charge(i)*Prop%rx(i)
  Prop%dipole_moment(2) = Prop%dipole_moment(2) + Prop%charge(i)*Prop%ry(i)
  Prop%dipole_moment(3) = Prop%dipole_moment(3) + Prop%charge(i)*Prop%rz(i)
 end do
!
 Prop%dmx(:)=0d0
 Prop%dmy(:)=0d0
 Prop%dmz(:)=0d0
!
 midn=0
 do i=1,Prop%Np
  if(Prop%atomType(i)==1)then
   mid=Prop%molID(i)
   midn=midn+1
   m=0
!
   Prop%molx(midn)=Prop%rx(i)
   Prop%moly(midn)=Prop%ry(i)
   Prop%molz(midn)=Prop%rz(i)
!
   do j=1,Prop%Np
      if(Prop%molID(j)==mid)then
       if(Prop%atomType(j)==2)then
        Prop%dmx(midn)=Prop%dmx(midn)+(Prop%rx(j)-Prop%rx(i))*Prop%charge(j)
        Prop%dmy(midn)=Prop%dmy(midn)+(Prop%ry(j)-Prop%ry(i))*Prop%charge(j)
        Prop%dmz(midn)=Prop%dmz(midn)+(Prop%rz(j)-Prop%rz(i))*Prop%charge(j)
        m=m+1
       end if
      end if
      if(m==2)exit
     end do
     end if
    end do
    Prop%mp=midn

 return
 end subroutine
!
!---------------------------------------------------------
! read configuration data form gz compressed file
 subroutine read_configuration_gz(cnfile,cxxx,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128, cxxx*128
 integer, intent(in) :: iunit
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
!
 integer :: i
 real*8 :: dum
 character :: dummy
!
 character :: cnfilexxx*128, cnfilexxxf*128
 character :: command1*10, command2*128, command3*128
 character :: command*128
!
! ******************************************************************
!
 cnfilexxx = trim(cnfile)//'.gz'
!
 command1 = 'gunzip -c '
 command2 = ' > ./'//trim(cxxx)
 command3 = 'rm -f ./'//trim(cxxx)
 command  = command1//trim(cnfilexxx)//trim(command2)
 command3 = trim(command3)
!
 call system(command)
!
 cnfilexxxf = './'//trim(cxxx)
 cnfilexxxf = trim(cnfilexxxf)
 open(unit=iunit,file=cnfilexxxf,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) Prop%np
  read(iunit,*) dummy
  read(iunit,*) Prop%xlow, Prop%xhigh
  read(iunit,*) Prop%ylow, Prop%yhigh
  read(iunit,*) Prop%zlow, Prop%zhigh
  read(iunit,*) dummy
!
  print*,(Prop%xhigh-Prop%xlow),(Prop%yhigh-Prop%ylow),(Prop%zhigh-Prop%zlow)
!
  do i = 1, Prop%np
   read(iunit,*) Prop%atomID(i), Prop%molID(i),Prop%atomType(i), Prop%mass(i), Prop%charge(i),&
                 Prop%rx(i), Prop%ry(i), Prop%rz(i), &
                 Prop%vx(i), Prop%vy(i), Prop%vz(i)
  end do
 close(iunit)
!
 call system(command3)
!
 Prop%lx = Prop%xhigh - Prop%xlow
 Prop%ly = Prop%yhigh - Prop%ylow
 Prop%lz = Prop%zhigh - Prop%zlow
!
 Prop%NpType(:)=0
! *** normalize position data to the range -0.5 to +0.5 ***
 do i = 1, Prop%np
  Prop%rx(i) = ( Prop%rx(i) - (Prop%xlow+Prop%xhigh)/2d0 )/Prop%lx
  Prop%ry(i) = ( Prop%ry(i) - (Prop%ylow+Prop%yhigh)/2d0 )/Prop%ly
  Prop%rz(i) = ( Prop%rz(i) - (Prop%zlow+Prop%zhigh)/2d0 )/Prop%lz
!
  if (Prop%rx(i).eq.0.5d0) Prop%rx(i) = -0.5d0
  if (Prop%ry(i).eq.0.5d0) Prop%ry(i) = -0.5d0
  if (Prop%rz(i).eq.0.5d0) Prop%rz(i) = -0.5d0
!
! *** count type
  Prop%NpType(Prop%AtomType(i))=Prop%NpType(Prop%AtomType(i)) +1
 end do
!
 return
 end subroutine
!
!---------------------------------------------------------! read configuration data form gz compressed file
 subroutine read_configuration_msd_gz(cnfile,cxxx,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128, cxxx*128
 integer, intent(in) :: iunit
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
!
 integer :: i
 real*8 :: dum
 character :: dummy
!
 character :: cnfilexxx*128, cnfilexxxf*128
 character :: command1*10, command2*128, command3*128
 character :: command*128
!
! ******************************************************************
!
 cnfilexxx = trim(cnfile)//'.gz'
!
 command1 = 'gunzip -c '
 command2 = ' > ./'//trim(cxxx)
 command3 = 'rm -f ./'//trim(cxxx)
 command  = command1//trim(cnfilexxx)//trim(command2)
 command3 = trim(command3)
!
 call system(command)
!
 cnfilexxxf = './'//trim(cxxx)
 cnfilexxxf = trim(cnfilexxxf)
 open(unit=iunit,file=cnfilexxxf,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) Prop%np
  read(iunit,*) dummy
  read(iunit,*) Prop%xlow, Prop%xhigh
  read(iunit,*) Prop%ylow, Prop%yhigh
  read(iunit,*) Prop%zlow, Prop%zhigh
  read(iunit,*) dummy
!
 ! print*,(Prop%xhigh-Prop%xlow),(Prop%yhigh-Prop%ylow),(Prop%zhigh-Prop%zlow)
!
  do i = 1, Prop%np
   read(iunit,*) Prop%atomID(i), Prop%molID(i),Prop%atomType(i), Prop%mass(i), Prop%rx(i), Prop%ry(i), Prop%rz(i)
  end do
 close(iunit)
!
 call system(command3)
!
 Prop%lx = Prop%xhigh - Prop%xlow
 Prop%ly = Prop%yhigh - Prop%ylow
 Prop%lz = Prop%zhigh - Prop%zlow
!
 return
 end subroutine
!
!---------------------------------------------------------! read configuration data form gz compressed file
 subroutine read_configuration_dipole_gz(cnfile,cxxx,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128, cxxx*128
 integer, intent(in) :: iunit
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
!
 integer :: i
 real*8 :: dum
 character :: dummy
!
 character :: cnfilexxx*128, cnfilexxxf*128
 character :: command1*10, command2*128, command3*128
 character :: command*128
!
! ******************************************************************
!
 cnfilexxx = trim(cnfile)//'.gz'
!
 command1 = 'gunzip -c '
 command2 = ' > ./'//trim(cxxx)
 command3 = 'rm -f ./'//trim(cxxx)
 command  = command1//trim(cnfilexxx)//trim(command2)
 command3 = trim(command3)
!
 call system(command)
!
 cnfilexxxf = './'//trim(cxxx)
 cnfilexxxf = trim(cnfilexxxf)
 open(unit=iunit,file=cnfilexxxf,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) Prop%np
  read(iunit,*) dummy
  read(iunit,*) Prop%xlow, Prop%xhigh
  read(iunit,*) Prop%ylow, Prop%yhigh
  read(iunit,*) Prop%zlow, Prop%zhigh
  read(iunit,*) dummy
!
 ! print*,(Prop%xhigh-Prop%xlow),(Prop%yhigh-Prop%ylow),(Prop%zhigh-Prop%zlow)
!
  do i = 1, Prop%np
   read(iunit,*) Prop%atomID(i), Prop%molID(i),Prop%atomType(i), Prop%mass(i), Prop%charge(i), Prop%rx(i), Prop%ry(i), Prop%rz(i)
  end do
 close(iunit)
!
 call system(command3)
!
 Prop%lx = Prop%xhigh - Prop%xlow
 Prop%ly = Prop%yhigh - Prop%ylow
 Prop%lz = Prop%zhigh - Prop%zlow
!
 return
 end subroutine
!
!
!---------------------------------------------------------! read configuration data form gz compressed file
 subroutine read_configuration_ripple_gz(cnfile,cxxx,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128, cxxx*128
 integer, intent(in) :: iunit
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
!
 integer :: i
 real*8 :: dum
 character :: dummy
!
 character :: cnfilexxx*128, cnfilexxxf*128
 character :: command1*10, command2*128, command3*128
 character :: command*128
!
! ******************************************************************
!
 cnfilexxx = trim(cnfile)//'.gz'
!
 command1 = 'gunzip -c '
 command2 = ' > ./'//trim(cxxx)
 command3 = 'rm -f ./'//trim(cxxx)
 command  = command1//trim(cnfilexxx)//trim(command2)
 command3 = trim(command3)
!
 call system(command)
!
 cnfilexxxf = './'//trim(cxxx)
 cnfilexxxf = trim(cnfilexxxf)
 open(unit=iunit,file=cnfilexxxf,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) Prop%np
  read(iunit,*) dummy
  read(iunit,*) Prop%xlow, Prop%xhigh
  read(iunit,*) Prop%ylow, Prop%yhigh
  read(iunit,*) Prop%zlow, Prop%zhigh
  read(iunit,*) dummy
!
 ! print*,(Prop%xhigh-Prop%xlow),(Prop%yhigh-Prop%ylow),(Prop%zhigh-Prop%zlow)
!
  do i = 1, Prop%np
   read(iunit,*) Prop%atomID(i), Prop%molID(i),Prop%atomType(i), Prop%mass(i), Prop%charge(i), Prop%rx(i), Prop%ry(i), Prop%rz(i)
  end do
 close(iunit)
!
 call system(command3)
!
 Prop%lx = Prop%xhigh - Prop%xlow
 Prop%ly = Prop%yhigh - Prop%ylow
 Prop%lz = Prop%zhigh - Prop%zlow
!

!print*,prop%lx
!stop
 return
 end subroutine
!
!
!---------------------------------------------------------! read configuration data form gz compressed file
 subroutine read_configuration_ca_gz(cnfile,cxxx,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128, cxxx*128
 integer, intent(in) :: iunit
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
!
 integer :: i
 real*8 :: dum
 character :: dummy
!
 character :: cnfilexxx*128, cnfilexxxf*128
 character :: command1*10, command2*128, command3*128
 character :: command*128
!
! ******************************************************************
!
 cnfilexxx = trim(cnfile)//'.gz'
!
 command1 = 'gunzip -c '
 command2 = ' > ./'//trim(cxxx)
 command3 = 'rm -f ./'//trim(cxxx)
 command  = command1//trim(cnfilexxx)//trim(command2)
 command3 = trim(command3)
!
 call system(command)
!
 cnfilexxxf = './'//trim(cxxx)
 cnfilexxxf = trim(cnfilexxxf)
 open(unit=iunit,file=cnfilexxxf,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) Prop%np
  read(iunit,*) dummy
  read(iunit,*) Prop%xlow, Prop%xhigh
  read(iunit,*) Prop%ylow, Prop%yhigh
  read(iunit,*) Prop%zlow, Prop%zhigh
  read(iunit,*) dummy
!
 ! print*,(Prop%xhigh-Prop%xlow),(Prop%yhigh-Prop%ylow),(Prop%zhigh-Prop%zlow)
!
  do i = 1, Prop%np
   read(iunit,*) Prop%atomID(i), Prop%atomType(i), Prop%rx(i), Prop%ry(i), Prop%rz(i)
  end do
 close(iunit)
!
 call system(command3)
!
 Prop%lx = Prop%xhigh - Prop%xlow
 Prop%ly = Prop%yhigh - Prop%ylow
 Prop%lz = Prop%zhigh - Prop%zlow

 return
 end subroutine

!-----------------------------------------------------------------
! read configuration data form gz compressed file
 subroutine read_configuration_gz_old(cnfile,cxxx,iunit,Prop)
!
 implicit none
!
 character, intent(in) :: cnfile*128, cxxx*128
 integer, intent(in) :: iunit
 type(ReadFilesProp),intent(inout) :: Prop !structure of type readfilesProp
!
 integer :: i
 real*8 :: dum
 character :: dummy
!
 character :: cnfilexxx*128, cnfilexxxf*128
 character :: command1*10, command2*128, command3*128
 character :: command*128
!
! ******************************************************************
!
 cnfilexxx = trim(cnfile)//'.gz'
!
 command1 = 'gunzip -c '
 command2 = ' > ./'//trim(cxxx)
 command3 = 'rm -f ./'//trim(cxxx)
 command  = command1//trim(cnfilexxx)//trim(command2)
 command3 = trim(command3)
!
 call system(command)
!
 cnfilexxxf = './'//trim(cxxx)
 cnfilexxxf = trim(cnfilexxxf)
 open(unit=iunit,file=cnfilexxxf,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) Prop%np
  read(iunit,*) dummy
  read(iunit,*) Prop%xlow, Prop%xhigh
  read(iunit,*) Prop%ylow, Prop%yhigh
  read(iunit,*) Prop%zlow, Prop%zhigh
  read(iunit,*) dummy
  do i = 1, Prop%np
   read(iunit,*) Prop%atomID(i), Prop%atomType(i), Prop%mass(i), Prop%rx(i), Prop%ry(i), Prop%rz(i), &
                 Prop%vx(i), Prop%vy(i), Prop%vz(i)
  end do
 close(iunit)
!
 call system(command3)
!
 Prop%lx = Prop%xhigh - Prop%xlow
 Prop%ly = Prop%yhigh - Prop%ylow
 Prop%lz = Prop%zhigh - Prop%zlow
!
 Prop%NpType(:)=0
! *** normalize position data to the range -0.5 to +0.5 ***
 do i = 1, Prop%np
  Prop%rx(i) = ( Prop%rx(i) - (Prop%xlow+Prop%xhigh)/2d0 )/Prop%lx
  Prop%ry(i) = ( Prop%ry(i) - (Prop%ylow+Prop%yhigh)/2d0 )/Prop%ly
  Prop%rz(i) = ( Prop%rz(i) - (Prop%zlow+Prop%zhigh)/2d0 )/Prop%lz
!
  if (Prop%rx(i).eq.0.5d0) Prop%rx(i) = -0.5d0
  if (Prop%ry(i).eq.0.5d0) Prop%ry(i) = -0.5d0
  if (Prop%rz(i).eq.0.5d0) Prop%rz(i) = -0.5d0
!
! *** count type
  Prop%NpType(Prop%AtomType(i))=Prop%NpType(Prop%AtomType(i)) +1
 end do
!
 return
 end subroutine
!

 end module
!
