 program poly
! 
!***************************************************************
 use class_arguments
 use class_num2char
 use class_random
 implicit none
!
    integer,parameter :: ntmax=5000000
    integer :: i,j,k,n
    integer :: count
    integer :: seed
    real(8) :: xlbox,ylbox,zlbox !<box size
    real(8) :: xx,yy,zz !<each position
    integer :: mm
    integer :: npp
    real(8) :: pl,rr,th,rx
!
    integer :: nw !number of water molecule
    real(8) :: dw0 !<length per one water
    real(8) :: dw1 !<length per one water
    real(8) :: sig1 !<sigma ratio
    real(8) :: er !<exclusion radius
    integer :: id !<charge distribution type
    integer :: nt !<total solvent particles
    real(8) :: fraction !<molar fraction
!
    character :: aID*2
    character :: ff*128
!
    type(argProp) :: arg !<structure for arguments
!
!**************************************************************
!command line argument
    call argInit(arg)
    call argGetReal("-lx",xlbox,arg,30d0)
    call argGetReal("-ly",ylbox,arg,30d0)
    call argGetReal("-lz",zlbox,arg,30d0)
    call argGetReal("-lw",dw0,arg,1.1d0)
    call argGetInt("-pn",mm,arg,10)
    call argGetInt("-nt",nt,arg,ntmax)
    call argGetInt("-npp",npp,arg,1)
    call argGetReal("-pl",pl,arg,1.0d0)
    call argGetReal("-s1",sig1,arg,1.0d0)
    call argGetReal("-er",er,arg,3.0d0)
    call argGetInt("-ch",id,arg,0)
    call argGetReal("-a",fraction,arg,0d0)
    call argGetInt("-se",seed,arg,345)

    if(arg%narg==0)then
     write(*,"(a)")" ----make polymer configuration----"
     write(*,"(a)")" [Usage]"
     write(*,"(a)")" -nt [int]: [int] total particle if specified"
     write(*,"(a)")" -lx [real]: [real] is box size in x (necessary)"
     write(*,"(a)")" -ly [real]: [real] is box size in y (necessary)"
     write(*,"(a)")" -lz [real]: [real] is box size in z (necessary)"
     write(*,"(a)")" -lw [real]: [real] is spacing for molecules (default 1.1)"
     write(*,"(a)")" -npp [int]: [int] number of polymers (default 1)"
     write(*,"(a)")" -pn [int]: [int] number of monomer (default 10)"
     write(*,"(a)")" -pl [real]: [real] polymer interval (default 1.0)"
     write(*,"(a)")" -er [real]: [real] exclution radius (for ch=2, default 3.0)"
     write(*,"(a)")" -s1 : sigma ratio"
     write(*,"(a)")" -a [real]: [real] is molar fraction of solute (default=0)"
     write(*,"(a)")" -ch [int]: charge distribution type [1-4] (default=1)"
     write(*,"(a)")"            0: graphene water-ethanol"
     write(*,"(a)")"            1: LJ with two volumes  pure A1 / A1+A2"
     write(*,"(a)")" -se [int]: [int] is a seed for randam number (default=345)"
     stop
    end if
!**************************************************************
!
    call sGenrand(seed,RandomProp)
!
!----position of wall atoms in 2d
!
    select case(id)
!
    case(1)
!print*,xlbox,ylbox,zlbox
!stop
!wall LJ mixture membrane
!
!-----make mt header

 print*,'write("Data Atoms") {'

 if(mm==1)then
   i=1
   xx=0d0
   yy=0d0
   zz=0d0
   ff='$atom:M'//trim(int2char(i,1))//' '//'$mol:.   @atom:M 0.0000 '//trim(dble2char(xx,1,3))//' '//trim(dble2char(yy,1,3))//' '//trim(dble2char(zz,1,3))
   write(*,"(a)")ff
 else
  do i=1,mm
   rx=pl/10d0
   rr=3d0
   th=dasin(pl/2d0/rr)
   xx=-dble(mm-1)*rx/2d0+dble(i)*rx
   yy=rr*dcos(th*dble(i))
   zz=rr*dsin(th*dble(i))
   ff='$atom:M'//trim(int2char(i,1))//' '//'$mol:.   @atom:M 0.0000 '//trim(dble2char(xx,1,3))//' '//trim(dble2char(yy,1,3))//' '//trim(dble2char(zz,1,3))
   write(*,"(a)")ff
  end do
 end if

 print*,'}'
 print*,' '
 print*,'write_once("Data Masses") {'
 print*,' # atomType mass'
 print*,'@atom:M    1.0'
 print*,'}'
 print*,' '
 print*,'write("Data Bonds") {'
 print*,'# bondID bondType atomID1 atomID2'

 do i=1,mm-1
  ff='$bond:MM'//trim(int2char(i,1))//' '//'@bond:MM $atom:M'//trim(int2char(i,1))//' $atom:M'//trim(int2char(i+1,1))
  write(*,"(a)")ff
 end do

 print*,'}'
 print*,' '

 print*,'# --- Force-field parameters go in the "In Settings" section: ---'
 print*,' '
 print*,'write_once("In Settings") {'
 print*,'# -- Non-bonded (Pair) interactions --'
 print*,'#  atomType1 atomType2 parameter-list (epsilon, sigma)'
 print*,'pair_coeff @atom:M @atom:M     1.0  1.0'
 print*,' '
 print*,'# -- Bonded interactions --'
 print*,'# bondType parameter list (k_bond, r0)'
 print*,'bond_coeff @bond:MM 1000.00 1.0'

 print*,'}'
 print*,'} # Pl'

     xx=xlbox/2d0
     yy=ylbox/2d0
     zz=zlbox/2d0

     aID="Pl"
     count=1
    do i=1,npp
     if(npp>1)then
      yy=ylbox/4d0+ylbox/2d0/dble(npp)*dble(i)
     end if
     write(*,'(a)')"Poly"//trim(int2char(i,1))//" = new "//trim(aID)//".move(" &
          //trim(dble2char(xx,1,3))//", " &
          //trim(dble2char(yy,1,3))//", " &
          //trim(dble2char(zz,1,3))//")"
    end do

     call place_A1A2
!
    case(2)
    
     xx=xlbox/2d0
     yy=ylbox/2d0
     zz=zlbox/2d0

     aID="Pl"
     count=1
     write(*,'(a)')"Poly"//trim(int2char(count,1))//" = new "//trim(aID)//".move(" &
          //trim(dble2char(xx,1,3))//", " &
          //trim(dble2char(yy,1,3))//", " &
          //trim(dble2char(zz,1,3))//")"
!
     call place_A1A2
!
!
    case(3)
!
     call place_porous
!
    end select 
!
     write(*,*) ""
     write(*,"(a)") 'write_once("Data Boundary") {'
     write(*,"(a)") "    "//trim(dble2char(0d0,1,1))//" " &
                          //trim(dble2char(xlbox,1,2))//" xlo xhi"   
     write(*,"(a)") "    "//trim(dble2char(0d0,1,1))//" " &
                          //trim(dble2char(ylbox,1,2))//" ylo yhi"   
     write(*,"(a)") "    "//trim(dble2char(0d0,1,1))//" " &
                          //trim(dble2char(zlbox,1,2))//" zlo zhi"   
     write(*,"(a)") "}"
!
!-----------------------------------------------------------------------------------     
 contains
!-----place mixture upper and lower regions
!
!-----place A1-A2 mixture in bulk
  subroutine place_A1A2
   implicit none
   integer :: n,nw,ns,n1
   integer :: count1,count2,count3
   real(8),allocatable :: xn(:),yn(:),zn(:)
   integer,allocatable :: sol(:)
   real(8) :: p
   real(8) :: dw
   real(8) :: wx0,wy0,wz0
   real(8) :: x0,y0,z0,rr,r0

    x0=xlbox/2d0
    y0=ylbox/2d0
    z0=zlbox/2d0

    dw=((1d0-fraction)*dw0**3+fraction*(dw0*sig1)**3)**(1d0/3d0)
    nw=2*int(xlbox/dw)*int(ylbox/dw)*int(zlbox/dw)
    allocate(xn(nw),yn(nw),zn(nw),sol(nw))
    wx0=dw/2d0
    wy0=dw/2d0
    wz0=dw/2d0

    xx=wx0
    yy=wy0
    zz=wz0
!
    count1=1
!
    do n=1,nw
!
     xn(count1)=xx
     yn(count1)=yy
     zn(count1)=zz

     xx=xx+dw
     if(xx+dw/2d0>xlbox)then
      yy=yy+dw
      xx=wx0
      if(yy+dw/2d0>ylbox)then
!
       zz=zz+dw
       yy=wy0

      if (zz>zlbox)then
       exit
      end if
     end if
    end if

    if(count1==nt)exit

     rr=dsqrt((xx-x0)**2+(yy-y0)**2+(zz-z0)**2)
     if(id==2 .and. rr<er)then
      cycle
     else
      count1=count1+1
     end if

   end do
!
   if(count1<nt .and. nt<ntmax)then
    print*,"<<error not enough particles; try larger density>>"
    stop
   end if
   ns=nint(dble(count1)*fraction)
   nw=count1

!
   sol(:)=0
   if(ns>0)then
    do n=1,ns
      do
        p=dgenrand(RandomProp)
        n1=dble(nw)*p+1
        if(n1>nw)n1=nw
        if(sol(n1)==0)then
         sol(n1)=1
         exit
        end if
       end do
     end do
   end if

!
    do n=1,nw
      if(sol(n)==1)then    
        write(*,'(a)')"Sol"//trim(int2char(n,1))//" = new A2.move(" &
            //trim(dble2char(xn(n),1,3))//", " &
            //trim(dble2char(yn(n),1,3))//", " &
            //trim(dble2char(zn(n),1,3))//")"
       else if(sol(n)==0)then
        write(*,'(a)')"Sol"//trim(int2char(n,1))//" = new A1.move(" &
            //trim(dble2char(xn(n),1,3))//", " &
            //trim(dble2char(yn(n),1,3))//", " &
            //trim(dble2char(zn(n),1,3))//")"
       end if
     end do
!    print*,count_eth,count_res

    deallocate(xn,yn,zn,sol)
   return
  end subroutine    
!
!
!-----place porous geometry
  subroutine place_porous
   implicit none
   integer :: n,nw,ns,n1
   integer :: count1,count2,count3
   real(8),allocatable :: xn(:),yn(:),zn(:)
   integer,allocatable :: sol(:)
   real(8) :: p
   real(8) :: dw
   real(8) :: wx0,wy0,wz0
   real(8) :: x0,y0,z0,rr,r0

    x0=xlbox/2d0
    y0=ylbox/2d0
    z0=zlbox/2d0

    dw=((1d0-fraction)*dw0**3+fraction*(dw0*sig1)**3)**(1d0/3d0)
    nw=2*int(xlbox/dw)*int(ylbox/dw)*int(zlbox/dw)
    allocate(xn(nw),yn(nw),zn(nw),sol(nw))
    wx0=dw/2d0
    wy0=dw/2d0
    wz0=dw/2d0

    xx=wx0
    yy=wy0
    zz=wz0
!
    count1=1
!
    do n=1,nw
!
     xn(count1)=xx
     yn(count1)=yy
     zn(count1)=zz

     xx=xx+dw
     if(xx+dw/2d0>xlbox)then
      yy=yy+dw
      xx=wx0
      if(yy+dw/2d0>ylbox)then
!
       zz=zz+dw
       yy=wy0

      if (zz>zlbox)then
       exit
      end if
     end if
    end if

    if(count1==nt)exit

     rr=dsqrt((xx-x0)**2+(yy-y0)**2+(zz-z0)**2)
     if(id==2 .and. rr<er)then
      cycle
     else
      count1=count1+1
     end if

   end do
!
   if(count1<nt .and. nt<ntmax)then
    print*,"<<error not enough particles; try larger density>>"
    stop
   end if
   ns=nint(dble(count1)*fraction)
   nw=count1

!
   sol(:)=0
   if(ns>0)then
    do n=1,ns
      do
        p=dgenrand(RandomProp)
        n1=dble(nw)*p+1
        if(n1>nw)n1=nw
!        if(sol(n1)==0 .and. xn(n1)>xlbox/2d0)then
        if(sol(n1)==0)then
         sol(n1)=1
         exit
        end if
       end do
     end do
   end if

!
    do n=1,nw
      if(sol(n)==1)then    
        write(*,'(a)')"Sol"//trim(int2char(n,1))//" = new A2.move(" &
            //trim(dble2char(xn(n),1,3))//", " &
            //trim(dble2char(yn(n),1,3))//", " &
            //trim(dble2char(zn(n),1,3))//")"
       else if(sol(n)==0)then
        write(*,'(a)')"Sol"//trim(int2char(n,1))//" = new A1.move(" &
            //trim(dble2char(xn(n),1,3))//", " &
            //trim(dble2char(yn(n),1,3))//", " &
            //trim(dble2char(zn(n),1,3))//")"
       end if
     end do
!    print*,count_eth,count_res

    deallocate(xn,yn,zn,sol)
   return
  end subroutine    
!
   
  end program
!
!
!

