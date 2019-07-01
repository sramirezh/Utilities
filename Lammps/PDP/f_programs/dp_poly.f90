program dp_poly
!***************************************************************
! analyze diffusio-phoresis polymers
!***************************************************************
 use class_arguments
 use class_num2char
 use class_quick_sort
    implicit none
! 
    real(8) :: pi,kboltzmann,unit_charge,epsilon0
    parameter ( pi  = 4d0*datan(1d0))
    parameter ( kboltzmann  =   1.38064852d-23)
    parameter ( unit_charge  =   1.6021766208d-19)
    parameter ( epsilon0  =   8.854187817d-12)

!**************************************************************
    integer :: np
    real(8) :: xlow,xhigh,ylow,yhigh,zlow,zhigh
    real(8) :: lx,ly,lz
    integer,allocatable :: aType(:), at(:), NpType(:)
    real(8),allocatable ::  aID(:), aID_sort(:), rx0(:), ry0(:), rz0(:), ifx0(:), ify0(:), ifz0(:), vx0(:), vy0(:), vz0(:)
    real(8),allocatable :: rx(:), ry(:), rz(:), ifx(:), ify(:), ifz(:), vx(:), vy(:), vz(:)
!**************************************************************

    integer :: nstep, dstep, rstep, sstep
    integer :: istepmax,rstepmax,sstepmax
    integer :: start
!
    integer :: i, j, k,kk

    integer :: c
    integer :: nc
    real(8) :: xc, yc, zc
    real(8) :: xu, yu, zu
    real(8) :: rr
    real(8) :: cad_av,conrg
    real(8) :: rg_av,rg2_av
    real(8) :: vsph
    real(8),allocatable :: rg2(:),rg(:),cad(:)
    real(8),allocatable :: xp(:), yp(:), zp(:)
    real(8),allocatable :: xcc(:), ycc(:), zcc(:)
    real(8),allocatable :: xca(:), yca(:), zca(:)
    real(8) :: timestep
!
    integer :: count
!
    character :: cnfile*128
    character :: file_dump*128
    character :: file_adsorption*128
!
    type(argProp) :: arg !<structure for arguments
    type(quickSortProp) :: QS !<structure for quick sort
!
!**************************************************************
!command line argument
    call argInit(arg)
    call argGetChar("-f",file_dump,arg,"./conf/dumpfile")
    call argGetInt("-n",nstep,arg,100000) ! total step
    call argGetInt("-d",dstep,arg,10000)   ! interval of step
 !   call argGetInt("-r",rstep,arg,10000)   ! interval of step
    call argGetInt("-s",sstep,arg,0)   ! start step
!    call argExist("-o2d",IsOut2D,arg,.false.)   ! to output 2d data
    call argGetReal("-dt",timestep,arg,0.005d0)   ! timestep (fs)
!
    if(arg%narg==0)then
     write(*,"(a)")" ----comput mean square displacement----"
     write(*,"(a)")" [Usage]"
     write(*,"(a)")" -f [str]: [str] is file name to be read without step number (default=./conf/dumpfile)"
     write(*,"(a)")" -n [int]: [int] is total step number (default 100000)"
     write(*,"(a)")" -d [int]: [int] is step number interval (default 10000)"
     write(*,"(a)")" -s [int]: [int] is start step number (default= 0)"
!     write(*,"(a)")" -o2d: output 2d data"
     write(*,"(a)")" -dt [real]: [real] is time step (default=0.005d0)"
     stop
    end if
!
    file_adsorption="adsorption.dat"
!
   istepmax = (nstep-sstep)/dstep 
   sstepmax = sstep/dstep
!   rstepmax = rstep/dstep
!
  allocate( xp(200), yp(200), zp(200) )
  allocate( cad(istepmax), rg(istepmax), rg2(istepmax))
  allocate( xcc(istepmax), ycc(istepmax), zcc(istepmax))
  allocate( xca(istepmax), yca(istepmax), zca(istepmax))
!
! *** start iteration ************************
!
 !
 do k = sstepmax,sstepmax+istepmax
    kk=k-sstepmax+1
     call readFromFile(k)

     c=0
     do i=1,np
       if(at(i)==3)then
        c=c+1
        xp(c)=rx(i)
        yp(c)=ry(i)
        zp(c)=rz(i)
        if(c>1)then
         if(xp(c)-xp(c-1)>lx/2d0) xp(c)=xp(c)-lx
         if(xp(c)-xp(c-1)<-lx/2d0) xp(c)=xp(c)+lx
         if(yp(c)-yp(c-1)>ly/2d0) yp(c)=yp(c)-ly
         if(yp(c)-yp(c-1)<-ly/2d0) yp(c)=yp(c)+ly
        end if
       end if
      end do
      nc=c

!******************************************************************************
!  compute center of mass of polymer
!******************************************************************************
      xc=0d0
      yc=0d0
      zc=0d0
      do c=1,nc
        xc=xc+xp(c)
        yc=yc+yp(c)
        zc=zc+zp(c)
      end do
      xc=xc/dble(nc)
      yc=yc/dble(nc)
      zc=zc/dble(nc)

      xcc(kk)=xc
      ycc(kk)=yc
      zcc(kk)=zc
      if(xcc(kk)>lx) xcc(kk)=xcc(kk)-lx
      if(xcc(kk)<0d0) xcc(kk)=xcc(kk)+lx
      if(ycc(kk)>ly) ycc(kk)=ycc(kk)-ly
      if(ycc(kk)<0d0) ycc(kk)=ycc(kk)+ly
      if(zcc(kk)>lz) zcc(kk)=zcc(kk)-lz
      if(zcc(kk)<0d0) zcc(kk)=zcc(kk)+lz


!******************************************************************************
!  compute absolute position of polymer
!******************************************************************************
     if(kk>1)then
      xca(kk)=xca(kk-1)+xcc(kk)-xcc(kk-1)
      yca(kk)=yca(kk-1)+ycc(kk)-ycc(kk-1)
      zca(kk)=zca(kk-1)+zcc(kk)-zcc(kk-1)
      yca(kk)=0d0
      zca(kk)=0d0
      if(xcc(kk)-xcc(kk-1)>lx/2d0) xca(kk)=xca(kk)-lx
      if(xcc(kk)-xcc(kk-1)<-lx/2d0) xca(kk)=xca(kk)+lx
      if(ycc(kk)-ycc(kk-1)>ly/2d0) yca(kk)=yca(kk)-ly
      if(ycc(kk)-ycc(kk-1)<-ly/2d0) yca(kk)=yca(kk)+ly
      if(zcc(kk)-zcc(kk-1)>lz/2d0) zca(kk)=zca(kk)-lz
      if(zcc(kk)-zcc(kk-1)<-lz/2d0) zca(kk)=zca(kk)+lz
     else
      xca(kk)=0d0
      yca(kk)=0d0
      zca(kk)=0d0
     end if

!******************************************************************************
!  compute Rg
!******************************************************************************
     rg2(kk)=0d0
     do c=1,nc
       rg2(kk)=rg2(kk)+(xp(c)-xc)**2+(yp(c)-yc)**2+(zp(c)-zc)**2
     end do
     rg2(kk)=rg2(kk)/dble(nc)
     rg(kk)=dsqrt(rg2(kk))

!******************************************************************************
!  comute adsorption 1
!******************************************************************************
     count=0
     do i=1,np
       if(at(i)==2)then
        xu=rx(i)
        yu=ry(i)
        zu=rz(i)
        if(xu-xcc(kk)>lx/2d0)xu=xu-lx
        if(xu-xcc(kk)<-lx/2d0)xu=xu+lx
        if(yu-ycc(kk)>ly/2d0)yu=yu-ly
        if(yu-ycc(kk)<-ly/2d0)yu=yu+ly
        if(zu-zcc(kk)>lz/2d0)zu=zu-lz
        if(zu-zcc(kk)<-lz/2d0)zu=zu+lz
        rr=dsqrt((xu-xcc(kk))**2+(yu-ycc(kk))**2+(zu-zcc(kk))**2)

        if(rr<rg(kk))then
         count=count+1
        end if
       end if
     end do
     cad(kk)=dble(count)


 end do
!
  cnfile = trim(file_adsorption)
  open(unit=83,file=cnfile,status='unknown')
   write(*,*)"write time series data in ", trim(cnfile)
    rewind(83)
!     write(83,"(a)") 'variables="time[fs]" "survrat" "lifetime" '
     do i=1,istepmax
       write(83,'(1p20e17.8)')  dble(i)*timestep*dble(dstep),xca(i),yca(i),zca(i),rg(i),cad(i)
     end do
    close(83)
!
  cad_av=0d0
  rg2_av=0d0
  do i=1,istepmax
    cad_av=cad_av+cad(i)
    rg2_av=rg2_av+rg2(i)
  end do
  cad_av=cad_av/dble(istepmax)
  rg2_av=rg2_av/dble(istepmax)
  rg_av=dsqrt(rg2_av)
  
  vsph=4d0*pi*rg_av**3/3d0
  conrg=cad_av/vsph

  print*,"solute_count =",cad_av
  print*,"rg_ave =",rg_av
  print*,"rg_concentration =",conrg

  cnfile = "average_info.dat"
  open(unit=83,file=cnfile,status='unknown')
   write(*,*)"write time series data in ", trim(cnfile)
    rewind(83)
!     write(83,"(a)") 'variables="time[fs]" "survrat" "lifetime" '
       write(83,*)  "soute count =", cad_av
       write(83,*)  "rg ave =", rg_av
       write(83,*)  "rg concentration =", conrg
    close(83)
!


!
 stop 'normal termination!' 
!    
!
! **************************************************************************
 contains 
!---------------------------------------------------------
! read configuration data form gz compressed file
 subroutine readFromFile(k)

 implicit none
!
 integer,intent(in) :: k

 character :: cnfile*128, cxxx*128
 integer :: iunit
!
 logical,save :: first=.true.
 integer :: i,id
 real*8 :: dum
 character :: dummy
!
 character :: cnfilexxx*128, cnfilexxxf*128
 character :: command1*10, command2*128, command3*128
 character :: command*128
!
! ******************************************************************
!
 cxxx = 'exttmp'
 cxxx = trim(cxxx)

 cnfile = trim(file_dump)//trim( int2char( dstep*k ,1 ) )
 print *, k, "/", sstepmax+istepmax, trim(cnfile)

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

 iunit=77
 open(unit=iunit,file=cnfilexxxf,status='unknown')
  rewind(iunit)
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) dummy
  read(iunit,*) np
  read(iunit,*) dummy
  read(iunit,*) xlow, xhigh
  read(iunit,*) ylow, yhigh
  read(iunit,*) zlow, zhigh
  read(iunit,*) dummy
!
  if(first)then
   allocate( aType(np), rx0(np), ry0(np), rz0(np), ifx0(np), ify0(np), ifz0(np), vx0(np), vy0(np), vz0(np) )
   allocate( aID(np), at(np), rx(np), ry(np), rz(np), ifx(np), ify(np), ifz(np), vx(np), vy(np), vz(np) )
   allocate( aID_sort(np), NpType(20) )
   first=.false.
  end if
!
  do i = 1, np
   read(iunit,*) aID(i), aType(i), rx0(i), ry0(i), rz0(i), ifx0(i), ify0(i), ifz0(i), vx0(i), vy0(i), vz0(i)
  end do
 close(iunit)
!
 call system(command3)
!
  lx = xhigh - xlow
  ly = yhigh - ylow
  lz = zhigh - zlow
!  print*,lx,ly,lz
!
 NpType(:)=0
! *** normalize position data to the range -0.5 to +0.5 ***
 do i = 1, np
! *** count type
  NpType(aType(i))=NpType(aType(i)) +1
 end do
!
   call quick_sort(aID,aID_sort,np,QS)  
   do i=1,np
    id=quick_sort_id(i,QS)
    at(i)=aType(id)
    rx(i)=rx0(id)-xlow; ry(i)=ry0(id)-ylow; rz(i)=rz0(id)-zlow
    ifx(i)=ifx0(id); ify(i)=ify0(id); ifz0(i)=ifz(id)
    vx(i)=vx0(id); vy(i)=vy0(id); vz(i)=vz0(id)
    do
       if(0d0<=rx(i) .and. rx(i)<lx .and. 0d0<=ry(i) .and. ry(i)<ly)exit
       if (rx(i)>=lx) rx(i)=rx(i)-lx
       if (ry(i)>=ly) ry(i)=ry(i)-ly
       if (rx(i)<0d0) rx(i)=rx(i)+lx
       if (ry(i)<0d0) ry(i)=ry(i)+ly
    end do
   end do


 return
 end subroutine readFromFile
!
 end program
!

