!************************************************************************************************
!>Lagged Fibonacci Sequence   FP100Q37 [FORTRAN9x VERSION]
!>\details
!>Author: \n
!>Lagged Fibonacci sequence was first introduced by 
!>Mitchell & Moore in 1958, where lags were (55,24). \n
!>\n
!>Taken from: \n
!>    Donald E. Knuth(1994), The art of computer programming,
!>    3rd edition, Vol. 2 p.602 (The answer of execise 3.6.11). \n
!>    (Translated to FORTRAN77 from original ANSI-C program) \n
!>Note: \n
!>    The original version was written by H.Sugimoto of Kyoto Univ. \n
!>    Inherited by H.Yoshida of Toyota central R&D labs.,inc. \n
!>\n
!>  Note: \n
!>    Lags (RANDOM_KK,RANDOM_LL) may be (100,37) or (127,30), (83,258), etc. \n
!>    The cycle of random sequence is 2^51*(2^RANDOM_KK-1). \n
!>    The "guaranteed separation time" is 2^70. \n
!>    The generator produces different random sequence in guaranteed separation time \n
!>    for every given seed, which must be in the range [0,2^30-3]. \n
!>  Note: \n
!>    Correct if and only if the operating system uses IEEE standard arithmetics. \n
!>  Release: \n
!>    Fri Jul  9 19:00:44 JST 1999  Version 0.01 \n
!>    Wed Jul 26 14:51:11 JST 2000  Version 0.02  RANDOM_ARRAYS2D added \n
!>    Fri Jul 11 20:43:08 JST 2003  Version 0.03  Fortran95 \n
!>    Tue Sep 14 13:00:47 JST 2004  Version 0.04  WhatRandom(), Sequence of Structure \n
!>    Thu Apr 24 17:45:30 JST 2008  Version 0.05  reformed, dBoxMuller added (by H.Yoshida) \n
!>    Thu Apr 24 17:45:30 JST 2008  Version 0.06  dBoxMuller added
!>
!>  This package provides following functions:
!
!*************************************************************************************************
!-----------------------------------------------------------
module class_random
  integer(4),parameter::      Random_KK=100 !<  the long lag  
  integer(4),parameter::      Random_LL=37  !<  the short lag 
  integer(4),parameter::      Random_TT=70  !<  guaranteed separation between systems 
!
!
!> property of random number generator
!>\details
!>all members are PRIVATE
  type  RandomProperty  
    private
    integer(4)::      index !<index for random generator  
    real(8)::     state(0:Random_KK-1)  !<state
  end type RandomProperty

  character(LEN=*),parameter::    Random_Kind="fp100q37" !< kind of random sequence
!
  type(RandomProperty) :: RandomProp
!
!-----------------------------------------------------------
  contains

!> comput residue
!>\details
!><<<PRIVATE>>>
  real(8) function  Random_Mod1(x)
    real(8):: x !input value
    Random_Mod1 =x-INT(x)
    return
  end function
!



!-----------------------------------------------------------
!>TEST PROGRAM TO GENERATE  0.274526263073941567**
!>\details
!>----------------------------------------------------------- \n
!>  ORIGINAL PROGRAM BY D.E.KNUTH \n
!>-----------------------------------------------------------
  subroutine  RandomTest
    implicit  none
    real(8):: A(0:2008)
    real(8):: randu(0:Random_KK-1)
    integer(4)::  m
!-----------------------------------------------------------
    call  RandomStart(310952,randu)
    do m=0,2008
      call  RandomArray(randu,a,1009)
    enddo
    write (*,*) randu(0),': Computed method-1.'
    call  RandomStart(310952,randu)
    do m=0,1008
      call  RandomArray(randu,a,2009)
    enddo
    write (*,*) randu(0),': Computed method-2.'
    write (*,*) 0.27452626307394156768d0,': Both should be this value.'
    return
  end subroutine



!-----------------------------------------------------------
!>  Random Number Generator
!>\details
!><<<PRIVATE>>> \n
!>  RandomArray(RSTAT,AA,N) obtains n random numbers in AA \n
!>    array size: real*8  RSTAT(0:Random_KK-1) \n
!>    array size: real*8  AA(0:*) \n
!>-----------------------------------------------------------\n
!>  ORIGINAL PROGRAM BY D.E.KNUTH \n
!>-----------------------------------------------------------
  subroutine    RandomArray(randu,aa,n)
    implicit    none
    integer,intent(in)::    n !< dimension of the array
    real(8),intent(out)::   aa(0:n-1) !<array in which random numbers are stored
    real(8),intent(inout)::   randu(0:Random_KK-1) !<state of random number
    integer(4)::    i,j
    do j=0,Random_KK-1
      aa(j) =randu(j)
    enddo
    do j=Random_KK,n-1
      aa(j) =Random_Mod1(aa(j-Random_KK)+aa(j-Random_LL))
    enddo
    do i=0,Random_LL-1
      randu(i)=Random_Mod1(aa(n+i-Random_KK)+aa(n+i-Random_LL))
    enddo
    do i=Random_LL,Random_KK-1
    randu(i)=Random_Mod1(aa(n+i-Random_KK)+randu(i-Random_LL))
    enddo
    return
  end subroutine



!-----------------------------------------------------------
!>  Initialize Random Generator
!>\details
!><<<PRIVATE>>> \n
!>  RandomStart(SEED,RSTAT) setup the state of random number generator \n
!>    array size: real*8  RSTAT(0:Random_KK-1) \n
!>----------------------------------------------------------- \n
!>    ORIGINAL PROGRAM BY D.E.KNUTH \n
!>----------------------------------------------------------- 
  subroutine  RandomStart(seed,randu)
    implicit  none
    integer(4),intent(in):: seed!<seed integer for randum generator
    real(8),intent(inout):: randu(0:Random_KK-1) !<state of random number
    integer(4)::  j,s,t
    real(8):: u(0:Random_KK+Random_KK-1)
    real(8):: ul(0:Random_KK+Random_KK-1)
    real(8),parameter:: ulp=(1d0/(2**30))/(2**22)
    real(8):: ss
    !------------------------------------------------------------
    ss  =2d0*ulp*(seed+2)
    do j=0,Random_KK-1
      ul(j) =0d0
      u (j) =ss
      ss  =ss+ss
      if  (ss.ge.1d0) ss=ss-(1d0-2*ulp)
    enddo
    do j=Random_KK,Random_KK+Random_KK-2
    ul(j) =0d0
    u (j) =0d0
    enddo
    u (1) =u(1)+ulp
    ul(1) =ulp
    s =seed
    t =Random_TT-1
    do while(t.ne.0)
      do j=Random_KK-1,1,-1
        ul(j+j) =ul(j)
        u (j+j) =u(j)
      enddo
      do j=Random_KK+Random_KK-2,Random_KK-Random_LL+1,-2
        ul(Random_KK+Random_KK-1-j) =0d0
        u (Random_KK+Random_KK-1-j) =u(j)-ul(j)
      enddo
      do j=Random_KK+Random_KK-2,Random_KK,-1
        if  (ul(j).ne.0d0)  then
          ul(j-(Random_KK-Random_LL)) =ulp-ul(j-(Random_KK-Random_LL))
          u (j-(Random_KK-Random_LL)) =Random_Mod1(u(j-(Random_KK-Random_LL))+u(j))
          ul(j-Random_KK)   =ulp-ul(j-Random_KK)
          u (j-Random_KK)   =Random_Mod1(u(j-Random_KK)+u(j))
        endif
      enddo
      if  (MOD(s,2).eq.1) then
        do j=Random_KK,1,-1
          ul(j) =ul(j-1)
          u (j) =u (j-1)
        enddo
        ul(0) =ul(Random_KK)
        u (0) =u (Random_KK)
        if  (ul(Random_KK).ne.0d0)  then
          ul(Random_LL) =ulp-ul(Random_LL)
          u (Random_LL) =Random_Mod1(u(Random_LL)+u(Random_KK))
        endif
      endif
      if  (s.ne.0)  then
        s =ishft(s,-1)
      else
        t =t-1
      endif
    enddo
    do j=0,Random_LL-1
      randu(j+Random_KK-Random_LL)  =u(j)
    enddo
    do j=Random_LL,Random_KK-1
      randu(j-Random_LL)    =u(j)
    enddo
  return
  end subroutine
!
!
!
!-----------------------------------------------------------
!
!       Routines to be used in user's program
!
!-----------------------------------------------------------
!
!
!-----------------------------------------------------------
!> (user's) Initialize the Random Generator with given seed
!-----------------------------------------------------------
subroutine  sGenrand(seed,data)
  implicit  none
  type(RandomProperty),intent(out)::    data !<structure of type 'RandomProperty'
  integer(4),intent(in)::   seed !<seed integer
  !-----------------------------------------------------------
  if  (seed.lt.0) then
    write (*,*) 'Seed out of range [0,2^30-3].'
    stop
  endif
  if  (seed.gt.2**30-3) then
    write (*,*) 'Seed out of range [0,2^30-3].'
    stop
  endif
  call  RandomStart(seed,data%state)
  data%index  =0
  return
end subroutine
!
!
!-----------------------------------------------------------
!> (user's) Get a random number
!-----------------------------------------------------------
real(8) function  dGenrand(random)
  implicit      none
  type(RandomProperty),intent(inout)::    random !<structure of type 'RandomProperty'
  !-----------------------------------------------------------
  dGenrand      =random%state(random%index)
  random%state(random%index)  =dGenrand+random%state(MOD(random%index+Random_KK-Random_LL,Random_KK))
  random%state(random%index)  =Random_Mod1(random%state(random%index))
  random%index      =MOD(random%index+1,Random_KK)
  return
end function
!
!-----------------------------------------------------------
!> (user's) Get an array of random number
!-----------------------------------------------------------
subroutine aGenrand(random,a,n)
  implicit      none
  type(RandomProperty),intent(inout)::    random!<structure of type 'RandomProperty'
  integer,intent(in) :: n !< dimension of array
  real(8),intent(out) :: a(n) !<array in which random numbers are stored
  real(8),allocatable :: b(:)
  integer :: nb
  !-----------------------------------------------------------
  nb=max(Random_KK,n)
  allocate(b(nb))
  call  RandomArray(random%state,b,nb)
  a(1:n)=b(1:n)
  deallocate(b)
  return
end subroutine
!
!-----------------------------------------------------------
!> (user's) Get one or two random numbers normally distributed
!-----------------------------------------------------------
subroutine  dBoxMuller(random,r1,r2)
  implicit      none
  type(RandomProperty),intent(inout)::    random!<structure of type 'RandomProperty'
  real(8),intent(out) :: r1 !<normally distributed random number 1
  real(8),intent(out),optional :: r2 !<(can be omitted) normally distributed random number 2
  real(8) :: a,k,c,t
  real(8),parameter :: PI=3.1415926535897932384626d0
  !-----------------------------------------------------------
        do
         a=dGenrand(random);k=dGenrand(random)
         if(a>0d0)exit
        end do
        c=dsqrt(-2d0*log(a))
        t=2d0*PI*k
        r1=c*dsin(t)
        if(present(r2))then
         r2=c*dcos(t)
        end if
        return
end subroutine
!-----------------------------------------------------------
!>  Put a random kind
!-----------------------------------------------------------
subroutine  WhatRandom(String)
  implicit      none
  character(LEN=*),intent(out)::    String!<kind of random number is put
  !-----------------------------------------------------------
  String  =""
  String  =Random_Kind(1:min(len(Random_Kind),len(String)))
 return
end subroutine
end module 
