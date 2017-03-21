program main
implicit none

character(len=79)buffer1
character(len=79)pdbfile,anglelist,outputfile,buffer
integer i,j,k,ios,nargs,N
real,allocatable :: x(:),y(:),z(:)
real :: temp(3)
integer :: temp1(3)
real dist,dist1,angle,inner
real :: pi=3.14159265358979

nargs=iargc()
 if(nargs==4)then
   call getarg(1,pdbfile)
   call getarg(2,anglelist)
   call getarg(3,outputfile)
   call getarg(4,buffer)
   read(buffer,"(I)") N
   allocate(x(N))
   allocate(y(N))
   allocate(z(N))
 else
   write(*,*) "Wrong Input!"
   write(*,*) "Usage: pdbfile anglelist outputfile beadnum"
  goto 1000
 end if

buffer1="                                                                               "
ios=0
open(10,file=pdbfile,status="old")
do i=1,N
 read(10,"(A)",iostat=ios) buffer
  if(ios/=0)exit
! read(buffer(32:54),"(2(F7.3,1X),F7.3)") (temp(k),k=1,3)
   read(buffer(31:38),"(F8.3)") x(i)
   read(buffer(39:46),"(F8.3)") y(i)
   read(buffer(47:54),"(F8.3)") z(i)
!    x(i)=temp(1)
!    y(i)=temp(2)
!    z(i)=temp(3)
end do
close(10)

ios=0
open(20,file=anglelist,status="old")
open(30,file=outputfile,status="replace")
 do while(.true.)
 read(20,"(A)",iostat=ios)buffer
  if(ios/=0)exit
 read(buffer,*) (temp1(j),j=1,3)
 dist=sqrt((x(temp1(1))-x(temp1(2)))**2+(y(temp1(1))-y(temp1(2)))**2+(z(temp1(1))-z(temp1(2)))**2)
 dist1=sqrt((x(temp1(3))-x(temp1(2)))**2+(y(temp1(3))-y(temp1(2)))**2+(z(temp1(3))-z(temp1(2)))**2)
 inner=(x(temp1(1))-x(temp1(2)))*(x(temp1(3))-x(temp1(2)))+(y(temp1(1))-y(temp1(2)))*(y(temp1(3))-y(temp1(2)))+(z(temp1(1))-z(temp1(2)))*(z(temp1(3))-z(temp1(2)))
 angle=acos(inner/(dist*dist1))
 angle=(angle*180)/pi
 write(buffer1(3:6),"(I4)") temp1(1)
 write(buffer1(10:13),"(I4)") temp1(2)
 write(buffer1(17:20),"(I4)") temp1(3)
 buffer1(22:22)="1"
 write(buffer1(26:40),"(E15.9)") angle
 buffer1(44:58)="0.400000000E+02"
 write(30,"(A)") buffer1
end do
close(20)
close(30)

1000 end program
