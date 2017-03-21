program main
implicit none

 character(len=79)buffer
integer i,j,k,ios
integer,parameter :: first=230
integer :: c(first)
integer :: c1(first)
integer :: c2(first)
integer :: c3(first)
integer :: c4(first)
real :: x(first)

ios=0
open(10,file="anoutput",status="old")
do i=1,first
 read(10,"(A)",iostat=ios)buffer
 if(ios/=0)exit
 read(buffer(8:10),"(I3)") c(i)
 read(buffer(24:30),"(F7.3)") x(i)
end do
close(10)

ios=0
open(20,file="DIH.list",status="old")
do i=1,first
 read(20,"(A)",iostat=ios)buffer
 if(ios/=0)exit
 read(buffer(10:12),"(I3)") c1(i)
 read(buffer(22:24),"(I3)") c2(i)
 read(buffer(34:36),"(I3)") c3(i)
 read(buffer(46:48),"(I3)") c4(i)
end do
close(20)

open(30,file="dihs.dat",status="replace")
do i=1,first
 write(30,"(3X,I3,4X,I3,4X,I3,4X,I3,1X,I1,3X,E15.9,3X,E15.9,1X,I1/3X,I3,4X,I3,4X,I3,4X,I3,1X,I1,3X,E15.9,3X,E15.9,1X,I1)") c1(i),c2(i),c3(i),c4(i),1,x(i)+180.0,0.100000000E+01,1,c1(i),c2(i),c3(i),c4(i),1,3*x(i)+540.0,0.500000000E+00,3
end do
close(30)

!open(40,file="dih_apo-improper.dat",status="replace")
!do i=1,first
! write(40,"(3X,I3,4X,I3,4X,I3,4X,I3,1X,I1,3X,E15.9,3X,E15.9)") c1(i),c2(i),c3(i),c4(i),2,x(i),0.100000000E+02
!end do
!close(40)

end program
