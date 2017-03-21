program WriteBonds
implicit none

character(len=79)ARGUMENT
integer i,k,j,ios
integer,parameter :: first=196
integer c1(first),c2(first),bondNum
real x(first),y(first),z(first),dist

!open(1,file="bonds.list",status="replace")
!do i=1,first-1
!  write(1,"(2X,I4,3X,I4)") i,i+1
!enddo
!close(1)
ios=0
i=0
open(10,file="BOND.list",status="old")
do while(.true.)
read(10,"(A)",iostat=ios)ARGUMENT
if(ios/=0)exit
i=i+1
read(ARGUMENT(3:6),"(I4)") c1(i)
read(ARGUMENT(10:13),"(I4)") c2(i)
enddo
bondNum=i
close(10)
write(*,*) bondNum

ios=0
open(2,file="ca.pdb",status="old")
open(3,file="Bonds.dat",status="replace")
do i=1,first
  read(2,"(A)",iostat=ios)ARGUMENT
  if(ios/=0)exit
  read(ARGUMENT(31:38),"(F8.3)") x(i)
  read(ARGUMENT(39:46),"(F8.3)") y(i)
  read(ARGUMENT(47:54),"(F8.3)") z(i)
enddo
do j=1,bondNum
  dist=sqrt((x(c1(j))-x(c2(j)))**2+(y(c1(j))-y(c2(j)))**2+(z(c1(j))-z(c2(j)))**2)
  write(3,"(2X,I4,2X,I4,2X,I1,3X,E15.9,3X,E15.9)") c1(j),c2(j),1,dist/10,0.200000000E+05
enddo
close(2)
close(3)


end program
