program MJpotential
implicit none

character(len=119)ARGUMENT,tmp(12)
integer i,j,k,ios,nargs,m,n,pairsNumCACA,pairsNum,p1,p2
integer resNum(182)
character(len=20) :: resType(182)
real :: tmp1(20),matrix(20,20),e10,e12,e_total
character(len=20) :: Res(20),Input,Output
real,allocatable :: e_MJ(:)

nargs=iargc()
 if(nargs==1)then 
 call getarg(1,Input)
   Output=Input
else
  write(*,*) "Wrong Input!"
  write(*,*) "./MJ.x pairs.list"
 goto 1000
endif
 do i=1,11
 if(Output(i:i+4)==".list")then
  Output(i:i+6)=".listmj"
 goto 200
 endif
 enddo


200 ios=0
i=0
open(11,file="proteinCA.pdb",status="old")
do while(.true.)
 read(11,"(A)",iostat=ios)ARGUMENT
 if(ios/=0)exit
 i=i+1
 read(ARGUMENT(9:11),"(I3)") resNum(i)
 read(ARGUMENT(18:20),"(A3)") resType(i)
enddo
close(11)
  write(*,*) "Residue Number is :",i

ios=0
i=0
open(12,file="index.txt",status="old")
do while(.true.)
 read(12,"(A)",iostat=ios)ARGUMENT
 if(ios/=0)exit
 i=i+1
 read(ARGUMENT(1:3),"(A3)") Res(i)
enddo
close(12)
  write(*,*) "Index Number is :",i

ios=0
open(13,file="mj1.txt",status="old")
do i=1,20
 read(13,"(A)",iostat=ios)ARGUMENT
 if(ios/=0)exit
 read(ARGUMENT,*) (tmp1(k),k=1,20)
 do j=1,20
    matrix(i,j)=tmp1(j)
!  read(tmp1(j),"(F5.2)") matrix(i,j)
 enddo
enddo
close(13)
  write(*,*) matrix(20,1)

ios=0
i=0
j=0
open(14,file=Input,status="old")
do while(.true.)
 read(14,"(A)",iostat=ios)ARGUMENT
 if(ios/=0)exit
  i=i+1
 read(ARGUMENT,*) (tmp(k),k=1,3)
! if(tmp(8)=="CA-CA")then
  j=j+1
! endif
enddo
allocate(e_MJ(j))
    pairsNum=i
    pairsNumCACA=j
  write(*,*) "There are ",pairsNum,"pairs!"
  write(*,*) "There are ",pairsNumCACA,"CA-CA pairs!"
rewind(14)
ios=0
i=0
do while(.true.)
 read(14,"(A)",iostat=ios)ARGUMENT
 if(ios/=0)exit
 read(ARGUMENT,*) (tmp(k),k=1,3)
! if(tmp(8)=="CA-CA")then
  i=i+1
 read(tmp(1),"(I4)") p1
 read(tmp(2),"(I4)") p2
 do m=1,20
  do n=1,20
  if(Res(m)==resType(p1) .and. Res(n)==resType(p2))then
    e_MJ(i)=matrix(m,n)
    e_total=e_total+e_MJ(i)
  endif
  enddo
 enddo
! endif
enddo

  write(*,*) "e_total is",e_total
rewind(14)
open(15,file=Output,status="replace")
ios=0
i=0
do while(.true.)
 read(14,"(A)",iostat=ios)ARGUMENT
 if(ios/=0)exit
 read(ARGUMENT,*) (tmp(k),k=1,3)
! if(tmp(8)=="CA-CA")then
  i=i+1
! read(tmp(4),"(E15.9)") e10
! read(tmp(5),"(E15.9)") e12
! write(ARGUMENT(15:29),"(E15.9)") e10*(e_MJ(i)*pairsNumCACA/e_total)
! write(ARGUMENT(33:47),"(E15.9)") e12*(e_MJ(i)*pairsNumCACA/e_total)
 write(ARGUMENT(29:38),"(F10.3)") e_MJ(i)*pairsNumCACA/e_total
! endif
 write(15,"(A50)") ARGUMENT
enddo
close(14)
close(15)
1000 end program
