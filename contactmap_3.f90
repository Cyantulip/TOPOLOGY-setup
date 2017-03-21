program Contactmap
implicit none

character(len=79)ARGUMENT,ARGUMENT1,PDB,PDB1
integer,parameter :: first=40000!!arrayDimension for atomCoordinates(big enough)
integer,parameter :: second=1000!!arrayDimension for residueCooidinates
character(len=5) atomtype(first),restype(first),chain(first)
integer i,j,k,ios,nargs,m,resNum(first),atomNum,ResidueNumber,tmp
real :: x_Protein(first),y_Protein(first),z_Protein(first),x_Residue(second),y_Residue(second),z_Residue(second),Rcut,AngleCutoff
real dist_ji,dist_ki,dist_ij,dist_kj,dist_Res,angle_jik,angle_ijk,inner_jik,inner_ijk
real :: pi=3.14159265358979

nargs=iargc()
  if(nargs==4)then
     call getarg(1,PDB)
     call getarg(2,PDB1)
     call getarg(3,ARGUMENT)
       read(ARGUMENT,"(F)") Rcut
     call getarg(4,ARGUMENT1)
       read(ARGUMENT1,"(F)") AngleCutoff
  else
    write(*,*) "Wrong Input!"
    write(*,*) "./Contactmap.x    allatomPDBfile  CGPDBfile    Rcut     AngleCutoff"
    write(*,*) "./Contactmap.x    ../complexs.pdb   ca.pdb      9.0             30.0 "
    goto 100
  endif

ios=0
i=0
open(10,file=PDB,status="old")
do while(.true.)
 read(10,"(A)",iostat=ios) ARGUMENT
 if(ios/=0)exit
 if(ARGUMENT(1:4)=="ATOM")then
  i=i+1
  read(ARGUMENT(14:16),"(A3)") atomtype(i)
  read(ARGUMENT(18:20),"(A3)") restype(i)
  read(ARGUMENT(22:22),"(A1)") chain(i)
  read(ARGUMENT(23:26),"(I4)") resNum(i)
!  resNum(i)=tmp-2
  read(ARGUMENT(31:38),"(F8.3)") x_Protein(i)
  read(ARGUMENT(39:46),"(F8.3)") y_Protein(i)
  read(ARGUMENT(47:54),"(F8.3)") z_Protein(i)
 endif
enddo
  atomNum=i
   write(*,*)"There are",atomNum,"  atoms!"
close(10)
	
ios=0
i=0
open(100,file=PDB1,status="old")
do while(.true.)
 read(100,"(A)",iostat=ios) ARGUMENT
 if(ios/=0)exit
 if(ARGUMENT(14:14)=="C")then
  i=i+1
  read(ARGUMENT(31:38),"(F8.3)") x_Residue(i)
  read(ARGUMENT(39:46),"(F8.3)") y_Residue(i)
  read(ARGUMENT(47:54),"(F8.3)") z_Residue(i)
 endif
enddo
  ResidueNumber=i
  write(*,*)"There are",i,"  Residues!"
close(100)
	
m=0
open(30,file="output",status="replace")
do j=1,atomNum
  do i=j+1,atomNum
	if(resNum(i)>resNum(j)+3)then!!!!!!!!!!!i>j+3 or i>j+2!!!!!!!!!!!!!!!!!!!!!!
	dist_ji=sqrt((x_Protein(j)-x_Protein(i))**2+(y_Protein(j)-y_Protein(i))**2+(z_Protein(j)-z_Protein(i))**2)
	if(dist_ji<Rcut .and. chain(i)/=chain(j))then
       do k=1,atomNum
		if(k/=i .and. k/=j)then
!		if(k/=i .and. k/=j .and. chain(k)==chain(j))then
            dist_ki=sqrt((x_Protein(k)-x_Protein(i))**2+(y_Protein(k)-y_Protein(i))**2+(z_Protein(k)-z_Protein(i))**2)
            inner_jik=(x_Protein(j)-x_Protein(i))*(x_Protein(k)-x_Protein(i))+(y_Protein(j)-y_Protein(i))*(y_Protein(k)-y_Protein(i))+(z_Protein(j)-z_Protein(i))*(z_Protein(k)-z_Protein(i))
            angle_jik=acos(inner_jik/(dist_ji*dist_ki))
            angle_jik=(angle_jik*180)/pi
 						
            dist_kj=sqrt((x_Protein(k)-x_Protein(j))**2+(y_Protein(k)-y_Protein(j))**2+(z_Protein(k)-z_Protein(j))**2)
            inner_ijk=(x_Protein(i)-x_Protein(j))*(x_Protein(k)-x_Protein(j))+(y_Protein(i)-y_Protein(j))*(y_Protein(k)-y_Protein(j))+(z_Protein(i)-z_Protein(j))*(z_Protein(k)-z_Protein(j))
            angle_ijk=acos(inner_ijk/(dist_ji*dist_kj))
            angle_ijk=(angle_ijk*180)/pi
 						

			if((dist_ki<=dist_ji .and. angle_jik<=AngleCutoff) .or. (dist_kj<=dist_ji .and. angle_ijk<=AngleCutoff))then
			m=1
			endif
		endif
       enddo
		if(m==0)then
            dist_Res=sqrt((x_Residue(resNum(j))-x_Residue(resNum(i)))**2+(y_Residue(resNum(j))-y_Residue(resNum(i)))**2+(z_Residue(resNum(j))-z_Residue(resNum(i)))**2)
                !  write(30,"(I5,2X,I5,3X,I5,2X,I5,3X,F10.7)")j,i,resNum(j),resNum(i),dist_ji
                  write(30,"(I5,2X,I5,3X,F10.7)")resNum(j),resNum(i),dist_Res
                !  write(*,*) angle_jik,angle_ijk
		else
		m=0
		endif
	endif
	endif
   enddo
enddo
close(30)

100 end program
