program WriteTOP
implicit none

character(len=79)ARGUMENT,TEMPLATE,OUTPUT,PAIRLIST,PDB,BONDLIST,ANGLELIST,DIHLIST
character(len=9) tmp1(5),atomtype(10000),restype(1000)
integer i,j,k,ios,ios1,ios2,nargs
real R_repulsive,tmp(9),e_LJ,c6,c12,e_Repul,e_P,e_R,e_i
real e_DH,epsilon,K_Coulomb,C_Salt,Kai
parameter(epsilon=80.0,K_Coulomb=138.94,C_Salt=0.15,Kai=1.2394)

nargs=iargc()
   if(nargs==13)then
       call getarg(1,TEMPLATE)
       call getarg(2,OUTPUT)
       call getarg(3,PAIRLIST)
       call getarg(4,PDB)
       call getarg(5,BONDLIST)
       call getarg(6,ANGLELIST)
       call getarg(7,DIHLIST)
       call getarg(8,ARGUMENT)
          read(ARGUMENT,"(F)")R_repulsive
       call getarg(9,ARGUMENT)
          read(ARGUMENT,"(F)")e_LJ
 !         e_DH=e_LJ/(K_Coulomb*exp(-Kai*0.5)/(80*0.5))
          e_DH=1.0
       call getarg(10,ARGUMENT)
          read(ARGUMENT,"(F)")e_Repul
       call getarg(11,ARGUMENT)
          read(ARGUMENT,"(F)")e_P
       call getarg(12,ARGUMENT)
          read(ARGUMENT,"(F)")e_R
       call getarg(13,ARGUMENT)
          read(ARGUMENT,"(F)")e_i
   else
       write(*,*) ""
       write(*,*) "  Wrong Input!"
       write(*,*) "./WriteTOP.x TEMPLATE  OUTPUT PAIRLIST  CGPDB BONDLIST ANGLELIST DIHLIST  R_repulsive(angstrom) e_LJ  e_Repulsive e_Protein  e_RNA  e_inter"
       write(*,*) "./writeTOP.x 1025_22.top  1019_x.top PAIR-3.list  forTOPca.pdb Bonds.dat Angles.dat dihs.dat  4.0  1.0  1.0  1.0 2.0 0.5"
       write(*,*) ""
       goto 100
   endif

ios=0
open(1,file=TEMPLATE,status="old")
open(2,file=OUTPUT,status="replace")
open(3,file=PAIRLIST,status="old")
do while(.true.)
	read(1,"(A)",iostat=ios)ARGUMENT
	if(ios/=0)exit
	if(ARGUMENT(2:14)=="[ atomtypes ]")then
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(ARGUMENT(37:51),"(E15.9)") (R_repulsive/10)**12*e_Repul
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
	elseif(ARGUMENT(2:10)=="[ atoms ]")then
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
open(4,file=PDB,status="old")
		ios2=0
		i=0
		do while(.true.)
			read(4,"(A)",iostat=ios2)ARGUMENT
			if(ios2/=0)exit
			i=i+1
			read(ARGUMENT(9:11),*) tmp1(1)
			read(ARGUMENT(14:16),"(A3)") tmp1(2)
			atomtype(i)=tmp1(2)
			read(ARGUMENT(18:20),"(A3)") tmp1(3)
			restype(i)=tmp1(3)
			read(ARGUMENT(24:26),*) tmp1(4)
			if(tmp1(3)=="ARG" .or. tmp1(3)=="LYS")then
			write(2,"(3X,A3,2X,A3,3X,A3,1X,A3,2X,A3,5X,A3,3X,F5.3,3X,F5.3)") tmp1(1),tmp1(2),tmp1(4),tmp1(3),tmp1(2),tmp1(1),sqrt(e_DH),1.000
			elseif(tmp1(3)=="GLU" .or. tmp1(3)=="ASP" .or. tmp1(2)=="CP")then
			write(2,"(3X,A3,2X,A3,3X,A3,1X,A3,2X,A3,5X,A3,2X,F6.3,3X,F5.3)") tmp1(1),tmp1(2),tmp1(4),tmp1(3),tmp1(2),tmp1(1),-sqrt(e_DH),1.000
			else
			write(2,"(3X,A3,2X,A3,3X,A3,1X,A3,2X,A3,5X,A3,3X,F5.3,3X,F5.3)") tmp1(1),tmp1(2),tmp1(4),tmp1(3),tmp1(2),tmp1(1),0.000,1.000
			endif
		enddo
close(4)
		do while(.true.)
			read(1,"(A)")ARGUMENT
			if(ARGUMENT=="")then
			write(2,"(A)")
			exit
			endif
		enddo
	elseif(ARGUMENT(2:10)=="[ pairs ]")then
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
!!!!!!! Contact Map!!!!!!!!!!!
		ios1=0
		do while(.true.)
			read(3,*,iostat=ios1)(tmp(k),k=1,4)
			if(ios1/=0)exit
			if(tmp(1)<72 .and. tmp(2)<72)then
			c6=6*(tmp(3)/10)**10*e_LJ*e_P*tmp(4)
			c12=5*(tmp(3)/10)**12*e_LJ*e_P*tmp(4)
			elseif(tmp(1)>71 .and. tmp(2)>71)then
			c6=6*(tmp(3)/10)**10*e_LJ*e_R
			c12=5*(tmp(3)/10)**12*e_LJ*e_R
			else
			c6=6*(tmp(3)/10)**10*e_LJ*e_i
			c12=5*(tmp(3)/10)**12*e_LJ*e_i
			endif
			select case(restype(INT(ANINT(tmp(1)))))
			case("ARG","LYS")
			 select case(restype(INT(ANINT(tmp(2)))))
			  case("ASP","GLU")
			  c6=0.1*c6
			  c12=0.1*c12
			 end select
			 select case(atomtype(INT(ANINT(tmp(2)))))
			  case("CP")
			  c6=0.1*c6
			  c12=0.1*c12
			 end select
			case("ASP","GLU")
			 select case(restype(INT(ANINT(tmp(2)))))
			  case("ARG","LYS")
			  c6=0.1*c6
			  c12=0.1*c12
			 end select
			end select
			write(2,"(2X,I5,2X,I5,2X,I1,3X,E15.9,3X,E15.9)") INT(ANINT(tmp(1))),INT(ANINT(tmp(2))),1,c6,c12
		enddo
close(3)
		do while(.true.)
			read(1,"(A)")ARGUMENT
			if(ARGUMENT=="")then
			write(2,"(A)")
			exit
			endif
		enddo
	elseif(ARGUMENT(2:10)=="[ bonds ]")then
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
!		do while(.true.)
!			read(1,"(A)")ARGUMENT
!			if(ARGUMENT=="")exit
!			write(2,"(A)")ARGUMENT
!		enddo
!		write(2,"(A)")
open(5,file=BONDLIST,status="old")
		ios2=0
		do while(.true.)
			read(5,*,iostat=ios2)(tmp(k),k=1,5)
			if(ios2/=0)exit
			write(2,"(2X,I5,3X,I5,2X,I1,3X,E15.9,3X,E15.9)") INT(ANINT(tmp(1))),INT(ANINT(tmp(2))),1,tmp(4),tmp(5)
		enddo
close(5)
		do while(.true.)
			read(1,"(A)")ARGUMENT
			if(ARGUMENT=="")then
			write(2,"(A)")
			exit
			endif
		enddo
	elseif(ARGUMENT(2:15)=="[ exclusions ]")then
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
!!!!!!!!!!exclusions!!!!!!!!!!!!
open(3,file=PAIRLIST,status="old")
		ios2=0
		do while(.true.)
			read(3,*,iostat=ios2)(tmp(k),k=1,4)
			if(ios2/=0)exit
			write(2,"(2X,I5,3X,I5)") INT(ANINT(tmp(1))),INT(ANINT(tmp(2)))
		enddo
close(3)
		do while(.true.)
			read(1,"(A)")ARGUMENT
			if(ARGUMENT=="")then
			write(2,"(A)")
			exit
			endif
		enddo
!	elseif(ARGUMENT(2:11)=="[ angles ]")then
!		write(2,"(A)")ARGUMENT
!		read(1,"(A)")ARGUMENT
!		write(2,"(A)")ARGUMENT
!		do while(.true.)
!			read(1,"(A)")ARGUMENT
!			if(ARGUMENT=="")then
!			write(2,"(A)")
!			exit
!			endif
!		enddo
	elseif(ARGUMENT(2:11)=="[ angles ]")then
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
!		do while(.true.)
!			read(1,"(A)")ARGUMENT
!			if(ARGUMENT=="")exit
!			write(2,"(A)")ARGUMENT
!		enddo
!		write(2,"(A)")
open(6,file=ANGLELIST,status="old")
		ios2=0
		do while(.true.)
			read(6,*,iostat=ios2)(tmp(k),k=1,6)
			if(ios2/=0)exit
			write(2,"(2X,I5,3X,I5,3X,I5,2X,I1,3X,E15.9,3X,E15.9)") INT(ANINT(tmp(1))),INT(ANINT(tmp(2))),INT(ANINT(tmp(3))),1,tmp(5),tmp(6)
		enddo
close(6)
		do while(.true.)
			read(1,"(A)")ARGUMENT
			if(ARGUMENT=="")then
			write(2,"(A)")
			exit
			endif
		enddo
	elseif(ARGUMENT(2:14)=="[ dihedrals ]")then
		write(2,"(A)")ARGUMENT
		read(1,"(A)")ARGUMENT
		write(2,"(A)")ARGUMENT
!		do while(.true.)
!			read(1,"(A)")ARGUMENT
!			if(ARGUMENT=="")exit
!			write(2,"(A)")ARGUMENT
!		enddo
!		write(2,"(A)")
open(7,file=DIHLIST,status="old")
		ios2=0
		do while(.true.)
			read(7,*,iostat=ios2)(tmp(k),k=1,8)
			if(ios2/=0)exit
			write(2,"(2X,I5,3X,I5,3X,I5,3X,I5,2X,I1,3X,E15.9,3X,E15.9,1X,I1)") INT(ANINT(tmp(1))),INT(ANINT(tmp(2))),INT(ANINT(tmp(3))),INT(ANINT(tmp(4))),1,tmp(6),tmp(7),INT(ANINT(tmp(8)))
		enddo
close(7)
		do while(.true.)
			read(1,"(A)")ARGUMENT
			if(ARGUMENT=="")then
			write(2,"(A)")
			exit
			endif
		enddo
	else
	write(2,"(A)")ARGUMENT
	endif
enddo

close(1)
close(2)
close(3)
100 end program
