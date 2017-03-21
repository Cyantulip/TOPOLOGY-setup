!======================================================
! Yong Wang, wyongciac@gmail.com
! Reference document:
! https://groups.google.com/forum/#!topic/bionet.xtallography/LIBWXBm8hy4
! http://math.stackexchange.com/questions/47059/how-do-i-calculate-a-dihedral-angle-given-cartesian-coordinates
!"Given the coordinates of the four points, obtain the vectors b1, b2, and b3 by vector subtraction.
!Let me use the nonstandard notation ⟨v⟩ to denote v/∥v∥, the unit vector in the direction of the vector v. 
!Compute n1=⟨b1×b2⟩ and n2=⟨b2×b3⟩, the normal vectors to the planes containing b1 and b2, and b2 and b3 respectively. 
!The angle we seek is the same as the angle between n1 and n2.
!The three vectors n1, ⟨b2⟩, and m1:=n1×⟨b2⟩ form an orthonormal frame. 
!Compute the coordinates of n2 in this frame: x=n1⋅n2 and y=m1⋅n2. 
!(You don't need to compute ⟨b2⟩⋅n2 as it should always be zero.)
!The dihedral angle, with the correct sign, is atan2(y,x).
!(The reason I recommend the two-argument atan2 function to the traditional cos−1 in this case 
!is both because it naturally produces an angle over a range of 2π, 
!and because cos−1 is poorly conditioned when the angle is close to 0 or ±π.)"
!======================================================

        program Dihedral_Angle

        implicit none
        integer :: i,j,resi,atomi,nargs,ndih,ios, I3, J3, K3, L3, JN, COL1, COL2,COL3,COL4,COL5
        character(80) :: line,ctmp(10),arg,pdbfile,listfile,HEADER, PDB, PDBp
	character(200) :: fmat1,fmat2
        real,dimension(8) :: tmp
        character(5) :: atomtype,restype, atomt(5000), rest(5000)
        real,dimension(5000,3) :: Coor
        integer,dimension(5000) :: dihi, dihj, dihk, dihl
        real :: pi=3.141592653589793,zero=0.d0,one=1.d0,tenm3=1.0d-3,TM24=1.0d-24
        real :: AP1,AP0,Z10,Z20,Z12,S,CT0,CT1,CT2
        real,dimension(5000) :: XT,YT,ZT,Z1,Z2,FXI,FYI,XIJ,YIJ,ZIJ,XKJ,YKJ,ZKJ,XKL,YKL,ZKL,DX,DY,DZ,GX,GY,GZ,CT

        ! find out how many command line arguments there are
        nargs = iargc()
        if(nargs /= 2)then
                write(6,*) "the number of argument is incorrect"
                write(6,*) "usage: Dihedral_Angle.x CACB.pdb dih.list"
                goto 1000
        else if (nargs==2) then
                call getarg(1,arg)
                pdbfile = arg
                call getarg(2,arg)
                listfile=arg
        endif

        !==========================================
        ! read the coordinate from PDB file
        !==========================================
        PDB="(A6,I5,TR1,A4,TR1,A3,TR2,I4,TR4,3F8.3)"
        PDBp="(A6,I5,TR1,A4,TR1,A3,TR1,A1,I4,TR4,3F8.3)"
	fmat1="(I6,1X,I4,1X,A4,1X,A3,3(F7.3,1X))"
        open(20,file=pdbfile,status='old')
        readpdb: do
                read(20,('(A70)'),iostat=ios) line
                if(ios .ne. 0) exit readpdb
                read(line,*) HEADER

                if ( HEADER == 'ATOM' ) then
			!ATOM      7  SD  MET A   1      23.930  23.959   5.904  1.00 17.17           S
	                read(line,PDB) Col1,atomi,atomtype,restype,resi,(Coor(atomi,i),i=1,3)
			atomt(atomi)=atomtype
			rest(atomi)=restype
                        !write(*,fmat1) atomi,resi,atomtype,restype,(Coor(atomi,i),i=1,3)
                endif
        enddo readpdb

        !==========================================
        ! read the definition of dihedral angle
        !==========================================
        open(30,file=listfile,status='old')
        ndih=0
        readlist: do
                read(30,('(A70)'),iostat=ios) line
                if(ios .ne. 0) exit readlist
                read(line,*) (tmp(i),i=1,4)
                ndih=ndih+1
                dihi(ndih)=tmp(1)
                dihj(ndih)=tmp(2)
                dihk(ndih)=tmp(3)
                dihl(ndih)=tmp(4)
                write(*,*) ndih, dihi(ndih),dihj(ndih),dihk(ndih),dihl(ndih)
        enddo readlist

        !==========================================
        ! calculate the dihedral angle
        !==========================================
        do i=1,ndih
                !----- CALCULATION OF ij, kj, kl VECTORS -----
                I3=dihi(i)
                J3=dihj(i)
                K3=dihk(i)
                L3=dihl(i)
                XT(I3)=Coor(I3,1)
                YT(I3)=Coor(I3,2)
                ZT(I3)=Coor(I3,3)
                XT(J3)=Coor(J3,1)
                YT(J3)=Coor(J3,2)
                ZT(J3)=Coor(J3,3)
                XT(K3)=Coor(K3,1)
                YT(K3)=Coor(K3,2)
                ZT(K3)=Coor(K3,3)
                XT(L3)=Coor(L3,1)
                YT(L3)=Coor(L3,2)
                ZT(L3)=Coor(L3,3)

                XIJ(i) = XT(I3)-XT(J3)
                YIJ(i) = YT(I3)-YT(J3)
                ZIJ(i) = ZT(I3)-ZT(J3)
                XKJ(i) = XT(K3)-XT(J3)
                YKJ(i) = YT(K3)-YT(J3)
                ZKJ(i) = ZT(K3)-ZT(J3)
                XKL(i) = XT(K3)-XT(L3)
                YKL(i) = YT(K3)-YT(L3)
                ZKL(i) = ZT(K3)-ZT(L3)
        enddo

        !----- GET THE NORMAL VECTOR -----
        open(300,file="anoutput",status="replace")
        do i=1,ndih
                DX(i) = YIJ(i)*ZKJ(i)-ZIJ(i)*YKJ(i)
                DY(i) = ZIJ(i)*XKJ(i)-XIJ(i)*ZKJ(i)
                DZ(i) = XIJ(i)*YKJ(i)-YIJ(i)*XKJ(i)
                GX(i) = ZKJ(i)*YKL(i)-YKJ(i)*ZKL(i)
                GY(i) = XKJ(i)*ZKL(i)-ZKJ(i)*XKL(i)
                GZ(i) = YKJ(i)*XKL(i)-XKJ(i)*YKL(i)
                FXI(i) = SQRT(DX(i)*DX(i)+DY(i)*DY(i)+DZ(i)*DZ(i)+TM24)
                FYI(i) = SQRT(GX(i)*GX(i)+GY(i)*GY(i)+GZ(i)*GZ(i)+TM24)
                CT(i) = DX(i)*GX(i)+DY(i)*GY(i)+DZ(i)*GZ(i)
                Z10 = one/FXI(i)
                Z20 = one/FYI(i)
                if (tenm3 .gt. FXI(i)) Z10 = zero
                if (tenm3 .gt. FYI(i)) Z20 = zero

                Z12 = Z10*Z20
                Z1(i) = Z10
                Z2(i) = Z20

                CT0 = MIN(one,CT(i)*Z12)
                CT1 = MAX(-one,CT0)
                S = XKJ(i)*(DZ(i)*GY(i)-DY(i)*GZ(i))+YKJ(i)*(DX(i)*GZ(i)-DZ(i)*GX(i))+ZKJ(i)*(DY(i)*GX(i)-DX(i)*GY(i))
                AP0 = ACOS(CT1)
                AP1 = PI-SIGN(AP0,S)

		! output the dihedral angles
                !write(*,'(6(F9.3,2X),E15.9,1X,E15.9)') AP0,CT0,CT1,S,AP1,AP1*180/PI,180+AP1*180/PI,3*(180+(AP1*180/PI))
                I3=dihi(i)
                J3=dihj(i)
                K3=dihk(i)
                L3=dihl(i)
		fmat2="(A5,1X,I4,1X,F9.3,1X,F9.3,4(1X,A4,1X,A3))"
         write(*,fmat2) "Dih: ", i, AP1, AP1*180/PI, atomt(I3),rest(I3),atomt(J3),rest(J3), atomt(K3),rest(K3), atomt(L3),rest(L3)
         write(300,fmat2) "Dih: ", i, AP1, AP1*180/PI, atomt(I3),rest(I3),atomt(J3),rest(J3), atomt(K3),rest(K3), atomt(L3),rest(L3)
         enddo

close(300)

1000    end program Dihedral_Angle
