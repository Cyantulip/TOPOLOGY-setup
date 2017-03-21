program RNABonds
implicit none

character(len=79)ARGUMENT
integer i,j,k,ios

open(1,file="bondRNA.list",status="replace")
do i=71,196
       if(mod(i,3)==2)then
     write(1,"(2X,I4,3X,I4/2X,I4,3X,I4/2X,I4,3X,I4)") i,i+1,i+1,i+2,i+1,i+3
       endif
enddo
close(1)
end program
