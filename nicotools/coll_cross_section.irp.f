program sigma
use gnufor2
implicit none

double precision, parameter :: pi = dacos(-1d0)

integer :: nsta
double precision, dimension(:), allocatable :: esta

integer :: nbproj
double precision, dimension(:), allocatable :: bproj, sig
double precision, dimension(:,:), allocatable :: prob
double precision :: norm

integer :: istate

character(40) :: finput
character(1) :: again

integer :: ilen, jlen
logical :: file_e
integer :: i, j, k, opt

call getarg(1,finput)

ilen=index(finput,' ')
inquire( file="./"//finput(1:ilen-1), exist=file_e )
if ( file_e .eqv. .false. ) then
 write(*,*) finput(1:ilen-1), " does not exist"
 stop
endif

open(unit=10,file=finput)
 read(10,*)nbproj,nsta
 allocate(esta(nsta),sig(nsta))
 do i = 1, nsta
   read(10,*)esta(i)
 enddo
 allocate(bproj(nbproj),prob(nsta,nbproj))
 do i = 1, nbproj
   read(10,*)bproj(i),(prob(j,i),j=1,nsta),norm
   if(norm<0.9d0 .or. norm> 1.1d0) print*, "WARNING: NORM IS NOT CONSERVED ==>",norm
 enddo
close(10)

!! computes the cross sections

 print*, "Cross sections in 10^-16 cm^2"
 sig(:) = 0d0
do i = 1, nsta
   sig(i) = 0.5d0*bproj(1)*bproj(1)*prob(i,1)
 do j = 1, nbproj-1
   sig(i) =  sig(i) + 0.5d0*(bproj(j+1)-bproj(j))*(bproj(j)*prob(i,j)+bproj(j+1)*prob(i,j+1))
 enddo
  sig(i) =  sig(i)*2d0*pi
  print*, i, esta(i), sig(i)/3.57d0
enddo

 print*, sum(sig(6:nsta))/3.57d0

deallocate(esta,bproj,prob,sig)
end


