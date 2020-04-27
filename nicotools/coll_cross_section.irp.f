program sigma
use gnufor2
implicit none

double precision, parameter :: pi = dacos(-1d0)

integer :: nsta
double precision, dimension(:), allocatable :: esta

integer :: nbproj
double precision, dimension(:), allocatable :: bproj, sig, pion
double precision, dimension(:,:), allocatable :: prob
double precision :: norm, sig_ion, sig_si, sig_di, sig_exc, sip, dip

integer :: istate

character(40) :: finput, csip, cdip
character(1) :: again

integer :: ilen, jlen
logical :: file_e
integer :: i, j, k, opt

call getarg(1,finput)
call getarg(2,csip)
call getarg(3,cdip)
read(csip,*)sip
read(cdip,*)dip

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
 allocate(bproj(nbproj),prob(nsta,nbproj),pion(nbproj))
 do i = 1, nbproj
   read(10,*)bproj(i),(prob(j,i),j=1,nsta),norm
   pion(i) = 1d0 - norm
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
!   write(10+i,*)bproj(j),prob(i,j)
 enddo
  sig(i) =  sig(i)*2d0*pi
!  print*, i, esta(i), sig(i)/3.57d0
enddo

sig_exc = 0d0
sig_si = 0d0
sig_di = 0d0
do i = 2, nsta
 if(esta(i)<sip) then
   sig_exc = sig_exc + sig(i)
 elseif(esta(i)<dip) then
   sig_si = sig_si + sig(i)
 else
   sig_di = sig_di + sig(i)
 endif
enddo

sig_ion= 0.5d0*bproj(1)*bproj(1)*pion(1)
do j = 1, nbproj-1
  sig_ion =  sig_ion + 0.5d0*(bproj(j+1)-bproj(j))*(bproj(j)*pion(j)+bproj(j+1)*pion(j+1))
enddo
  sig_ion =  sig_ion*2d0*pi
do i = 1, nsta
 write(*,'(i5,1X,10(f16.9,1X))'), i, esta(i), sig(i)/3.57d0
enddo
 write(*,*), 
! write(*,*)'Ionization = ', sig_ion/3.57d0
 write(*,*)'Exc, SI, DI = ', sig_exc/3.57d0, sig_si/3.57d0, sig_di/3.57d0

deallocate(esta,bproj,prob,sig)
end


