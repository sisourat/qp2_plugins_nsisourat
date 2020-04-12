program stieltjes
use interpolation
implicit none

integer :: npt
integer, parameter :: QR_K = selected_real_kind (32)
real (kind=QR_K), allocatable, dimension(:) :: e, g

integer :: nmax = 30 ! according to Mueller-Plathe and Diercksen Stieltjes is inaccurate for n>=15
!real (kind=QR_K), dimension(0:nmax) :: sk
!real (kind=QR_K), dimension(nmax,nmax) ::  e1, g1, e2, g2
real (kind=QR_K), allocatable, dimension(:) :: sk
real (kind=QR_K), allocatable, dimension(:,:) ::  e1, g1, e2, g2
real (kind=QR_K) :: shift1 = 0.4d0, shift2 = 0.8d0
integer :: imax1, imax2
integer :: inmax
integer :: i, j
integer :: exit_cycle

! Tsveta
real (kind=QR_K) :: g_
real (kind=QR_K), allocatable, dimension(:) :: g1_3o, g2_3o
real (kind=QR_K), allocatable, dimension(:) :: e1_3o, e2_3o
integer :: ind
real (kind=QR_K) :: temp, pi

pi = dacos(-1d0)

! reads and sorts energy and matrix elements

read(*,*)npt
allocate(e(npt),g(npt))

do i = 1, npt
  read(*,*)e(i),g(i)
enddo
 
if(sum(g)==0d0) then
 write(*,*)"All matrix elements are zero, I stop"
 stop
endif

call sort2(npt,e,g)

open(237,file="gamma.sch2.sh0.4.nmax.dat")

do inmax = 8, nmax
   allocate(sk(0:inmax))
   allocate(e1(inmax,inmax), g1(inmax,inmax))
!   allocate(e2(inmax,inmax), g2(inmax,inmax))
   shift1 = 0.4d0   
!   shift2 = 0.8d0   
   
   ! computes the nmax moments of the matrix elements
    
    open(unit=10,file='moments.txt')
    sk(:) = 0q0
    do i = 0, inmax
       do j = 1, npt
          sk(i) = sk(i) + e(j)**i*g(j)
       enddo
     write(10,*)i,abs(sk(i))
    enddo
    close(10)
   
   ! "images" the matrix elements for two different shift values, may be used to evaluate the highest accurate order
   
    shift1 = shift1 + abs(e(1))
    e(:) = e(:) + shift1
    call imaging(npt,e,g,inmax,imax1,e1,g1)
    e(:) = e(:) - shift1
   
    do i = 1, imax1
     do  j = 1, i-1
      write(100+i,'(2(f20.16,1X))')e1(j,i)-shift1,g1(j,i)
     enddo
    enddo
   
!    shift2 = shift2 + abs(e(1))
!    e(:) = e(:) + shift2
!    call imaging(npt,e,g,inmax,imax2,e2,g2)
!    e(:) = e(:) - shift2
!   
!    do i = 1, imax2
!     do  j = 1, i-1
!      write(200+i,'(2(f20.16,1X))')e2(j,i)-shift2,g2(j,i)
!     enddo
!    enddo
   
   ! TSVETA
   ! Strategy 1 - take the 3 highest orders and interpolate
   ! linearly each of them
   !            - Gamma is the average of the 3 values
   
!    write(*,*)"The average decay width after interpolation of the three highest orders separately:"
!   
!    g_ = 0.0d0
!    do i = imax1-2,imax1
!        call interp(e1(:,i),g1(:,i),i-1,shift1,temp)
!        g_ = g_ + temp
!    end do
!    g_ = 2.0d0*pi*g_/3.0d0
!    write(*,'(1X,A9,f8.3)')"Shift 1: ",shift1-abs(e(1))
!    write(*,'(E23.15)')g_
!   
!    g_ = 0.0d0
!    do i = imax2-2,imax2
!        call interp(e2(:,i),g2(:,i),i-1,shift2,temp)
!        g_ = g_ + temp
!    end do
!    g_ = 2.0d0*pi*g_/3.0d0
!    write(*,'(1X,A9,f8.3)')"Shift 2: ",shift2-abs(e(1))
!    write(*,'(E23.15)')g_
   
   ! Strategy 2 - take the 3 highest orders and interpolate
   ! them together linearly
   !            - Gamma is the obtained value
    allocate(e1_3o(imax1-2))
    allocate(g1_3o(imax1-2))
   
    write(*,*)
    write(*,*)
    write(*,*)"The decay width after interpolation of the three highest orders together:"
    ind = 1
    do i = imax1-2,imax1
       do j = 1, i-1
         e1_3o(ind) = e1(j,i)
         g1_3o(ind) = g1(j,i)
         ind = ind + 1
       end do
    end do
   
    call interp(e1(:,imax1),g1(:,imax1),imax1-3,shift1,g_)
    g_ = g_*2.0d0*pi
    write(*,'(1X,A9,f8.3)')"Shift 1: ",shift1-abs(e(1))
    write(*,'(E23.15)')g_
    ! March 03, 2016 Tsveta
    write(237,'(I3,A3,E23.15)')inmax," ",g_
   
   
   
    deallocate(e1_3o)
    deallocate(g1_3o)
    deallocate(sk,e1,g1)
end do
close(237)
deallocate(e,g)

end program stieltjes
