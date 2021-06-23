program cippres_stieltjes
  use general
  use interpolation
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
! cippres_stieltjes computes the decay widths/photoionization cross sections using Stieltjes imaging and the Fano couplings/dipole matrix elements obtained from cippres_fano/cippres_dip
  END_DOC

  integer :: npt
  integer, parameter :: QR_K = 8 !selected_real_kind (32)
  real (kind=QR_K), allocatable, dimension(:) :: e, g

  integer :: nmin, nmax ! according to Mueller-Plathe and Diercksen Stieltjes is inaccurate for n>=15
  real (kind=QR_K), allocatable, dimension(:) :: gord
  real (kind=QR_K), allocatable, dimension(:,:) ::  e1, g1
  real (kind=QR_K), allocatable, dimension(:) :: eallord,gallord
  real (kind=QR_K) :: shift1 
  integer :: imax1, imax2
  integer :: inmax, ishift
  integer :: i, j, k, ichan, kord, n_ord_av
  integer :: exit_cycle
  character(len=60) :: fname

! Tsveta
  real (kind=QR_K) :: g_
  real (kind=QR_K) :: temp, gav, stadev, gtot

  character(len=lenmax) :: finput
  integer :: ilen, jlen
  logical :: file_e

  finput='stj_input'
  
  if(i_stj_job==1) then ! Fano case

      if(ifanosta==0) then
        print*, "Please set ifanosta (the initial state)"
        print*, "qp set cippres ifanosta X "
        stop
      endif

      open(unit=10,file=finput)   
      write(10,*) i_stj_job
      write(10,*) n_csf_cippres(ici2)
      do i = 1, n_csf_cippres(ici2) 
        write(10,'(100(e24.16,1X))') e_couplings_cippres(i,ifanosta), twoe_couplings_cippres(i,ifanosta)
      enddo
      close(10)

  elseif(i_stj_job==2) then ! Dipole case

      if(idipsta==0) then
        print*, "Please set idipsta (the initial state)"
        print*, "qp set cippres idipsta X "
        stop
      endif

!      open(unit=10,file=finput)   
!      write(10,*) i_stj_job
!      write(10,*) n_csf_cippres(ici2)
!      do i = 1, n_csf_cippres(ici2) 
!        write(10,*), edip_couplings_cippres(i,idipsta), dip_couplings_cippres(i,idipsta)
!      enddo
!      close(10)


  else
   print*, "Please set up i_stj_job correctly (1==Fano and 2==Dipole)"
   stop
  endif
  
  call system('$QP_ROOT/plugins/qp2_plugins_nsisourat/cippres/libs/stieltjes/stieltjes < '//trim(finput))

  open(unit=10,file='stieltjes.info.txt')
   read(10,*)nmin,nmax,shift1
   nmin = nmin + 2
  close(10)

  allocate(e1(nmax,nmax), g1(nmax,nmax),gord(nmax))       
  allocate(eallord(nmax**2),gallord(nmax**2))
  e1(:,:) = 0d0
  g1(:,:) = 0d0
  do i = nmin, nmax
    write(fname, '(A16,I2,A4)')"stieltjes.order.",i,".txt"
    open(238,file=fname)
    k = 0
    do  j = nmin-1, i
       k = k + 1
       read(238,*)e1(k,i),g1(k,i)
    enddo
    close(238)
  enddo

  if(i_stj_job==1) then

       k = 1
       gord(:) = 0d0
       n_ord_av = 0
       do imax1 = nmin, nmax
         k = k + 1
         write(*,*)"#Interp at order", imax1
         if(e1(1,imax1)<0d0) then 
           call interp(e1(1:k,imax1),g1(1:k,imax1),k,0d0,g_)  
           n_ord_av += 1
           g_ = 2.0d0*pi*g_
           gord(imax1) = g_
           write(*,'(I3,A3,F23.15,A)')imax1," ",g_*27211,' in meV'
         else
           write(*,*)"Wrong energy range, Stieltjes order might be too low"
         endif
       
       end do

       gav = sum(gord(:))/n_ord_av
       stadev = sqrt(sum( (gord(:)-gav)**2) )/n_ord_av

       write(*,'(a)')""
       write(*,'(a)')"TOTAL RATE FROM SCHEME 1"
       write(*,'(2(f20.12,1X),a)')gav,stadev, 'in au'
       write(*,'(2(f20.12,1X),a)')gav*27210,stadev*27210, 'in meV'
       
       kord = 0
       k = 1
       do i = nmin, nmax
         k = k + 1
         do  j = 1, k
            kord = kord + 1
            eallord(kord) = e1(j,i)
            gallord(kord) = g1(j,i)
!            write(*,*)i,j,e1(j,i),kord,eallord(kord)
         enddo
       enddo

       call sort2(kord,eallord,gallord)

       open(unit=222,file='c.allorder.txt')
       do i = 1, kord
         write(222,'(2(f20.12,1X))')eallord(i),gallord(i)*27211*2d0*pi
       enddo
       close(222)
       call interp(eallord,gallord,kord,0d0,g_)  
       write(*,'(a)')"TOTAL RATE FROM SCHEME 2"
       write(*,'(1(f20.12,1X),a)')g_*2d0*pi, 'in au'
       write(*,'(1(f20.12,1X),a)')g_*27210*2d0*pi,'in meV'

   elseif(i_stj_job==2) then

  do i = nmin, nmax
    write(fname, '(A22,I2,A4)')"photoionization.order.",i,".txt"
    open(238,file=fname)
    write(238,*)"#Photoionization cross sections in Mb, Energy in eV"
    k = 0
    do  j = nmin-1, i
       k = k + 1
       write(238,*)e1(k,i)*27.211d0,0.428*g1(k,i)*e1(k,i)
    enddo
    close(238)
  enddo


   endif

!  deallocate(eallord,gallord)
!  deallocate(gord,e,g)
!  deallocate(e1,g1)

end program cippres_stieltjes

subroutine imaging(npt,e_point,g_point,nmax,maxord,eorder,gorder)
implicit none

 integer :: nmax

 integer, parameter :: QR_K = selected_real_kind (32) ! use quadruple precision
 double precision, parameter :: overmax=1d0 ! set arbitrary to 1 to determine the maximum order of pol.

 integer, intent(in) :: npt
 real (kind=QR_K), dimension(npt), intent(in) :: e_point, g_point

 real (kind=QR_K), dimension(0:nmax,npt) :: qpol
 real (kind=QR_K), dimension(nmax) :: acoef
 real (kind=QR_K), dimension(0:nmax) :: bcoef

 real (kind=QR_K), dimension(nmax) :: diag, offdiag
 real (kind=QR_K), dimension(nmax,nmax) :: abvec
 real (kind=QR_K) :: asum, bprod, qnorm, qoverlap
 integer :: iord, ierr, maxord, min, max

 real (kind=QR_K), dimension(nmax) :: enew, gnew
 real (kind=QR_K), dimension(nmax,nmax) :: eorder, gorder

 integer :: i, j

! initiate the recursive computation of the a,b coefficients and the orthogonal 
! polynomials according to (3.3.20-23) of Mueller-Plathe & Dierksen (1990)
       bcoef(0)=0.q0
       acoef(1)=0.q0
       do i=1,npt
          bcoef(0)=bcoef(0)+g_point(i)
          acoef(1)=acoef(1)+g_point(i)/e_point(i)
       end do
       acoef(1)=acoef(1)/bcoef(0)

       do i=1,npt
          qpol(0,i)=1.q0
          qpol(1,i)=1.q0/e_point(i)-acoef(1)
       end do

       bcoef(1)=0.q0
       acoef(2)=0.q0
       do i=1,npt
          bcoef(1)=bcoef(1)+qpol(1,i)*g_point(i)/e_point(i)
          acoef(2)=acoef(2)+qpol(1,i)*g_point(i)/(e_point(i)**2)
       end do
       bcoef(1)= bcoef(1)/bcoef(0)
       acoef(2)=acoef(2)/(bcoef(0)*bcoef(1))-acoef(1)

! calculate the higher-order coefficients and polynomials recursively
! up to the (NMAX-1)th order (total of NMAX polynomials)

       asum=acoef(1)
       do i=3,nmax

          asum=asum+acoef(i-1)

          do j=1,npt
             qpol(i-1,j)=(1.q0/e_point(j)-acoef(i-1))*qpol(i-2,j)-bcoef(i-2)*qpol(i-3,j)
          end do

          bprod=bcoef(0)
          do j=1,i-2
             bprod=bprod*bcoef(j)
          end do

          bcoef(i-1)=0.q0
          do j=1,npt
             bcoef(i-1)=bcoef(i-1)+qpol(i-1,j)*g_point(j)/(e_point(j)**(i-1))
          end do
          bcoef(i-1)=bcoef(i-1)/bprod

          bprod=bprod*bcoef(i-1)

          acoef(i)=0.q0
          do j=1,npt
             acoef(i)=acoef(i)+qpol(i-1,j)*g_point(j)/(e_point(j)**i)
          end do
          acoef(i)=acoef(i)/bprod-asum

       end do

! calculate the nmax-th order polynomial just for the orthogonality check 
       do j=1,npt
          qpol(nmax,j)=(1.q0/e_point(j)-acoef(nmax))*qpol(nmax-1,j)-bcoef(nmax-1)*qpol(nmax-2,j)
       end do

! check the orthogonality of the polynomials to define the maximal approximation order 
! if the orthogonality is preserved for all orders, MAXORD is set to NMAX
       maxord=nmax
       qnorm=bcoef(0)
       do i=1,nmax
          qnorm=0.q0
          qoverlap=0.q0
          do j=1,npt
             qnorm=qnorm+qpol(i,j)**2*g_point(j)
             qoverlap=qoverlap+qpol(i,j)*qpol(i-1,j)*g_point(j)
          end do
          if (qabs(qoverlap).lt.1.q-50) qoverlap=1.q-50
          
          if (qnorm/qabs(qoverlap).le.overmax) then
! MAXORD=I-1 is appropriate since the polynomial failing 
! the orthogonality check should not be used
             maxord=i-1
             go to 10
          end if
       end do
 10    continue
!!       maxord=nmax

! look how many Stieltjes orders are available
       if (maxord.lt.5) then
          min=maxord
          max=maxord
          print*, '***WARNING*** Stieltjes:'
          print*, ' only very low-order approximation is available'
          print*, ' MAXORD=',maxord
       else
          min=5
          max=maxord
          print*, ' MAXORD=',maxord
       end if

! perform the gamma calculation using the successive approximations 
! n=5,...,nmax

   do iord=5,maxord

     write(*,*)"Performs Stieltjes at order",iord 

! fill the coefficients matrix
       do i=1,iord
          diag(i)=acoef(i)
!          write(*,*)"diag",i,diag(i)
       end do
       do i=2,iord
          offdiag(i)=-qsqrt(bcoef(i-1))
!          write(*,*)"offdiag",i,offdiag(i)
       end do

! diagonalize the coefficients matrix
! initialize the arrays
       do i=1,nmax
          do j=1,nmax
             abvec(i,j)=0.q0
          end do
          abvec(i,i)=1.q0
       end do
       call tql2(nmax,iord,diag,offdiag,abvec,ierr)
       if (ierr.ne.0) then
          print*, '***WARNING*** Stieltjes:'
          print*, ' the eigenvalue no. ',ierr,' failed to converge'
       end if

! fill the Stieltjes energy and gamma arrays
! note that the eigenvalues are inverse energies and are given in ascending order 
       do i=1,iord
!          print*,diag(iord+1-i),abvec(1,iord+1-i)**2
          enew(i)=1.q0/diag(iord+1-i)
          gnew(i)=bcoef(0)*abvec(1,iord+1-i)**2
       end do

       call eigsrtnico(enew,gnew,iord,nmax)

! calculate the gamma's by simple numerical differentiation at the middle 
! point of each [ENEW(I),ENEW(I+1)] interval
       do i=1,iord-1
          eorder(i,iord)=0.5d0*(enew(i)+enew(i+1))
          gorder(i,iord)=0.5d0*(gnew(i+1)+gnew(i))/(enew(i+1)-enew(i))
       end do

   enddo ! loop over iord  


end subroutine

