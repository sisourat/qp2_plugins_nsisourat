program cippres_flighttime
  use general
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! cippres_timedep_entropy computes the von Neumann entropy of a time dependent wavefunction
  END_DOC

  character(len=lenmax) :: finput
  integer :: ilen, jlen
  logical :: file_e

  integer :: i, j, ista, ib, it, k, l
  complex (kind=8), dimension(:,:), allocatable :: cmattmp
  double precision, dimension(:,:), allocatable :: dys_norm
  complex (kind=8), dimension(:), allocatable :: tdcoef
  double precision, dimension(:), allocatable :: rtdcoef, ctdcoef
  double precision, dimension(:), allocatable :: esta, norm, flighttime, probcation, tdprob

  double precision :: b, time, told, rt, ct, tinit
  double precision :: sip, pb
 
  integer :: irun, nb, nsta, ntime, nsta_t, ncation, i1, i2, i3, npseudo
   irun = 1
   print*,'irun=',irun
   
   open(unit=20,file='Dyson_norms.txt') 
      read(20,*)nsta,ncation,sip
     allocate(dys_norm(ncation,nsta),norm(nsta),flighttime(ncation),probcation(ncation),tdprob(ncation))
     do i = 1, nsta
      read(20,*)ista
      read(20,*)(dys_norm(j,i),j=1,ncation),norm(i)
      dys_norm(:,i)=dys_norm(:,i)/norm(i)
     enddo
   close(20)

   allocate(cmattmp(nmo,nmo))

   open(unit=10,file='Psit_collision.out') 
     read(10,*)nb,nsta_t,ntime

     if(nsta/=nsta_t) then
      print*,"error in nb of nsta",nsta, nsta_t
      stop
     endif
     allocate(tdcoef(nsta),rtdcoef(nsta),ctdcoef(nsta))
     allocate(esta(nsta))
     do i = 1, nsta
      read(10,*)esta(i)
     enddo

     do ib = 1,  nb 
      flighttime(:)=0d0
      probcation(:)=0d0

        read(10,*)b,told,(rtdcoef(ista), ctdcoef(ista),ista=1,nsta)
        tinit=told
        told=told-tinit
      do it = 2, ntime-1
        read(10,*)b,time,(rtdcoef(ista), ctdcoef(ista),ista=1,nsta)
        time=time-tinit
        tdcoef(:) = dcmplx(rtdcoef,ctdcoef)+0d-15
       
        
        npseudo = 0 
        tdprob(:) = 0d0
        do i = 1, nsta
          if(esta(i)>sip) then
           npseudo = npseudo + 1
           do j = 1, ncation
            pb = conjg(tdcoef(i))*tdcoef(i)*dys_norm(j,i)
            if(pb>1d-15) then
              flighttime(j) = flighttime(j) + pb*abs(time-told)*(told+time)*0.5d0
              probcation(j) = probcation(j) + pb*abs(time-told)*0.5d0
            endif
              tdprob(j) = tdprob(j) + pb
           enddo
          endif
        enddo
        told=time

         do j = 1, ncation
           write(ib*100+j,'(100(f15.9,1X))')time,tdprob(j)
         enddo

      enddo
        read(10,*)b,time,(rtdcoef(ista), ctdcoef(ista),ista=1,nsta)
        tdcoef(:) = dcmplx(rtdcoef,ctdcoef)+0d-15

!          write(*,'(100(f25.20,1X))')probcation(:)
          flighttime(:) = flighttime(:)/probcation(:)
          write(13,'(100(f15.9,1X))')b,flighttime(:)

!          write(*,'(100(f15.5,1X))')b,(flighttime(:)-told/1d0)*240d0
          write(14,'(100(f15.5,1X))')b,(flighttime(:)-(flighttime(1)+flighttime(2)+flighttime(3))/3d0)*0.024d0
!          write(*,'(100(f25.20,1X))')
     enddo
         
   close(10)

  deallocate(esta,tdcoef,ctdcoef,rtdcoef,dys_norm,norm,flighttime,probcation,tdprob)

end program cippres_flighttime
