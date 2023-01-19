program cippres_timedep_entropy
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
  double precision, dimension(:,:,:,:), allocatable :: densmatsta
  complex (kind=8), dimension(:,:), allocatable :: densmatvec
  double precision, dimension(:,:), allocatable :: densmat
  complex (kind=8), dimension(:,:), allocatable :: cmattmp
  complex (kind=8), dimension(:), allocatable :: tdcoef
  double precision, dimension(:), allocatable :: rtdcoef, ctdcoef
  double precision, dimension(:), allocatable :: densmatval, esta

  double precision :: lin_entrop, vNentrop, b, time, rt, ct
  double precision :: lin_entrop_approx1, vNentrop_approx1
  double precision :: lin_entrop_approx2, vNentrop_approx2, dens
 
   integer :: irun, nb, nsta, ntime, nsta_t, i1, i2, i3
   irun = 1
   print*,'irun=',irun

   open(unit=10,file='Density_matrices.txt') 
     read(10,*)nmo, nsta
     write(*,*)nmo,nsta
     allocate(densmatsta(nmo,nmo,nsta,nsta))
     do i = 1, nsta
       do j = 1, nsta
            read(10,*)ib,it
        do k = 1, nmo
          do l = 1, nmo
            write(*,*)l,k,j,i
            read(10,*)densmatsta(l,k,j,i)
          enddo
        enddo
      enddo
     enddo
   close(10)

   allocate(cmattmp(nmo,nmo))
   allocate(densmat(nmo,nmo))
   allocate(densmatvec(nmo,nmo))
   allocate(densmatval(nmo))

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
     do ib = 1, 1!, nb 
      do it = 1, ntime
        read(10,*)b,time,(rtdcoef(ista), ctdcoef(ista),ista=1,nsta)
        tdcoef(:) = dcmplx(rtdcoef,ctdcoef)

          cmattmp(:,:) = 0d0
        do i = 1, nsta
         do j = 1, nsta
                cmattmp(:,:) = cmattmp(:,:) + densmatsta(:,:,j,i)*dconjg(tdcoef(j))*tdcoef(i)*exp(-dcmplx(0d0,esta(i)-esta(j))*time)
         enddo
!                write(*,*)'nico',i,i,dconjg(tdcoef(i))*tdcoef(i)
        enddo

        vNentrop_approx1 = 0d0
        lin_entrop_approx1 = 1d0
        do i = 1, nmo
          dens = cdabs(cmattmp(i,i)) + 1d-15
!          write(*,*)time,i,dens
          vNentrop_approx1 = vNentrop_approx1 - abs(dens) * log(dens)/log(2d0)
          lin_entrop_approx1 = lin_entrop_approx1 - abs(dens) **2
        enddo
!          write(40,'(100(f15.9,1X))')time,cmattmp(1,1),cmattmp(2,2),cmattmp(1,2),cmattmp(2,1)

!        densmat(:,:) = real(cmattmp(:,:))!**2
        call lapack_diagc(densmatval,densmatvec,cmattmp,nmo) 
        densmatval(:) = densmatval(:) + 1d-15
        vNentrop = 0d0
        lin_entrop = 1d0
        do i = 1, nmo
!          write(*,*)time,i,densmatval(i)
          vNentrop = vNentrop - abs(densmatval(i)) * log(abs(densmatval(i)))/log(2d0)
          lin_entrop = lin_entrop - abs(densmatval(i)) **2
        enddo

          cmattmp(:,:) = 0d0
          densmatval(:) = 0d0
        do i = 1, nsta
                cmattmp(:,:) = cmattmp(:,:) + densmatsta(:,:,i,i)*dconjg(tdcoef(i))*tdcoef(i)
        enddo
          call lapack_diagc(densmatval,densmatvec,cmattmp,nmo) 
          densmatval(:) = densmatval(:) + 1d-15
          vNentrop_approx2 = 0d0
          lin_entrop_approx2 = 1d0
          do i = 1, nmo
            vNentrop_approx2 = vNentrop_approx2 - abs(densmatval(i)) * log(abs(densmatval(i)))/log(2d0)
            lin_entrop_approx2 = lin_entrop_approx2 - abs(densmatval(i)) **2
          enddo

          write(*,'(100(f15.9,1X))')time,vNentrop,vNentrop_approx1,vNentrop_approx2,lin_entrop,lin_entrop_approx1,lin_entrop_approx2
!          write(*,'(100(f15.9,1X))')time,vNentrop,lin_entrop
    enddo
   enddo
         
   close(10)

!  do ista = 1, N_states
!   print*,""
!  print*,"DENSITY MATRIX", ista
!  print*,""
!  write(30,*)ista
!  do i = 1, mo_num
!   do j = 1, mo_num
!       print*,j,i,one_e_dm_mo_beta(j,i,ista)+one_e_dm_mo_alpha(j,i,ista)
!      densmat(j,i) = 0.5d0*(one_e_dm_mo_beta(j,i,ista)+one_e_dm_mo_alpha(j,i,ista))
!    enddo
!   enddo

!   call lapack_diagd(densmatval,densmatvec,densmat,mo_num,mo_num) 
!      vNentrop = 0d0
!    do i = 1, mo_num
!      write(30,*)densmatval(i)
!      vNentrop = vNentrop - abs(densmatval(i)) * log(abs(densmatval(i)))/log(2d0)
!    enddo
!   write(20,*)eigvalues_cippres(ista,1),vNentrop
!  enddo
!  deallocate(densmat,densmatvec,densmatval)
!

  deallocate(densmat,densmatvec,densmatval,densmatsta,esta,tdcoef)

end program cippres_timedep_entropy

subroutine lapack_diagc(eigvalues,eigvectors,H,n)
  implicit none
  BEGIN_DOC
  ! Diagonalize complex Hermitian matrix H
  !
  ! H is untouched between input and ouptut
  !
  ! eigevalues(i) = ith lowest eigenvalue of the H matrix
  !
  ! eigvectors(i,j) = <i|psi_j> where i is the basis function and psi_j is the j th eigenvector
  !
  END_DOC
  integer, intent(in)            :: n
  COMPLEX*16, intent(in)   :: H(n,n)

  double precision, intent(out)  :: eigvectors(n,n)
  double precision, intent(out)  :: eigvalues(n)

  DOUBLE PRECISION :: W( N ), RWORK( 3*N-2 )
  integer          :: lda, lwmax, info, i,j,l,k, lwork
  parameter ( LWMAX = 1000 )
  COMPLEX*16       :: A( N, N ), WORK( LWMAX )

  A=H
  lda = n
  LWORK = -1
  CALL ZHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,INFO)
  LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )

  CALL ZHEEV( 'Vectors', 'Lower', N, A, LDA, W, WORK, LWORK, RWORK,INFO)
!
!     Check for convergence.
!
      IF( INFO.GT.0 ) THEN
         WRITE(*,*)'The algorithm failed to compute eigenvalues.'
         STOP
      END IF
!
!     Print eigenvalues.
!
!      CALL PRINT_RMATRIX( 'Eigenvalues', 1, N, W, 1 )

  eigvectors = 0.d0
  eigvalues = 0.d0
  do j = 1, n
    eigvalues(j) = w(j)
!    write(*,*)w(j),A(j,j),H(j,j)
    do i = 1, n
      eigvectors(i,j) = A(i,j)
    enddo
  enddo
end

