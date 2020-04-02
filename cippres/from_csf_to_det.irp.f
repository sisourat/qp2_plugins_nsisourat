subroutine myprint_det(string,Nint,output) 
  use bitmasks
  implicit none
  BEGIN_DOC
  ! Subroutine to print the content of a determinant using the '+-' notation
  END_DOC
  integer, intent(in)            :: Nint
  integer(bit_kind), intent(in)  :: string(Nint,2)
  character*(2048)                :: output(2)

  call bitstring_to_str( output(1), string(1,1), Nint )
  call bitstring_to_str( output(2), string(1,2), Nint )
!  print *,  trim(output(1))
!  print *,  trim(output(2))

end
 

subroutine from_csf_to_det(irun)
 implicit none
 use bitmasks ! you need to include the bitmasks_module.f90 features
 integer, intent(in) :: irun
 integer :: i,j,k,idet_tmp,istate

 character*(2048)   :: output(2)

 print*,'*************************************'
 print*,'n_det_max_csf, n_csf_max',n_det_max_csf, n_csf_max
 print*,'*************************************'
 idet_tmp = 0
 do i = 1, n_csf_cippres(irun) ! first loop on the csf of the space ispace 
   do k = 1, n_det_csf_cippres(i,irun) ! then on the determinants belonging to the ith CSF of space ispace
    idet_tmp += 1
   enddo 
  enddo

 N_det = idet_tmp 
 N_states = n_csf_max
 touch N_det N_states
 print*,'N_det,psi_det_size',N_det,psi_det_size
 idet_tmp = 0
 do i = 1, n_csf_cippres(irun) ! first loop on the csf of the space ispace 
   do k = 1, n_det_csf_cippres(i,irun) ! then on the determinants belonging to the ith CSF of space ispace
    idet_tmp += 1
    psi_det(:,:,idet_tmp) = csf_basis(:,:,k,i,irun) 
    do istate = 1, n_csf_max
     psi_coef(idet_tmp,istate) = coef_det_csf_basis(k,i,irun) * eigvectors_cippres(i,istate,irun)
    enddo
   enddo 
  enddo
 print*,'idet_tmp          = ',idet_tmp
 print*,'N_det,N_states    = ',N_det,N_states
 touch psi_det psi_coef 
 print*,'N_det_alpha_unique',N_det_alpha_unique
 print*,'N_det_beta_unique ',N_det_beta_unique
 print*,'*N_det,N_states    = ',N_det,N_states

 double precision, allocatable :: coef_alpha_beta(:,:,:)
 integer(bit_kind), allocatable :: psi_alpha_uniq_tmp(:,:), psi_beta_uniq_tmp(:,:)
 allocate(psi_alpha_uniq_tmp(N_int,N_det_alpha_unique),psi_beta_uniq_tmp(N_int,N_det_beta_unique))

 psi_alpha_uniq_tmp(:,1:N_det_alpha_unique) = psi_det_alpha_unique(:,1:N_det_alpha_unique)
 psi_beta_uniq_tmp(:,1:N_det_beta_unique) = psi_det_beta_unique(:,1:N_det_beta_unique)

 allocate ( coef_alpha_beta(N_det_alpha_unique,N_det_beta_unique,N_states)  )
 
 integer :: get_index_in_psi_det_alpha_unique,get_index_in_psi_det_beta_unique
 integer :: n_alpha_tmp,n_beta_tmp
 n_alpha_tmp = N_det_alpha_unique
 n_beta_tmp = N_det_beta_unique
 coef_alpha_beta = 0.d0
 do k=1,N_det
   i = get_index_in_psi_det_alpha_unique(psi_det(1,1,k),N_int)
   ASSERT (i>0)
   ASSERT (i<=N_det_alpha_unique)

   j = get_index_in_psi_det_beta_unique (psi_det(1,2,k),N_int)
   ASSERT (j>0)
   ASSERT (j<=N_det_beta_unique)
   do istate = 1, N_states
    coef_alpha_beta(i,j,istate) += psi_coef(k,istate)
   enddo
 enddo
 print*,'N_det_alpha_unique',N_det_alpha_unique
 print*,'N_det_beta_unique ',N_det_beta_unique

 N_det = N_det_alpha_unique * N_det_beta_unique 
 touch N_det 
 idet_tmp = 0
 do j = 1, n_alpha_tmp 
  do i = 1, n_beta_tmp 
   idet_tmp += 1
   psi_det(:,1,idet_tmp) = psi_alpha_uniq_tmp(:,j)
   psi_det(:,2,idet_tmp) = psi_beta_uniq_tmp(:,i)
!   print*,'idet = ', idet_tmp
!   call print_det(psi_det(1,1,idet_tmp),N_int)
!   print*,''
   do istate = 1, N_states
    psi_coef(idet_tmp,istate) = coef_alpha_beta(j,i,istate)
   enddo
  enddo
 enddo

 touch psi_det psi_coef 

  open(unit=10,file='cistates_det.txt')
  write(10,*),mo_num,N_det,N_states
  do i = 1, N_det
   write(10,*),i
    call myprint_det(psi_det(1,1,i),N_int,output)
    write(10,*),  trim(output(1))
    write(10,*),  trim(output(2))
  enddo
  do istate = 1, N_states
     write(10,*),istate, eigvalues_cippres(istate,irun), N_states
    do i = 1, N_det
     write(10,*), psi_coef(i,istate)
    enddo  
  enddo
  close(10)
  stop
! do istate = 1, N_states
!  double precision :: accu,hij
!  accu = 0.d0
!  do i = 1, N_det
!   do j = 1, N_det
!!    print*,'i,j',i,j
!!    call print_det(psi_det(1,1,i),N_int)
!!    call print_det(psi_det(1,1,j),N_int)
!    call i_H_j(psi_det(1,1,i),psi_det(1,1,j),N_int,hij)
!    accu += hij * psi_coef(j,istate) * psi_coef(i,istate) 
!   enddo
!  enddo
!  print*,'accu              = ',accu+nuclear_repulsion
!  print*,'eigvalues_cippres = ',eigvalues_cippres(istate,irun)
! enddo
 
end

